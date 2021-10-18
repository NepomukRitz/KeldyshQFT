// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <cmath>
#include <complex>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp_stub.h"
#endif

#include "config.hpp"
#include "domain.hpp"
#include "invoker.hpp"
#include "queue.hpp"
#include "rule.hpp"
#include "types.hpp"

namespace paid {

template <typename T, typename Key>
struct PAIDOutput {
 public:
  T operator[](std::size_t key) const { // gives "value" of a certain "key"
    auto search = valmap_.find(key);
    if (search == valmap_.cend()) throw new std::out_of_range("key not found");
    return search->second;
  };

  std::size_t iterations_;
  std::size_t function_evaluations_;
  error_type error_estimate_;
  std::unordered_map<Key, T> valmap_;
  TerminationReason termination_reason_;
};

template <std::size_t N, typename F, typename Key>
struct PAIDInput {
 public:
  Domain<N> d;
  F f;
  Key idx;
};

template <std::size_t N, typename F, typename T, typename Key,
          typename... Args>
class PAID {
 public:
  PAID() = delete;
  explicit PAID(PAIDConfig& config) noexcept
      : q_(config.min_size,
           config.min_error),  // too_small_area and too_small_error
        rule_(config.order),
        config_(config),
        valid_(true) {
    static_assert(N == sizeof...(Args) || sizeof...(Args) == 1,
                  "Wrong number of Arguments for F");
    using T2 = typename std::result_of<F(Args...)>::type;
    static_assert(std::is_same<T2, T>::value,
                  "F(Args...) does not return T type");
  }

  PAID(const PAID&) = default;               // Copy Constructor
  PAID(PAID&&) = default;                    // Move constructor
  PAID& operator=(const PAID&) & = default;  // Copy assignment
  PAID& operator=(PAID&&) & = default;       // Move assignment
  ~PAID() {}

  PAIDOutput<T, Key> solve(const std::vector<PAIDInput<N, F, Key>>& inputs) {
    assert(valid_ == true);  // TODO
    valid_ = false;

    const int num_threads = omp_get_max_threads();
    std::vector<error_type> cur_error_in(num_threads, 0.0);
    error_type error_estimate = 0;
    std::size_t tasks_in_flight = 0;
    bool done = false;
    TerminationReason termination_reason = TerminationReason::Unknown;
    std::size_t counter = 0;
    std::size_t correct_error_counter = 1;
    std::size_t fevals = 0;

#pragma omp parallel
    {
      // Private variables
      std::vector<Task<N, F, T>> in_queue;
      std::vector<Task<N, F, T>> out_queue;
      const int thread_id = omp_get_thread_num();
      std::size_t work_size = 0;
      error_type error_out = 0;

#pragma omp for schedule(static), nowait
      for (std::size_t i = 0; i < inputs.size(); ++i) {
        const auto& task = inputs[i];
        const auto af = task.d.getTransform();
        const auto res = rule_.apply<N, F, T, Args...>(task.f, af);

        out_queue.push_back(
            Task<N, F, T>{task.d, task.f, res.value, res.error, task.idx});

#pragma omp atomic
        fevals += res.fevals_;
      }

        // we ensure that all the inital tasks are in the queue before entering
        // the while loop
#pragma omp critical(tasks)
      {
        error_out = 0;
        for (auto& task : out_queue) {
          error_out += task.error;
          q_.push(task);
        }
        error_estimate += error_out;
      }
      out_queue.clear();
#pragma omp barrier

      while (!done) {
        assert(counter >= 0);
        assert(fevals >= 0);

#pragma omp critical(tasks)
        {
          for (auto& task : out_queue) {
            q_.push(task);
          }

          // instead of dealing with errors for tasks in flight---those that
          // are out of the queue for processing---we over-approximate:
          //
          // We do not subtract the error of the tasks in flight immediately
          // Instead we first add the errors for the newly generated tasks,
          // and then subtract the error for the old tasks. To put in another
          // way, each time we are done with a set of tasks, we update only
          // the difference

          if (config_.correct_error) {
            error_out = 0;
            for (auto& task : out_queue) {
              error_out += task.error;
            }

            error_estimate += error_out - cur_error_in[thread_id];
#pragma omp flush(error_estimate)
          }

          out_queue.clear();
          in_queue.clear();

          // now get new items
#pragma omp flush(tasks_in_flight)
          tasks_in_flight -= work_size;
          work_size = std::min(q_.size() / num_threads + 1,
                               config_.ntasks_per_iteration);
          work_size = std::min(q_.size(), work_size);
          for (std::size_t i = 0; i < work_size; ++i) {
            in_queue.push_back(q_.top());
            q_.pop();
          }
          tasks_in_flight += work_size;
#pragma omp flush(tasks_in_flight)

          if (config_.correct_error) {
            cur_error_in[thread_id] = 0;
            for (auto& task : in_queue) {
              cur_error_in[thread_id] += task.error;
            }
          }

            // check if done
#pragma omp master
          {
            if (work_size > 0) counter++;

            if (counter % config_.check_every_iteration == 0 ||
                counter < config_.check_below_iteration) {
              // Correct the error from the queue
              // We space the corrections in exponential increasing intervals
              if (config_.correct_error &&
                  static_cast<double>(counter) >
                      std::exp(static_cast<double>(correct_error_counter))) {
                correct_error_counter++;

                // To correct we sum over all errors in the queue
                auto res = q_.sum();
                error_type error_update =
                    res.first + std::accumulate(cur_error_in.begin(),
                                                cur_error_in.end(), 0.0);
                error_estimate = error_update;
#pragma omp flush(error_estimate)
              }
            }

            if (fevals >= config_.max_f_evals)
              termination_reason = TerminationReason::Max_Evals_Reached;

            if (config_.check_error(error_estimate))
              termination_reason = TerminationReason::Tolerance_Reached;

            if (counter >= config_.max_iterations)
              termination_reason = TerminationReason::Max_Iter_Reached;

            if (tasks_in_flight == 0 && q_.size() == 0) {
              termination_reason = TerminationReason::No_Work_Left;
            }

            if (termination_reason != TerminationReason::Unknown) {
            //  std::cout << "termination reason: " << termination_reason << "\n";
              done = true;
#pragma omp flush(done)
            }

          }  // omp master

        }  // omp critical(tasks)

        // Process the tasks
        for (const auto& task : in_queue) {
          const auto doms = task.d.split();
          for (const auto& dom : doms) {
            const auto af = dom.getTransform();
            const auto res = rule_.apply<N, F, T, Args...>(task.f, af);
            out_queue.push_back(
                Task<N, F, T>{dom, task.f, res.value, res.error, task.idx});
#pragma omp atomic
            fevals += res.fevals_;
          }
        }

      }  // End while

#pragma omp critical(tasks)
      {
        // We have still have work_size items, so we store them!
        for (auto& task : out_queue) {
          q_.push(task);
        }
        out_queue.clear();
      }

    }  // end parallel

    auto res = q_.sum();

    return PAIDOutput<T, Key>{counter, fevals, static_cast<double>(res.first),
                              res.second, termination_reason};
  }

 private:
  Queue<N, F, T, Key> q_;
  const class ClenshawCurtis rule_;
  const PAIDConfig config_;
  bool valid_;
};

}  // namespace paid
