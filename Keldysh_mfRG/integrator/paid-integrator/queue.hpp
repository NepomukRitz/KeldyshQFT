// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:
# ifndef  PAID_INTEGRATOR_QUEUE_HPP
# define PAID_INTEGRATOR_QUEUE_HPP

#pragma once

#include <deque>
#include <queue>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp_stub.h"
#endif

#include "task.hpp"
#include "types.hpp"

namespace paid {

// template as any T, with T::value_type?
template <std::size_t Dim, typename F, typename T, typename Key>
class Queue {
 public:
  Queue(double min_size, double min_error)
      : min_size_(min_size), min_error_(min_error), tasks_(), small_tasks_() {}

  bool empty() const { return tasks_.empty(); }
  std::size_t size() const { return tasks_.size(); } // how many "tasks"

  void pop() noexcept { // removes first element from "tasks_" and keeps heap structure
    assert(!tasks_.empty());
    std::pop_heap(tasks_.begin(), tasks_.end());
    tasks_.pop_back();
  }

  Task<Dim, F, T, Key>& top() noexcept { return tasks_.front(); } // first element of "tasks"

  // if domain size or error small: "task" is moved to "small_tasks", else: "task" is moved to "tasks" and keeps heap structure
  void push(Task<Dim, F, T, Key>& task) noexcept {
    if (task.d.size() < min_size_ ||
        task.error < min_error_) {  // TODO min value?
      push_small(task);
    } else {
      tasks_.push_back(task);
      std::push_heap(tasks_.begin(), tasks_.end());
    }
  }

  // sums up errors and value maps of "tasks" and "small_tasks", returns "{error, valmap}"
  std::pair<error_type, std::unordered_map<Key, T>> sum() const noexcept {
    std::unordered_map<Key, T> valmap;
    error_type error = 0;

    for (auto task : small_tasks_) {
      error += task.error;
      valmap[task.idx] += task.val;
    }

    for (auto task : tasks_) {
      error += task.error;
      valmap[task.idx] += task.val;
    }
    return {error, valmap};
  }

 private:
  void push_small(Task<Dim, F, T, Key>& task) { small_tasks_.push_back(task); }

  const double min_size_;
  const double min_error_;

  // std::priority_queue does not have .cbegin()
  std::deque<Task<Dim, F, T, Key>> tasks_; // "deque" is a container with fast insertion and deletion of first and last element
  std::deque<Task<Dim, F, T, Key>> small_tasks_;
};

}  // namespace paid

#endif