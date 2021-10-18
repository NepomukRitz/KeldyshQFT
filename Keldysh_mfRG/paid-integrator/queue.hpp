// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

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
template <std::size_t N, typename F, typename T, typename Key>
class Queue {
 public:
  Queue(double min_size, double min_error)
      : min_size_(min_size), min_error_(min_error), tasks_(), small_tasks_() {}

  bool empty() const { return tasks_.empty(); }
  std::size_t size() const { return tasks_.size(); } // how many tasks

  void pop() noexcept { // removes first element from "tasks_"
    assert(!tasks_.empty());
    std::pop_heap(tasks_.begin(), tasks_.end());
    tasks_.pop_back();
  }

  Task<N, F, T>& top() noexcept { return tasks_.front(); } // first element of "tasks"

  // "task" is moved to "small_tasks_" for small size/error or moved to the beginning of "tasks"
  void push(Task<N, F, T>&& task) noexcept {
    if (task.d.size() < min_size_ ||
        task.error < min_error_) {  // TODO min value?
      push_small(std::move(task));
    } else {
      tasks_.push_back(task);
      std::push_heap(tasks_.begin(), tasks_.end());
    }
  }

  // sums up errors and value maps of tasks and small_tasks, returns {error, valmap}
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
  void push_small(Task<N, F, T>&& task) { small_tasks_.push_back(task); }

  const double min_size_;
  const double min_error_;

  // std::priority_queue does not have .cbegin()
  std::deque<Task<N, F, T>> tasks_; // "deque" is a container with fast insertion and deletion of first and last element
  std::deque<Task<N, F, T>> small_tasks_;
};

}  // namespace paid
