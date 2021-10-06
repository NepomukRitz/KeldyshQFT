#ifndef QUEUE_HPP
#define QUEUE_HPP

#pragma once
#include <map>

template <class T>
struct Task {
  T d;
  typename T::f_type f;
  typename T::value_type val;
  Base<typename T::value_type> error;
  std::size_t idx;

  bool operator<(const Task& rhs) { return error > rhs.error; }
};

template <typename T>
std::ostream& operator<<(std::ostream& oss_, const Task<T>& rhs) {
  std::ostringstream oss;
  oss << "(" << rhs.d.size() << "," << rhs.val << "," << rhs.error << ")";
  oss_ << oss.str();
  return oss_;
}

template <class T>
class Queue {
 public:
  bool empty() { return tasks.empty(); }

  // pop version that takes iterators to fill
  // put version that takes iterators to fill

  Task<T> pop() {
    std::pop_heap(tasks.begin(), tasks.end());
    auto& task = tasks.back();
    tasks.pop_back();
    return task;
  }

  void put(Task<T> task) {
    if (task.d.size() < 1e-16 || task.error / task.d.size() < 1e-14) {
      // std::cout << "task too small" << task << "\n";
      small_tasks.push_back(task);
      return;
    }
    tasks.push_back(task);
    std::push_heap(tasks.begin(), tasks.end());
  }

  std::pair<Base<typename T::value_type>,
            std::map<std::size_t, typename T::value_type>>
  sum() {
    // std::cout << "size queue: " << tasks.size() << " size small queue " <<
    // small_tasks.size() << "\n";
    std::map<std::size_t, typename T::value_type> valmap;
    Base<typename T::value_type> error = 0;
    typename T::value_type value = 0;

    for (auto task : small_tasks) {
      error += task.error / task.d.size();
      valmap[task.idx] += task.val;
    }

    for (auto task : tasks) {
      error += task.error / task.d.size();
      valmap[task.idx] += task.val;
    }
    return {error, valmap};
  }

 private:
  std::vector<Task<T>> tasks;
  std::vector<Task<T>> small_tasks;
  // static const array_size = 5;
};

# endif //QUEUE_HPP