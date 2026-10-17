#pragma once

#include <chrono>
#include <cstdint>

namespace adsim {

  // Wall-clock breakdown of a solve, used by the benchmark runner to see where
  // time goes (residual vs IDA internals vs logging vs CSS). Collection is off by
  // default; when `enabled` is false every timer below is a no-op. Nothing here
  // feeds back into the numerics.
  struct SolverProfile {
    // Error-test failures are counted over this many internal steps after the
    // start of each cycle step, to expose the cost of crossing BC discontinuities.
    static constexpr long kBoundaryWindow = 20;

    bool enabled = false;

    std::int64_t residual_ns = 0;   // all residual calls, incl. DQ Jacobian columns
    long residual_calls = 0;
    std::int64_t one_step_ns = 0;   // total time inside IDASolve(ONE_STEP)
    std::int64_t log_ns = 0;        // EquationLog/BoundaryLog recording during the solve
    std::int64_t css_ns = 0;        // CSS check + acceleration

    long cycle_steps = 0;              // cycle steps integrated
    long steps_at_max_step = 0;        // accepted internal steps limited by max_step
    long boundary_err_test_fails = 0;  // see kBoundaryWindow

    void Reset() {
      const bool keep = enabled;
      *this = SolverProfile{};
      enabled = keep;
    }
  };

  // Adds the elapsed time of its scope to `target` when `on` is true.
  class ProfileTimer {
  public:
    ProfileTimer(bool on, std::int64_t& target) : on_(on), target_(target) {
      if (on_) start_ = std::chrono::steady_clock::now();
    }
    ~ProfileTimer() {
      if (on_) {
        target_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                       std::chrono::steady_clock::now() - start_).count();
      }
    }
    ProfileTimer(const ProfileTimer&) = delete;
    ProfileTimer& operator=(const ProfileTimer&) = delete;

  private:
    bool on_;
    std::int64_t& target_;
    std::chrono::steady_clock::time_point start_{};
  };

} // end adsim namespace
