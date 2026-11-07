/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/state_cache.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace adsim {

  CycleStateCache::CycleStateCache(size_t capacity) : capacity_(capacity == 0 ? 1 : capacity) {}

  CycleStateCache::CycleStateCache(std::vector<double> lower, std::vector<double> upper, size_t capacity)
      : lower_(std::move(lower)), upper_(std::move(upper)), capacity_(capacity == 0 ? 1 : capacity) {
    if (lower_.size() != upper_.size()) {
      throw std::invalid_argument("CycleStateCache: " + std::to_string(lower_.size()) + " lower bounds but " +
                                  std::to_string(upper_.size()) + " upper bounds");
    }
    for (size_t i = 0; i < lower_.size(); ++i) {
      if (!std::isfinite(lower_[i]) || !std::isfinite(upper_[i]) || !(lower_[i] < upper_[i])) {
        throw std::invalid_argument("CycleStateCache: variable " + std::to_string(i) +
                                    " needs finite bounds with lower < upper");
      }
    }
  }

  size_t CycleStateCache::Dimension() const {
    if (!lower_.empty()) return lower_.size();
    return entries_.empty() ? 0 : entries_.front().key.size();
  }

  void CycleStateCache::CheckSize(const std::vector<double>& x) const {
    const size_t n = Dimension();
    if (n != 0 && x.size() != n) {
      throw std::invalid_argument("CycleStateCache: a key needs " + std::to_string(n) + " values, got " +
                                  std::to_string(x.size()));
    }
  }

  std::vector<double> CycleStateCache::Normalise(const std::vector<double>& x) const {
    CheckSize(x);
    if (lower_.empty()) return x;
    std::vector<double> key(x.size());
    for (size_t i = 0; i < x.size(); ++i) key[i] = (x[i] - lower_[i]) / (upper_[i] - lower_[i]);
    return key;
  }

  void CycleStateCache::Add(const std::vector<double>& x, CycleState state) {
    if (state.Empty()) throw std::invalid_argument("CycleStateCache: cannot cache an empty state");
    std::vector<double> key = Normalise(x);
    for (double v : key) {
      if (!std::isfinite(v)) throw std::invalid_argument("CycleStateCache: a key must be finite");
    }
    for (auto& e : entries_) {
      if (e.key == key) {
        e.state = std::move(state);
        e.serial = next_serial_++;
        return;
      }
    }
    if (entries_.size() >= capacity_) {
      size_t oldest = 0;
      for (size_t i = 1; i < entries_.size(); ++i)
        if (entries_[i].serial < entries_[oldest].serial) oldest = i;
      entries_.erase(entries_.begin() + static_cast<std::ptrdiff_t>(oldest));
    }
    entries_.push_back(Entry{std::move(key), std::move(state), next_serial_++});
  }

  const CycleState* CycleStateCache::Nearest(const std::vector<double>& x, double* distance) const {
    if (entries_.empty()) return nullptr;
    const std::vector<double> key = Normalise(x);
    const Entry* best = nullptr;
    double best_d2 = 0.0;
    for (const auto& e : entries_) {
      double d2 = 0.0;
      for (size_t i = 0; i < key.size(); ++i) {
        const double d = key[i] - e.key[i];
        d2 += d * d;
      }
      if (!best || d2 < best_d2 || (d2 == best_d2 && e.serial > best->serial)) {
        best = &e;
        best_d2 = d2;
      }
    }
    if (distance) *distance = std::sqrt(best_d2);
    return &best->state;
  }

  RunResult RunWithStateCache(Simulation& sim, CycleStateCache* cache, const std::vector<double>& x,
                              RunOptions options, CachedRunInfo* info) {
    if (options.initial_state) {
      throw std::invalid_argument("RunWithStateCache chooses the initial state: options.initial_state must be null");
    }
    CachedRunInfo local;
    CachedRunInfo& out = info ? *info : local;
    out = CachedRunInfo{};
    auto Note = [&](const std::string& text) {
      if (!out.note.empty()) out.note += "; ";
      out.note += text;
      if (options.on_message) options.on_message(LogLevel::LOW, text);
    };

    const CycleState* start = nullptr;
    double distance = std::numeric_limits<double>::quiet_NaN();
    if (cache && !cache->Empty()) {
      start = cache->Nearest(x, &distance);
      std::string reason;
      try {
        reason = sim.CheckState(*start);
      } catch (const std::exception& e) {
        reason = e.what();
      }
      if (!reason.empty()) {
        Note("warm start skipped, the cached state does not fit the deck: " + reason);
        out.skipped = true;
        start = nullptr;
      }
    }

    RunResult r;
    if (start) {
      // The cache may change in Add below; the run needs the state only
      // during the call.
      options.initial_state = start;
      out.warm = true;
      out.distance = distance;
      out.start_cycles = start->cycles;
      r = sim.Run(options);
      if (r.status == RunStatus::Failed) {
        Note("the warm-started run failed (" + r.message + "); repeated from the initial conditions");
        out.retried_cold = true;
        out.warm = false;
        out.distance = std::numeric_limits<double>::quiet_NaN();
        out.start_cycles = 0;
        options.initial_state = nullptr;
        const double failed_time = r.wall_time_s;
        r = sim.Run(options);
        r.wall_time_s += failed_time;
      }
    } else {
      r = sim.Run(options);
    }

    // Only a cyclic steady state is cached: the state at the end of a run that
    // stopped at the cycle limit depends on where that run started, so
    // starting others from it would change their results, not just their
    // cost. A deck without a CSS check therefore always runs cold.
    if (cache && r.Ok() && r.css_reached && sim.HasFinalState()) {
      try {
        cache->Add(x, sim.FinalState());
        out.cached = true;
      } catch (const std::exception& e) {
        Note(std::string("the final state was not cached: ") + e.what());
      }
    }
    return r;
  }

}  // namespace adsim
