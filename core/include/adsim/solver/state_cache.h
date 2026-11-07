/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// CycleStateCache: end-of-run cycle states of earlier evaluations, keyed by
// their decision vector, so that each new run of an optimisation or a study
// can warm-start from the state of the NEAREST point evaluated so far.
//
//   adsim::CycleStateCache cache(lower, upper);       // the variables' bounds
//   for (each candidate x) {
//     adsim::Simulation sim = adsim::Simulation::FromJson(deck);
//     ... set the decision variables to x ...
//     adsim::RunResult r = adsim::RunWithStateCache(sim, &cache, x);
//     ... read the KPIs ...                            // the state of x is now cached
//   }
//
// Distances are Euclidean in the variables scaled to [0, 1] by their bounds
// (or in the keys as given, for a cache made without bounds: a sensitivity
// study keys its runs by their probability coordinates, already in (0, 1)).
// Nearest, not the most recent: a global optimiser (MaxLIPO) jumps across the
// whole box early on, so the previous candidate is often far away.
//
// The cache is an object the caller owns and fills; nothing in it is global,
// so an ask/tell loop (a front end proposing candidates one at a time) can
// keep one per optimisation. It holds whole CycleStates (DAE vector and
// endpoint histories), so it keeps at most `capacity` of them and forgets the
// oldest first. It is not thread-safe: share one between threads only with
// your own lock, or give each worker its own.
//
// RunWithStateCache does the run: it starts from the nearest cached state
// through the ordinary warm-start path (RunOptions::initial_state, which
// remaps the state onto another grid and stretches its histories onto changed
// step times), skips the warm start with a note when the state does not fit
// the deck (a structural change: components, layers, steps or endpoints),
// falls back to a cold run when the warm-started run fails, and caches the
// final state of a run that reached cyclic steady state. Only CSS states are
// cached: a run that stops at the cycle limit ends in a state that depends on
// where it started, and a fixed-cycle study must give every run the same
// start. So a deck without a CSS check, or a study at a fixed number of cycles
// that never reaches CSS, runs every candidate cold, exactly as before.
//
// A warm start changes where the run stops, not the cyclic steady state it
// stops at: warm and cold results agree within the CSS tolerance's bound
// (docs/cyclic-steady-state.md). An optimiser should still re-run its optimum
// from the deck's initial conditions and report that value (OptimizeSimulation
// in optimisation.h does).

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"

namespace adsim {

  class CycleStateCache {
  public:
    static constexpr size_t kDefaultCapacity = 64;

    // A cache whose keys are used as given (already scaled).
    explicit CycleStateCache(size_t capacity = kDefaultCapacity);
    // A cache that scales each key to [0, 1] with these bounds. Throws
    // std::invalid_argument when the sizes differ or a bound pair is not
    // finite with lower < upper.
    CycleStateCache(std::vector<double> lower, std::vector<double> upper,
                    size_t capacity = kDefaultCapacity);

    // The key of `x`: (x - lower) / (upper - lower) per variable, or `x`
    // itself without bounds. Throws std::invalid_argument for a wrong size.
    std::vector<double> Normalise(const std::vector<double>& x) const;

    // Stores `state` as the state of `x` (in the variables' own units). An
    // entry with the same key is replaced; when full, the oldest entry goes.
    // Throws std::invalid_argument for a wrong size, a non-finite value or an
    // empty state.
    void Add(const std::vector<double>& x, CycleState state);

    // The cached state nearest to `x`, or nullptr when the cache is empty.
    // `distance` (optional) receives the scaled Euclidean distance to it. Ties
    // go to the most recent entry. The pointer is valid until the next Add or
    // Clear.
    const CycleState* Nearest(const std::vector<double>& x, double* distance = nullptr) const;

    size_t Size() const { return entries_.size(); }
    bool Empty() const { return entries_.empty(); }
    size_t Capacity() const { return capacity_; }
    // The number of variables: the bounds' size, or the first key's (0 for an
    // empty cache without bounds).
    size_t Dimension() const;
    void Clear() { entries_.clear(); }

    const std::vector<double>& Lower() const { return lower_; }
    const std::vector<double>& Upper() const { return upper_; }

  private:
    struct Entry {
      std::vector<double> key;  // scaled
      CycleState state;
      std::uint64_t serial = 0;
    };
    void CheckSize(const std::vector<double>& x) const;

    std::vector<double> lower_;
    std::vector<double> upper_;
    size_t capacity_ = kDefaultCapacity;
    std::vector<Entry> entries_;
    std::uint64_t next_serial_ = 0;
  };

  // What RunWithStateCache did.
  struct CachedRunInfo {
    bool warm = false;           // the run that produced the result started from a cached state
    double distance = std::numeric_limits<double>::quiet_NaN();  // scaled distance to that state
    int start_cycles = 0;        // CycleState::cycles of that state
    bool skipped = false;        // a cached state was found but did not fit the deck
    bool retried_cold = false;   // the warm-started run failed and was repeated cold
    bool cached = false;         // the run reached CSS and its final state was added to the cache
    std::string note;            // why the warm start was skipped or retried (empty otherwise)
  };

  // Runs `sim` from the cached state nearest to `x` (a cold run when `cache`
  // is null or empty), then caches the final state of a run that reached
  // cyclic steady state under `x`. `options.initial_state` must be
  // null: the cache chooses the start. When the nearest state does not fit
  // the deck (Simulation::CheckState) the run is cold and the reason is in
  // info->note; when a warm-started run fails it is repeated cold. Notes are
  // also passed to options.on_message (LogLevel::LOW). A cancelled run is not
  // repeated.
  RunResult RunWithStateCache(Simulation& sim, CycleStateCache* cache, const std::vector<double>& x,
                              RunOptions options = {}, CachedRunInfo* info = nullptr);

}  // namespace adsim
