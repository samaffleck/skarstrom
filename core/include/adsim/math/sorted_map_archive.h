/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Deterministic saving of unordered maps.
//
// cereal writes a std::unordered_map in iteration order, which depends on the
// hash function and the bucket count. Saving the same deck twice, or on two
// platforms, can then give files that differ only in the order of the entries
// (components, isotherm parameters, adsorbents). SerializeMap saves the
// entries in key order and loads them as cereal does. The archive format is the
// same as the one cereal writes, so older files load unchanged and the
// containers (and thus the order of any sum over them in the solver) are not
// changed.

#include <algorithm>
#include <cstddef>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/types/concepts/pair_associative_container.hpp>
#include <cereal/types/unordered_map.hpp>

namespace adsim {

namespace detail {

// The entries of a map, sorted by key. Only used to save.
template <class Map>
struct SortedMapView {
  std::vector<const typename Map::value_type*> items;
};

template <class Archive, class Map>
void CEREAL_SAVE_FUNCTION_NAME(Archive& ar, SortedMapView<Map> const& view) {
  ar(cereal::make_size_tag(static_cast<cereal::size_type>(view.items.size())));
  for (const auto* item : view.items)
    ar(cereal::make_map_item(item->first, item->second));
}

}  // namespace detail

// Saves `map` under `name` with its entries sorted by key; loads it as cereal
// does. `Map` is a std::unordered_map (or any map with `value_type` and `<`
// on its keys).
template <class Archive, class Map>
void SerializeMap(Archive& ar, const char* name, Map& map) {
  if constexpr (Archive::is_saving::value) {
    detail::SortedMapView<Map> view;
    view.items.reserve(map.size());
    for (const auto& item : map) view.items.push_back(&item);
    std::sort(view.items.begin(), view.items.end(),
              [](const auto* a, const auto* b) { return a->first < b->first; });
    ar(cereal::make_nvp(name, view));
  } else {
    ar(cereal::make_nvp(name, map));
  }
}

}  // namespace adsim
