/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/experimental_data.h"

#include <charconv>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include <cctype>
#include <string_view>

namespace adsim {

  static inline void chomp_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
  }

  static inline std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c) { return !is_space(c); }).base(), s.end());
    return s;
  }

static bool parse_double(std::string_view sv, double& out) {
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) sv.remove_prefix(1);
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back())))  sv.remove_suffix(1);
    if (sv.empty()) return false;

    const char* b = sv.data();
    char* end;
    out = std::strtod(b, &end);

    return end == b + sv.size();
}

  static void split_row_2cols(const std::string& line, std::string& c0, std::string& c1) {
    // Allow ',', ';', or '\t'
    size_t pos = std::string::npos;
    for (char d : {',', ';', '\t'}) {
      size_t p = line.find(d);
      if (p != std::string::npos) { pos = p; break; }
    }
    if (pos == std::string::npos) {
      c0 = trim(line);
      c1.clear();
    }
    else {
      c0 = trim(line.substr(0, pos));
      c1 = trim(line.substr(pos + 1));
    }
  }

  bool ExperimentalData::LoadData(const std::filesystem::path& path, bool headers_included) {
    std::ifstream in(path);
    if (!in) return false;

    XYData d;
    d.name = path.stem().string();

    std::string line;
    bool first_line = true;

    while (std::getline(in, line)) {
      chomp_cr(line);
      if (line.empty()) continue;

      if (first_line && headers_included) {
        first_line = false;
        std::string h0, h1;
        split_row_2cols(line, h0, h1);
        if (!h0.empty()) d.x_label = h0;
        if (!h1.empty()) d.y_label = h1;
        continue;
      }
      first_line = false;

      std::string s0, s1;
      split_row_2cols(line, s0, s1);
      if (s0.empty() || s1.empty()) continue;

      double xv, yv;
      if (!parse_double(s0, xv)) continue;
      if (!parse_double(s1, yv)) continue;

      d.x.push_back(xv);
      d.y.push_back(yv);
    }

    if (d.x.empty() || d.y.empty() || d.x.size() != d.y.size()) 
      return false;

    data_.push_back(std::move(d));
    return true;
  }
  
  bool ExperimentalData::AddData(const XYData& data) {
    data_.push_back(data);
    data_.back().name = "Experiment " + std::to_string(data_.size());
    return true;
  }

  bool ExperimentalData::DeleteDataByIndex(size_t idx) {
    if (idx >= data_.size()) return false;
    data_.erase(data_.begin() + static_cast<std::ptrdiff_t>(idx));
    return true;
  }

} // end adsim namespace


CEREAL_CLASS_VERSION(adsim::ExperimentalData, 1);
