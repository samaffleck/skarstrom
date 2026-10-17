/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Curated adsorbent library catalog (data/adsorbents/library.json).
//
// The library is curated in-house: every sourced entry carries a citation (DOI
// + table/figure references) and per-parameter provenance. Each entry's
// physics lives in its own cereal single-Adsorbent blob next to the catalog
// (`file`), loadable unchanged by LoadAdsorbentFromJSON; the catalog only adds
// the metadata (no C++ class changes, so every existing adsorbent file and
// deck still loads).
//
// This header is a small read-only view of the catalog (list + citation) for
// front ends and the library tests. scripts/docs/gen_adsorbent_library.py
// validates the catalog against the blobs and renders the docs page from it.

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cereal/external/rapidjson/document.h>

namespace adsim {

  struct AdsorbentLibrarySource {
    std::string id;
    std::string citation;  // full reference, human readable
    std::string doi;       // bare DOI ("10.xxxx/…"); empty only for legacy leads
  };

  struct AdsorbentLibrarySourceRef {
    std::string source_id;  // key into AdsorbentLibrary::sources
    std::string refs;       // tables / figures / sections used
  };

  struct AdsorbentLibraryEntry {
    std::string key;       // stable id (also the web add_example_adsorbent key)
    std::string file;      // cereal blob, relative to the catalog directory
    std::string label;     // display name (== the blob's Name)
    std::string material;  // material name without the source suffix
    std::string family;    // zeolite / carbon / MOF / silica / amine / ...
    std::string status;    // "sourced" or "legacy"
    std::string legacy_note;
    std::vector<std::string> applications;
    std::vector<std::string> components;
    std::vector<std::string> example_decks;
    std::vector<AdsorbentLibrarySourceRef> sources;

    bool IsLegacy() const { return status == "legacy"; }
  };

  struct AdsorbentLibrary {
    int schema_version = 0;
    std::vector<AdsorbentLibrarySource> sources;
    std::vector<AdsorbentLibraryEntry> entries;

    const AdsorbentLibrarySource* FindSource(const std::string& id) const {
      for (const auto& s : sources)
        if (s.id == id) return &s;
      return nullptr;
    }

    const AdsorbentLibraryEntry* Find(const std::string& key) const {
      for (const auto& e : entries)
        if (e.key == key) return &e;
      return nullptr;
    }

    // One line per cited source: "<citation> doi:<doi> (<refs>)".
    std::string CitationText(const AdsorbentLibraryEntry& e) const {
      if (e.IsLegacy()) return "Source unknown - legacy entry. " + e.legacy_note;
      std::string out;
      for (const auto& r : e.sources) {
        const AdsorbentLibrarySource* s = FindSource(r.source_id);
        if (!s) continue;
        if (!out.empty()) out += "\n";
        out += s->citation;
        if (!s->doi.empty()) out += " doi:" + s->doi;
        if (!r.refs.empty()) out += " (" + r.refs + ")";
      }
      return out;
    }
  };

  namespace adsorbent_library_detail {
    inline std::string Str(const rapidjson::Value& o, const char* k) {
      if (!o.IsObject() || !o.HasMember(k) || !o[k].IsString()) return {};
      return o[k].GetString();
    }
    inline std::vector<std::string> StrList(const rapidjson::Value& o, const char* k) {
      std::vector<std::string> out;
      if (!o.IsObject() || !o.HasMember(k) || !o[k].IsArray()) return out;
      for (const auto& v : o[k].GetArray())
        if (v.IsString()) out.emplace_back(v.GetString());
      return out;
    }
  }  // namespace adsorbent_library_detail

  // Parses the catalog JSON text. Throws std::runtime_error on malformed input.
  inline AdsorbentLibrary ParseAdsorbentLibrary(const std::string& text) {
    namespace d = adsorbent_library_detail;
    rapidjson::Document doc;
    doc.Parse(text.c_str());
    if (doc.HasParseError() || !doc.IsObject())
      throw std::runtime_error("adsorbent library: catalog is not a JSON object");
    AdsorbentLibrary lib;
    if (doc.HasMember("schema_version") && doc["schema_version"].IsInt())
      lib.schema_version = doc["schema_version"].GetInt();
    if (doc.HasMember("sources") && doc["sources"].IsObject()) {
      for (auto it = doc["sources"].MemberBegin(); it != doc["sources"].MemberEnd(); ++it) {
        AdsorbentLibrarySource s;
        s.id = it->name.GetString();
        s.citation = d::Str(it->value, "citation");
        s.doi = d::Str(it->value, "doi");
        lib.sources.push_back(std::move(s));
      }
    }
    if (!doc.HasMember("entries") || !doc["entries"].IsArray())
      throw std::runtime_error("adsorbent library: catalog has no 'entries' array");
    for (const auto& v : doc["entries"].GetArray()) {
      AdsorbentLibraryEntry e;
      e.key = d::Str(v, "key");
      e.file = d::Str(v, "file");
      e.label = d::Str(v, "label");
      e.material = d::Str(v, "material");
      e.family = d::Str(v, "family");
      e.status = d::Str(v, "status");
      e.legacy_note = d::Str(v, "legacy_note");
      e.applications = d::StrList(v, "applications");
      e.components = d::StrList(v, "components");
      e.example_decks = d::StrList(v, "example_decks");
      if (v.HasMember("sources") && v["sources"].IsArray()) {
        for (const auto& r : v["sources"].GetArray())
          e.sources.push_back({d::Str(r, "id"), d::Str(r, "refs")});
      }
      lib.entries.push_back(std::move(e));
    }
    return lib;
  }

  inline AdsorbentLibrary LoadAdsorbentLibraryFile(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) throw std::runtime_error("adsorbent library: cannot open " + path);
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ParseAdsorbentLibrary(ss.str());
  }

}  // namespace adsim
