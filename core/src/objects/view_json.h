/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Internal helpers for the JSON the deck view and the manifests emit
// (rapidjson DOM, written compact or indented). Not installed.

#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/prettywriter.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>

namespace adsim::view_json {

  using Value = rapidjson::Value;
  using Document = rapidjson::Document;
  using Allocator = rapidjson::Document::AllocatorType;

  inline Value Str(std::string_view s, Allocator& a) {
    return Value(s.data(), static_cast<rapidjson::SizeType>(s.size()), a);
  }

  // A finite number, or null.
  inline Value Num(double x) {
    Value v;
    if (std::isfinite(x)) v.SetDouble(x);
    return v;
  }

  inline Value StrList(const std::vector<std::string>& items, Allocator& a) {
    Value v(rapidjson::kArrayType);
    for (const auto& s : items) v.PushBack(Str(s, a), a);
    return v;
  }

  // obj[key] = value (key copied).
  inline void Set(Value& obj, std::string_view key, Value value, Allocator& a) {
    obj.AddMember(Str(key, a), value, a);
  }
  inline void Set(Value& obj, std::string_view key, std::string_view s, Allocator& a) {
    obj.AddMember(Str(key, a), Str(s, a), a);
  }
  inline void Set(Value& obj, std::string_view key, const char* s, Allocator& a) {
    obj.AddMember(Str(key, a), Str(s ? s : "", a), a);
  }
  inline void Set(Value& obj, std::string_view key, const std::string& s, Allocator& a) {
    obj.AddMember(Str(key, a), Str(s, a), a);
  }
  inline void Set(Value& obj, std::string_view key, bool b, Allocator& a) {
    obj.AddMember(Str(key, a), Value(b), a);
  }
  inline void SetInt(Value& obj, std::string_view key, long long n, Allocator& a) {
    obj.AddMember(Str(key, a), Value(static_cast<int64_t>(n)), a);
  }
  inline void SetNull(Value& obj, std::string_view key, Allocator& a) {
    obj.AddMember(Str(key, a), Value(), a);
  }

  inline std::string Write(const Value& v, bool pretty) {
    rapidjson::StringBuffer buffer;
    if (pretty) {
      rapidjson::PrettyWriter<rapidjson::StringBuffer> w(buffer);
      w.SetIndent(' ', 2);
      v.Accept(w);
    } else {
      rapidjson::Writer<rapidjson::StringBuffer> w(buffer);
      v.Accept(w);
    }
    return std::string(buffer.GetString(), buffer.GetSize());
  }

}  // namespace adsim::view_json
