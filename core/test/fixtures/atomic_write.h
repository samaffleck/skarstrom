#pragma once

// Atomic file replacement for files the test suite rewrites while other tests
// may be reading them (the example decks in examples/decks/, re-saved by the
// paper tests). Under `ctest --parallel` a reader that opens a deck while a
// plain std::ofstream is rewriting it sees truncated JSON. Writing to a
// temporary file in the same directory and renaming it over the target means a
// reader opens either the old or the new complete file, never a partial one.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>

#include <cereal/archives/json.hpp>

#include "adsim/objects/user_data.h"

namespace adsim::test {

  namespace detail {
    // A temporary name next to `target`: ".<name>.tmp-<16 hex digits>". It is
    // shell-safe (cmake/check_data_names.cmake) and does not end in ".json",
    // so tools that glob the decks never pick it up.
    inline std::filesystem::path AtomicTempPath(const std::filesystem::path& target) {
      static std::atomic<std::uint64_t> counter{0};
      std::random_device rd;
      std::uint64_t r = (static_cast<std::uint64_t>(rd()) << 32) ^ rd();
      r ^= std::hash<std::thread::id>{}(std::this_thread::get_id());
      r ^= static_cast<std::uint64_t>(
          std::chrono::steady_clock::now().time_since_epoch().count());
      r += counter.fetch_add(1) * 0x9E3779B97F4A7C15ull;
      std::ostringstream name;
      name << '.' << target.filename().string() << ".tmp-" << std::hex;
      name.width(16);
      name.fill('0');
      name << r;
      return target.parent_path() / name.str();
    }
  }  // namespace detail

  // Writes `target` atomically: `write` fills a temporary file in the same
  // directory, which is flushed, closed and renamed over `target`.
  //
  // Throws std::runtime_error (temporary file cannot be created or written,
  // e.g. the directory is missing or not writable), std::filesystem_error (the
  // rename keeps failing) or whatever `write` throws. On any failure the
  // temporary file is removed and `target` is left as it was.
  //
  // POSIX rename() replaces the target atomically even while readers hold it
  // open. On Windows the replace (MoveFileExW) fails with a sharing violation
  // while another process has the target open without FILE_SHARE_DELETE, which
  // is how std::ifstream opens it; that window is short (a test reading a
  // deck), so the rename is retried with a growing back-off for about 5 s
  // before giving up with the last error.
  inline void WriteFileAtomically(const std::filesystem::path& target,
                                  const std::function<void(std::ostream&)>& write) {
    namespace fs = std::filesystem;
    const fs::path tmp = detail::AtomicTempPath(target);
    auto remove_tmp = [&tmp] {
      std::error_code ignored;
      fs::remove(tmp, ignored);
    };

    try {
      {
        // Text mode, like SaveToJson, so the bytes match it on every platform.
        std::ofstream out(tmp, std::ios::out | std::ios::trunc);
        if (!out.is_open())
          throw std::runtime_error("cannot create temporary file '" + tmp.string() +
                                   "' to write '" + target.string() + "'");
        write(out);
        out.flush();
        if (!out)
          throw std::runtime_error("failed writing temporary file '" + tmp.string() + "'");
        out.close();
        if (out.fail())
          throw std::runtime_error("failed closing temporary file '" + tmp.string() + "'");
      }

      std::error_code ec;
      auto delay = std::chrono::milliseconds(5);
      const auto give_up = std::chrono::steady_clock::now() + std::chrono::seconds(5);
      for (;;) {
        fs::rename(tmp, target, ec);
        if (!ec) break;
        if (std::chrono::steady_clock::now() >= give_up)
          throw fs::filesystem_error("cannot replace file (still failing after retries)",
                                     tmp, target, ec);
        std::this_thread::sleep_for(delay);
        delay = std::min(delay * 2, std::chrono::milliseconds(200));
      }
    } catch (...) {
      remove_tmp();
      throw;
    }
  }

  // SaveToJson (adsim/objects/user_data.h) with an atomic write: the same
  // bytes (cereal JSON, root "value0") and the same sim_version bump, but a
  // concurrent reader never sees a partly written file, and failures throw
  // instead of being ignored.
  inline void SaveToJsonAtomically(UserData& data, const std::filesystem::path& file_name) {
    data.info.sim_version++;
    WriteFileAtomically(file_name, [&data](std::ostream& out) {
      cereal::JSONOutputArchive archive(out);  // flushes the closing braces on destruction
      archive(data);
    });
  }

}  // namespace adsim::test
