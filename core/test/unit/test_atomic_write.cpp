// WriteFileAtomically / SaveToJsonAtomically (test/fixtures/atomic_write.h):
// the helper the suite uses to re-save example decks while other tests may be
// reading them. A reader must only ever see a complete old or new file, the
// temporary file must not be left behind, and failures must throw.

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "../src/pch.h"
#include "../fixtures/atomic_write.h"
#include "adsim/objects/user_data.h"

namespace adsim {
  namespace {

    namespace fs = std::filesystem;

    std::string ReadAll(const fs::path& p) {
      std::ifstream in(p, std::ios::binary);
      return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }

    std::vector<std::string> Listing(const fs::path& dir) {
      std::vector<std::string> names;
      for (const auto& e : fs::directory_iterator(dir)) names.push_back(e.path().filename().string());
      std::sort(names.begin(), names.end());
      return names;
    }

    class AtomicWrite : public ::testing::Test {
    protected:
      void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        dir_ = fs::temp_directory_path() / (std::string("adsim_atomic_write_") + info->name());
        fs::remove_all(dir_);
        fs::create_directories(dir_);
      }
      void TearDown() override {
        std::error_code ec;
        fs::permissions(dir_, fs::perms::owner_all, fs::perm_options::add, ec);
        fs::remove_all(dir_, ec);
      }
      fs::path dir_;
    };

    TEST_F(AtomicWrite, WritesTheFileAndLeavesNoTemporary) {
      const fs::path target = dir_ / "deck.json";
      test::WriteFileAtomically(target, [](std::ostream& out) { out << "{\"a\": 1}"; });
      EXPECT_EQ(ReadAll(target), "{\"a\": 1}");
      EXPECT_EQ(Listing(dir_), std::vector<std::string>{"deck.json"});

      test::WriteFileAtomically(target, [](std::ostream& out) { out << "{\"b\": 2}"; });
      EXPECT_EQ(ReadAll(target), "{\"b\": 2}");
      EXPECT_EQ(Listing(dir_), std::vector<std::string>{"deck.json"});
    }

    TEST_F(AtomicWrite, ReaderNeverSeesAPartialFile) {
      // A writer alternates two large contents while a reader keeps reading
      // the target: every read must be exactly one of them (with a plain
      // ofstream rewrite, reads of an empty or truncated file show up).
      const fs::path target = dir_ / "deck.json";
      const std::string a(512 * 1024, 'a');
      const std::string b(384 * 1024, 'b');
      test::WriteFileAtomically(target, [&](std::ostream& out) { out << a; });

      std::atomic<bool> done{false};
      std::atomic<int> reads{0}, bad{0}, unopened{0};
      std::thread reader([&] {
        while (!done.load()) {
          std::ifstream in(target, std::ios::binary);
          if (!in.is_open()) { ++unopened; continue; }
          const std::string s{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
          if (s != a && s != b) ++bad;
          ++reads;
        }
      });
      for (int i = 0; i < 60; ++i) {
        const std::string& next = (i % 2 == 0) ? b : a;
        // On Windows the rename can wait for the reader to close the file;
        // the helper retries, so this must not throw.
        ASSERT_NO_THROW(test::WriteFileAtomically(target, [&](std::ostream& out) { out << next; }));
      }
      done = true;
      reader.join();

      EXPECT_EQ(bad.load(), 0) << "of " << reads.load() << " reads";
      EXPECT_GT(reads.load(), 0);
#ifndef _WIN32
      // POSIX rename() swaps the directory entry atomically: the target can
      // always be opened. (Windows may refuse an open during the replace.)
      EXPECT_EQ(unopened.load(), 0);
#endif
      EXPECT_EQ(Listing(dir_), std::vector<std::string>{"deck.json"});
    }

    TEST_F(AtomicWrite, FailingWriterLeavesTargetAndRemovesTemporary) {
      const fs::path target = dir_ / "deck.json";
      test::WriteFileAtomically(target, [](std::ostream& out) { out << "old"; });
      EXPECT_THROW(test::WriteFileAtomically(target,
                                             [](std::ostream& out) {
                                               out << "half";
                                               throw std::runtime_error("writer failed");
                                             }),
                   std::runtime_error);
      EXPECT_EQ(ReadAll(target), "old");
      EXPECT_EQ(Listing(dir_), std::vector<std::string>{"deck.json"});
    }

    TEST_F(AtomicWrite, MissingDirectoryThrows) {
      const fs::path target = dir_ / "no_such_dir" / "deck.json";
      EXPECT_THROW(test::WriteFileAtomically(target, [](std::ostream& out) { out << "x"; }),
                   std::runtime_error);
      EXPECT_FALSE(fs::exists(dir_ / "no_such_dir"));
      EXPECT_TRUE(Listing(dir_).empty());
    }

    TEST_F(AtomicWrite, ReadOnlyDirectoryThrowsAndLeavesTargetUnchanged) {
      const fs::path target = dir_ / "deck.json";
      test::WriteFileAtomically(target, [](std::ostream& out) { out << "old"; });
      fs::permissions(dir_, fs::perms::owner_read | fs::perms::owner_exec |
                                fs::perms::group_read | fs::perms::group_exec |
                                fs::perms::others_read | fs::perms::others_exec);
      {
        // Root and Windows ignore POSIX directory permissions: nothing to test.
        std::ofstream probe(dir_ / "probe");
        if (probe.is_open()) {
          probe.close();
          fs::permissions(dir_, fs::perms::owner_all, fs::perm_options::add);
          fs::remove(dir_ / "probe");
          GTEST_SKIP() << "directory permissions are not enforced for this user";
        }
      }
      EXPECT_THROW(test::WriteFileAtomically(target, [](std::ostream& out) { out << "new"; }),
                   std::runtime_error);
      fs::permissions(dir_, fs::perms::owner_all, fs::perm_options::add);
      EXPECT_EQ(ReadAll(target), "old");
      EXPECT_EQ(Listing(dir_), std::vector<std::string>{"deck.json"});
    }

    TEST_F(AtomicWrite, SaveToJsonAtomicallyMatchesSaveToJson) {
      // Same bytes and the same sim_version bump as SaveToJson, so re-saved
      // example decks do not change because of the atomic write.
      UserData data;
      const auto version = data.info.sim_version;
      SaveToJson(data, (dir_ / "plain.json").string());
      EXPECT_EQ(data.info.sim_version, version + 1);

      data.info.sim_version = version;
      test::SaveToJsonAtomically(data, dir_ / "atomic.json");
      EXPECT_EQ(data.info.sim_version, version + 1);

      const std::string plain = ReadAll(dir_ / "plain.json");
      ASSERT_FALSE(plain.empty());
      EXPECT_EQ(ReadAll(dir_ / "atomic.json"), plain);
      EXPECT_EQ(Listing(dir_), (std::vector<std::string>{"atomic.json", "plain.json"}));
    }

  }  // namespace
}  // namespace adsim
