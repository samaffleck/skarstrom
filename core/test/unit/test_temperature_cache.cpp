// The per-cell temperature memo (TemperatureCache) and the isotherm loops that
// read it: the cached loadings must be bit-identical to the plain per-cell
// functions, through temperature and parameter changes.

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/binary_diffusion_model.h"
#include "adsim/equations/equation.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/temperature_cache.h"
#include "adsim/math/constants.h"

namespace adsim {

  namespace {
    bool SameBits(double a, double b) {
      return std::memcmp(&a, &b, sizeof a) == 0;
    }

    // A two-component state for the isotherm loops: y, q* maps and T, P vectors.
    struct IsothermState {
      ComponentEquation yi;
      ComponentIndependentEquation qi_sat;
      std::vector<double> T, P;
      std::vector<std::string> components{ "CarbonDioxide", "Nitrogen" };

      explicit IsothermState(size_t n) : T(n), P(n) {
        for (size_t i = 0; i < n; ++i) {
          T[i] = 290.0 + 7.5 * i;
          P[i] = 1.0e5 + 3.0e4 * i;
        }
        for (size_t k = 0; k < components.size(); ++k) {
          const auto& c = components[k];
          auto [it, ok] = yi.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                                     std::forward_as_tuple(c, CellType::CellCentre, EquationType::Differential,
                                                           [](UserData*) {}, 1e-8));
          it->second.var.resize(n);
          for (size_t i = 0; i < n; ++i) it->second.var[i] = k == 0 ? 0.1 + 0.05 * i : 0.9 - 0.05 * i;
          auto [jt, ok2] = qi_sat.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                                          std::forward_as_tuple(c, CellType::CellCentre, [](UserData*) {}));
          jt->second.var.assign(n, 0.0);
        }
      }
      IsothermUpdateParams Params(const std::string& component) {
        return IsothermUpdateParams(yi, qi_sat, T, P, component, 0, static_cast<int>(T.size()) - 1);
      }
    };
  } // namespace

  TEST(TemperatureCache, HitsOnSameBitsAndRecomputesOnAnyKeyChange) {
    TemperatureCache cache;
    int calls = 0;
    auto f = [&] { ++calls; return 42.0; };
    EXPECT_EQ(cache.Get(3, 300.0, 1.0, 2.0, f), 42.0);
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(cache.Get(3, 300.0, 1.0, 2.0, f), 42.0);
    EXPECT_EQ(calls, 1);                      // same cell, same keys: no call
    cache.Get(3, 300.0000001, 1.0, 2.0, f);
    EXPECT_EQ(calls, 2);                      // T changed
    cache.Get(3, 300.0000001, 1.5, 2.0, f);
    EXPECT_EQ(calls, 3);                      // first parameter changed
    cache.Get(3, 300.0000001, 1.5, 2.5, f);
    EXPECT_EQ(calls, 4);                      // second parameter changed
    cache.Get(4, 300.0000001, 1.5, 2.5, f);
    EXPECT_EQ(calls, 5);                      // another cell
    cache.Get(3, -0.0, 1.5, 2.5, f);
    cache.Get(3, 0.0, 1.5, 2.5, f);
    EXPECT_EQ(calls, 7);                      // +0 and -0 are different keys
    EXPECT_EQ(cache.Size(), 5u);
  }

  TEST(TemperatureCache, DualSiteLangmuirLoopMatchesPerCellFunctionBitForBit) {
    DualSiteLangmuirIsotherm iso;
    iso.params_["CarbonDioxide"] = DualSiteLangmuirIsothermParameters{ 3.0, 1e-12, 4.7e4, 2.0, 1e-13, 4.2e4 };
    iso.params_["Nitrogen"] = DualSiteLangmuirIsothermParameters{ 1.0, 2e-12, 1.8e4, 0.5, 3e-13, 1.5e4 };
    IsothermState s(6);

    auto check = [&] {
      for (const auto& comp : s.components) {
        auto p = s.Params(comp);
        iso.Update(p);
        for (size_t n = 0; n < s.T.size(); ++n) {
          double expected = 0.0;
          IsothermFunction::DualSiteLangmuir(expected, s.T[n], s.P[n], s.yi.at(comp).var[n], iso.params_.at(comp));
          EXPECT_TRUE(SameBits(s.qi_sat.at(comp).var[n], expected)) << comp << " cell " << n;
        }
      }
    };
    check();
    check();                       // every cell a cache hit
    s.T[2] = 317.3;                // one cell's temperature moves
    check();
    iso.params_["CarbonDioxide"].q1 = 4.9e4;   // a parameter edit between runs
    check();
    s.T.assign(s.T.size(), 305.0); // all cells at one temperature
    check();
  }

  TEST(TemperatureCache, ExtendedDualSiteConcentrationUpdateAllMatchesTheFormula) {
    ExtendedDualSiteLangmuirIsotherm_Concentration iso;
    iso.params_["CarbonDioxide"] = DualSiteLangmuirIsothermParametersConcentration{ 3.0, 1e-8, 4.0e4, 2.0, 1e-7, 3.5e4 };
    iso.params_["Nitrogen"] = DualSiteLangmuirIsothermParametersConcentration{ 1.0, 2e-8, 1.8e4, 0.5, 3e-7, 1.5e4 };
    IsothermState s(5);
    const double R = constants::kIdealGasConstant;

    // The competitive dual-site Langmuir loading written out as the former
    // per-component loop did it (denominators summed in params_ map order).
    auto reference = [&](const std::string& comp, size_t n) {
      double d1 = 0.0, d2 = 0.0;
      for (const auto& [c, p] : iso.params_) {
        double ci = s.P[n] * s.yi.at(c).var[n] / (R * s.T[n]);
        d1 += p.b1 * std::exp(p.q1 / (R * s.T[n])) * ci;
        d2 += p.b2 * std::exp(p.q2 / (R * s.T[n])) * ci;
      }
      const auto& p = iso.params_.at(comp);
      double ci = s.P[n] * s.yi.at(comp).var[n] / (R * s.T[n]);
      return p.m1 * p.b1 * std::exp(p.q1 / (R * s.T[n])) * ci / (1 + d1) +
             p.m2 * p.b2 * std::exp(p.q2 / (R * s.T[n])) * ci / (1 + d2);
    };

    auto check_all = [&] {
      auto p = s.Params(s.components.front());
      iso.UpdateAll(p, s.components);
      for (const auto& comp : s.components)
        for (size_t n = 0; n < s.T.size(); ++n)
          EXPECT_DOUBLE_EQ(s.qi_sat.at(comp).var[n], reference(comp, n)) << comp << " cell " << n;
    };
    check_all();
    check_all();
    s.T[1] = 299.9;
    s.yi.at("Nitrogen").var[3] = 0.4;
    check_all();

    // The per-component Update is the same pass for one component.
    std::vector<double> from_all = s.qi_sat.at("Nitrogen").var;
    auto p = s.Params(s.components[1]);  // "Nitrogen"
    iso.Update(p);
    for (size_t n = 0; n < s.T.size(); ++n)
      EXPECT_TRUE(SameBits(s.qi_sat.at("Nitrogen").var[n], from_all[n])) << n;
  }

  TEST(TemperatureCache, ExtendedSips2AndArrheniusSipsLoopsAreStableAcrossRepeats) {
    // A second pass over an unchanged state (all cache hits) reproduces the
    // first pass bit for bit, and a parameter edit changes the result.
    {
      ExtendedSips2Isotherm iso;
      iso.params_["CarbonDioxide"];
      iso.params_["Nitrogen"];
      iso.params_["Nitrogen"].B2 = 9.0e3;
      IsothermState s(5);
      auto p = s.Params(s.components.front());
      iso.UpdateAll(p, s.components);
      std::vector<double> first = s.qi_sat.at("CarbonDioxide").var;
      iso.UpdateAll(p, s.components);
      for (size_t n = 0; n < first.size(); ++n) EXPECT_TRUE(SameBits(first[n], s.qi_sat.at("CarbonDioxide").var[n]));
      iso.params_["CarbonDioxide"].B2 = 16.0e3;
      iso.UpdateAll(p, s.components);
      EXPECT_NE(first[2], s.qi_sat.at("CarbonDioxide").var[2]);
    }
    {
      ArrheniusSipsIsotherm iso;
      iso.params_["CarbonDioxide"];
      iso.params_["Nitrogen"];
      iso.params_["CarbonDioxide"].c_amp = 0.2;
      iso.params_["CarbonDioxide"].c_beta = 0.01;
      IsothermState s(5);
      auto p = s.Params(s.components.front());
      iso.UpdateAll(p, s.components);
      std::vector<double> first = s.qi_sat.at("CarbonDioxide").var;
      iso.UpdateAll(p, s.components);
      for (size_t n = 0; n < first.size(); ++n) EXPECT_TRUE(SameBits(first[n], s.qi_sat.at("CarbonDioxide").var[n]));
      iso.params_["CarbonDioxide"].k_B = 2.0e3;
      iso.UpdateAll(p, s.components);
      EXPECT_NE(first[2], s.qi_sat.at("CarbonDioxide").var[2]);
    }
  }

  TEST(TemperatureCache, BinaryDiffusionTableIsBitIdenticalAcrossRepeatsAndTemperatureEdits) {
    BinaryDiffusionCoefficientModel model;
    const std::vector<std::string> comps{ "CarbonDioxide", "Nitrogen" };
    const int N = 6;
    model.Initialise(N, comps);
    std::vector<double> T(N), P(N);
    for (int i = 0; i < N; ++i) { T[i] = 280.0 + 11.0 * i; P[i] = 1.2e5 + 2.0e4 * i; }

    // Reference: a fresh model with no memo history, evaluated once.
    auto reference = [&](const std::vector<double>& Tv, const std::vector<double>& Pv) {
      BinaryDiffusionCoefficientModel fresh;
      fresh.Initialise(N, comps);
      fresh.Update(Tv, Pv);
      return fresh.m_DiffusionCoefficients.at({ comps[0], comps[1] }).Dm_ij;
    };

    model.Update(T, P);
    model.Update(T, P);                                  // all memo hits
    auto ref = reference(T, P);
    const auto& got = model.m_DiffusionCoefficients.at({ comps[0], comps[1] }).Dm_ij;
    for (int i = 0; i < N; ++i) EXPECT_TRUE(SameBits(got[i], ref[i])) << i;

    T[2] = 333.3;                                        // one cell's temperature moves
    P[4] = 9.0e4;                                        // pressure is not memoised
    model.Update(T, P);
    ref = reference(T, P);
    for (int i = 0; i < N; ++i) EXPECT_TRUE(SameBits(got[i], ref[i])) << i;
  }

} // namespace adsim
