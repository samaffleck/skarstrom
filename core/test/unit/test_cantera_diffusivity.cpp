// Molecular diffusivities against Cantera (values recorded in
// cantera_diffusion_reference.h by scripts/gen_cantera_diffusion_reference.py,
// with the engine's Lennard-Jones pairs). This tests the Chapman-Enskog and
// Neufeld implementation (BinaryDiffusionCoefficientModel) and the
// mixture-averaged rule (WilkeMixtureDiffusivity), not the parameter choice.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/binary_diffusion_model.h"
#include "adsim/objects/fluid.h"
#include "cantera_diffusion_reference.h"

namespace adsim {

  namespace {

    namespace ref = cantera_reference;

    // The engine's D_ij at the reference temperatures and pressure. A model of
    // two components holds the pair; a model of one holds the self pair.
    std::vector<double> EngineBinary(const std::string& a, const std::string& b) {
      const std::vector<double> T(ref::kTemperatures.begin(), ref::kTemperatures.end());
      const std::vector<double> P(T.size(), ref::kPressure);
      BinaryDiffusionCoefficientModel model;
      model.Initialise(static_cast<int>(T.size()), a == b ? std::vector<std::string>{ a }
                                                          : std::vector<std::string>{ a, b });
      model.Update(T, P);
      return model.m_DiffusionCoefficients.at({ a, b }).Dm_ij;
    }

  } // namespace

  // The recorded values hold only for the Lennard-Jones table they were
  // recorded with. After a change of the table, run the script again.
  TEST(CanteraDiffusivity, RecordedWithTheEngineLennardJonesTable) {
    for (const auto& in : ref::kInputs) {
      const auto& d = Components::GetComponentData(in.name);
      EXPECT_EQ(d.eps, in.eps) << in.name << ": run scripts/gen_cantera_diffusion_reference.py";
      EXPECT_EQ(d.sigma, in.sigma) << in.name << ": run scripts/gen_cantera_diffusion_reference.py";
    }
  }

  // Every pair, self pairs included, at 200, 300 and 650 K (0.2 < T* < 70).
  // Measured: worst 0.26 %, mean -0.05 %: the Neufeld fit against Cantera's
  // collision-integral tables. (With the rounded textbook constant 2.66e-7 of
  // earlier versions the mean was -0.19 %.)
  TEST(CanteraDiffusivity, BinaryDiffusivitiesMatch) {
    double worst = 0.0, sum = 0.0;
    std::string worst_case;
    int count = 0;
    for (const auto& r : ref::kBinary) {
      const auto D = EngineBinary(r.a, r.b);
      for (size_t k = 0; k < ref::kTemperatures.size(); ++k) {
        const double dev = D[k] / r.D[k] - 1.0;
        EXPECT_LT(std::abs(dev), 0.005) << r.a << "-" << r.b << " at " << ref::kTemperatures[k] << " K";
        sum += dev;
        ++count;
        if (std::abs(dev) > worst) {
          worst = std::abs(dev);
          worst_case = r.a + "-" + r.b + " at " + std::to_string(ref::kTemperatures[k]) + " K";
        }
      }
    }
    ASSERT_GT(count, 1000);
    const double mean = sum / count;
    EXPECT_GT(mean, -0.0015);
    EXPECT_LT(mean, 0.001);
    std::cout << "D_ij vs Cantera " << ref::kCanteraVersion << ": " << count << " values, mean "
              << 100 * mean << " %, worst " << 100 * worst << " % (" << worst_case << ")\n";
  }

  // The mixture-averaged D_m,i (Wilke's equation) with the engine's D_ij, against
  // Cantera's mix_diff_coeffs_mole. Measured: within 0.2 %.
  TEST(CanteraDiffusivity, MixtureAveragedDiffusivitiesMatch) {
    for (const auto& m : ref::kMixtures) {
      for (size_t i = 0; i < m.components.size(); ++i) {
        std::vector<double> y, D;
        for (size_t j = 0; j < m.components.size(); ++j) {
          if (j == i) continue;
          y.push_back(m.y[j]);
          D.push_back(EngineBinary(m.components[i], m.components[j])[1]);   // 300 K
        }
        const double Dm = WilkeMixtureDiffusivity(y.data(), D.data(), y.size());
        EXPECT_NEAR(Dm, m.D[i], 0.005 * m.D[i]) << m.title << ": " << m.components[i];
      }
    }
  }

} // namespace adsim
