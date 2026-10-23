// adsim includes
#include "../src/pch.h"
#include "adsim/objects/compressor.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/objects/fluid.h"

// STL includes
#include <cmath>
#include <sstream>
#include <string>

// cereal includes
#include <cereal/archives/json.hpp>

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  TEST(PumpEfficiencyTest, ConstantModelDefault) {
    ConstantPumpEfficiencyModel model;
    EXPECT_DOUBLE_EQ(model.efficiency_, 0.7);
    EXPECT_DOUBLE_EQ(model.CalculateEfficiency(1e5), 0.7);
    EXPECT_DOUBLE_EQ(model.CalculateEfficiency(5e3), 0.7);
    EXPECT_EQ(model.GetType(), PumpEfficiencyType::CONSTANT);
  }

  TEST(PumpEfficiencyTest, VacuumPumpModelAtAtmospheric) {
    VacuumPumpEfficiencyModel model;
    // At atmospheric: eta = 0.8 * 19.55 * 1.0 / (1 + 19.55 * 1.0)
    double expected = (0.8 * 19.55 * 1.0) / (1.0 + 19.55 * 1.0);
    EXPECT_NEAR(model.CalculateEfficiency(1e5), expected, 1e-10);
    EXPECT_EQ(model.GetType(), PumpEfficiencyType::VACUUM_PUMP);
  }

  TEST(PumpEfficiencyTest, VacuumPumpModelAtLowPressure) {
    VacuumPumpEfficiencyModel model;
    // At 5000 Pa (0.05 bar): eta = 0.8 * 19.55 * 0.05 / (1 + 19.55 * 0.05)
    double p_ratio = 5000.0 / 1e5;
    double expected = (0.8 * 19.55 * p_ratio) / (1.0 + 19.55 * p_ratio);
    EXPECT_NEAR(model.CalculateEfficiency(5000.0), expected, 1e-10);
    // Should be much lower than atmospheric efficiency
    EXPECT_LT(model.CalculateEfficiency(5000.0), model.CalculateEfficiency(1e5));
  }

  TEST(PumpEfficiencyTest, VacuumPumpModelApproachesZeroAtVacuum) {
    VacuumPumpEfficiencyModel model;
    // At very low pressure, efficiency should approach 0
    EXPECT_NEAR(model.CalculateEfficiency(0.0), 0.0, 1e-10);
  }

  TEST(PumpEfficiencyTest, VacuumPumpModelCustomParams) {
    VacuumPumpEfficiencyModel model;
    model.eta_max_ = 0.9;
    model.alpha_ = 10.0;

    double p_ratio = 0.5; // 50000 Pa
    double expected = (0.9 * 10.0 * p_ratio) / (1.0 + 10.0 * p_ratio);
    EXPECT_NEAR(model.CalculateEfficiency(50000.0), expected, 1e-10);
  }

  TEST(PumpEfficiencyTest, PolyModelDefault) {
    PolyPumpEfficiencyModel model;
    // All coefficients zero by default
    EXPECT_DOUBLE_EQ(model.CalculateEfficiency(1e5), 0.0);
    EXPECT_DOUBLE_EQ(model.CalculateEfficiency(0.0), 0.0);
    EXPECT_EQ(model.GetType(), PumpEfficiencyType::POLYNOMIAL);
  }

  TEST(PumpEfficiencyTest, PolyModelEvaluation) {
    PolyPumpEfficiencyModel model;
    model.a_ = 0.1;
    model.b_ = 2e-6;
    model.c_ = 0.0;
    model.d_ = 0.0;
    model.e_ = 0.0;

    double P = 50000.0;
    // eta = 0.1 + 2e-6 * 50000 = 0.1 + 0.1 = 0.2
    EXPECT_NEAR(model.CalculateEfficiency(P), 0.2, 1e-10);
  }

  TEST(PumpEfficiencyTest, PolyModelAllTerms) {
    PolyPumpEfficiencyModel model;
    model.a_ = 1.0;
    model.b_ = 1.0;
    model.c_ = 1.0;
    model.d_ = 1.0;
    model.e_ = 1.0;
    model.f_ = 1.0;

    double P = 2.0;
    // eta = 1 + 2 + 4 + 8 + 16 + 32 = 63
    EXPECT_NEAR(model.CalculateEfficiency(P), 63.0, 1e-10);
  }

  TEST(PumpEfficiencyTest, PolyModelEachCoefficientMultipliesItsPower) {
    const double P = 3.0;
    for (int k = 0; k <= 5; ++k) {
      SCOPED_TRACE("coefficient of P^" + std::to_string(k));
      PolyPumpEfficiencyModel model;
      ExpressionParam* coefficient[] = {&model.a_, &model.b_, &model.c_, &model.d_, &model.e_, &model.f_};
      *coefficient[k] = 1.0;
      EXPECT_DOUBLE_EQ(model.CalculateEfficiency(P), std::pow(P, k));
    }
  }

  // PolyPumpEfficiencyModel v3 added the P^2 term. A v2 archive's c, d, e were
  // the P^3, P^4, P^5 coefficients; they load as d, e, f with c = 0, so an old
  // deck's efficiency is unchanged.
  TEST(PumpEfficiencyTest, PolyModelVersion2LoadsWithTheSameEfficiency) {
    PolyPumpEfficiencyModel saved;
    saved.a_ = 0.1;
    saved.b_ = 2e-6;
    saved.c_ = 3e-17;  // v2: x P^3
    saved.d_ = -4e-23; // v2: x P^4
    saved.e_ = 5e-29;  // v2: x P^5
    std::ostringstream out;
    {
      cereal::JSONOutputArchive archive(out);
      archive(cereal::make_nvp("model", saved));
    }
    std::string json = out.str();
    // The first class version in the document is the model's own (the
    // leftover "f" is not read by the v2 layout).
    const std::string v3 = "\"cereal_class_version\": 3";
    const auto version_pos = json.find(v3);
    ASSERT_NE(version_pos, std::string::npos) << json;
    json.replace(version_pos, v3.size(), "\"cereal_class_version\": 2");

    PolyPumpEfficiencyModel loaded;
    {
      std::istringstream in(json);
      cereal::JSONInputArchive archive(in);
      archive(cereal::make_nvp("model", loaded));
    }
    EXPECT_EQ(loaded.c_.Get(), 0.0);
    EXPECT_EQ(loaded.d_.Get(), 3e-17);
    EXPECT_EQ(loaded.e_.Get(), -4e-23);
    EXPECT_EQ(loaded.f_.Get(), 5e-29);
    for (double P : {1e4, 1e5, 3e5}) {
      const double v2_eta = 0.1 + 2e-6 * P + 3e-17 * std::pow(P, 3) - 4e-23 * std::pow(P, 4) + 5e-29 * std::pow(P, 5);
      EXPECT_NEAR(loaded.CalculateEfficiency(P), v2_eta, 1e-12 * std::abs(v2_eta)) << "P = " << P;
    }
  }

  TEST(PumpEfficiencyTest, FactoryCreation) {
    auto constant = CreatePumpEfficiencyModel(PumpEfficiencyType::CONSTANT);
    ASSERT_NE(constant, nullptr);
    EXPECT_EQ(constant->GetType(), PumpEfficiencyType::CONSTANT);

    auto vac = CreatePumpEfficiencyModel(PumpEfficiencyType::VACUUM_PUMP);
    ASSERT_NE(vac, nullptr);
    EXPECT_EQ(vac->GetType(), PumpEfficiencyType::VACUUM_PUMP);

    auto poly = CreatePumpEfficiencyModel(PumpEfficiencyType::POLYNOMIAL);
    ASSERT_NE(poly, nullptr);
    EXPECT_EQ(poly->GetType(), PumpEfficiencyType::POLYNOMIAL);
  }

  TEST(PumpEfficiencyTest, IPumpDefaultEfficiencyModel) {
    auto pump = CreatePump(PumpType::VACUUM_PUMP);
    ASSERT_NE(pump->efficiency_model_, nullptr);
    EXPECT_EQ(pump->efficiency_model_->GetType(), PumpEfficiencyType::CONSTANT);
  }

  // E-19: n/(n-1) in the polytropic work is infinite at n = 1, so a pump in
  // use must reject poly_index_ <= 1 at check time.
  TEST(PumpChecksTest, PolytropicIndexAtOrBelowOneIsRejected) {
    for (double n : {1.0, 0.9, -1.0}) {
      SCOPED_TRACE("poly_index = " + std::to_string(n));
      auto pump = CreatePump(PumpType::COMPRESSOR);
      pump->use_pump_ = true;
      pump->poly_index_ = n;
      std::string err;
      EXPECT_FALSE(pump->PerformChecks(err));
      EXPECT_NE(err.find("polytropic index must be greater than 1"), std::string::npos) << err;

      // An unused pump is never evaluated, so it is not rejected.
      pump->use_pump_ = false;
      std::string err_unused;
      EXPECT_TRUE(pump->PerformChecks(err_unused)) << err_unused;
    }

    auto pump = CreatePump(PumpType::VACUUM_PUMP);
    pump->use_pump_ = true;
    pump->poly_index_ = 1.001;
    std::string err;
    EXPECT_TRUE(pump->PerformChecks(err)) << err;
  }

  TEST(PumpChecksTest, BoundaryChecksRejectPolytropicIndexOfOne) {
    Fluid fluid;
    fluid.components_ = {"N2"};

    // Pump-driven stream (the pump sets the flow).
    PumpBoundaryStream pump_stream;
    pump_stream.SetDirection(BoundaryDirection::INLET);
    pump_stream.mole_fraction["N2"] = 1.0;
    std::string ok_err;
    EXPECT_TRUE(pump_stream.PerformChecks(ok_err, fluid)) << ok_err;
    pump_stream.pump_->poly_index_ = 1.0;
    std::string err;
    EXPECT_FALSE(pump_stream.PerformChecks(err, fluid));
    EXPECT_NE(err.find("polytropic index"), std::string::npos) << err;

    // Specified-flow stream with the optional power calculation enabled.
    BoundaryStream flow_stream;
    flow_stream.SetDirection(BoundaryDirection::INLET);
    flow_stream.mole_fraction["N2"] = 1.0;
    flow_stream.pump_ = CreatePump(PumpType::COMPRESSOR);
    flow_stream.pump_->use_pump_ = true;
    flow_stream.pump_->poly_index_ = 1.0;
    std::string flow_err;
    EXPECT_FALSE(flow_stream.PerformChecks(flow_err, fluid));
    EXPECT_NE(flow_err.find("polytropic index"), std::string::npos) << flow_err;

    // Pressure boundary with the optional power calculation enabled.
    PressureBoundaryCondition pressure_bc;
    pressure_bc.direction_ = BoundaryDirection::OUTLET;
    pressure_bc.pump_ = CreatePump(PumpType::VACUUM_PUMP);
    pressure_bc.pump_->use_pump_ = true;
    std::string p_ok_err;
    EXPECT_TRUE(pressure_bc.PerformChecks(p_ok_err, fluid)) << p_ok_err;
    pressure_bc.pump_->poly_index_ = 1.0;
    std::string p_err;
    EXPECT_FALSE(pressure_bc.PerformChecks(p_err, fluid));
    EXPECT_NE(p_err.find("polytropic index"), std::string::npos) << p_err;
  }

  // F-04: efficiencies are stored as fractions; a pump in use must have
  // 0 < eta <= 1 (> 0 and <= 100 %). eta <= 0 made the power infinite or
  // negative with no diagnostic.
  namespace {
    const char* kEfficiencyRangeMsg = "Pump efficiency must be > 0 and <= 100 %";

    std::unique_ptr<IPump> PumpWithConstantEfficiency(PumpType type, double eta) {
      auto pump = CreatePump(type);
      pump->use_pump_ = true;
      auto model = std::make_unique<ConstantPumpEfficiencyModel>();
      model->efficiency_ = eta;
      pump->efficiency_model_ = std::move(model);
      return pump;
    }
  } // namespace

  TEST(PumpChecksTest, ConstantEfficiencyOutsideZeroToHundredPercentIsRejected) {
    for (auto type : {PumpType::COMPRESSOR, PumpType::VACUUM_PUMP}) {
      for (double eta : {0.0, -0.5, 1.2, std::nan("")}) {
        SCOPED_TRACE("pump " + std::to_string(static_cast<int>(type)) + ", eta = " + std::to_string(eta));
        auto pump = PumpWithConstantEfficiency(type, eta);
        std::string err;
        EXPECT_FALSE(pump->PerformChecks(err));
        EXPECT_NE(err.find(kEfficiencyRangeMsg), std::string::npos) << err;
        EXPECT_NE(err.find("Constant efficiency model"), std::string::npos) << err;

        // An unused pump is never evaluated, so it is not rejected.
        pump->use_pump_ = false;
        std::string err_unused;
        EXPECT_TRUE(pump->PerformChecks(err_unused)) << err_unused;
      }
      // The message reports the value in percent too.
      {
        auto pump = PumpWithConstantEfficiency(type, 1.2);
        std::string err;
        EXPECT_FALSE(pump->PerformChecks(err));
        EXPECT_NE(err.find("1.2 (120 %)"), std::string::npos) << err;
      }
      for (double eta : {1e-3, 0.7, 1.0}) {
        SCOPED_TRACE("valid eta = " + std::to_string(eta));
        auto pump = PumpWithConstantEfficiency(type, eta);
        std::string err;
        EXPECT_TRUE(pump->PerformChecks(err)) << err;
        EXPECT_TRUE(err.empty()) << err;
      }
    }
  }

  TEST(PumpChecksTest, DefaultPumpsPassTheEfficiencyCheck) {
    for (auto type : {PumpType::COMPRESSOR, PumpType::VACUUM_PUMP}) {
      auto pump = CreatePump(type);
      pump->use_pump_ = true;
      std::string err;
      EXPECT_TRUE(pump->PerformChecks(err)) << err;
      // Default vacuum-pump correlation as well.
      pump->efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::VACUUM_PUMP);
      EXPECT_TRUE(pump->PerformChecks(err)) << err;
    }
  }

  TEST(PumpChecksTest, PolynomialEfficiencyWithDefaultZeroCoefficientsIsRejected) {
    for (auto type : {PumpType::COMPRESSOR, PumpType::VACUUM_PUMP}) {
      SCOPED_TRACE("pump " + std::to_string(static_cast<int>(type)));
      auto pump = CreatePump(type);
      pump->use_pump_ = true;
      pump->efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::POLYNOMIAL);
      std::string err;
      EXPECT_FALSE(pump->PerformChecks(err));
      EXPECT_NE(err.find(kEfficiencyRangeMsg), std::string::npos) << err;
      EXPECT_NE(err.find("Polynomial efficiency model gives 0 (0 %) at a suction pressure of 100000 Pa"),
                std::string::npos) << err;
    }
  }

  // The polynomial is checked at the pump's reference suction pressure: the
  // compressor's suction pressure, the vacuum pump's discharge pressure.
  TEST(PumpChecksTest, PolynomialEfficiencyIsCheckedAtTheReferencePressure) {
    auto make_poly = [] {
      auto poly = std::make_unique<PolyPumpEfficiencyModel>();
      poly->a_ = 0.4;
      poly->b_ = 1e-6; // eta = 0.4 + 1e-6 P: 0.5 at 1 bar, 1.1 at 7 bar
      return poly;
    };

    Compressor compressor;
    compressor.use_pump_ = true;
    compressor.efficiency_model_ = make_poly();
    compressor.suction_pressure_ = 1e5;
    std::string err;
    EXPECT_TRUE(compressor.PerformChecks(err)) << err;
    compressor.suction_pressure_ = 7e5;
    EXPECT_FALSE(compressor.PerformChecks(err));
    EXPECT_NE(err.find("1.1 (110 %) at a suction pressure of 700000 Pa"), std::string::npos) << err;

    VacuumPump vacuum;
    vacuum.use_pump_ = true;
    vacuum.efficiency_model_ = make_poly();
    vacuum.discharge_pressure_ = 1e5;
    std::string vac_err;
    EXPECT_TRUE(vacuum.PerformChecks(vac_err)) << vac_err;
    vacuum.discharge_pressure_ = 7e5;
    EXPECT_FALSE(vacuum.PerformChecks(vac_err));
    EXPECT_NE(vac_err.find(kEfficiencyRangeMsg), std::string::npos) << vac_err;
  }

  // The vacuum-pump correlation lies in (0, eta_max) for every p > 0, so its
  // parameters are checked rather than one pressure.
  TEST(PumpChecksTest, VacuumPumpEfficiencyParametersAreChecked) {
    struct Case { double eta_max; double alpha; bool valid; };
    for (auto c : {Case{0.8, 19.55, true}, Case{1.0, 5.0, true}, Case{0.0, 19.55, false},
                   Case{1.2, 19.55, false}, Case{-0.5, 19.55, false}, Case{0.8, 0.0, false},
                   Case{0.8, -1.0, false}}) {
      SCOPED_TRACE("eta_max = " + std::to_string(c.eta_max) + ", alpha = " + std::to_string(c.alpha));
      auto pump = CreatePump(PumpType::VACUUM_PUMP);
      pump->use_pump_ = true;
      auto model = std::make_unique<VacuumPumpEfficiencyModel>();
      model->eta_max_ = c.eta_max;
      model->alpha_ = c.alpha;
      pump->efficiency_model_ = std::move(model);
      std::string err;
      EXPECT_EQ(pump->PerformChecks(err), c.valid) << err;
      if (!c.valid) {
        EXPECT_NE(err.find(kEfficiencyRangeMsg), std::string::npos) << err;
        EXPECT_NE(err.find("Vacuum Pump efficiency model"), std::string::npos) << err;
      }
    }
  }

  TEST(PumpChecksTest, BoundaryChecksRejectZeroEfficiencyAndNameTheStream) {
    Fluid fluid;
    fluid.components_ = {"N2"};

    PumpBoundaryStream pump_stream;
    pump_stream.name_ = "Compressor feed";
    pump_stream.SetDirection(BoundaryDirection::INLET);
    pump_stream.mole_fraction["N2"] = 1.0;
    pump_stream.pump_->efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::POLYNOMIAL);
    std::string err;
    EXPECT_FALSE(pump_stream.PerformChecks(err, fluid));
    EXPECT_NE(err.find(kEfficiencyRangeMsg), std::string::npos) << err;
    EXPECT_NE(err.find("stream 'Compressor feed'"), std::string::npos) << err;

    BoundaryStream flow_stream;
    flow_stream.name_ = "Feed";
    flow_stream.SetDirection(BoundaryDirection::INLET);
    flow_stream.mole_fraction["N2"] = 1.0;
    flow_stream.pump_ = PumpWithConstantEfficiency(PumpType::COMPRESSOR, 0.0);
    std::string flow_err;
    EXPECT_FALSE(flow_stream.PerformChecks(flow_err, fluid));
    EXPECT_NE(flow_err.find(kEfficiencyRangeMsg), std::string::npos) << flow_err;
    EXPECT_NE(flow_err.find("stream 'Feed'"), std::string::npos) << flow_err;

    PressureBoundaryCondition pressure_bc;
    pressure_bc.direction_ = BoundaryDirection::OUTLET;
    pressure_bc.pump_ = PumpWithConstantEfficiency(PumpType::VACUUM_PUMP, 0.0);
    std::string p_err;
    EXPECT_FALSE(pressure_bc.PerformChecks(p_err, fluid));
    EXPECT_NE(p_err.find(kEfficiencyRangeMsg), std::string::npos) << p_err;
  }

  // The checks are the user-facing guard; the maths must still stay finite
  // and non-negative if an efficiency <= 0 reaches it.
  TEST(PumpChecksTest, PowerAndFlowStayFiniteForNonPositiveEfficiency) {
    for (double eta : {0.0, -0.5, std::nan("")}) {
      SCOPED_TRACE("eta = " + std::to_string(eta));
      auto compressor = PumpWithConstantEfficiency(PumpType::COMPRESSOR, eta);
      auto* comp = dynamic_cast<Compressor*>(compressor.get());
      comp->suction_pressure_ = 1e5;
      const double power = comp->CalculatePower(3e5, 0.1, 300.0);
      EXPECT_TRUE(std::isfinite(power));
      EXPECT_GE(power, 0.0);
      const double flow = comp->CalculateMolarFlow(3e5, 1000.0, 300.0);
      EXPECT_TRUE(std::isfinite(flow));
      EXPECT_GE(flow, 0.0);

      auto vacuum = PumpWithConstantEfficiency(PumpType::VACUUM_PUMP, eta);
      const double vac_power = vacuum->CalculatePower(0.2e5, 0.1, 300.0);
      EXPECT_TRUE(std::isfinite(vac_power));
      EXPECT_GE(vac_power, 0.0);
      const double vac_flow = vacuum->CalculateMolarFlow(0.2e5, 1000.0, 300.0);
      EXPECT_TRUE(std::isfinite(vac_flow));
      EXPECT_GE(vac_flow, 0.0);
    }
  }

} // end adsim namespace
