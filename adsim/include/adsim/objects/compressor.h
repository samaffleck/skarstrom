/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <memory>
#include <string>

// Cereal includes
#include <cereal/access.hpp>

// Adsim includes
#include "adsim/math/general_functions.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  class Compressor;

  static const char* CompressorEnergyList[] = {
		"Isothermal",	// 0 
		"Adiabatic"	// 1
	};

	enum class CompressorEnergyType {
		ISOTHERMAL = 0,
		ADIABATIC = 1
	};

  static const char* PumpList[] = {
		"Compressor",	// 0 
		"Vacuum Pump"	// 1
	};

	enum class PumpType {
		COMPRESSOR = 0,
		VACUUM_PUMP = 1
	};

  static const char* PumpPowerTypeList[] = {
		"Constant",	// 0 
		"Polynomial"	// 1
	};

	enum class PumpPowerType {
		CONSTANT = 0,
		POLYNOMIAL = 1
	};

  struct PumpPowerModelParams {
    double pressure;
  };

  class IPumpPowerModel {
  public:
    virtual ~IPumpPowerModel() = default;
    // Outlet pressure for a compressor, inlet pressure for a vacuum pump
    virtual double CalculatePower(const PumpPowerModelParams& params) const = 0;
    virtual PumpPowerType GetType() const = 0;
    virtual void EvaluateExpressions(const EvalContext& ctx) {}

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class ConstantPumpPowerModel : public IPumpPowerModel {
  public:
    double CalculatePower(const PumpPowerModelParams& params) const override;
    PumpPowerType GetType() const override { return PumpPowerType::CONSTANT; }
    void EvaluateExpressions(const EvalContext& ctx) override { power_.Evaluate(ctx); }

    ExpressionParam power_{0.0};

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IPumpPowerModel>(this),
        cereal::make_nvp("Power", power_)
			);
		}
  };

  class PolyPumpPowerModel : public IPumpPowerModel {
  public:
    double CalculatePower(const PumpPowerModelParams& params) const override;
    PumpPowerType GetType() const override { return PumpPowerType::POLYNOMIAL; }
    void EvaluateExpressions(const EvalContext& ctx) override { poly_.EvaluateExpressions(ctx); }

    Polynomial5 poly_{};

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IPumpPowerModel>(this),
        cereal::make_nvp("Polynomial_Function", poly_)
			);
		}
  };

  // ---------------------------------------------------------------------------
  // Pump efficiency models
  // ---------------------------------------------------------------------------
  static const char* PumpEfficiencyTypeList[] = {
    "Constant",       // 0
    "Vacuum Pump",    // 1
    "Polynomial"      // 2
  };

  enum class PumpEfficiencyType {
    CONSTANT = 0,
    VACUUM_PUMP = 1,
    POLYNOMIAL = 2
  };

  // Efficiencies are stored as FRACTIONS (1.0 = 100 %). A valid efficiency is
  // 0 < eta <= 1, i.e. > 0 and <= 100 %: eta <= 0 makes the pump power
  // infinite or negative, so IPump::PerformChecks rejects it.
  class IPumpEfficiencyModel {
  public:
    virtual ~IPumpEfficiencyModel() = default;
    // `pressure` is the pump's suction (inlet) pressure P1 for both pump kinds.
    virtual double CalculateEfficiency(double pressure) const = 0;
    virtual PumpEfficiencyType GetType() const = 0;
    virtual void EvaluateExpressions(const EvalContext& ctx) {}
    // Checks 0 < eta <= 1 (> 0 and <= 100 %). The default evaluates the model
    // at `reference_pressure` (see IPump::EfficiencyReferencePressure); a model
    // whose range is known for every pressure overrides it. Appends a message
    // to error_msg and returns false on failure.
    virtual bool PerformChecks(double reference_pressure, std::string& error_msg) const;

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {}
  };

  class ConstantPumpEfficiencyModel : public IPumpEfficiencyModel {
  public:
    double CalculateEfficiency(double pressure) const override;
    PumpEfficiencyType GetType() const override { return PumpEfficiencyType::CONSTANT; }
    void EvaluateExpressions(const EvalContext& ctx) override { efficiency_.Evaluate(ctx); }

    ExpressionParam efficiency_{0.7};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IPumpEfficiencyModel>(this),
        cereal::make_nvp("Efficiency", efficiency_)
      );
    }
  };

  // eta = (eta_max * alpha * p/1e5) / (1 + alpha * p/1e5)
  // where p is the pump's inlet (suction-side) pressure
  // (eta_max_ and alpha_ are not user-editable in the GUI today, but expose
  //  them via ExpressionParam anyway for consistency / future use)
  class VacuumPumpEfficiencyModel : public IPumpEfficiencyModel {
  public:
    double CalculateEfficiency(double pressure) const override;
    PumpEfficiencyType GetType() const override { return PumpEfficiencyType::VACUUM_PUMP; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      eta_max_.Evaluate(ctx); alpha_.Evaluate(ctx);
    }

    // 0 < eta < eta_max for every p > 0 iff 0 < eta_max <= 1 and alpha > 0,
    // so the parameters are checked instead of sampling a pressure.
    bool PerformChecks(double reference_pressure, std::string& error_msg) const override;

    ExpressionParam eta_max_{0.8};
    ExpressionParam alpha_{19.55};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IPumpEfficiencyModel>(this),
        cereal::make_nvp("Eta_Max", eta_max_),
        cereal::make_nvp("Alpha", alpha_)
      );
    }
  };

  // eta = a + b*P + c*P^3 + d*P^4 + e*P^5
  class PolyPumpEfficiencyModel : public IPumpEfficiencyModel {
  public:
    double CalculateEfficiency(double pressure) const override;
    PumpEfficiencyType GetType() const override { return PumpEfficiencyType::POLYNOMIAL; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      a_.Evaluate(ctx); b_.Evaluate(ctx); c_.Evaluate(ctx);
      d_.Evaluate(ctx); e_.Evaluate(ctx);
    }

    ExpressionParam a_{0.0};
    ExpressionParam b_{0.0};
    ExpressionParam c_{0.0};
    ExpressionParam d_{0.0};
    ExpressionParam e_{0.0};

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::virtual_base_class<IPumpEfficiencyModel>(this),
        cereal::make_nvp("a", a_),
        cereal::make_nvp("b", b_),
        cereal::make_nvp("c", c_),
        cereal::make_nvp("d", d_),
        cereal::make_nvp("e", e_)
      );
    }
  };

  struct EnergyModelParams {
    double P1{};
    double P2{};
    double n{};
    double T1{};
  };

  class ICompressorEnergyModel {
  public:
    virtual ~ICompressorEnergyModel() = default;
    virtual double CalculateOutletTemperature(const EnergyModelParams& params) const = 0;
    virtual CompressorEnergyType GetType() const = 0;
  
  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
		}
  };
  
  class IsothermalCompressorEnergyModel : public ICompressorEnergyModel {
  public:
    double CalculateOutletTemperature(const EnergyModelParams& params) const override;
    CompressorEnergyType GetType() const override { return CompressorEnergyType::ISOTHERMAL; }

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<ICompressorEnergyModel>(this)
			);
		}
  };

  class AdiabaticCompressorEnergyModel : public ICompressorEnergyModel {
  public:
    double CalculateOutletTemperature(const EnergyModelParams& params) const override;
    CompressorEnergyType GetType() const override { return CompressorEnergyType::ADIABATIC; }

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<ICompressorEnergyModel>(this)
			);
		}
  };

  class IPump {
  public:
    IPump();
    virtual ~IPump() = default;
    virtual double CalculatePower(double pressure, double mol_flow, double temperature) const = 0;
    virtual double CalculateMolarFlow(double pressure, double power, double temperature) const = 0;
    virtual double CalculateOutletTemperature(double pressure, double temperature) const = 0;
    virtual PumpType GetType() const = 0;
    // Validates a pump that is in use (use_pump_): poly_index_ must be > 1,
    // since the polytropic work n/(n-1) is infinite at n = 1, and the
    // efficiency must be > 0 and <= 100 % (CheckEfficiency). Appends a
    // message to error_msg and returns false on failure.
    bool PerformChecks(std::string& error_msg) const;
    // The efficiency part of PerformChecks (regardless of use_pump_); also
    // used read-only by the desktop GUI to warn next to the efficiency inputs.
    bool CheckEfficiency(std::string& error_msg) const;
    // Suction pressure P1 at which a pressure-dependent efficiency model is
    // checked: the compressor's suction pressure (exactly the P1 it uses,
    // unless its stream is connected), or the vacuum pump's discharge
    // pressure (the top of its suction range; the bed-side P1 varies during
    // the step and is only known at run time).
    virtual double EfficiencyReferencePressure() const = 0;
    virtual void EvaluateExpressions(const EvalContext& ctx) {
      poly_index_.Evaluate(ctx);
      if (efficiency_model_) efficiency_model_->EvaluateExpressions(ctx);
    }

    bool use_pump_ = false;
    ExpressionParam poly_index_{1.4};
    std::unique_ptr<ICompressorEnergyModel> energy_eq_;
    std::unique_ptr<IPumpEfficiencyModel> efficiency_model_;

  private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Use_Pump", use_pump_),
        cereal::make_nvp("Polytropic_Index", poly_index_),
        cereal::make_nvp("Energy_Equation", energy_eq_),
        cereal::make_nvp("Efficiency_Model", efficiency_model_)
      );
    }
  };

  class Compressor : public IPump {
  public:
    double CalculatePower(double outlet_pressure, double mol_flow, double inlet_temperature) const override;
    double CalculateMolarFlow(double outlet_pressure, double power, double inlet_temperature) const override;
    double CalculateOutletTemperature(double outlet_pressure, double inlet_temperature) const override;
    // Same as above but with an explicit suction pressure in place of
    // suction_pressure_ (e.g. a connected upstream stream's pressure). Used so
    // callers never have to overwrite the user's ExpressionParam.
    double CalculateMolarFlow(double outlet_pressure, double power, double inlet_temperature, double suction_pressure) const;
    double CalculateOutletTemperature(double outlet_pressure, double inlet_temperature, double suction_pressure) const;
    PumpType GetType() const override { return PumpType::COMPRESSOR; }
    double EfficiencyReferencePressure() const override { return suction_pressure_; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      IPump::EvaluateExpressions(ctx);
      suction_pressure_.Evaluate(ctx);
    }

    ExpressionParam suction_pressure_{1e5};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::virtual_base_class<IPump>(this),
				cereal::make_nvp("Suction_Pressure", suction_pressure_)
			);
		}
  };

  class VacuumPump : public IPump {
  public:
    double CalculatePower(double inlet_pressure, double mol_flow, double inlet_temperature) const override;
    double CalculateMolarFlow(double inlet_pressure, double power, double inlet_temperature) const override;
    double CalculateOutletTemperature(double inlet_pressure, double inlet_temperature) const override;
    // Same as above but with an explicit discharge pressure in place of
    // discharge_pressure_ (e.g. a connected downstream stream's pressure).
    double CalculateMolarFlow(double inlet_pressure, double power, double inlet_temperature, double discharge_pressure) const;
    PumpType GetType() const override { return PumpType::VACUUM_PUMP; }
    double EfficiencyReferencePressure() const override { return discharge_pressure_; }
    void EvaluateExpressions(const EvalContext& ctx) override {
      IPump::EvaluateExpressions(ctx);
      discharge_pressure_.Evaluate(ctx);
    }

    ExpressionParam discharge_pressure_{1e5};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
        cereal::virtual_base_class<IPump>(this),
				cereal::make_nvp("Discharge_Pressure", discharge_pressure_)
			);
		}
  };

  /* Factory function */
	inline std::unique_ptr<ICompressorEnergyModel> CreateCompressorEnergyModel(CompressorEnergyType type) {
		switch (type) {
		case CompressorEnergyType::ISOTHERMAL: return std::make_unique<IsothermalCompressorEnergyModel>();
		case CompressorEnergyType::ADIABATIC: return std::make_unique<AdiabaticCompressorEnergyModel>();
		default: return nullptr;
		}
	}

	inline std::unique_ptr<IPump> CreatePump(PumpType type) {
		switch (type) {
		case PumpType::COMPRESSOR: return std::make_unique<Compressor>();
		case PumpType::VACUUM_PUMP: return std::make_unique<VacuumPump>();
		default: return nullptr;
		}
	}

	inline std::unique_ptr<IPumpPowerModel> CreatePumpPowerModel(PumpPowerType type) {
		switch (type) {
      case PumpPowerType::CONSTANT: return std::make_unique<ConstantPumpPowerModel>();
      case PumpPowerType::POLYNOMIAL: return std::make_unique<PolyPumpPowerModel>();
      default: return nullptr;
		}
	}

	inline std::unique_ptr<IPumpEfficiencyModel> CreatePumpEfficiencyModel(PumpEfficiencyType type) {
		switch (type) {
      case PumpEfficiencyType::CONSTANT: return std::make_unique<ConstantPumpEfficiencyModel>();
      case PumpEfficiencyType::VACUUM_PUMP: return std::make_unique<VacuumPumpEfficiencyModel>();
      case PumpEfficiencyType::POLYNOMIAL: return std::make_unique<PolyPumpEfficiencyModel>();
      default: return nullptr;
		}
	}

} // end adsim namespace
