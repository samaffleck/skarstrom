/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/expression_param.h"

#include <cmath>
#include <optional>
#include <string>


namespace adsim {

	static const char* BoundaryValueTypeList[] = {
		"Constant",		// 0
		"Linear",		// 1
		"Exponential",	// 2
		"Log-linear",	// 3
		"Double exponential"	// 4
	};

	enum class BoundaryValueType {
		CONSTANT = 0,
		LINEAR = 1,
		EXPONENTIAL = 2,
		LOG_LINEAR = 3,
		DOUBLE_EXPONENTIAL = 4
	};

	class IBoundaryValue {
	public:
		virtual ~IBoundaryValue() = default;

		virtual BoundaryValueType GetType() = 0;
		virtual double GetValue(double current_time) = 0;
		// Called at every step start by pressure-driven parents (pressure BC,
		// valve stream): the start value is the current bed-end pressure and,
		// unless the profile has its own ramp time (use_step_time_ off), the
		// ramp spans the whole step. These are RUNTIME values held in
		// transient members (never serialised) — they must not overwrite the
		// user's ExpressionParams, which are what the deck saves.
		virtual void InitialiseStartValue(double start_value) = 0;
		virtual void InitialiseStepTime(double step_time) {};

		// Re-evaluate any ExpressionParam fields using the provided context.
		// Call once before starting the solver so user-variable references are fresh.
		virtual void EvaluateExpressions(const EvalContext& ctx) {}

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar) {}
	};

	class ConstantBoundaryValue : public IBoundaryValue {
	public:
		ConstantBoundaryValue() = default;
		ConstantBoundaryValue(double value) : constant_value_(value) {}

		BoundaryValueType GetType() override { return BoundaryValueType::CONSTANT; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override {}
		void EvaluateExpressions(const EvalContext& ctx) override { constant_value_.Evaluate(ctx); }

		ExpressionParam constant_value_{1e5};

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(cereal::virtual_base_class<IBoundaryValue>(this));
			if (version >= 2) {
				ar(cereal::make_nvp("Constant_Value", constant_value_));
			} else {
				// v1: plain double on disk — read into a temp, then set ExpressionParam
				double v = 1e5;
				ar(cereal::make_nvp("Constant_Value", v));
				constant_value_ = v;
			}
		}
	};

	class LinearBoundaryValue : public IBoundaryValue {
	public:
		BoundaryValueType GetType() override { return BoundaryValueType::LINEAR; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override { runtime_start_value_ = start_value; }
		void InitialiseStepTime(double step_time) override {
			if (use_step_time_) runtime_ramp_time_ = step_time;
			else runtime_ramp_time_.reset();
		}
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			ramp_time_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		ExpressionParam ramp_time_{1};
		// Pressure-driven parents only (pressure BC, valve stream): true = the
		// ramp lasts the whole step (ramp_time_ is ignored), false = it lasts
		// ramp_time_, which must not exceed the step time (CheckPressureProfile).
		// Flow and temperature profiles always use ramp_time_.
		bool use_step_time_ = true;

	private:
		// Transient runtime overrides set by InitialiseStartValue / InitialiseStepTime
		// (pressure/valve parents only); when present GetValue uses them instead of
		// start_value_ / ramp_time_. Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_ramp_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Ramp_Time", ramp_time_)
			);
			// v4: the ramp time of a pressure profile can differ from the step
			// time; older decks always ramped over the whole step.
			if (version >= 4) ar(cereal::make_nvp("Use_Step_Time", use_step_time_));
			else use_step_time_ = true;
		}
	};

	// A first-order approach from the start value v0 to the end value v1 whose
	// remaining change has shrunk by the factor F after the ramp time t_r:
	//   v(t) = v1 + (v0 - v1) F^(-t / t_r) = v1 + (v0 - v1) e^(-ln(F) t / t_r).
	// The default F = 100 has covered 99 % of the change at t_r (the fixed
	// factor before class version 4). F must be greater than 1 (CheckProfileFactor).
	class ExponentialBoundaryValue : public IBoundaryValue {
	public:
		BoundaryValueType GetType() override { return BoundaryValueType::EXPONENTIAL; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override { runtime_start_value_ = start_value; }
		void InitialiseStepTime(double step_time) override {
			if (use_step_time_) runtime_step_time_ = step_time;
			else runtime_step_time_.reset();
		}
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			step_time_.Evaluate(ctx);
			factor_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		double gain_ = 0.95;        // unused; kept so the deck field round-trips
		ExpressionParam step_time_{1};   // [s] the ramp time t_r
		ExpressionParam factor_{100.0};  // [-] F: the change left at t_r is 1/F of the initial
		// See LinearBoundaryValue::use_step_time_.
		bool use_step_time_ = true;

	private:
		// Transient runtime overrides (see LinearBoundaryValue). Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_step_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Proportional_Gain", gain_),
				cereal::make_nvp("Step_Time", step_time_)
			);
			// v4: a user factor (was fixed at 100) and a pressure ramp time that
			// can differ from the step time.
			if (version >= 4) {
				ar(cereal::make_nvp("Factor", factor_),
				   cereal::make_nvp("Use_Step_Time", use_step_time_));
			} else {
				factor_ = 100.0;
				use_step_time_ = true;
			}
		}
	};

	// A two-term exponential decay towards the end value: the initial excess
	// over the end value splits into a share w that decays at rate_1 and the
	// rest at rate_2,
	//   v(t) = v1 + (v0 - v1) [w e^(-r1 t) + (1 - w) e^(-r2 t)].
	// The form Estupinan Perez et al. (2019) fitted to the evacuation pressure
	// of their VSA rig (Sect. 4.3, SI Table S2): a fast term for the column
	// gas and a slow one for the pump and piping. The rates are given directly
	// in 1/s (not a 99 % time like the exponential profile) so a fitted value
	// goes in as printed; a rate of zero freezes its term at the start value.
	class DoubleExponentialBoundaryValue : public IBoundaryValue {
	public:
		BoundaryValueType GetType() override { return BoundaryValueType::DOUBLE_EXPONENTIAL; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override { runtime_start_value_ = start_value; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			rate_1_.Evaluate(ctx);
			rate_2_.Evaluate(ctx);
			weight_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		ExpressionParam rate_1_{0.02};   // [1/s] rate of the slow term
		ExpressionParam rate_2_{0.1};    // [1/s] rate of the fast term
		ExpressionParam weight_{0.5};    // [-] share of the initial excess in the slow term

	private:
		std::optional<double> runtime_start_value_;  // not serialised (see LinearBoundaryValue)

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Rate_1", rate_1_),
				cereal::make_nvp("Rate_2", rate_2_),
				cereal::make_nvp("Weight", weight_)
			);
		}
	};

	// A geometric ramp: the value is multiplied by the same factor every second,
	// v(t) = v0 * (v1 / v0)^(t / ramp), so ln(v) rises linearly and the rate of
	// rise relative to the value, d(ln v)/dt = ln(v1 / v0) / ramp, is constant.
	// For a pressurisation from vacuum this spreads the relative pressure rise
	// (which sets the inlet pressure gradient) evenly over the step, where a
	// linear ramp concentrates it at the start. Holds the end value after the
	// ramp, like LinearBoundaryValue. The ratio needs both values positive: when
	// either is zero or negative the ramp is linear instead.
	class LogLinearBoundaryValue : public IBoundaryValue {
	public:
		BoundaryValueType GetType() override { return BoundaryValueType::LOG_LINEAR; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override { runtime_start_value_ = start_value; }
		void InitialiseStepTime(double step_time) override {
			if (use_step_time_) runtime_ramp_time_ = step_time;
			else runtime_ramp_time_.reset();
		}
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			ramp_time_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		ExpressionParam ramp_time_{1};
		// See LinearBoundaryValue::use_step_time_.
		bool use_step_time_ = true;

	private:
		// Transient runtime overrides (see LinearBoundaryValue). Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_ramp_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Ramp_Time", ramp_time_)
			);
			// v2: see LinearBoundaryValue (v4).
			if (version >= 2) ar(cereal::make_nvp("Use_Step_Time", use_step_time_));
			else use_step_time_ = true;
		}
	};

	// ---- checks -----------------------------------------------------------------
	// The ramp time a linear, exponential or log-linear profile holds itself
	// (ramp_time_ / step_time_), or nullopt for a constant or double exponential.
	inline std::optional<double> ProfileRampTime(const IBoundaryValue* v) {
		if (auto* l = dynamic_cast<const LinearBoundaryValue*>(v)) return l->ramp_time_.Get();
		if (auto* e = dynamic_cast<const ExponentialBoundaryValue*>(v)) return e->step_time_.Get();
		if (auto* g = dynamic_cast<const LogLinearBoundaryValue*>(v)) return g->ramp_time_.Get();
		return std::nullopt;
	}

	// Whether a pressure-driven profile ramps over the whole step (true for the
	// kinds without a ramp time, which never use one).
	inline bool ProfileUsesStepTime(const IBoundaryValue* v) {
		if (auto* l = dynamic_cast<const LinearBoundaryValue*>(v)) return l->use_step_time_;
		if (auto* e = dynamic_cast<const ExponentialBoundaryValue*>(v)) return e->use_step_time_;
		if (auto* g = dynamic_cast<const LogLinearBoundaryValue*>(v)) return g->use_step_time_;
		return true;
	}

	// Set the switch on a linear, exponential or log-linear profile (no-op for
	// the other kinds). Takes effect at the next step start.
	inline void SetProfileUsesStepTime(IBoundaryValue* v, bool use_step_time) {
		if (auto* l = dynamic_cast<LinearBoundaryValue*>(v)) l->use_step_time_ = use_step_time;
		else if (auto* e = dynamic_cast<ExponentialBoundaryValue*>(v)) e->use_step_time_ = use_step_time;
		else if (auto* g = dynamic_cast<LogLinearBoundaryValue*>(v)) g->use_step_time_ = use_step_time;
	}

	// Whether the kind has a ramp time (and so the switch): linear, exponential, log-linear.
	inline bool ProfileHasRampTime(const IBoundaryValue* v) { return ProfileRampTime(v).has_value(); }

	// An exponential profile's factor F must exceed 1 (F = 1 never moves off the
	// start value, F < 1 grows without bound). Any parent: pressure, valve,
	// flow or inlet temperature. `what` names the profile in the message.
	inline bool CheckProfileFactor(const IBoundaryValue* v, const std::string& what,
	                               std::string& error_message) {
		auto* e = dynamic_cast<const ExponentialBoundaryValue*>(v);
		if (!e) return true;
		const double f = e->factor_.Get();
		if (std::isfinite(f) && f > 1.0) return true;
		error_message += "The exponential factor of the " + what + " is " + std::to_string(f)
		                 + "; it must be greater than 1 (100 covers 99 % of the change by the ramp time). ";
		return false;
	}

	// A pressure-driven profile (pressure boundary, valve stream) with its own
	// ramp time: the ramp must fit in the step (and not be negative). A small
	// tolerance absorbs the rounding of a FRACTION-mode step time.
	inline bool CheckPressureProfile(const IBoundaryValue* v, double step_time, const std::string& what,
	                                 std::string& error_message) {
		if (!v) return true;
		if (!CheckProfileFactor(v, what, error_message)) return false;
		if (ProfileUsesStepTime(v)) return true;
		const auto ramp = ProfileRampTime(v);
		if (!ramp) return true;
		if (!std::isfinite(*ramp) || *ramp < 0.0) {
			error_message += "The ramp time of the " + what + " is " + std::to_string(*ramp)
			                 + " s; it must be zero or positive. ";
			return false;
		}
		if (*ramp > step_time * (1.0 + 1e-9) + 1e-9) {
			error_message += "The ramp time of the " + what + " (" + std::to_string(*ramp)
			                 + " s) is longer than the step (" + std::to_string(step_time)
			                 + " s); shorten it or ramp over the step time. ";
			return false;
		}
		return true;
	}

	/* Factory function */
	inline std::unique_ptr<IBoundaryValue> CreateBoundaryValue(BoundaryValueType type) {
		switch (type) {
		case BoundaryValueType::CONSTANT: return std::make_unique<ConstantBoundaryValue>();
		case BoundaryValueType::LINEAR: return std::make_unique<LinearBoundaryValue>();
		case BoundaryValueType::EXPONENTIAL: return std::make_unique<ExponentialBoundaryValue>();
		case BoundaryValueType::LOG_LINEAR: return std::make_unique<LogLinearBoundaryValue>();
		case BoundaryValueType::DOUBLE_EXPONENTIAL: return std::make_unique<DoubleExponentialBoundaryValue>();
		default: return nullptr;
		}
	}

	// ---- inlet-temperature profiles -------------------------------------------
	// An inlet temperature is a constant ExpressionParam plus an optional time
	// profile (the same Linear / Exponential / Log-linear classes the flow boundaries use,
	// with the user's own start value and ramp time — no runtime overrides).
	// nullptr or a ConstantBoundaryValue means "constant": the ExpressionParam.
	// Linear holds its end value after the ramp (no extrapolation, cf. bug E-05)
	// and Exponential saturates on its own.
	inline double EvaluateInletTemperatureProfile(const ExpressionParam& constant,
	                                              IBoundaryValue* profile, double time) {
		if (!profile || profile->GetType() == BoundaryValueType::CONSTANT) return constant.Get();
		return profile->GetValue(time);
	}

	inline BoundaryValueType InletTemperatureProfileType(IBoundaryValue* profile) {
		return profile ? profile->GetType() : BoundaryValueType::CONSTANT;
	}

	// A fresh inlet-temperature profile of `type`: nullptr for CONSTANT, else a
	// ramp that starts and ends at `temperature` [K] over `ramp_time` [s] (so
	// switching the type changes nothing until the user edits the values).
	inline std::unique_ptr<IBoundaryValue> CreateInletTemperatureProfile(
	    BoundaryValueType type, double temperature, double ramp_time) {
		if (type == BoundaryValueType::LINEAR) {
			auto p = std::make_unique<LinearBoundaryValue>();
			p->start_value_ = temperature;
			p->end_value_ = temperature;
			p->ramp_time_ = ramp_time;
			return p;
		}
		if (type == BoundaryValueType::EXPONENTIAL) {
			auto p = std::make_unique<ExponentialBoundaryValue>();
			p->start_value_ = temperature;
			p->end_value_ = temperature;
			p->step_time_ = ramp_time;
			return p;
		}
		if (type == BoundaryValueType::LOG_LINEAR) {
			auto p = std::make_unique<LogLinearBoundaryValue>();
			p->start_value_ = temperature;
			p->end_value_ = temperature;
			p->ramp_time_ = ramp_time;
			return p;
		}
		if (type == BoundaryValueType::DOUBLE_EXPONENTIAL) {
			// Rates such that the slow term has covered 99 % of its change by
			// ramp_time and the fast one has settled well before.
			auto p = std::make_unique<DoubleExponentialBoundaryValue>();
			p->start_value_ = temperature;
			p->end_value_ = temperature;
			const double slow = ramp_time > 0.0 ? std::log(100.0) / ramp_time : 0.0;
			p->rate_1_ = slow;
			p->rate_2_ = 5.0 * slow;
			p->weight_ = 0.5;
			return p;
		}
		return nullptr;
	}

} // end adsim namespace
