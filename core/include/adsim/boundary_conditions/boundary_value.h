/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/expression_param.h"

#include <optional>


namespace adsim {

	static const char* BoundaryValueTypeList[] = {
		"Constant",		// 0
		"Linear",		// 1
		"Exponential",	// 2
		"Log-linear"	// 3
	};

	enum class BoundaryValueType {
		CONSTANT = 0,
		LINEAR = 1,
		EXPONENTIAL = 2,
		LOG_LINEAR = 3
	};

	class IBoundaryValue {
	public:
		virtual ~IBoundaryValue() = default;

		virtual BoundaryValueType GetType() = 0;
		virtual double GetValue(double current_time) = 0;
		// Called at every step start by pressure-driven parents (pressure BC,
		// valve stream): the start value is the current bed-end pressure and
		// the ramp spans the whole step. These are RUNTIME values held in
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
		void InitialiseStepTime(double step_time) override { runtime_ramp_time_ = step_time; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			ramp_time_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		ExpressionParam ramp_time_{1};

	private:
		// Transient runtime overrides set by InitialiseStartValue / InitialiseStepTime
		// (pressure/valve parents only); when present GetValue uses them instead of
		// start_value_ / ramp_time_. Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_ramp_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Ramp_Time", ramp_time_)
			);
		}
	};

	class ExponentialBoundaryValue : public IBoundaryValue {
	public:
		BoundaryValueType GetType() override { return BoundaryValueType::EXPONENTIAL; }
		double GetValue(double time) override;
		void InitialiseStartValue(double start_value) override { runtime_start_value_ = start_value; }
		void InitialiseStepTime(double step_time) override { runtime_step_time_ = step_time; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			step_time_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		double gain_ = 0.95;        // internal — not user-editable
		ExpressionParam step_time_{1};

	private:
		// Transient runtime overrides (see LinearBoundaryValue). Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_step_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Proportional_Gain", gain_),
				cereal::make_nvp("Step_Time", step_time_)
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
		void InitialiseStepTime(double step_time) override { runtime_ramp_time_ = step_time; }
		void EvaluateExpressions(const EvalContext& ctx) override {
			end_value_.Evaluate(ctx);
			start_value_.Evaluate(ctx);
			ramp_time_.Evaluate(ctx);
		}

		ExpressionParam end_value_{1e5};
		ExpressionParam start_value_{1};
		ExpressionParam ramp_time_{1};

	private:
		// Transient runtime overrides (see LinearBoundaryValue). Not serialised.
		std::optional<double> runtime_start_value_;
		std::optional<double> runtime_ramp_time_;

		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const /*version*/) {
			ar(
				cereal::virtual_base_class<IBoundaryValue>(this),
				cereal::make_nvp("End_Value", end_value_),
				cereal::make_nvp("Start_Value", start_value_),
				cereal::make_nvp("Ramp_Time", ramp_time_)
			);
		}
	};

	/* Factory function */
	inline std::unique_ptr<IBoundaryValue> CreateBoundaryValue(BoundaryValueType type) {
		switch (type) {
		case BoundaryValueType::CONSTANT: return std::make_unique<ConstantBoundaryValue>();
		case BoundaryValueType::LINEAR: return std::make_unique<LinearBoundaryValue>();
		case BoundaryValueType::EXPONENTIAL: return std::make_unique<ExponentialBoundaryValue>();
		case BoundaryValueType::LOG_LINEAR: return std::make_unique<LogLinearBoundaryValue>();
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
		return nullptr;
	}

} // end adsim namespace
