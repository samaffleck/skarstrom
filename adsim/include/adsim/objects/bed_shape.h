/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <cassert>
#include <variant>

// Cereal includes
#include <cereal/access.hpp>

// Adsim includes
#include "adsim/math/constants.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  struct CylinderBed {
    ExpressionParam diameter{0.1}; // [m]

    double GetArea() const {
      assert(diameter.Get() > 0.0);
      return ((constants::kPI / 4.0) * diameter.Get() * diameter.Get());
    }

    double GetHydraulicDiameter() const {
      return diameter.Get();
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      diameter.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("Diameter", diameter));
    }
  };

  struct SquareBed {
    ExpressionParam width{0.1}; // [m]

    double GetArea() const {
      assert(width.Get() > 0.0);
      return (width.Get() * width.Get());
    }

    double GetHydraulicDiameter() const {
      return width.Get(); // 4*A/P = 4*w^2/(4w) = w
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      width.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(cereal::make_nvp("Width", width));
    }
  };

  struct RectangularBed {
    ExpressionParam width{0.1};  // [m]
    ExpressionParam length{0.1}; // [m]

    double GetArea() const {
      assert(width.Get() > 0.0);
      assert(length.Get() > 0.0);
      return (width.Get() * length.Get());
    }

    double GetHydraulicDiameter() const {
      const double w = width.Get();
      const double l = length.Get();
      return 2.0 * w * l / (w + l); // 4*A/P
    }

    void EvaluateExpressions(const EvalContext& ctx) {
      width.Evaluate(ctx);
      length.Evaluate(ctx);
    }

  private:
    friend class cereal::access;
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Width", width),
        cereal::make_nvp("Length", length)
      );
    }
  };

  using BedShape = std::variant<CylinderBed, SquareBed, RectangularBed>;

  inline double GetArea(const BedShape& shape) {
    return std::visit([](const auto& s) { return s.GetArea(); }, shape);
  }

  inline double GetHydraulicDiameter(const BedShape& shape) {
    return std::visit([](const auto& s) { return s.GetHydraulicDiameter(); }, shape);
  }

  inline void EvaluateBedShapeExpressions(BedShape& shape, const EvalContext& ctx) {
    std::visit([&](auto& s) { s.EvaluateExpressions(ctx); }, shape);
  }

} // end adsim namespace
