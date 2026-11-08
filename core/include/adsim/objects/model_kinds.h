/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// One definition per model kind.
//
// Every selectable model (an isotherm, a rate coefficient, a dispersion or
// heat-transfer correlation, a boundary profile, a pump model, ...) belongs to
// a family and is described once here: its stable tag (the snake_case name the
// builder, the registry's choice lists, the Python package and front ends use),
// a display label, and its ordered parameters with name, unit, default, search
// bounds and whether the optimiser may vary them. The descriptors are listed in
// the family's enum order, so a kind's index is its stored enum value (which
// never changes).
//
// What is built from them:
//   * the builder's kind lists (deck_builder.h: RateCoefficientKinds(), ...),
//   * the parameter registry's entries for each kind (units, bounds, tunable),
//   * the isotherm fitter's kind list and the IAST family table (isotherms),
//   * the parameter carry-over of the kind switches (deck_edit.h): a parameter
//     keeps its value across a switch when the old and new kinds both have it
//     under the same name and unit.
// ModelFields() binds a live model's parameters to these names.

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "adsim/equations/isotherm.h"
#include "adsim/math/expression_param.h"

namespace adsim {

class IMassTransferModel;
class IMassTransferCoefficientModel;
class IDispersionModel;
class IHeatDispersionModel;
class IInternalHeatTransferModel;
class IExternalHeatTransferModel;
class IParticleHeatTransferModel;
class IBoundaryValue;
class IPump;
class IPumpEfficiencyModel;
class IPumpPowerModel;
class ValveBoundaryStream;

// The unit text of a profile parameter that takes the unit of the profiled
// quantity (a pressure in Pa, a temperature in K, a flow in its kind's unit).
inline constexpr const char* kProfileQuantityUnit = "@";

struct KindParameter {
  std::string name;   // registry leaf / field name, e.g. "tortuosity"
  std::string unit;   // SI unit text as the registry shows it; "-" dimensionless
  double default_value = 0.0;
  double lower = 0.0;  // suggested search bounds
  double upper = 0.0;
  bool tunable = true;
  std::string label;  // display name; empty: the name
  const std::string& Label() const { return label.empty() ? name : label; }
};

struct KindDescriptor {
  std::string tag;
  std::string label;
  std::vector<KindParameter> parameters;
  const KindParameter* Find(const std::string& name) const;
};

enum class ModelFamily {
  RateModel,             // MassTransferModelType: ldf, ldf2, dual_rate_ldf
  RateCoefficient,       // MassTransferCoefficientType
  Momentum,              // MomentumModel
  Dispersion,            // AxialDispersionType
  HeatDispersion,        // AxialHeatDispersionType
  InternalHeatTransfer,  // InternalHeatTransferType
  ExternalHeatTransfer,  // ExternalHeatTransferType
  ParticleHeatTransfer,  // ParticleHeatTransferType
  Boundary,              // BoundaryType: pressure, flow, closed
  Stream,                // BoundaryStreamType: flow, valve, pump
  Profile,               // BoundaryValueType
  Flow,                  // FlowType
  Pump,                  // "none" then PumpType + 1
  PumpEnergy,            // CompressorEnergyType
  PumpEfficiency,        // PumpEfficiencyType
  PumpPower,             // PumpPowerType
};

// Every family, in enum order.
const std::vector<ModelFamily>& AllModelFamilies();
// snake_case family name ("rate_coefficient", ...).
const char* ModelFamilyName(ModelFamily family);

// The kinds of a family in enum order (index = stored enum value; for Pump,
// index 0 is "none" and index k is PumpType k - 1).
const std::vector<KindDescriptor>& KindDescriptors(ModelFamily family);
// The tags of a family, in enum order.
const std::vector<std::string>& KindTags(ModelFamily family);
// The descriptor with this tag, or nullptr.
const KindDescriptor* FindKind(ModelFamily family, const std::string& tag);
// The descriptor at enum value `index`, or nullptr when out of range.
const KindDescriptor* KindAt(ModelFamily family, std::size_t index);

// The rate constant of an LDF2 model is in kg/(mol s), not 1/s, and its
// suggested bounds are a decade wider on each side: |q* - q| is typically 0.1
// to 10 mol/kg. The rate-coefficient parameters this applies to are those whose
// descriptor unit is "1/s".
inline constexpr const char* kLdf2RateConstantUnit = "kg/(mol-s)";
KindParameter Ldf2RateConstant(const KindParameter& p);

// The slow site of a dual-rate LDF model has its own coefficient model, of the
// same kind as the fast site's. Its parameters are those of the
// rate-coefficient descriptor under the name "slow_<name>" (registry leaf
// adsorbent[A].kinetics[C].slow_<name>), labelled "slow-site <label>".
inline constexpr const char* kSlowSitePrefix = "slow_";
KindParameter SlowSiteParameter(const KindParameter& p);

// ---- isotherms -------------------------------------------------------------------

// Pure-component form IAST integrates (iast.cpp): a kind maps onto one of
// these, or has none (the water-coupled kinds).
enum class IastFamily { None, Inert, Henry, Langmuir, DualSiteLangmuir, Freundlich, Sips, DualSiteSips,
                        AntiLangmuir, RedlichPeterson, Toth, Gab };

struct IsothermKindDescriptor : KindDescriptor {
  IsothermType type = IsothermType::INERT;
  IastFamily iast = IastFamily::None;
  bool competitive = false;     // IsCompetitiveIsothermKind
  bool water_coupled = false;   // IsWaterCoupledIsothermKind
  bool fittable = false;        // offered by the isotherm fitter (isotherm_fit.h)
  int fit_order = 0;            // position in the fitter's kind list (FittableKinds), from 1
  bool temperature_dependent = false;  // the fitter can identify a T dependence
  bool has_pair = false;        // PairedIsothermKind
  IsothermType pair = IsothermType::INERT;
};

// One per IsothermType, in enum order (index = stored enum value).
const std::vector<IsothermKindDescriptor>& IsothermKindDescriptors();
const IsothermKindDescriptor& IsothermKind(IsothermType type);

// ---- binding a live model's parameters to the descriptor names ----------------------

using ModelFieldList = std::vector<std::pair<std::string, ExpressionParam*>>;

// The model's parameters by descriptor name, in descriptor order (empty for a
// null model or a kind with no parameters). Profile names are the plain
// descriptor names ("start", "end", "time", "value", ...).
ModelFieldList ModelFields(IMassTransferModel* model);  // the rate model's own (dual_rate_ldf: split)
ModelFieldList ModelFields(IMassTransferCoefficientModel* model);
ModelFieldList ModelFields(IDispersionModel* model);
ModelFieldList ModelFields(IHeatDispersionModel* model);
ModelFieldList ModelFields(IInternalHeatTransferModel* model);
ModelFieldList ModelFields(IExternalHeatTransferModel* model);
ModelFieldList ModelFields(IParticleHeatTransferModel* model);
ModelFieldList ModelFields(IBoundaryValue* model);
ModelFieldList ModelFields(IPumpEfficiencyModel* model);
ModelFieldList ModelFields(IPumpPowerModel* model);
// A pump's own fields by kind (suction or discharge pressure); the polytropic
// index is common to every pump and not part of a kind.
ModelFieldList ModelFields(IPump* pump);

// The slow site's descriptor (the coefficient kind's, every parameter through
// SlowSiteParameter) and its live fields under those names; nullptr / empty
// when the rate model has no slow-site coefficient model.
const KindDescriptor* SlowSiteKind(IMassTransferModel* model);
ModelFieldList SlowSiteFields(IMassTransferModel* model);

// Copy every parameter `to` shares with `from` (same name and same unit in
// their descriptors), value and expression text. Returns the names copied.
std::vector<std::string> CarryOverParameters(const KindDescriptor& from_kind, const ModelFieldList& from,
                                             const KindDescriptor& to_kind, const ModelFieldList& to);

}  // namespace adsim
