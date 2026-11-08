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
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/math/expression_param.h"

namespace adsim {

class IMassTransferModel;
class IMassTransferCoefficientModel;
class MassTransferCoefficientSeriesResistance;
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

// The particle shape a diffusion resistance is written for: packed spheres
// (packed layers), a slab of half-thickness l (structured layers), a cylinder
// (no fragments yet).
enum class DiffusionShape { Sphere, Slab, Cylinder };

// One resistance of the generic mass-transfer coefficient model
// (series_resistance; MassTransferCoefficientSeriesResistance): its switch tag,
// label, the term it refines, whether the capacity factor K multiplies it, its
// parameters (names in the kind's `parameters`, display order) and the LaTeX
// fragments of its term per particle shape. Placeholders in the fragments:
// {K} becomes "K\," when the capacity-factor mode is not none and the term is
// scaled by capacity (else nothing); {D} becomes D_{m,i}, or D_{B,i} with
// Knudsen diffusion on. MassTransferEquationLatex composes them.
struct ResistanceDescriptor {
  std::string tag;                      // "film", "macropore", "knudsen", "micropore", "lumped"
  std::string label;                    // "external film", ...
  std::string parent;                   // "" or the term it refines ("knudsen" -> "macropore")
  bool scaled_by_capacity = false;      // multiplied by K (film, macropore)
  std::vector<std::string> parameters;  // names in the kind's `parameters`, in display order
  struct Latex {
    std::string term;   // the term of 1/k (empty: refines its parent's term)
    std::string where;  // the definitions it uses (may be empty)
    // The definitions when the term's activation energy is the number 0
    // (lumped: k_L = k; micropore: D_c = D_c0); empty: `where` always.
    std::string where_isothermal;
  };
  Latex sphere, slab, cylinder;  // per DiffusionShape (cylinder empty for now)
  std::string source;            // literature, for docs and tooltips
  const Latex& Fragments(DiffusionShape shape) const {
    return shape == DiffusionShape::Slab ? slab : shape == DiffusionShape::Cylinder ? cylinder : sphere;
  }
};

// A named switch set of the generic model (the former coefficient kinds): the
// resistances on, the capacity-factor mode (CapacityFactorMode value) and the
// parameter values the preset fixes.
struct PresetDescriptor {
  std::string tag;
  std::vector<std::string> on;
  int capacity_mode = 0;
  std::vector<std::pair<std::string, double>> values;
};

// One capacity-factor mode (CapacityFactorMode order): tag, label and the LaTeX
// definition of K (empty for none).
struct CapacityModeDescriptor {
  std::string tag, label, latex;
};

struct KindDescriptor {
  std::string tag;
  std::string label;
  std::vector<KindParameter> parameters;
  // A kind kept for loading old decks and as a preset tag only: not offered
  // in the builder's kind lists (the six former rate-coefficient kinds).
  bool legacy = false;
  std::vector<ResistanceDescriptor> resistances;       // series_resistance only, in series order
  std::vector<CapacityModeDescriptor> capacity_modes;  // series_resistance only, CapacityFactorMode order
  std::vector<PresetDescriptor> presets;               // series_resistance only
  const KindParameter* Find(const std::string& name) const;
  const ResistanceDescriptor* FindResistance(const std::string& tag) const;
  const PresetDescriptor* FindPreset(const std::string& tag) const;
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
// The generic coefficient model (MassTransferCoefficientSeriesResistance)
// binds every parameter under the series_resistance names (ResistanceFields);
// an older coefficient class (load-only) binds its own kind's.
ModelFieldList ModelFields(IMassTransferCoefficientModel* model);
ModelFieldList ResistanceFields(MassTransferCoefficientSeriesResistance* model);
// The preset a coefficient model matches (MassTransferCoefficientSeriesResistance::
// PresetType; an older class: its kind), and the parameters a preset shows
// bound under that preset kind's names (arrhenius: pre_exp_factor is the
// lumped k; SERIES_RESISTANCE: every parameter). For the parameter carry-over
// of a preset switch (deck_edit.h). An older class gives ModelFields.
MassTransferCoefficientType CoefficientPreset(IMassTransferCoefficientModel* model);
ModelFieldList PresetFields(IMassTransferCoefficientModel* model, MassTransferCoefficientType preset);
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

// The descriptor of the generic coefficient model (series_resistance), which
// carries the resistances, capacity modes and presets.
const KindDescriptor& ResistanceKind();

// The display equation of a coefficient model for one particle shape, built
// only from the descriptor's fragments: "\frac{1}{k_i} = " and the enabled
// terms joined by " + " (series order); then, when any definitions apply,
// " \\ " and the definitions of the enabled terms and of the capacity mode
// joined by ",\qquad ", the whole inside \begin{aligned} ... \end{aligned}.
// Empty for a shape without fragments (cylinder).
std::string MassTransferEquationLatex(const MassTransferCoefficientSeriesResistance& model, DiffusionShape shape);

// Copy every parameter `to` shares with `from` (same name and same unit in
// their descriptors), value and expression text. Returns the names copied.
std::vector<std::string> CarryOverParameters(const KindDescriptor& from_kind, const ModelFieldList& from,
                                             const KindDescriptor& to_kind, const ModelFieldList& to);

}  // namespace adsim
