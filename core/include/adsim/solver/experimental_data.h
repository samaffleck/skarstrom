/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes 
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

// adsim includes
#include "adsim/objects/reactor.h"


namespace adsim {

	// ---------------------------------------------------------------------------
	// What a dataset measures (its observable binding)
	//
	// A dataset stored in the deck can name the simulated quantity it measures,
	// so the engine can extract the simulated counterpart from a run's logs and
	// compare the two (SimulatedSeries and CompareWithSimulation in
	// parameterisation.h) without the user describing it again. One dataset
	// measures one observable.
	//
	//   * Bulk: a bed variable (`variable` is one of BulkVariableNames(), the
	//     names Simulation::Series takes: "temperature", "pressure", ...) at the
	//     top, the bottom or a position z [m] from the bottom of the column, over
	//     one cycle (`cycle`, 0-based; -1 = the last logged cycle). x is time from
	//     the start of that cycle.
	//   * Component: the same for a component variable (ComponentVariableNames():
	//     "mole_fraction", ...) of `component` (a deck component name).
	//   * Kpi: a per-cycle KPI (sensitivity::KPINames(): "purity", "recovery",
	//     ...) of the tracked `component`. x is the cycle number: the solver's
	//     0-based cycle index is x + cycle_offset (-1 for 1-based data).
	//   * None: not bound (the default; decks saved before bindings existed).
	//
	// Units. The points are stored as given, in the units the data is in:
	//   * x: for a time series, `x_unit` is any time unit the engine's unit
	//     parser reads ("s", "ms", "min", "h"...); `x_scale` is the matching
	//     multiplier to seconds (t[s] = x * x_scale). SetExperimentalDataBinding
	//     (deck_edit.h) keeps the two consistent. For a KPI, x is a cycle number
	//     ("cycle", scale 1; the default "s" is replaced by it).
	//   * y: `y_unit` is a unit expression (ParseUnitExpression, unit_parser.h)
	//     of the observable's dimension ("K", "degC", "%", "kPa", "ccm"...); empty
	//     means the observable's SI unit. The text is kept rather than a factor
	//     because some units have an offset (degC, degF): SI = y * factor +
	//     offset. The data stays exactly as measured, so changing the unit never
	//     rescales the points.
	//
	// Weighting (for fits): `uncertainty` is an optional per-point standard
	// deviation in the y unit (empty, or the same length as y, every value > 0);
	// `weight` scales the dataset's share of a multi-dataset objective (> 0).
	// ---------------------------------------------------------------------------

	enum class ObservableKind : int { None = 0, Bulk = 1, Component = 2, Kpi = 3 };
	enum class ObservableLocation : int { Top = 0, Bottom = 1, Position = 2 };
	enum class ObservableXAxis : int { Time = 0, Cycle = 1 };

	// The tags the deck, the deck view, the edit ops and Python use.
	inline const std::vector<std::string>& ObservableKindTags() {
		static const std::vector<std::string> tags = {"none", "bulk", "component", "kpi"};
		return tags;
	}
	inline const std::vector<std::string>& ObservableLocationTags() {
		static const std::vector<std::string> tags = {"top", "bottom", "position"};
		return tags;
	}

	struct ObservableBinding {
		ObservableKind kind = ObservableKind::None;
		std::string variable;              // series name (Bulk, Component) or KPI name (Kpi)
		std::string component;             // Component and Kpi
		ObservableLocation location = ObservableLocation::Top;  // Bulk and Component
		double z = 0.0;                    // [m] from the bottom, with Position
		int cycle = -1;                    // Bulk and Component: 0-based cycle, -1 = last logged
		int cycle_offset = 0;              // Kpi: solver cycle index = x + cycle_offset
		std::string x_unit = "s";          // time unit of x ("cycle" for a KPI)
		double x_scale = 1.0;              // [s] per x unit (1 for a KPI)
		std::string y_unit;                // unit of y; empty = the observable's SI unit
		std::vector<double> uncertainty;   // per-point standard deviation in the y unit (optional)
		double weight = 1.0;               // the dataset's weight in a multi-dataset objective

		bool Bound() const { return kind != ObservableKind::None; }
		ObservableXAxis XAxis() const {
			return kind == ObservableKind::Kpi ? ObservableXAxis::Cycle : ObservableXAxis::Time;
		}

	private:
		friend class cereal::access;
		template <class Archive>
		void save(Archive& ar, std::uint32_t const) const {
			const std::string kind_tag = ObservableKindTags()[static_cast<size_t>(kind)];
			const std::string location_tag = ObservableLocationTags()[static_cast<size_t>(location)];
			ar(
				cereal::make_nvp("Kind", kind_tag),
				cereal::make_nvp("Variable", variable),
				cereal::make_nvp("Component", component),
				cereal::make_nvp("Location", location_tag),
				cereal::make_nvp("Z", z),
				cereal::make_nvp("Cycle", cycle),
				cereal::make_nvp("Cycle_Offset", cycle_offset),
				cereal::make_nvp("X_Unit", x_unit),
				cereal::make_nvp("X_Scale", x_scale),
				cereal::make_nvp("Y_Unit", y_unit),
				cereal::make_nvp("Uncertainty", uncertainty),
				cereal::make_nvp("Weight", weight)
			);
		}
		template <class Archive>
		void load(Archive& ar, std::uint32_t const) {
			std::string kind_tag, location_tag;
			ar(
				cereal::make_nvp("Kind", kind_tag),
				cereal::make_nvp("Variable", variable),
				cereal::make_nvp("Component", component),
				cereal::make_nvp("Location", location_tag),
				cereal::make_nvp("Z", z),
				cereal::make_nvp("Cycle", cycle),
				cereal::make_nvp("Cycle_Offset", cycle_offset),
				cereal::make_nvp("X_Unit", x_unit),
				cereal::make_nvp("X_Scale", x_scale),
				cereal::make_nvp("Y_Unit", y_unit),
				cereal::make_nvp("Uncertainty", uncertainty),
				cereal::make_nvp("Weight", weight)
			);
			kind = ObservableKind::None;  // an unknown tag (a newer engine's kind) reads as unbound
			for (size_t i = 0; i < ObservableKindTags().size(); ++i)
				if (ObservableKindTags()[i] == kind_tag) kind = static_cast<ObservableKind>(i);
			location = ObservableLocation::Top;
			for (size_t i = 0; i < ObservableLocationTags().size(); ++i)
				if (ObservableLocationTags()[i] == location_tag) location = static_cast<ObservableLocation>(i);
		}
	};

	class XYData {
	public:
		std::string name;
		std::vector<double> x;
		std::vector<double> y;
		std::string x_label;
		std::string y_label;
		ObservableBinding observable;  // what the points measure (version 1; unbound before)

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Name", name),
				cereal::make_nvp("x_data", x),
				cereal::make_nvp("y_data", y),
				cereal::make_nvp("x_label", x_label),
				cereal::make_nvp("y_label", y_label)
			);
			if (version >= 1) ar(cereal::make_nvp("Observable", observable));
		}
	};

	// Checks a binding against a dataset and a deck's components: a known kind,
	// variable and location, a component the deck has (and none for a bulk
	// variable), z >= 0, a time x unit (Bulk, Component), a y unit of the
	// observable's dimension, an uncertainty of y's length with positive
	// values, a positive weight. On success it also sets x_scale from x_unit
	// (and x_unit "cycle", x_scale 1 for a KPI). An unbound binding is valid.
	// Returns false with a message in `error` otherwise.
	bool NormaliseObservableBinding(ObservableBinding& binding, const XYData& data,
	                                const std::vector<std::string>& components, std::string& error);

	// The SI unit of a bound observable ("K", "-", "mol/mol", "J/kg"...), or ""
	// when the binding is unbound or names no known variable.
	std::string ObservableSIUnit(const ObservableBinding& binding);

  class ExperimentalData {
  public:
		bool LoadData(const std::filesystem::path& path, bool headers_included);
		bool AddData(const XYData& data);
		bool DeleteDataByIndex(size_t idx);
		void DeleteAllData() { data_.clear(); }		
		std::vector<XYData>& GetAllData() { return data_; }
		const std::vector<XYData>& GetAllData() const { return data_; }
		size_t size() const { return data_.size(); }
		// The first dataset named `name`, or nullptr.
		XYData* Find(const std::string& name) {
			for (auto& d : data_)
				if (d.name == name) return &d;
			return nullptr;
		}

	private:
		std::vector<XYData> data_;
		
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Data", data_)
			);
		}

  };

} // end adsim namespace

// Class versions live next to the inline serialize templates: every
// translation unit that saves or loads a deck must see the same numbers.
CEREAL_CLASS_VERSION(adsim::ObservableBinding, 0);
CEREAL_CLASS_VERSION(adsim::XYData, 1);  // 1: Observable (what the dataset measures)
CEREAL_CLASS_VERSION(adsim::ExperimentalData, 1);
