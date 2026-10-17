/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <vector>
#include <string>
#include <unordered_map>


namespace adsim {

	// Defined once in outlet_conditions.cpp (not `inline`: an inline function
	// would need a definition in every translation unit that uses it).
	double InterpolateTimeSeries(const std::vector<double>& times, const std::vector<double>& values, double time);

	class OutletConditions {
	public:
		void Clear();
		void LogOutletConditions(double u, double T, double P, double t, const std::unordered_map<std::string, double>& yi);
		double GetP(double time) const;
		double GetT(double time) const;
		double GetU(double time) const;
		double GetYi(double time, const std::string& component) const;

	// private:
		std::unordered_map<std::string, std::vector<double>> yi_{};
		std::vector<double> u_{};
		std::vector<double> T_{};
		std::vector<double> P_{};
		std::vector<double> t_{};
	};

} // end adsim namespace
