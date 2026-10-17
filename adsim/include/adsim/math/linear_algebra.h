/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <functional>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace adsim {

	void TDMA(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V);
	void TDMA(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);
	double getError(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V);
	void Newton(const Eigen::VectorXd& Jp, const Eigen::VectorXd& Je, const Eigen::VectorXd& Jw, const Eigen::VectorXd& R, Eigen::VectorXd& V, double relaxFactor = 1.0);
	void NewtonStep(const Eigen::MatrixXd& J, const Eigen::VectorXd& b, Eigen::VectorXd& x, double relaxFactor = 1.0);
	void NewtonStepLU(const Eigen::MatrixXd& J, const Eigen::VectorXd& b, Eigen::VectorXd& x, double relaxFactor = 1.0);
	void NewtonStepSparse(const Eigen::SparseMatrix<double>& J_sparse, const Eigen::VectorXd& b, Eigen::VectorXd& x, double relaxFactor = 1.0);
	void ApproximateJacobianForward(const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx, Eigen::MatrixXd& J, const Eigen::VectorXd& x, const Eigen::VectorXd& y);
	void ApproximateJacobianCentral(const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx, Eigen::MatrixXd& J, const Eigen::VectorXd& x, const Eigen::VectorXd& y);
	void ApproximateJacobianSparse(const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx, Eigen::SparseMatrix<double>& J, const Eigen::VectorXd& x, const Eigen::VectorXd& y, int mu, int ml);

} // end adsim namespace
