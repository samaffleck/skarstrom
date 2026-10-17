/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/math/linear_algebra.h"


namespace adsim {

	void TDMA(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V) {
		// Forward pass
		auto n = int(V.size());
		std::vector<double> cPrime(n);
		std::vector<double> dPrime(n);

		cPrime[0] = Ae(0) / Ap(0);
		dPrime[0] = S[0] / Ap(0);
		for (int i = 1; i < n - 1; i++)
		{
			cPrime[i] = Ae(i) / (Ap(i) - Aw(i) * cPrime[i - 1]);
			dPrime[i] = (S[i] - Aw(i) * dPrime[i - 1]) / (Ap(i) - Aw(i) * cPrime[i - 1]);
		}
		dPrime[n - 1] = (S[n - 1] - Aw(n - 1) * dPrime[n - 2]) / (Ap(n - 1) - Aw(n - 1) * cPrime[n - 2]);

		// Back substitution
		V[n - 1] = dPrime[n - 1];
		for (int i = n - 2; i >= 0; i--)
		{
			V[i] = dPrime[i] - cPrime[i] * V[i + 1];
		}
	}

	void LU(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x) {
		x = A.lu().solve(b);
		//Eigen::SparseLU<Eigen::MatrixXd> linsolver; // Should this be static?
		//linsolver.compute(A);
		//V = linsolver.solve(S);
	}

	void TDMA(const Eigen::MatrixXd& A, const Eigen::VectorXd& S, Eigen::VectorXd& V) {
		int n = int(V.size());
		std::vector<double> cPrime(n);
		std::vector<double> dPrime(n);

		// Forward sweep
		cPrime[0] = A(0, 1) / A(0, 0);  // superdiagonal at (0,1)
		dPrime[0] = S[0] / A(0, 0);

		for (int i = 1; i < n - 1; ++i) {
			double denom = A(i, i) - A(i, i - 1) * cPrime[i - 1];  // subdiagonal at (i, i-1)
			cPrime[i] = A(i, i + 1) / denom;                       // superdiagonal at (i, i+1)
			dPrime[i] = (S[i] - A(i, i - 1) * dPrime[i - 1]) / denom;
		}

		dPrime[n - 1] = (S[n - 1] - A(n - 1, n - 2) * dPrime[n - 2]) /
			(A(n - 1, n - 1) - A(n - 1, n - 2) * cPrime[n - 2]);

		// Back substitution
		V[n - 1] = dPrime[n - 1];
		for (int i = n - 2; i >= 0; --i) {
			V[i] = dPrime[i] - cPrime[i] * V[i + 1];
		}
	}

	double getError(const Eigen::VectorXd& Ap, const Eigen::VectorXd& Ae, const Eigen::VectorXd& Aw, const Eigen::VectorXd& S, Eigen::VectorXd& V) {
		const auto N = (int)V.size();
		Eigen::VectorXd residual;
		residual.setConstant(N, 0);

		residual(0) = Ap(0) * V(0) + Ae(0) * V(1) - S(0);

		for (size_t n = 1; n < N - 1; ++n)
		{
			residual(n) = Aw(n) * V(n - 1) + Ap(n) * V(n) + Ae(n) * V(n + 1) - S(n);
		}

		residual(N - 1) = Aw(N - 1) * V(N - 2) + Ap(N - 1) * V(N - 1) - S(N - 1);

		return residual.norm();
	}

	void Newton(const Eigen::VectorXd& Jp, const Eigen::VectorXd& Je, const Eigen::VectorXd& Jw, const Eigen::VectorXd& R, Eigen::VectorXd& V, double relaxFactor) {
		Eigen::VectorXd dV = V; // Copy of vector of unknowns. This will store the change in variable.
		dV.setConstant(0.0);
		TDMA(Jp, Je, Jw, -R, dV);
		V += dV * relaxFactor;
	}

	void NewtonStep(const Eigen::MatrixXd& J, const Eigen::VectorXd& R, Eigen::VectorXd& x, double relaxFactor) {
		Eigen::VectorXd dx = x; // Copy of vector of unknowns. This will store the change in variable.
		TDMA(J, -R, dx); // We can use this because we are assuming J is tridiagonal
		x += dx * relaxFactor;
	}

	void NewtonStepLU(const Eigen::MatrixXd& J, const Eigen::VectorXd& b, Eigen::VectorXd& x, double relaxFactor) {
		Eigen::VectorXd dx = x; // Copy of vector of unknowns. This will store the change in variable.
		LU(J, -b, dx); 
		x += dx * relaxFactor;
	}

	void NewtonStepSparse(const Eigen::SparseMatrix<double>& J_sparse, const Eigen::VectorXd& b, Eigen::VectorXd& x, double relaxFactor) {
		// Solve J dx = -b
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(J_sparse);

		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("SparseLU factorization failed: " + solver.lastErrorMessage());
		}

		Eigen::VectorXd dx = solver.solve(-b);

		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("SparseLU solve failed");
		}

		x += dx * relaxFactor;
	}

	void ApproximateJacobianForward(const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx, Eigen::MatrixXd& J, const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
		// Adapted from fsolve (https://people.sc.fsu.edu/~jburkardt/cpp_src/fsolve/fsolve.cpp)
		// fx is the function which updates the system of N non-linear equations. 
		// y is the vector of the residual of this function update evaluated at x. 
		// x is the vector of N unknowns
		// J is the (N by N) jacobian matrix that gets approximated by this function
		
		auto N = (int)x.size();
		double h = 1e-6;
		auto delx = x; // Intentional copy
		auto dely = y; // Intentional copy
		auto x_org = x; // Intentional copy

		for (int k = 0; k < 3; ++k)
		{
			for (int j = k; j < N; j += 3)
			{
				delx(j) = x(j) + h;
			}
			fx(delx, dely);
			for (int j = k; j < N; j += 3)
			{
				//delx(j) = x(j); // reset the step change in value 
				delx(j) = x_org(j); // reset the step change in value 

				// Loop through all the equations
				for (int n = 0; n < N; ++n)
				{
					if (j - 1 <= n && n <= j + 1)
					{
						//J(j, n) = (dely(n) - y(n)) / h;
						J(n, j) = (dely(n) - y(n)) / h;
					}
				}
			}
		}
	}

	void ApproximateJacobianCentral(
		const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx,
		Eigen::MatrixXd& J,
		const Eigen::VectorXd& x,
		const Eigen::VectorXd& y)
	{
		int N = static_cast<int>(x.size());
		double h = 1e-6;

		Eigen::VectorXd x_plus = x;
		Eigen::VectorXd x_minus = x;
		Eigen::VectorXd y_plus(N);
		Eigen::VectorXd y_minus(N);

		// Only compute for tridiagonal structure: j-1 <= i <= j+1
		for (int k = 0; k < 3; ++k) // Loop over the 3 diagonals
		{
			for (int j = k; j < N; j += 3)
			{
				x_plus(j) = x(j) + h;
				x_minus(j) = x(j) - h;
			}

			fx(x_plus, y_plus);
			fx(x_minus, y_minus);

			for (int j = k; j < N; j += 3)
			{
				// Reset perturbed values
				x_plus(j) = x(j);
				x_minus(j) = x(j);

				for (int i = std::max(0, j - 1); i <= std::min(N - 1, j + 1); ++i)
				{
					J(i, j) = (y_plus(i) - y_minus(i)) / (2.0 * h);
				}
			}
		}
	}

	void ApproximateJacobianSparse(const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& fx, Eigen::SparseMatrix<double>& J, const Eigen::VectorXd& x, const Eigen::VectorXd& y, int mu, int ml) {
		// Adapted from fsolve (https://people.sc.fsu.edu/~jburkardt/cpp_src/fsolve/fsolve.cpp)
		// fx is the function which updates the system of N non-linear equations. 
		// y is the vector of the residual of this function update evaluated at x. 
		// x is the vector of N unknowns
		// J is the (N by N) jacobian matrix that gets approximated by this 
		// function. J is sparse, so must be filled by triplets
		// ml is the number of lower non-zero diagonals
		// mu is the number of upper non-zero diagonals
		// msum is the total number of non-zero diagonals, including the main diagaonal

		int msum = 1 + mu + ml;
		auto N = (int)x.size();
		std::vector<Eigen::Triplet<double>> triplets; // Order in the list does not matter
		triplets.reserve(msum * N);

		double h = 1e-6;
		auto delx = x; // Intentional copy
		auto dely = y; // Intentional copy

		for (int k = 0; k < msum; ++k) 
		{
			for (int j = k; j < N; j += msum)
			{
				delx(j) = x(j) + h;
			}
			fx(delx, dely);
			for (int j = k; j < N; j +=  msum)
			{
				delx(j) = x(j); // reset the step change in value 

				// Loop through all the equations
				for (int n = 0; n < N; ++n)
				{
					if (j - mu <= n && n <= j + ml)
					{
						double Jnj = (dely(n) - y(n)) / h;
						triplets.emplace_back(n, j, Jnj);
					}
				}
			}
		}

		J.setZero();         // Clear old values. TODO: IS this needed?
		J.setFromTriplets(triplets.begin(), triplets.end());
	}

} // end adsim namespace
