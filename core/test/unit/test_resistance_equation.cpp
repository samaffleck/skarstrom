/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The display equation of the generic mass-transfer coefficient model
// (MassTransferEquationLatex): every valid switch combination times every
// capacity-factor mode times {sphere, slab} against a golden list written from
// the fragment table (the composer may only use descriptor strings), plus the
// activation-energy forms, the cylinder (no fragments) and the descriptor shape.

#include <gtest/gtest.h>

#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "adsim/equations/mass_transfer_model.h"
#include "adsim/objects/model_kinds.h"

using namespace adsim;

namespace {

struct Golden {
  const char* on;     // comma-separated switches
  const char* mode;   // capacity-factor mode tag
  const char* shape;  // "sphere" or "slab"
  const char* latex;
};

// 23 switch sets (15 of film/macropore/micropore/lumped, 8 more with Knudsen)
// x 4 modes x 2 shapes; activation energies 0.
const std::vector<Golden>& GoldenList() {
  static const std::vector<Golden> g{
    {"film", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} \\ k_f = \frac{Sh\,D_{m,i}}{d_p}\end{aligned})"},
    {"film", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} \\ k_f = \frac{Sh\,D_{m,i}}{D_h}\end{aligned})"},
    {"film", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = K_{fix}\end{aligned})"},
    {"film", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = K_{fix}\end{aligned})"},
    {"macropore", "none", "sphere", R"(\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}})"},
    {"macropore", "none", "slab", R"(\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}})"},
    {"macropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ K = K_{fix}\end{aligned})"},
    {"macropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ K = K_{fix}\end{aligned})"},
    {"macropore,knudsen", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}}\end{aligned})"},
    {"macropore,knudsen", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}}\end{aligned})"},
    {"macropore,knudsen", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = K_{fix}\end{aligned})"},
    {"micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"film,macropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p}\end{aligned})"},
    {"film,macropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h}\end{aligned})"},
    {"film,macropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}}\end{aligned})"},
    {"film,macropore,knudsen", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}}\end{aligned})"},
    {"film,macropore,knudsen", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad K = K_{fix}\end{aligned})"},
    {"film,micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k\end{aligned})"},
    {"film,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k\end{aligned})"},
    {"film,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"macropore,micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\end{aligned})"},
    {"macropore,micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"macropore,micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"macropore,knudsen,micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"macropore,knudsen,micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"macropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"macropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k\end{aligned})"},
    {"macropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,knudsen,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,knudsen,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,macropore,micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,macropore,micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,micropore", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,macropore,knudsen,micropore", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0}\end{aligned})"},
    {"film,macropore,knudsen,micropore", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,micropore", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,knudsen,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,knudsen,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"macropore,knudsen,micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{m,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "none", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "none", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{\ell}{k_f} + \frac{\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "slope", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "slope", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\rho_p\,\frac{\partial q_i^{*}}{\partial c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "secant", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "secant", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "fixed", "sphere", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,r_p}{3\,k_f} + \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{d_p},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
    {"film,macropore,knudsen,micropore,lumped", "fixed", "slab", R"(\begin{aligned}\frac{1}{k_i} = f_{film}\,\frac{K\,\ell}{k_f} + \frac{K\,\tau\,\ell^{2}}{3\,\varepsilon_p\,D_{B,i}} + \frac{r_c^{2}}{15\,D_{c,i}} + \frac{1}{k_{L,i}} \\ k_f = \frac{Sh\,D_{m,i}}{D_h},\qquad \frac{1}{D_{B,i}} = \frac{1}{D_{K,i}} + \frac{1}{D_{m,i}},\quad D_{K,i} = \tfrac{2}{3}\,r_{pore}\sqrt{\tfrac{8RT}{\pi M_i}},\qquad D_{c,i} = D_{c0},\qquad k_{L,i} = k,\qquad K = K_{fix}\end{aligned})"},
  };
  return g;
}

MassTransferCoefficientSeriesResistance Model(const std::string& on, const std::string& mode) {
  MassTransferCoefficientSeriesResistance m;
  MassTransferResistances r;
  r.lumped = false;
  std::stringstream ss(on);
  for (std::string t; std::getline(ss, t, ',');) {
    if (t == "film") r.film = true;
    else if (t == "macropore") r.macropore = true;
    else if (t == "knudsen") r.knudsen = true;
    else if (t == "micropore") r.micropore = true;
    else if (t == "lumped") r.lumped = true;
    else ADD_FAILURE() << "unknown switch " << t;
  }
  m.resistances_ = r;
  int i = 0;
  for (const char* tag : CapacityFactorModeList) {
    if (mode == tag) m.capacity_mode_ = i;
    ++i;
  }
  return m;
}

}  // namespace

TEST(ResistanceEquation, Compose) {
  ASSERT_EQ(GoldenList().size(), 184u);
  for (const Golden& g : GoldenList()) {
    const auto m = Model(g.on, g.mode);
    const DiffusionShape shape = std::string(g.shape) == "slab" ? DiffusionShape::Slab : DiffusionShape::Sphere;
    EXPECT_EQ(MassTransferEquationLatex(m, shape), g.latex) << g.on << " / " << g.mode << " / " << g.shape;
  }
}

TEST(ResistanceEquation, ActivationEnergiesShowTheArrheniusForms) {
  auto m = Model("lumped", "none");
  m.activation_energy_ = 2e4;
  EXPECT_EQ(MassTransferEquationLatex(m, DiffusionShape::Sphere),
            R"(\begin{aligned}\frac{1}{k_i} = \frac{1}{k_{L,i}} \\ k_{L,i} = k_0\,e^{-E_a/RT}\end{aligned})");
  // An expression is not "the number 0": the general form.
  m.activation_energy_.SetExpr("0 * 1");
  EXPECT_NE(MassTransferEquationLatex(m, DiffusionShape::Sphere).find("e^{-E_a/RT}"), std::string::npos);

  auto c = Model("micropore", "none");
  c.micropore_activation_energy_ = 1e4;
  EXPECT_EQ(MassTransferEquationLatex(c, DiffusionShape::Slab),
            R"(\begin{aligned}\frac{1}{k_i} = \frac{r_c^{2}}{15\,D_{c,i}} \\ D_{c,i} = D_{c0}\,e^{-E_{a,c}/RT}\end{aligned})");
}

TEST(ResistanceEquation, CylinderHasNoFragmentsYet) {
  EXPECT_EQ(MassTransferEquationLatex(Model("film,macropore", "slope"), DiffusionShape::Cylinder), "");
}

TEST(ResistanceEquation, DescriptorShape) {
  const KindDescriptor& k = ResistanceKind();
  EXPECT_EQ(k.tag, "series_resistance");
  EXPECT_FALSE(k.legacy);
  ASSERT_EQ(k.resistances.size(), ResistanceTags().size());
  for (std::size_t i = 0; i < k.resistances.size(); ++i) {
    const ResistanceDescriptor& r = k.resistances[i];
    EXPECT_EQ(r.tag, ResistanceTags()[i]);
    EXPECT_FALSE(r.label.empty());
    EXPECT_FALSE(r.source.empty());
    EXPECT_EQ(r.parameters, MassTransferCoefficientSeriesResistance::ResistanceParameters(r.tag)) << r.tag;
    for (const auto& p : r.parameters) EXPECT_NE(k.Find(p), nullptr) << p;
    EXPECT_TRUE(r.cylinder.term.empty() && r.cylinder.where.empty()) << r.tag;
    if (r.parent.empty()) EXPECT_FALSE(r.sphere.term.empty() || r.slab.term.empty()) << r.tag;
  }
  EXPECT_EQ(k.FindResistance("knudsen")->parent, "macropore");
  EXPECT_TRUE(k.FindResistance("film")->scaled_by_capacity);
  EXPECT_TRUE(k.FindResistance("macropore")->scaled_by_capacity);
  EXPECT_FALSE(k.FindResistance("micropore")->scaled_by_capacity);
  ASSERT_EQ(k.capacity_modes.size(), 4u);
  for (std::size_t i = 0; i < 4; ++i) EXPECT_EQ(k.capacity_modes[i].tag, CapacityFactorModeList[i]);

  // Every former kind is legacy and a preset; the preset reproduces its
  // switches and mode.
  std::set<std::string> presets;
  for (const auto& p : k.presets) presets.insert(p.tag);
  for (std::size_t i = 0; i < static_cast<std::size_t>(MassTransferCoefficientType::SERIES_RESISTANCE); ++i) {
    const KindDescriptor* old = KindAt(ModelFamily::RateCoefficient, i);
    ASSERT_NE(old, nullptr);
    EXPECT_TRUE(old->legacy) << old->tag;
    EXPECT_EQ(presets.count(old->tag), 1u) << old->tag;
  }
  for (const auto& p : k.presets) {
    MassTransferCoefficientSeriesResistance m;
    m.resistances_ = MassTransferResistances{false, false, false, false, false};
    for (const auto& t : p.on) m.Switch(t) = true;
    m.capacity_mode_ = p.capacity_mode;
    for (const auto& [name, v] : p.values) {
      for (auto& [field, param] : ResistanceFields(&m))
        if (field == name) *param = v;
    }
    if (p.tag == "arrhenius") m.activation_energy_ = 1e4;
    EXPECT_EQ(m.PresetTag(), p.tag);
  }
}
