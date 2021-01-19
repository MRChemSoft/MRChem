/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "MRCPP/Gaussians"
#include "MRCPP/Plotter"
#include "MRCPP/Printer"

#include "Nucleus.h"
#include "chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "qmoperators/one_electron/H_E_dip.h"

namespace mrchem {

/** @brief computes the repulsion self energy of a set of nuclei
 *
 * @param[in] nucs the set of nuclei
 *
 */
double chemistry::compute_nuclear_repulsion(const Nuclei &nucs) {
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        const double Z_i = nuc_i.getCharge();
        const mrcpp::Coord<3> &R_i = nuc_i.getCoord();
        for (int j = i + 1; j < nNucs; j++) {
            const Nucleus &nuc_j = nucs[j];
            const double Z_j = nuc_j.getCharge();
            const mrcpp::Coord<3> &R_j = nuc_j.getCoord();
            double R_ij = math_utils::calc_distance(R_i, R_j);
            E_nuc += (Z_i * Z_j) / R_ij;
        }
    }
    return E_nuc;
}

/** @brief Returns the sum of atomic charges*/
double chemistry::get_total_charge(const Nuclei &nucs) {
    double charge = 0;
    for (const auto &nuc : nucs) charge += nuc.getCharge();
    return charge;
}

/** @breif computes the nuclear density as a sum of narrow Gaussians */
Density chemistry::compute_nuclear_density(double prec, const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / MATHCONST::pi, 3.0 / 2.0);
    auto gauss = mrcpp::GaussExp<3>();
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta * Z_i, R_i);
        gauss.append(gauss_f);
    }

    if ((*MRA).getWorldBox().isPeriodic()) {
        auto period = (*MRA).getWorldBox().getScalingFactor();
        for (auto &p : period) p *= 2.0;
        gauss.makePeriodic(period);
    }
    Density rho(false);
    density::compute(prec, rho, gauss);
    return rho;
}
Density chemistry::compute_nuclear_density_smeared(double prec, Nuclei nucs, double rc, double period) {
    Density rho(false);
    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    nucs = periodic::periodify_nuclei(nucs, period);

    auto b_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto g_i = 0.0;
            if (R <= rc and R >= 0) {
                g_i = -21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
                      (5.0 * mrcpp::pi * std::pow(rc, 8.0));
            }
            g_rc += g_i * nuc.getCharge();
        }
        return g_rc;
    };

    auto beta = 1.0e4;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    mrcpp::GaussExp<3> f_exp;
    for (auto &nuc : nucs) {
        auto g_func = mrcpp::GaussFunc<3>(beta, alpha, nuc.getCoord());
        f_exp.append(g_func);
    }

    mrcpp::build_grid(rho.real(), f_exp);
    mrcpp::project<3>(prec, rho.real(), b_smear);
    return rho;
}
/** @breif computes the nuclear density as a sum of narrow Gaussians */
double chemistry::compute_nuclear_self_repulsion(const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / MATHCONST::pi, 3.0 / 2.0);

    auto self_rep = 0.0;
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta, R_i);
        self_rep += Z_i * Z_i * gauss_f.calcCoulombEnergy(gauss_f);
    }
    return self_rep;
}


Density chemistry::hack_density(double prec, Nuclei nucs, double rc, double period, OrbitalVector Phi) {
    Density rho(false);
    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
    dip_oper.setup(prec);
    DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
    DoubleVector dip_el = dip_oper.trace(Phi).real();
    DoubleVector tot_dip = nuc_dip + dip_el;

    auto new_charge = tot_dip[2]/8.0;

    mrcpp::Coord<3> pos_coord{0.0, 0.0, 4.0};
    mrcpp::Coord<3> neg_coord{0.0, 0.0, -4.0};

    println(0, "tot_dip " << tot_dip[0] << " " << tot_dip[1] << " " << tot_dip[2])
    dip_oper.clear();

    std::vector<double> charges{-new_charge, new_charge};
    std::vector<mrcpp::Coord<3>> coords{pos_coord, neg_coord};

    auto b_smear = [charges, coords, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto i = 0; i < charges.size(); i++) {
            auto R = math_utils::calc_distance(r, coords[i]);
            auto g_i = 0.0;
            if (R <= rc and R >= 0) {
                g_i = -21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
                      (5.0 * mrcpp::pi * std::pow(rc, 8.0));
            }
            g_rc += g_i * charges[i];
        }
        return g_rc;
    };

    Density rho_tree(false);
    if (not rho_tree.hasReal()) rho_tree.alloc(NUMBER::Real);

    mrcpp::project<3>(prec, rho_tree.real(), b_smear);

    Density tmp_rho = chemistry::compute_nuclear_density_smeared(prec, nucs, rc, period);


    mrcpp::add(prec, rho.real(), 1.0, rho_tree.real(), 1.0, tmp_rho.real());

    PositionOperator r({0.0, 0.0, 0.0});
    r.setup(prec);

    DoubleVector rho_tree_mu = r.trace(rho_tree).real();
    println(0, "rho_tree_mu " << rho_tree_mu[0] << " " << rho_tree_mu[1] << " " << rho_tree_mu[2])

    DoubleVector tmp_rho_mu = r.trace(tmp_rho).real();
    println(0, "tmp_rho_mu " << tmp_rho_mu[0] << " " << tmp_rho_mu[1] << " " << tmp_rho_mu[2])

    DoubleVector rho_mu = r.trace(rho).real();
    println(0, "rho_mu " << rho_mu[0] << " " << rho_mu[1] << " " << rho_mu[2])

    r.clear();

    return rho;
}
} // namespace mrchem
