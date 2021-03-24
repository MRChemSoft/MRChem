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

#pragma once

#include "NuclearGradientOperator.h"
#include "PositionOperator.h"
#include "qmoperators/RankTwoTensorOperator.h"

namespace mrchem {

/** @class H_MB_dia
 *
 * @brief Diamagnetic NMR operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt the
 * nuclear magnetic moment M_K of nucleus K and the external magnetic field B:
 *
 * d^2H/dM_KdB = H_MB_dia
 *
 * H_MB_dia = \frac{alpha^2}{2}
 *            \sum_j \frac{(r_{jO} \cdot r_{jK})1 - r_{jO}r_{jK}^T}{r_jK}^3}
 *
 * This operator is the transpose of H_BM_dia.
 */

class H_MB_dia final : public RankTwoTensorOperator<3, 3> {
public:
    H_MB_dia(const mrcpp::Coord<3> &o, const mrcpp::Coord<3> &k, double c) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;

        PositionOperator r_o(o);
        NuclearGradientOperator r_rm3(1.0, k, c);

        RankZeroTensorOperator &o_x = r_o[0];
        RankZeroTensorOperator &o_y = r_o[1];
        RankZeroTensorOperator &o_z = r_o[2];
        RankZeroTensorOperator &k_x = r_rm3[0];
        RankZeroTensorOperator &k_y = r_rm3[1];
        RankZeroTensorOperator &k_z = r_rm3[2];

        // Invoke operator= to assign *this operator
        RankTwoTensorOperator<3, 3> &h = (*this);
        h[0][0] = -(alpha_2 / 2.0) * (o_y * k_y + o_z * k_z);
        h[0][1] = (alpha_2 / 2.0) * (o_y * k_x);
        h[0][2] = (alpha_2 / 2.0) * (o_z * k_x);
        h[1][0] = (alpha_2 / 2.0) * (o_x * k_y);
        h[1][1] = -(alpha_2 / 2.0) * (o_x * k_x + o_z * k_z);
        h[1][2] = (alpha_2 / 2.0) * (o_z * k_y);
        h[2][0] = (alpha_2 / 2.0) * (o_x * k_z);
        h[2][1] = (alpha_2 / 2.0) * (o_y * k_z);
        h[2][2] = -(alpha_2 / 2.0) * (o_x * k_x + o_y * k_y);
        h[0][0].name() = "h_MB_dia[x,x]";
        h[0][1].name() = "h_MB_dia[x,y]";
        h[0][2].name() = "h_MB_dia[x,z]";
        h[1][0].name() = "h_MB_dia[y,x]";
        h[1][1].name() = "h_MB_dia[y,y]";
        h[1][2].name() = "h_MB_dia[y,z]";
        h[2][0].name() = "h_MB_dia[z,x]";
        h[2][1].name() = "h_MB_dia[z,y]";
        h[2][2].name() = "h_MB_dia[z,z]";
    }
};

} // namespace mrchem
