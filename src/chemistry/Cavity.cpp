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

#include <array>
#include <cmath>
#include <vector>

#include "Cavity.h"
#include "utils/math_utils.h"

namespace mrchem {

Cavity::Cavity(std::vector<mrcpp::Coord<3>> &coord, std::vector<double> &R, double slope) {
    this->pos = coord;
    this->R = R;
    this->d = slope;
}

double Cavity::evalf(const mrcpp::Coord<3> &r) const {
    double C = 1.0;
    double s, O;

    for (int i = 0; i < pos.size(); i++) {
        s = math_utils::calc_distance(pos[i], r) - R[i];
        O = 0.5 * (1 + std::erf(s / d));
        C *= 1 - (1 - O);
    }
    C = 1 - C;
    return C;
}

} // namespace mrchem
