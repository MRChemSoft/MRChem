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

#include "CUBEfunction.h"
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace mrchem {

CUBEfunction::CUBEfunction() {
    normalize_basis();
}

// Do a quadrature of the file
double CUBEfunction::evalf(const mrcpp::Coord<3> &r) const {

    /** I'm assuming that we can describe any arbitrary point as r = (c_i + d_i)X + (c_j + dc_j)Y + (c_k + dc_k)Z
     * where X, Y, and Z are the three voxel basis vectors descibing the step size and direction of the three main vertexes of the CUBE.
     * Here I am assuming two things, first that the X, Y, and Z basis vectors are orthogonal to each other
     * (X \cdot Y = X \cdot Z=Y \cdot Z = 0),and that they are respectivelly orthogonal to the cartesian basis vectors, for example, with
     * the \hat{x} basis vector:\hat{x} \cdot Y = \hat{x} \cdot Z = 0 and \hat{x} \cdot X = X_x.
     * I believe that this might be generalizable to non-orthogonal basis, but I will not confirm that right now.
     * I also do not want to convert to cartesian coordinates either, as I think it is a waste of time.
     * I then use these d_i, d_j and d_k coefficients as parameters in a trilinear interpolation.
     **/

    // perform NX_j \cdot r to find the indices i, j and k of the cubefile.
    Eigen::Map<const Eigen::Vector3d> r_vec(&r[0]);
    Eigen::Vector3d coeff = normalized_basis * (r_vec - corner); // coefficients i, j and k in r = i*X + j*Y + k*Z
    Eigen::Vector3d d_index;
    Eigen::Vector3i index;
    for (auto i = 0; i < 3; i++) {
        index(i) = floor(coeff(i));       // nearest voxel point index that is lower than the point r
        d_index(i) = coeff(i) - index(i); // deviation from the nearest voxel point thats lower that the point r, it always in the interval [0,1)
    }

    // do the trilinear interpolation naively without loops or any logic (just plug in the equations)
    // TODO use the linear system form

    auto c000 = CUBE[0][(index(0)) * N_steps[1] * N_steps[2] + (index(1)) * N_steps[2] + (index(2))];
    auto c001 = CUBE[0][(index(0)) * N_steps[1] * N_steps[2] + (index(1)) * N_steps[2] + (1 + index(2))];
    auto c010 = CUBE[0][(index(0)) * N_steps[1] * N_steps[2] + (1 + index(1)) * N_steps[2] + (index(2))];
    auto c011 = CUBE[0][(index(0)) * N_steps[1] * N_steps[2] + (1 + index(1)) * N_steps[2] + (1 + index(2))];
    auto c100 = CUBE[0][(1 + index(0)) * N_steps[1] * N_steps[2] + (index(1)) * N_steps[2] + (index(2))];
    auto c101 = CUBE[0][(1 + index(0)) * N_steps[1] * N_steps[2] + (index(1)) * N_steps[2] + (1 + index(2))];
    auto c110 = CUBE[0][(1 + index(0)) * N_steps[1] * N_steps[2] + (1 + index(1)) * N_steps[2] + (index(2))];
    auto c111 = CUBE[0][(1 + index(0)) * N_steps[1] * N_steps[2] + (1 + index(1)) * N_steps[2] + (1 + index(2))];

    auto c00 = c000 * (1 - d_index(0)) + c100 * d_index(0);
    auto c01 = c001 * (1 - d_index(0)) + c101 * d_index(0);
    auto c10 = c010 * (1 - d_index(0)) + c110 * d_index(0);
    auto c11 = c011 * (1 - d_index(0)) + c111 * d_index(0);

    auto c0 = c00 * (1 - d_index(1)) + c10 * d_index(1);
    auto c1 = c01 * (1 - d_index(1)) + c11 * d_index(1);

    double c = c0 * (1 - d_index(2)) + c1 * d_index(2);

    // Alternativelly i can solve this as a linear system

    return c;
}

void CUBEfunction::normalize_basis() {
    // normalize the basis as X_j/(X_j\cdot X_j) = NX_j
    for (int i = 0; i < 3; i++) {
        normalized_basis.row(i) = voxel_axes.row(i) / voxel_axes.row(i).squaredNorm(); // should set the new normalized matrix with normalized vectors.
    }
}

} // namespace mrchem
