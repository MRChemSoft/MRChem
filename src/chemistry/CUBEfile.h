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
#include <MRCPP/MWFunctions>

namespace mrchem {

class CUBEfile final : public mrcpp::RepresentableFunction<3> {
public:
    CUBEfile(std::string file_path);
    double evalf(const mrcpp::Coord<3> &r) const override;

protected:
    void readFile(std::string file_path);
    void normalize_basis();
    std::string comments;
    Eigen::Vector3d corner;
    int N_atoms;
    int N_val = 1;
    std::array<int, 3> N_steps; // size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
    Eigen::Matrix3d voxel_axes; // size 3x3 matrix of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate
    std::vector<int> atom_numbers;
    std::vector<double> atom_charges;
    std::vector<std::vector<double>> atom_coords;
    std::vector<int> DSET_IDS;             // vector containing important information about data stored in each voxel point.
    std::vector<std::vector<double>> CUBE; // indexing here works as  [value_number-1][x_step + y_step + z_step].
    Eigen::Matrix3d normalized_basis;      //"normalized"", just multiply each row by its 1/norm^2
};
} // namespace mrchem
