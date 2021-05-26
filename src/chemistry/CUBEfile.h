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

    std::string comments;
    std::vector<double> corner;
    int N_atoms;
    int N_val = 1;
    int N_steps[3];          // size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
    double voxel_axes[3][3]; // size 3x3 array of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate
    std::vector<int> atom_numbers;
    std::vector<double> atom_charges;
    std::vector<std::vector<double>> atom_coords;
    std::vector<int> DSET_IDS; // vector containing important information about data stored in each voxel point.
    std::vector<std::vector<std::vector<std::vector<double>>>>
        CUBE; // indexing here works as [x_step][y_step][z_step][value_number-1]. If there is only one value per voxel, the vector will still be 4 dimensional.
};
} // namespace mrchem
