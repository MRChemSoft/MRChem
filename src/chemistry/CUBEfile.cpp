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

#include "CUBEfile.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace mrchem {

CUBEfile::CUBEfile(std::string file_path) {
    readFile(file_path);
}

// Do a quadrature of the file
double CUBEfile::evalf(const mrcpp::Coord<3> &r) const {}

// first draft on reading a cubefile
void CUBEfile::readFile(std::string file_path) {
    // open cubefile and read its contents into a vector

    std::ifstream raw_cube_file;
    std::vector<std::string> raw_cube_vector;

    raw_cube_file.open(file_path, std::ios::in);

    if (!raw_cube_file) {
        std::cout << "No such file!";
    } else {
        std::string line;
        while (std::getline(raw_cube_file, line, '\n')) { raw_cube_vector.push_back(line); }
    }
    raw_cube_file.close();

    // get out important info from vector

    // get comments
    std::string comments = raw_cube_vector[0] + "\n" + raw_cube_vector[1];

    // get corner of cubeplot
    std::vector<double> corner;
    int N_atoms;
    int N_val = 1; // NUmber of values per voxel.
    double value;
    std::stringstream origin(raw_cube_vector[2]);
    int i = 0;
    while (origin >> value) {
        i++;
        if (i == 1) {
            N_atoms = int(value);
            continue;
        }
        if (i == 5) {
            N_val = int(value);
            continue;
        }
        corner.push_back(value);
    }

    // get all voxel axis values all of this might be better to do while reading the file
    int N_steps[3];          // size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
    double voxel_axes[3][3]; // size 3x3 array of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate

    for (int j = 0; j != 3; j++) {
        std::stringstream datastream(raw_cube_vector[3 + j]);
        i = 0;
        while (datastream >> value) {
            if (i == 0) {
                N_steps[j] = int(value);
            } else {
                voxel_axes[j][i - 1] = value;
            }
            i++;
        }
    }
    std::cout << N_steps[0] << "  " << voxel_axes[0][0] << "  " << voxel_axes[0][1] << "  " << voxel_axes[0][2] << "\n";
    std::cout << N_steps[1] << "  " << voxel_axes[1][0] << "  " << voxel_axes[1][1] << "  " << voxel_axes[1][2] << "\n";
    std::cout << N_steps[2] << "  " << voxel_axes[2][0] << "  " << voxel_axes[2][1] << "  " << voxel_axes[2][2] << "\n";

    // get all atom coordinates
    int atom_numbers[std::abs(N_atoms)];
    double atom_charges[std::abs(N_atoms)];
    double atom_coords[std::abs(N_atoms)][3]; // 3D coordinates of each atom
    for (int j = 0; j < std::abs(N_atoms); j++) {
        std::stringstream datastream(raw_cube_vector[6 + j]);
        i = 0;
        while (datastream >> value) {
            if (i == 0) {
                atom_numbers[j] = value;
            } else if (i == 1) {
                atom_charges[j] = value;
            } else {
                atom_coords[j][i - 2] = value;
            }
            i++;
        }
    }
    std::cout << atom_numbers[0] << "  " << atom_charges[0] << "  " << atom_coords[0][0] << "  " << atom_coords[0][1] << "  " << atom_coords[0][2] << "\n";
    std::cout << atom_numbers[1] << "  " << atom_charges[1] << "  " << atom_coords[1][0] << "  " << atom_coords[1][1] << "  " << atom_coords[1][2] << "\n";
    std::cout << atom_numbers[2] << "  " << atom_charges[2] << "  " << atom_coords[2][0] << "  " << atom_coords[2][1] << "  " << atom_coords[2][2] << "\n";

    // get the DSET_IDS if there are any
    std::vector<int> DSET_IDS; // vector containing important information about data stored in each voxel point.
    if (N_atoms < 0) {
        std::stringstream DSET_IDS_stream(raw_cube_vector[6 + std::abs(N_atoms)]);
        i = 0;
        while (DSET_IDS_stream >> value) { DSET_IDS.push_back(int(value)); }
    } else {
        DSET_IDS.push_back(1);
    }
    for (auto iter = DSET_IDS.begin(); iter != DSET_IDS.end(); iter++) { std::cout << *iter << "    "; }
    std::cout << "\n";

    // get the voxel data
    std::vector<double> cube_data;
    int start = (N_atoms < 0) ? (6 + std::abs(N_atoms) + 1) : (6 + N_atoms);
    // first extract the space separated voxel values into a one dimensional vector
    for (auto iter = raw_cube_vector.begin() + start; iter != raw_cube_vector.end(); iter++) {
        std::stringstream linestream(*iter);
        while (linestream >> value) { cube_data.push_back(value); }
    }
    // now extract into a 3-dimensional array of size N_steps[0]*N_steps[1]*N_steps[2] must include the amount of DSET_IDS and N_vals afterwards
    double CUBE_array[N_steps[0]][N_steps[1]][N_steps[2]];
    for (auto i = 0; i < N_steps[0]; i++) { // i goes from 0 to N_X
        std::vector<std::vector<double>> j_vec;
        for (auto j = 0; j < N_steps[1]; j++) { // j goes from 0 to N_Y
            std::vector<double> k_vec;
            for (auto k = 0; k < N_steps[2]; k++) { // k goes from 0 to N_Z
                k_vec.push_back(cube_data[i + j + k]);
                CUBE_array[i][j][k] = cube_data[i + j + k]; // this should place the correct value into the correct place in the array
            }
            j_vec.push_back(k_vec);
        }
        (*this).CUBE.push_back(j_vec);
    }
}

} // namespace mrchem
