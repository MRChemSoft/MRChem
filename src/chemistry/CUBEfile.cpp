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
#include <Eigen/Core>
#include <Eigen/StdVector>
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
double CUBEfile::evalf(const mrcpp::Coord<3> &r) const {
    // normalize the basis as X_j/(X_j\cdot X_j) = NX_j
    Eigen::Matrix3d normalized_basis; //"normalized"", just multiply each row by its 1/norm^2
    for (int i = 0; i < 3; i++) {
        normalized_basis.row(i) = voxel_axes.row(i) / voxel_axes.row(i).squaredNorm(); // should set the new normalized matrix with normalized vectors.
    }

    // perform NX_j \cdot r to find the indices i, j and k of the cubefile.
    Eigen::Map<const Eigen::Vector3d> r_vec(&r[0]);
    std::cout << normalized_basis << "  *  " << r_vec << "\n";
    Eigen::Vector3d coeff_vec = normalized_basis * (r_vec - corner);
    std::cout << coeff_vec << "\n";

    double wop = 2.0;
    return wop;
}

// first draft on reading a cubefile
// for the second try I might want to go for a single while loop to extract everything into the respective variables with an index i and if else blocks
// if that does nothing on the performance then i can instead try to reuse some variables.
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
    comments = raw_cube_vector[0] + "\n" + raw_cube_vector[1];

    // get corner of cubeplot
    double value;

    std::stringstream origin(raw_cube_vector[2]); // might either use the same stringstream for all parameters or do it all at the start loop
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
        corner(i - 2) = value;
    }
    std::cout << corner << "\n";
    // get all voxel axis values all of this might be better to do while reading the file

    for (int j = 0; j != 3; j++) {
        std::stringstream datastream(raw_cube_vector[3 + j]);
        i = 0;
        while (datastream >> value) {
            if (i == 0) {
                N_steps[j] = int(value);
            } else {
                voxel_axes(j, (i - 1)) = value;
            }
            i++;
        }
    }

    // get all atom coordinates
    for (int j = 0; j < std::abs(N_atoms); j++) {
        std::stringstream datastream(raw_cube_vector[6 + j]);
        i = 0;
        std::vector<double> i_vec;
        while (datastream >> value) {
            if (i == 0) {
                atom_numbers.push_back(value);
            } else if (i == 1) {
                atom_charges.push_back(value);
            } else {
                i_vec.push_back(value);
            }
            i++;
        }
        atom_coords.push_back(i_vec);
    }

    // get the DSET_IDS if there are any.
    // TODO take into account the case where too many values are included that take more than just one line.
    if (N_atoms < 0) {
        std::stringstream DSET_IDS_stream(raw_cube_vector[6 + std::abs(N_atoms)]);
        i = 0;
        while (DSET_IDS_stream >> value) { DSET_IDS.push_back(int(value)); }
    } else {
        DSET_IDS.push_back(1);
        DSET_IDS.push_back(1);
    }

    // set number of vals dependent on DSET_IDS

    if (N_atoms < 0) { N_val = DSET_IDS[0]; }

    // get the voxel data
    std::vector<double> cube_data;
    int start = (N_atoms < 0) ? (6 + std::abs(N_atoms) + 1) : (6 + N_atoms);
    // first extract the space separated voxel values into a one dimensional vector
    for (auto iter = raw_cube_vector.begin() + start; iter != raw_cube_vector.end(); iter++) {
        std::stringstream linestream(*iter);
        while (linestream >> value) { cube_data.push_back(value); }
    }

    // now extract into a 4-dimensional vector of size N_steps[0]*N_steps[1]*N_steps[2]*N_val must include the amount of DSET_IDS and N_vals afterwards
    for (auto l = 0; l < N_val; l++) { //  l goes from 0 to N_val
        std::vector<double> mini_cube;
        for (auto i = 0; i < N_steps[0]; i++) {         // i goes from 0 to N_X
            for (auto j = 0; j < N_steps[1]; j++) {     // j goes from 0 to N_Y
                for (auto k = 0; k < N_steps[2]; k++) { // k goes from 0 to N_Z
                    mini_cube.push_back(cube_data[i + j + k + l]);
                }
            }
        }
        CUBE.push_back(mini_cube);
    }
}

} // namespace mrchem
