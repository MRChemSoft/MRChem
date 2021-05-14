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

void CUBEfile::readFile(std::string file_path) {
    double N_atoms;
    std::vector<double> corner;
    std::ifstream raw_cube_file;
    std::stringstream comments;
    bool DsetIDs_exists = false;
    double vals_per_grid_pt;
    double vol_coords[80][80][80];

    raw_cube_file.open(file_path, std::ios::in);

    if (!raw_cube_file) {
        std::cout << "No such file!!!! >:(";
    } else {
        std::string line;
        int i = 0;
        while (std::getline(raw_cube_file, line, '\n')) {
            if (i <= 1) comments << line << '\n';
            if (i == 2) {
                std::vector<double> origindata;
                std::stringstream origin(line);
                double value;
                while (origin >> value) { origindata.push_back(value); }
                if (origindata[0] < 0) DsetIDs_exists = true;
                N_atoms = std::fabs(origindata[0]);
                vals_per_grid_pt = origindata[4];
                corner = std::vector<double>(origindata.begin() + 1, origindata.end() - 1);
            }
            if ((3 <= i) and (i <= 5)) axes << line << "\n";
            if ((6 <= i) and (i <= (5 + N_atoms))) atoms << line << "\n";
            if ((i == (6 + N_atoms)) and (DSET_IDS_exists)) DSESET_IDS << line;
            if (i >= (7 + N_atoms)) i++;
        }
    }
    raw_cube_file.close();
    std::cout << comments.str();
    std::cout << corner[0] << " " << corner[1] << " " << corner[2] << "\n";
}

} // namespace mrchem
