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
    std::ifstream raw_cube_file;
    raw_cube_file.open(file_path, std::ios::in);
    std::string raw_cube_string;
    std::stringstream raw_str_stream;
    if (!raw_cube_file) {
        throw "No such file";
    } else {

        while (std::getline(raw_cube_file, raw_cube_string)) {
            raw_str_stream << raw_cube_string + "\n";
            std::cout << raw_cube_string;
        }
    }
    raw_cube_file.close();
    raw_cube_string.empty();
    raw_cube_string = raw_str_stream.str();
    std::cout << raw_cube_string;
}

} // namespace mrchem
