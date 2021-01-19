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

#include <array>
#include <vector>

#include "MRCPP/MWFunctions"

#include "chemistry/Nucleus.h"

namespace mrchem {

class NuclearFunction final : public mrcpp::RepresentableFunction<3> {
public:
    NuclearFunction() {}

    double evalf(const mrcpp::Coord<3> &r) const override;

    void push_back(const Nucleus &nuc, double S);

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;

protected:
    Nuclei nuclei;
    std::vector<double> smooth;
};

} // namespace mrchem
