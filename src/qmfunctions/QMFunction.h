/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include <memory>

#include "ComplexFunction.h"
#include "qmfunction_fwd.h"

namespace mrchem {

class QMFunction {
public:
    explicit QMFunction(bool share = false);
    QMFunction(const QMFunction &func);
    QMFunction &operator=(const QMFunction &func);
    QMFunction dagger();
    virtual ~QMFunction() = default;

    void release() { this->func_ptr.reset(); }
    void freeFunctions() { function().free(NUMBER::Total); }
    void clearFunctions();

    bool conjugate() const { return this->conj; }
    ComplexFunction &function() { return *this->func_ptr.get(); }
    const ComplexFunction &function() const { return *this->func_ptr.get(); }

    double norm() const;
    double squaredNorm() const;
    ComplexDouble integrate() const;

    void crop(double prec);
    void rescale(double c);
    void rescale(ComplexDouble c);
    void add(ComplexDouble c, QMFunction inp);

protected:
    bool conj;
    std::shared_ptr<ComplexFunction> func_ptr;
};

} // namespace mrchem
