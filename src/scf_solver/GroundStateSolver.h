/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <nlohmann/json.hpp>

#include "SCFSolver.h"
#include "properties/SCFEnergy.h"

/** @class GroundStateSolver
 *
 * @brief Ground state SCF optimization with kinetic energy operator
 *
 * This ground state SCF solver computes the Fock matrix explicitly by application of
 * the kinetic energy operator. This simplifies the algorithm significanly when a KAIN
 * iterative accelerator is present, or when diagonalization/localization occurs.
 * The derivative operator in the kinetic energy MIGHT affect the quadratic precision
 * in the energy, but this point seems to no longer be critical.
 *
 * NOTE: The old algorithm which completely avoided the kinetic energy operator has
 *       been removed, but can be recovered from the git history if needed in the
 *       future (search for EnergyOptimizer).
 */

namespace mrchem {

class Molecule;
class FockOperator;
class RankZeroTensorOperator;
template <int I> class RankOneTensorOperator;

class GroundStateSolver : public SCFSolver {
public:
    GroundStateSolver() = default;
    virtual ~GroundStateSolver() override = default;

    void setRotation(int iter) { this->rotation = iter; }
    void setLocalize(bool loc) { this->localize = loc; }
    void setCheckpointFile(const std::string &file) { this->chkFile = file; }
    void setZora(bool zora) { this->isZora = true; }
    void setKappaInv(std::shared_ptr<RankZeroTensorOperator> kappa_inv) { kappaInv = kappa_inv; }
    void setGradLnKappa(std::shared_ptr<RankOneTensorOperator<3>> grad_ln_kappa) { gradLnKappa = grad_ln_kappa; }

    nlohmann::json optimize(Molecule &mol, FockOperator &F);

protected:
    int rotation{0};      ///< Number of iterations between localization/diagonalization
    bool localize{false}; ///< Use localized or canonical orbitals
    bool isZora{false};   ///< ZORA scalar relativistic calculation
    std::string chkFile;  ///< Name of checkpoint file
    std::vector<SCFEnergy> energy;

    void reset() override;
    double calcPropertyError() const;
    void printProperty() const;
    void printParameters(const std::string &method) const;

    std::shared_ptr<RankZeroTensorOperator> kappaInv;
    std::shared_ptr<RankOneTensorOperator<3>> gradLnKappa;

    bool needLocalization(int nIter, bool converged) const;
    bool needDiagonalization(int nIter, bool converged) const;
};

} // namespace mrchem
