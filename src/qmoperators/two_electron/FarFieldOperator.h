#pragma once

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "FarFieldPotential.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/RankZeroTensorOperator.h"

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"

/** @class getFarFieldOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class FarFieldOperator final : public RankZeroTensorOperator {
public:
    FarFieldOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     const Nuclei &nucs,
                     double exp_prec,
                     bool mpi_share = false) {
        potential = std::make_shared<FarFieldPotential>(P, Phi, nucs, exp_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ~FarFieldOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }
    const Nuclei &getNuclei() const { return this->potential->getNuclei(); }

    double getNucPrec() { return this->potential->getNucPrec(); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }
    using RankZeroTensorOperator::trace;

private:
    std::shared_ptr<FarFieldPotential> potential{nullptr};
};

} // namespace mrchem
