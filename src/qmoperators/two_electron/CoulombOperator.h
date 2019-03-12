#pragma once

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class CoulombOperator final : public RankZeroTensorOperator {
public:
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P) {
        potential = std::make_shared<CoulombPotential>(P);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi) {
        potential = std::make_shared<CoulombPotentialD1>(P, Phi);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                    std::shared_ptr<OrbitalVector> Phi,
                    std::shared_ptr<OrbitalVector> X,
                    std::shared_ptr<OrbitalVector> Y) {
        potential = std::make_shared<CoulombPotentialD2>(P, Phi, X, Y);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    ~CoulombOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }

private:
    std::shared_ptr<CoulombPotential> potential{nullptr};
};

} // namespace mrchem
