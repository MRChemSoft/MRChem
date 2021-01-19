#pragma once

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "HartreePotential.h"
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
    explicit CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, bool mpi_share = false) {
        potential = std::make_shared<CoulombPotential>(P, nullptr, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                    std::shared_ptr<OrbitalVector> Phi,
                    bool mpi_share = false) {
        potential = std::make_shared<CoulombPotentialD1>(P, Phi, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                    std::shared_ptr<OrbitalVector> Phi,
                    std::shared_ptr<OrbitalVector> X,
                    std::shared_ptr<OrbitalVector> Y,
                    bool mpi_share = false) {
        potential = std::make_shared<CoulombPotentialD2>(P, Phi, X, Y, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }

    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                    std::shared_ptr<OrbitalVector> Phi,
                    const Nuclei &nucs,
                    const double &rc) {

        potential = std::make_shared<HartreePotential>(P, Phi, nucs, rc);
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ~CoulombOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }
    auto &getPotential() { return this->potential; }

private:
    std::shared_ptr<CoulombPotential> potential{nullptr};
};

} // namespace mrchem
