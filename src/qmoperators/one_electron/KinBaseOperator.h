#pragma once

#include "MomentumOperator.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class KineticOperator
 *
 * @brief Operator for kinetic energy
 *
 * This operator is constructed as the square of the more fundamental
 * MomentumOperator. The general base class functions for calculation of
 * expectation values are overwritten, as they can be improved due to
 * symmetry.
 *
 */

namespace mrchem {

class KinBaseOperator : public RankZeroTensorOperator {
public:

    virtual ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket) = 0;
    virtual ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket) = 0;

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

};

} // namespace mrchem
