#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "PositionOperator.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

#include "qmfunctions/Density.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

PositionPotential::PositionPotential(int d, const mrcpp::Coord<3> &o)
        : QMPotential(1) {
    auto f = [d, o](const mrcpp::Coord<3> &r) -> double { return (r[d] - o[d]); };
    this->func.set(f);
}

void PositionPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;
    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
}

ComplexVector PositionOperator::trace(Density rho) {
    ComplexVector mu = ComplexVector::Zero(3);
    mu[0] = qmfunction::dot(*this->r_x, rho);
    mu[1] = qmfunction::dot(*this->r_y, rho);
    mu[2] = qmfunction::dot(*this->r_z, rho);
    return mu;
}

void PositionPotential::clear() {
    free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();    // apply_prec = -1
}

} // namespace mrchem
