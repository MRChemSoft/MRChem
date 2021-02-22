#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "GradLnKappaOperator.h"
#include "ZoraPotential.h"

#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void GradLnKappaOperator::computeGradComponent(std::shared_ptr<QMPotential> component,
                                               int dir,
                                               std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                                               std::shared_ptr<ZoraPotential> Z) {

    if (D == nullptr) MSG_ERROR("No derivative operator");
    // Calc real part
    if (Z->hasReal()) {
        component->alloc(NUMBER::Real);
        mrcpp::apply(component->real(), *D, Z->real(), dir);
    }
    // Calc imag part
    if (Z->hasImag()) {
        component->alloc(NUMBER::Imag);
        mrcpp::apply(component->imag(), *D, Z->imag(), dir);
        if (Z->conjugate()) component->imag().rescale(-1.0);
    }
}
} // namespace mrchem
