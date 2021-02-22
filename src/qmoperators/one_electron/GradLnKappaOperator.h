#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/ZoraPotential.h"

namespace mrchem {

class GradLnKappaOperator final : public RankOneTensorOperator<3> {
public:
    GradLnKappaOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, std::shared_ptr<ZoraPotential> Z) {

        dkdx = std::make_shared<QMPotential>(-1, false);
        dkdy = std::make_shared<QMPotential>(-1, false);
        dkdz = std::make_shared<QMPotential>(-1, false);

        computeGradComponent(dkdx, 0, D, Z);
        computeGradComponent(dkdy, 1, D, Z);
        computeGradComponent(dkdz, 2, D, Z);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &d = (*this);
        //        p_x = std::make_shared<QMMomentum>(0, D);
        //        p_y = std::make_shared<QMMomentum>(1, D);
        //        p_z = std::make_shared<QMMomentum>(2, D);
        d[0] = dkdx;
        d[1] = dkdy;
        d[2] = dkdz;
        d[0].name() = "grad_x_kappa";
        d[1].name() = "grad_y_kappa";
        d[2].name() = "grad_z_kappa";
    }

private:
    std::shared_ptr<ZoraPotential> zora;
    std::shared_ptr<QMPotential> dkdx{nullptr};
    std::shared_ptr<QMPotential> dkdy{nullptr};
    std::shared_ptr<QMPotential> dkdz{nullptr};

    void computeGradComponent(std::shared_ptr<QMPotential> component,
                              int dir,
                              std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                              std::shared_ptr<ZoraPotential> Z);
};

} // namespace mrchem
