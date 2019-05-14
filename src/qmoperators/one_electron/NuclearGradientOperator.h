#pragma once

#include "QMPotential.h"
#include "analyticfunctions/NuclearGradientFunction.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class NuclearGradientPotential final : public QMPotential {
public:
    NuclearGradientPotential(int d, const Nucleus &nuc, double c)
            : QMPotential(1)
            , func(d, nuc, c) {}

    void setup(double prec) override;
    void clear() override;

private:
    NuclearGradientFunction func;
};

class NuclearGradientOperator final : public RankOneTensorOperator<3> {
public:
    NuclearGradientOperator(const Nucleus &nuc, double c) {
        x_rm3 = std::make_shared<NuclearGradientPotential>(0, nuc, c);
        y_rm3 = std::make_shared<NuclearGradientPotential>(1, nuc, c);
        z_rm3 = std::make_shared<NuclearGradientPotential>(2, nuc, c);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator &v = (*this);
        v[0] = x_rm3;
        v[1] = y_rm3;
        v[2] = z_rm3;
        v[0].name() = "x/r^3";
        v[1].name() = "y/r^3";
        v[2].name() = "z/r^3";
    }

private:
    std::shared_ptr<NuclearGradientPotential> x_rm3{nullptr};
    std::shared_ptr<NuclearGradientPotential> y_rm3{nullptr};
    std::shared_ptr<NuclearGradientPotential> z_rm3{nullptr};
};

} // namespace mrchem
