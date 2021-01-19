#pragma once

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class PositionPotential final : public QMPotential {
public:
    PositionPotential(int d, const mrcpp::Coord<3> &o);

protected:
    mrcpp::AnalyticFunction<3> func;

    void setup(double prec) override;
    void clear() override;
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const mrcpp::Coord<3> &o = {0.0, 0.0, 0.0}) {
        r_x = std::make_shared<PositionPotential>(0, o);
        r_y = std::make_shared<PositionPotential>(1, o);
        r_z = std::make_shared<PositionPotential>(2, o);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator &r = (*this);
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
        r[0].name() = "r[x]";
        r[1].name() = "r[y]";
        r[2].name() = "r[z]";
    }

    ComplexVector trace(Density rho);

protected:
    std::shared_ptr<PositionPotential> r_x{nullptr};
    std::shared_ptr<PositionPotential> r_y{nullptr};
    std::shared_ptr<PositionPotential> r_z{nullptr};
};

} // namespace mrchem
