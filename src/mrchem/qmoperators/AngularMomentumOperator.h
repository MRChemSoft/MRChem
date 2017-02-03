#ifndef ANGULARMOMENTUMOPERATOR_H
#define ANGULARMOMENTUMOPERATOR_H

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "MomentumOperator.h"
#include "DerivativeOperator.h"

class AngularMomentumOperator : public RankOneTensorOperator<3> {
public:
    AngularMomentumOperator(DerivativeOperator<3> &d, const double *o = 0)
            : r(o), p(d) {
        initializeTensorOperator();
    }
    virtual ~AngularMomentumOperator() { }

    void setup(double prec) {
        this->r.setup(prec);
        this->p.setup(prec);
    }

    void clear() {
        this->r.clear();
        this->p.clear();
    }

protected:
    PositionOperator r;
    MomentumOperator p;

    void initializeTensorOperator() {
        RankZeroTensorOperator &r_x = this->r[0];
        RankZeroTensorOperator &r_y = this->r[1];
        RankZeroTensorOperator &r_z = this->r[2];
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (r_y*p_z - r_z*p_y);
        h[1] = (r_z*p_x - r_x*p_z);
        h[2] = (r_x*p_y - r_y*p_x);
    }
};

#endif // ANGULARMOMENTUMOPERATOR_H
