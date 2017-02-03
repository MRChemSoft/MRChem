#ifndef KINETICOPERATOR_H
#define KINETICOPERATOR_H

#include "QMTensorOperator.h"
#include "MomentumOperator.h"

class KineticOperator : public RankZeroTensorOperator {
public:
    KineticOperator(DerivativeOperator<3> &d) : p(d) {
        initializeTensorOperator();
    }
    virtual ~KineticOperator() { }

    virtual void setup(double prec) { this->p.setup(prec); }
    virtual void clear() { this->p.clear(); }

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::adjoint;

protected:
    MomentumOperator p;

    void initializeTensorOperator() {
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankZeroTensorOperator &h = *this;
        h = -0.5*(p_x*p_x + p_y*p_y + p_z*p_z);
    }
};

#endif // KINETICOPERATOR_H
