#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class NuclearPotential : public QMPotential {
public:

    NuclearPotential(int refinement, bool mpi_share) : QMPotential(refinement, mpi_share) { }
    NuclearPotential(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false);
    ~NuclearPotential() override { free(NUMBER::Total); }

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

protected:
    NuclearFunction func;

    void allreducePotential(double prec, QMFunction &V_loc);

    Nuclei nucs;
    double proj_prec;
    double smooth_prec;
};

} // namespace mrchem
