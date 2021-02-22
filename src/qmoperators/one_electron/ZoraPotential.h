#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class ZoraPotential final : public QMPotential {
public:
    ZoraPotential(const Nuclei &nucs,
                  double zora_factor,
                  double proj_prec,
                  double smooth_prec = -1.0,
                  bool mpi_share = false,
                  int func_flag = 0);
    ~ZoraPotential() override { free(NUMBER::Total); }

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

private:
    double zoraFactor;
    void computeVKappaInv(double prec = -1.0);
    void computeKappaInv(double prec = -1.0);
    void computeKappa(double prec = -1.0);
    void computeLnKappa(double prec = -1.0);
    NuclearFunction func;
    void allreducePotential(double prec, QMFunction &V_loc);
};

} // namespace mrchem
