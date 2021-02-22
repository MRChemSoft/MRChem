#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/one_electron/ZoraPotential.h"

namespace mrchem {

class ZoraOperator final : public RankZeroTensorOperator {
public:
    ZoraOperator(const Nuclei &nucs,
                 double zora_factor,
                 double proj_prec,
                 double smooth_prec = -1.0,
                 bool mpi_share = false,
                 int func_flag = 0) {
        r_m1 = std::make_shared<ZoraPotential>(nucs, zora_factor, proj_prec, smooth_prec, mpi_share, func_flag);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
        v.name() = "ZORA";
    }

private:
    std::shared_ptr<ZoraPotential> r_m1{nullptr};
};

} // namespace mrchem
