#pragma once

#include "NuclearPotential.h"
#include "SmearedNuclearPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class NuclearOperator final : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false) {
        r_m1 = std::make_shared<NuclearPotential>(nucs, proj_prec, smooth_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
        v.name() = "V_nuc";
    }

    NuclearOperator(const Nuclei &nucs,
                    double proj_prec,
                    double smooth_prec = -1.0,
                    double rc = 1.0,
                    bool mpi_share = false,
                    std::shared_ptr<mrchem::OrbitalVector> Phi = nullptr) {
        r_m1 = std::make_shared<SmearedNuclearPotential>(nucs, proj_prec, smooth_prec, rc, mpi_share, Phi);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
    }

private:
    std::shared_ptr<NuclearPotential> r_m1{nullptr};
};

} // namespace mrchem
