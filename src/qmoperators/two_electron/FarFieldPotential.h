#pragma once

#include "CoulombPotential.h"
#include "chemistry/Nucleus.h"

namespace mrchem {

class FarFieldPotential final : public QMPotential {
public:
    FarFieldPotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                      std::shared_ptr<OrbitalVector> Phi,
                      const Nuclei &nucs,
                      double exp_prec,
                      bool mpi_share = false);
    ~FarFieldPotential() override = default;

    friend class FarFieldOperator;

private:
    double exp_prec{1.0};
    Nuclei nuclei{};
    Density density;

    std::shared_ptr<OrbitalVector> orbitals; ///< Unperturbed orbitals defining the ground-state electron density
    std::shared_ptr<mrcpp::PoissonOperator> poisson; ///< Operator used to compute the potential

    auto &getPoisson() { return this->poisson; }
    auto &getDensity() { return this->density; }
    const Nuclei &getNuclei() { return this->nuclei; }
    bool hasDensity() const { return (this->density.squaredNorm() < 0.0) ? false : true; }

    double getNucPrec() { return this->exp_prec; }

    void setup(double prec);

    void clear() override;

    void setupDensity(double prec);
    void setupPotential(double prec);
};

} // namespace mrchem
