#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "HartreePotential.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/H_E_dip.h"

using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */

HartreePotential::HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                                   std::shared_ptr<OrbitalVector> Phi,
                                   const Nuclei &nucs,
                                   const double &rc)

        : CoulombPotential(P, Phi)
        , nuclei(nucs)
        , rc(rc) {}

void HartreePotential::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    Density rho_el(false);
    density::compute(prec, rho_el, Phi, DensityType::Total);
    println(0, "do we enter each iter? Hartree")
    auto period = (*MRA).getWorldBox().getScalingFactor(0);
    this->b_smeared =
        chemistry::hack_density(this->apply_prec, this->nuclei, this->rc, period * 2.0, Phi);
        // chemistry::compute_nuclear_density_smeared(this->apply_prec, this->nuclei, this->rc, period * 2.0);
    qmfunction::add(rho, 1.0, rho_el, -1.0, this->b_smeared, -1.0);
    PositionOperator r({0.0, 0.0, 0.0});
    r.setup(prec);
    DoubleVector mu_el = r.trace(rho_el).real();
    println(0, "mu_el " << mu_el[0] << " " << mu_el[1] << " " << mu_el[2])

    DoubleVector mu_tot = r.trace(rho).real();
    println(0, "mu_tot " << mu_tot[0] << " " << mu_tot[1] << " " << mu_tot[2])
    r.clear();
    if (std::norm(rho.integrate()) > this->apply_prec) MSG_ABORT("Non-zero net charge on unit cell");

    timer.stop();
    // double t = timer.getWallTime();
    // int n = rho.getNNodes(NUMBER::Total);
    // Printer::printTree(0, "Hartree density", n, t);
}

} // namespace mrchem
