#include "MRCPP/Gaussians"
#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include <cmath>
#include <functional>

#include "ReactionPotential.h"
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "scf_solver/KAIN.h"

using mrcpp::ABGVOperator;
using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ReactionPotential::ReactionPotential(mrcpp::PoissonOperator *P,
                                     mrcpp::DerivativeOperator<3> *D,
                                     Cavity *C,
                                     const Nuclei &nucs,
                                     OrbitalVector *Phi,
                                     int hist)
        : QMPotential(1, false)
        , cavity(C)
        , nuclei(nucs)
        , orbitals(Phi)
        , poisson(P)
        , derivative(D)
        , rho_tot(false)
        , rho_el(false)
        , rho_nuc(false)
        , history(hist) {
    this->electronicEnergy = 0.0;
    this->nuclearEnergy = 0.0;
    this->totalEnergy = 0.0;
}

void ReactionPotential::setEpsilon(bool is_inv, QMFunction &cavity_func) {
    cavity->implementEpsilon(is_inv, cavity->isLinear());
    mrcpp::build_grid(cavity_func.real(), *cavity);
    qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec / 100);
    cavity->implementEpsilon(false, cavity->isLinear());
    if (not is_inv) this->d_coefficient = cavity->getCoefficient();
}

void ReactionPotential::setRhoEff(QMFunction const &inv_eps_func, QMFunction &rho_eff_func) {

    rho_nuc = chemistry::compute_nuclear_density(this->apply_prec, this->nuclei, 1000);
    density::compute(this->apply_prec, rho_el, *orbitals, DENSITY::Total);
    rho_el.rescale(-1.0);
    qmfunction::add(rho_tot, 1.0, rho_el, 1.0, rho_nuc, -1.0);

    auto onesf = [](const mrcpp::Coord<3> &r) { return 1.0; };

    QMFunction ones;
    QMFunction tmp; // can still be made shorter, possibly
    qmfunction::project(ones, onesf, NUMBER::Real, this->apply_prec);
    qmfunction::add(tmp, -1.0, ones, 1.0, inv_eps_func, -1.0);
    qmfunction::multiply(rho_eff_func, rho_tot, tmp, this->apply_prec);
}

void ReactionPotential::setGamma(QMFunction const &inv_eps_func,
                                 QMFunction &gamma_func,
                                 QMFunction &temp_func1,
                                 mrcpp::FunctionTreeVector<3> &d_cavity) {

    auto d_V = mrcpp::gradient(*derivative, temp_func1.real());
    if (cavity->isLinear()) {
        QMFunction temp_func2;
        temp_func2.alloc(NUMBER::Real);
        mrcpp::dot(this->apply_prec, temp_func2.real(), d_V, d_cavity);
        qmfunction::multiply(gamma_func, temp_func2, inv_eps_func, this->apply_prec);
    } else if (not cavity->isLinear()) {
        mrcpp::dot(this->apply_prec, gamma_func.real(), d_V, d_cavity);
    }
    gamma_func.rescale(1.0 / (4.0 * MATHCONST::pi));
    mrcpp::clear(d_V, true);
}

void ReactionPotential::grad_G(QMFunction &gamma_func,
                               QMFunction &cavity_func,
                               QMFunction &rho_tot,
                               QMFunction &grad_G_func) {
    // fix this
    QMFunction temp_func;
    qmfunction::multiply(temp_func, gamma_func, cavity_func, this->apply_prec);
    qmfunction::add(grad_G_func, 1.0, temp_func, -1.0, rho_tot, -1.0);
}

void ReactionPotential::accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain) {
    OrbitalVector phi_n(0);
    OrbitalVector dPhi_n(0);
    phi_n.push_back(Orbital(SPIN::Paired));
    dPhi_n.push_back(Orbital(SPIN::Paired));

    phi_n[0].QMFunction::operator=(temp);
    dPhi_n[0].QMFunction::operator=(diff_func);

    // double plevel = TelePrompter::setPrintLevel(-1);
    kain.accelerate(this->apply_prec, phi_n, dPhi_n);
    // TelePrompter::setPrintLevel(plevel);

    temp.QMFunction::operator=(phi_n[0]);
    diff_func.QMFunction::operator=(dPhi_n[0]);

    phi_n.clear();
    dPhi_n.clear();
}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);

    QMFunction V_0_func;
    QMFunction cavity_func;
    QMFunction inv_eps_func;
    QMFunction rho_eff_func;
    QMFunction gamma_func;
    QMFunction &temp = *this;
    mrcpp::FunctionTreeVector<3> d_cavity;

    cavity_func.alloc(NUMBER::Real);
    setEpsilon(false, cavity_func);

    inv_eps_func.alloc(NUMBER::Real);
    setEpsilon(true, inv_eps_func);

    setRhoEff(inv_eps_func, rho_eff_func);

    cavity_func.rescale(this->d_coefficient);
    d_cavity = mrcpp::gradient(*derivative, cavity_func.real());

    V_0_func.alloc(NUMBER::Real);
    mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

    if (not temp.hasReal()) {
        QMFunction tmp_poisson;
        tmp_poisson.alloc(NUMBER::Real);
        gamma_func.alloc(NUMBER::Real);
        setGamma(inv_eps_func, gamma_func, V_0_func, d_cavity);
        qmfunction::add(tmp_poisson, 1.0, gamma_func, 1.0, rho_eff_func, -1.0);
        QMFunction test_func;
        test_func.alloc(NUMBER::Real);
        mrcpp::apply(prec, test_func.real(), *poisson, tmp_poisson.real());
        temp = test_func;
    }
    // double A[3] = {0, 0, -7};
    // double B[3] = {0, 0, 7};
    // mrcpp::Plotter<3> plt(280000, A, B);
    // plt.linePlot(rho_eff_func.real(), "rho_eff");
    KAIN kain(this->history);
    auto error = 1.00;
    int iter = 1;
    while (error >= prec) {
        gamma_func.free(NUMBER::Real);
        QMFunction temp_func;
        QMFunction V_np1_func;
        QMFunction diff_func;
        QMFunction temp_func1;

        temp_func1.alloc(NUMBER::Real);
        qmfunction::add(temp_func1, 1.0, temp, 1.0, V_0_func, -1.0);
        gamma_func.alloc(NUMBER::Real);
        setGamma(inv_eps_func, gamma_func, temp_func1, d_cavity);
        V_np1_func.alloc(NUMBER::Real);
        qmfunction::add(temp_func, 1.0, rho_eff_func, 1.0, gamma_func, -1.0);
        mrcpp::apply(prec, V_np1_func.real(), *poisson, temp_func.real());
        qmfunction::add(diff_func, -1.0, temp, 1.0, V_np1_func, -1.0);
        error = diff_func.norm();

        // std::string gamma_file = "gamma";
        // gamma_file +=  std::to_string(iter);
        // plt.linePlot(gamma_func.real(), gamma_file);
        // plt.linePlot(temp.real(), "V_r");
        // plt.linePlot(diff_func.real(), "difference");
        if (iter > 1 and this->history > 0) accelerateConvergence(diff_func, temp, kain);
        V_np1_func.free(NUMBER::Real);
        qmfunction::add(V_np1_func, 1.0, temp, 1.0, diff_func, -1.0);
        // plt.linePlot(V_np1_func.real(), "V_np1(2)");
        temp = V_np1_func;

        std::cout << "iter.:\t" << iter << "\n"
                  << "error:\t" << error << std::endl;
        iter++;
    }

    mrcpp::clear(d_cavity, true);
}

double &ReactionPotential::getTotalEnergy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_tot, *this, this->apply_prec);
    totalEnergy = temp_prod_func.integrate().real();
    return totalEnergy;
}

double &ReactionPotential::getElectronicEnergy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_el, *this, this->apply_prec);
    electronicEnergy = temp_prod_func.integrate().real();
    return electronicEnergy;
}

double &ReactionPotential::getNuclearEnergy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_nuc, *this, this->apply_prec);
    nuclearEnergy = temp_prod_func.integrate().real();
    return nuclearEnergy;
}

void ReactionPotential::clear() {
    clearApplyPrec();
    rho_tot.free(NUMBER::Real);
    rho_el.free(NUMBER::Real);
    rho_nuc.free(NUMBER::Real);
    // QMFunction::free(NUMBER::Total);
}

} // namespace mrchem
