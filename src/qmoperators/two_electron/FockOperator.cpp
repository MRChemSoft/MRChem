#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "FarFieldOperator.h"
#include "FockOperator.h"
#include "XCOperator.h"
#include "chemistry/chemistry_utils.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using KineticOperator_p = std::shared_ptr<mrchem::KineticOperator>;
using NuclearOperator_p = std::shared_ptr<mrchem::NuclearOperator>;
using CoulombOperator_p = std::shared_ptr<mrchem::CoulombOperator>;
using FarFieldOperator_p = std::shared_ptr<mrchem::FarFieldOperator>;
using ExchangeOperator_p = std::shared_ptr<mrchem::ExchangeOperator>;
using XCOperator_p = std::shared_ptr<mrchem::XCOperator>;
using ElectricFieldOperator_p = std::shared_ptr<mrchem::ElectricFieldOperator>;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param t:   Kinetic operator
 * @param v:   Nuclear potential operator
 * @param j:   Coulomb operator
 * @param k:   HF exchange operator
 * @param xc:  Exchange-Correlation operator
 * @param ext: External Field operator
 *
 * Each of the arguments can be NULL, so this operators includes both core Hamiltonian,
 * the Hartree(-Fock) method and (pure/hybrid) Density Functional Theory.
 */
FockOperator::FockOperator(KineticOperator_p t,
                           NuclearOperator_p v,
                           CoulombOperator_p j,
                           FarFieldOperator_p f,
                           ExchangeOperator_p k,
                           XCOperator_p xc,
                           ElectricFieldOperator_p ext)
        : kin(t)
        , nuc(v)
        , coul(j)
        , ff(f)
        , ex(k)
        , xc(xc)
        , ext(ext) {}

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockOperator::build(double exx) {
    this->exact_exchange = exx;
    this->T = RankZeroTensorOperator();
    if (this->kin != nullptr) this->T += (*this->kin);

    this->V = RankZeroTensorOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->coul != nullptr) this->V += (*this->coul);
    if (this->ff != nullptr) this->V += (*this->ff);
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);

    RankZeroTensorOperator &F = (*this);
    F = this->kinetic() + this->potential();
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockOperator::setup(double prec) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Building Fock operator");
    mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
    mrcpp::print::separator(2, '-');
    this->kinetic().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec);
    t_tot.stop();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_tot);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockOperator::clear() {
    this->kinetic().clear();
    this->potential().clear();
    this->perturbation().clear();
}

/** @brief rotate orbitals of two-electron operators
 *
 * @param U: unitary transformation matrix
 *
 * This function should be used in case the orbitals are rotated *after* the FockOperator
 * has been setup. In particular the ExchangeOperator needs to rotate the precomputed
 * internal exchange potentials.
 */
void FockOperator::rotate(const ComplexMatrix &U) {
    if (this->ex != nullptr) this->ex->rotate(U);
}

/** @brief compute the SCF energy
 *
 * @param Phi: orbitals
 * @param F: Fock matrix
 *
 * This function will compute the total energy for a given OrbitalVector and
 * the corresponding Fock matrix. Tracing the kinetic energy operator is avoided
 * by tracing the Fock matrix and subtracting all other contributions.
 */
SCFEnergy FockOperator::trace(OrbitalVector &Phi, const Nuclei &nucs) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing molecular energy");

    double E_kin = 0.0;  // Kinetic energy
    double E_nn = 0.0;   // Nuclear repulsion
    double E_en = 0.0;   // Nuclear-electronic interaction
    double E_ee = 0.0;   // Electronic repulsion
    double E_x = 0.0;    // Exact Exchange
    double E_xc = 0.0;   // Exchange and Correlation
    double E_eext = 0.0; // External field contribution to the electronic energy
    double E_next = 0.0; // External field contribution to the nuclear energy
    double E_eeff = 0.0; // Far field contribution to the electronic energy
    double E_nff = 0.0;  // Far field contribution to the nuclear energy

    //// PERIODIC ENERGIES ////
    double E_se = 0.0;
    double E_ncor = 0.0;

    // Nuclear part
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();
    if (this->ff != nullptr) E_nff = -0.5 * this->ff->trace(nucs).real();

    // Electronic part
    if (this->kin != nullptr) E_kin = this->kin->trace(Phi).real();
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    if (this->ff != nullptr) E_eeff = this->ff->trace(Phi).real();

    if (this->coul->getPotential()->getRc() > 0.0) {
        E_nn = 0.0;
        E_se = -0.5 * mrcpp::dot<3>(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
        auto rc = this->coul->getPotential()->getRc();
        if (rc <= 0.0) MSG_ABORT("RC has to be positive");
        auto Ig = 10976.0 / (17875.0 * rc);
        for (auto &nuc : nucs) { E_ncor += 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig - 12.0 / (5.0 * rc)); }

        auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
        dip_oper.setup(1.0e-6);
        DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
        DoubleVector dip_el = dip_oper.trace(Phi).real();
        DoubleVector tot_dip = nuc_dip + dip_el;
        println(0, "dip_el nuc" << dip_el[0] << " " << dip_el[1] << " " << dip_el[2])
        dip_oper.clear();
        auto new_charge = tot_dip[2]/8.0;
        // Check RC TMP
	auto E_ncor_tmp = 0.0;
        E_ncor_tmp += new_charge * new_charge * (Ig - 12.0 / (5.0 * rc));
        E_ncor_tmp += new_charge * new_charge * (Ig - 12.0 / (5.0 * rc));
	println(0, "E_ncor_tmp " <<  E_ncor_tmp);  


    }

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, E_nff, E_eeff, E_se, E_ncor};
}

ComplexMatrix FockOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    auto t = this->getKineticOperator();
    auto v = this->potential();

    ComplexMatrix T = ComplexMatrix::Zero(bra.size(), ket.size());
    if (t != nullptr) T += (*t)(bra, ket);

    ComplexMatrix V = ComplexMatrix::Zero(bra.size(), ket.size());
    if (v.size() > 0) V += v(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T + V;
}

ComplexMatrix FockOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
