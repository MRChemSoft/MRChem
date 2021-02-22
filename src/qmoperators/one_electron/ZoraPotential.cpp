#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ZoraPotential.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief projects an analytic expression for a smoothed Zora potential
 *
 * @param[in] nucs collection of nuclei that defines the potential
 * @param[in] prec precision used both in smoothing and projection
 *
 * We create two different analytic functions for the Zora potential:
 *
 * 1) this->func: The total potential from all nuclei of the system.
 *                This is needed later for analytic calculations.
 *
 * 2) loc_func: Temporary function that contains only some of the
 *              nuclei, which are distributed among the available
 *              MPIs. This is used only for the projection below.
 */
ZoraPotential::ZoraPotential(const Nuclei &nucs,
                             double zora_factor,
                             double proj_prec,
                             double smooth_prec,
                             bool mpi_share,
                             int func_flag)
        : QMPotential(1, mpi_share)
        , zoraFactor(zora_factor) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    int pprec = Printer::getPrecision();
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 8;
    int w3 = 2 * w0 / 9;
    int w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4, ' ');
    o_head << std::setw(w3) << "Charge";
    o_head << std::setw(w3) << "Precision";
    o_head << std::setw(w3) << "Smoothing";

    mrcpp::print::header(2, "Projecting Zora potential");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    NuclearFunction loc_func;

    double c = 0.00435 * smooth_prec;
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double Z_5 = std::pow(Z, 5.0);
        double smooth = std::pow(c / Z_5, 1.0 / 3.0);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;

        this->func.push_back(nuc, smooth);
        if (mpi::orb_rank == proj_rank) loc_func.push_back(nuc, smooth);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth;
        println(2, o_row.str());
    }

    Timer t_tot;

    // Scale precision by system size
    int Z_tot = chemistry::get_total_charge(nucs);
    double abs_prec = proj_prec / Z_tot;

    QMFunction V_loc(false);

    Timer t_loc;
    qmfunction::project(V_loc, loc_func, NUMBER::Real, abs_prec);
    t_loc.stop();

    Timer t_com;
    allreducePotential(abs_prec, V_loc);
    t_com.stop();

    std::cout << "Zora Factor " << zoraFactor << std::endl;

    Timer t_map;
    switch (func_flag) {
        case 0: // required for kinetic energy and SCF
            computeKappa(proj_prec);
            break;
        case 1: // maybe required for SCF
            computeLnKappa(proj_prec);
            break;
        case 2: // required for SCF
            computeKappaInv(proj_prec);
            break;
        case 3: // maybe required for SCF
            computeVKappaInv(proj_prec);
            break;
        default:
            MSG_ABORT("Invalid func_flag");
    }
    t_map.stop();

    t_tot.stop();
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local potential", V_loc, t_loc);
    print_utils::qmfunction(2, "Allreduce potential", V_loc, t_com);
    print_utils::qmfunction(2, "Map potential", V_loc, t_map);
    mrcpp::print::footer(2, t_tot, 2);
}

void ZoraPotential::allreducePotential(double prec, QMFunction &V_loc) {
    QMFunction &V_tot = *this;

    // Add up local contributions into the grand master
    mpi::reduce_function(prec, V_loc, mpi::comm_orb);
    if (mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mpi::numerically_exact) V_loc.crop(prec);
    }

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mpi::broadcast_function(V_loc, mpi::comm_sh_group);
        if (mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mpi::share_function(V_tot, 0, tag, mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mpi::broadcast_function(V_loc, mpi::comm_orb);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
}

/** @brief maps the potential function into the kappa function
 *
 * @param[in] prec precision
 *
 *  \f \kappa = \frac{1}{1-V/2c^2} \f
 */
void ZoraPotential::computeKappa(double prec) {
    auto zf = this->zoraFactor;
    auto kappamap = [zf](double val) {
        double temp = 1.0 - val / zf;
        return 1.0 / temp;
    };
    QMFunction &ZoraFunction = *this;
    mrcpp::FunctionTree<3> &zora_real = ZoraFunction.real();
    zora_real.map(kappamap);
}

/** @brief maps the potential function into the logarithm of kappa
 *
 * @param[in] prec precision
 *
 *  \f \ln(k) = -\ln({1-V/2c^2}) \f
 */
void ZoraPotential::computeLnKappa(double prec) {
    auto zf = this->zoraFactor;
    auto kappamap = [zf](double val) {
        double temp = 1.0 - val / zf;
        return -std::log(temp);
    };
    QMFunction &ZoraFunction = *this;
    mrcpp::FunctionTree<3> &zora_real = ZoraFunction.real();
    zora_real.map(kappamap);
}

/** @brief maps the potential function into a scaled potential
 *
 * @param[in] prec precision
 *
 *  \f z = 1-V/2c^2 = 1/\kappa\f
 *
 */
void ZoraPotential::computeKappaInv(double prec) {
    auto zf = this->zoraFactor;
    auto zoramap = [zf](double val) { return 1.0 - val / zf; };
    QMFunction &ZoraFunction = *this;
    mrcpp::FunctionTree<3> &zora_real = ZoraFunction.real();
    zora_real.map(zoramap);
}

/** @brief maps the potential function into a scaled potential times the potential
 *
 * @param[in] prec precision
 *
 *  \f z = V(1-V/2c^2) = V/\kappa \f
 *
 */
void ZoraPotential::computeVKappaInv(double prec) {
    auto zf = this->zoraFactor;
    auto zoramap = [zf](double val) {
        auto temp = 1.0 - val / zf;
        return val * temp;
    };
    QMFunction &ZoraFunction = *this;
    mrcpp::FunctionTree<3> &zora_real = ZoraFunction.real();
    zora_real.map(zoramap);
}

} // namespace mrchem
