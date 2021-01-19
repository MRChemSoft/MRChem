#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SmearedNuclearPotential.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/print_utils.h"

#include "qmoperators/one_electron/H_E_dip.h"

using mrcpp::Printer;
using mrcpp::Timer;

using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;
namespace mrchem {

/** @brief projects an analytic expression for a smoothed nuclear potential
 *
 * @param[in] nucs collection of nuclei that defines the potential
 * @param[in] prec precision used both in smoothing and projection
 *
 * We create two different analytic functions for the nuclear potential:
 *
 * 1) this->func: The total potential from all nuclei of the system.
 *                This is needed later for analytic calculations.
 *
 * 2) loc_func: Temporary function that contains only some of the
 *              nuclei, which are distributed among the available
 *              MPIs. This is used only for the projection below.
 */
SmearedNuclearPotential::SmearedNuclearPotential(Nuclei nucs,
                                                 double proj_prec,
                                                 double smooth_prec,
                                                 double rc,
                                                 bool mpi_share,
                                                 OrbitalVector_p Phi)
        : NuclearPotential(1, mpi_share), rc(rc), orbitals(Phi) {
            this->nucs = nucs;
            this->proj_prec = proj_prec;
            this->smooth_prec = smooth_prec;
}

void SmearedNuclearPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    auto period = (*MRA).getWorldBox().getScalingFactor(0);
    auto nucs = periodic::periodify_nuclei(this->nucs, period * 2.0);
    auto rc = this->rc;
    println(0, "here?")
    OrbitalVector &Phi = *this->orbitals;
    println(0, "here?")

    auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
    dip_oper.setup(prec);
    DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
    DoubleVector dip_el = dip_oper.trace(Phi).real();
    DoubleVector tot_dip = nuc_dip + dip_el;
    println(0, "dip_el nuc" << dip_el[0] << " " << dip_el[1] << " " << dip_el[2])
    dip_oper.clear();
    auto new_charge = tot_dip[2]/8.0;
    mrcpp::Coord<3> pos_coord{0.0, 0.0, 4.0};
    mrcpp::Coord<3> neg_coord{0.0, 0.0, -4.0};
    std::vector<double> charges{new_charge, -new_charge};
    std::vector<mrcpp::Coord<3>> coords{pos_coord, neg_coord};

    std::vector<mrcpp::Coord<3>> pos_coords{};
    std::vector<mrcpp::Coord<3>> neg_coords{};

    pos_coords.push_back(pos_coord);
    neg_coords.push_back(neg_coord);

    // for (auto x  = -5; x < 6; x ++) {
    //     for (auto y = -5; y < 6; y++) {
    //         for (auto z = -5; z < 6; z++) {
    //             mrcpp::Coord<3> tmp_pos{pos_coord[0] + x*period*2.0, pos_coord[1] + y*period*2.0, pos_coord[2] + z*period*2.0 };
    //             mrcpp::Coord<3> tmp_neg{neg_coord[0] + x*period*2.0, neg_coord[1] + y*period*2.0, neg_coord[2] + z*period*2.0 };
    //             pos_coords.push_back(tmp_pos);
    //             neg_coords.push_back(tmp_neg);
    //         }
    //     }
    // }



    NuclearFunction loc_func;
    double c = 0.00435 * this->smooth_prec;
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double Z_5 = std::pow(Z, 5.0);
        double smooth = std::pow(c / Z_5, 1.0 / 3.0);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;

        this->func.push_back(nuc, smooth);
        if (mpi::orb_rank == proj_rank) loc_func.push_back(nuc, smooth);
    }

    auto V_smear = [nucs, rc, charges, pos_coords, neg_coords](const mrcpp::Coord<3> &r) -> double {
        auto v_g = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * nuc.getCharge();
        }
        for (auto i = 0; i < pos_coords.size(); i++) {
            auto R = math_utils::calc_distance(r, pos_coords[i]);
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * charges[1];
        }
        
        for (auto i = 0; i < neg_coords.size(); i++) {
            auto R = math_utils::calc_distance(r, neg_coords[i]);
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * charges[0];
        }

        return v_g;
    };

    // Scale precision by system size
    int Z_tot = chemistry::get_total_charge(nucs);
    double abs_prec = this->proj_prec / Z_tot;

    QMFunction V_loc(false);
    QMFunction V_smear_loc(false);

    qmfunction::project(V_loc, loc_func, NUMBER::Real, abs_prec);
    qmfunction::project(V_smear_loc, V_smear, NUMBER::Real, abs_prec);

    if (not this->hasReal()) this->alloc(NUMBER::Real);
    qmfunction::add(*this, 1.0, V_loc, -1.0, V_smear_loc, -1);
    println(0, "here? ")
}

} // namespace mrchem
