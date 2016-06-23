#include "SCFEnergy.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

SCFEnergy::SCFEnergy()
    : E_nuc(0.0),
      E_el(0.0),
      E_orb(0.0),
      E_kin(0.0),
      E_en(0.0),
      E_ee(0.0),
      E_xc(0.0),
      E_x(0.0) {
}

void SCFEnergy::clearNuclear() {
    this->E_nuc = 0.0;
}

void SCFEnergy::clearElectronic() {
    this->E_el = 0.0;

    this->E_orb = 0.0;
    this->E_kin = 0.0;
    this->E_en = 0.0;
    this->E_ee = 0.0;
    this->E_xc = 0.0;
    this->E_x = 0.0;
}

void SCFEnergy::compute(const Nuclei &nuclei) {
    Timer timer;
    timer.restart();

    TelePrompter::printHeader(0, "Calculating Nuclear Energy");

    clearNuclear();

    int nNucs = nuclei.size();
    for (int i = 0; i < nNucs; i++) {
        for (int j = i+1; j < nNucs; j++) {
            const Nucleus &nuc_i = nuclei[i];
            const Nucleus &nuc_j = nuclei[j];
            double Z_i = nuc_i.getCharge();
            double Z_j = nuc_j.getCharge();
            const double *R_i = nuc_i.getCoord();
            const double *R_j = nuc_j.getCoord();
            double r_ij = MathUtils::calcDistance(3, R_i, R_j);
            this->E_nuc += (Z_i*Z_j)/r_ij;
        }
    }
    int oldPrec = TelePrompter::setPrecision(15);
    println(0, " Nuclear energy              " << setw(30) << this->E_nuc    );
    TelePrompter::setPrecision(oldPrec);
    TelePrompter::printFooter(0, timer, 2);
}

void SCFEnergy::compute(FockOperator &f_oper, MatrixXd &f_mat, OrbitalVector &phi) {
    Timer timer;
    timer.restart();

    TelePrompter::printHeader(0, "Calculating Electronic Energy");

    clearElectronic();
    double E_xc2 = 0.0;

    NuclearPotential *V = f_oper.getNuclearPotential();
    CoulombOperator *J = f_oper.getCoulombOperator();
    ExchangeOperator *K = f_oper.getExchangeOperator();
    XCOperator *XC = f_oper.getXCOperator();

    if (XC != 0) this->E_xc = XC->getEnergy();
    for (int i = 0; i < phi.size(); i++) {
        Orbital &phi_i = phi.getOrbital(i);
        double occ = (double) phi_i.getOccupancy();
        double e_i = occ*f_mat(i,i);
        this->E_orb += e_i;

        if (V != 0) {
            println(2, "\n<" << i << "|V_nuc|" << i << ">");
            this->E_en += occ*(*V)(phi_i,phi_i);
        }
        if (J != 0) {
            println(2, "\n<" << i << "|J|" << i << ">");
            this->E_ee += 0.5*occ*(*J)(phi_i,phi_i);
        }
        if (K != 0) {
            println(2, "\n<" << i << "|K|" << i << ">");
            this->E_x += (*K)(phi_i,phi_i);
        }
        if (XC != 0) {
            println(2, "\n<" << i << "|V_xc|" << i << ">");
            E_xc2 += occ*(*XC)(phi_i,phi_i);
        }
    }
    double E_eex = this->E_ee + this->E_x;
    double E_orbxc2 = this->E_orb - E_xc2;
    this->E_kin = E_orbxc2 - 2.0*E_eex - this->E_en;
    this->E_el = E_orbxc2 - E_eex + this->E_xc;

    int oldPrec = TelePrompter::setPrecision(15);
    println(2, "                                                            ");
    println(0, " Electronic energy           " << setw(30) << this->E_el     );
    println(2, "                                                            ");
    TelePrompter::setPrecision(oldPrec);
    TelePrompter::printFooter(0, timer, 2);
}
