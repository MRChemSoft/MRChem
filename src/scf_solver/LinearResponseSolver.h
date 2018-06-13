#pragma once

#include "SCF.h"

/** @class LinearResponseSolver
 *
 */

namespace mrchem {

class LinearResponseSolver final : public SCF {
public:
    LinearResponseSolver(HelmholtzVector &h, Accelerator *k_x = 0, Accelerator *k_y = 0);

    void setupUnperturbed(double prec, FockOperator *fock, OrbitalVector *Phi, ComplexMatrix *F);
    void clearUnperturbed();

    void setup(FockOperator *fock, OrbitalVector *X);
    void setup(double omega, FockOperator *fock, OrbitalVector *X, OrbitalVector *Y);
    void clear();

    bool optimize();

protected:
    bool dynamic;
    double frequency;

    FockOperator *fOper_0;
    FockOperator *fOper_1;

    ComplexMatrix *fMat_0;
    ComplexMatrix *fMat_x;
    ComplexMatrix *fMat_y;

    OrbitalVector *orbitals_0;
    OrbitalVector *orbitals_x;
    OrbitalVector *orbitals_y;

    Accelerator *kain_x;
    Accelerator *kain_y;

    OrbitalVector setupHelmholtzArguments(OrbitalVector &dPhi,
                                          const ComplexMatrix &M,
                                          bool adjoint);

    double calcProperty();
    double calcPropertyError() const;
    void printProperty() const;
};

} //namespace mrchem
