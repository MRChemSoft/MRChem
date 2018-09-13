/**
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#pragma once

#include <string>
#include <vector>

#include "Nucleus.h"
#include "properties/SCFEnergy.h"

/** @class Molecule
 *
 * @brief Collection of properties related to a molecule.
 *
 * A molecule is basically a collection of nuclei and a collection of properties.
 * Mainly used for bookkeeping and printing of properties.
 *
 */

namespace mrchem {

class DipoleMoment;
class QuadrupoleMoment;
class GeometryDerivatives;
class Magnetizability;
class NMRShielding;
class HyperFineCoupling;
class SpinSpinCoupling;
class Polarizability;
class OpticalRotation;

class Molecule final {
public:
    Molecule(const Nuclei &nucs, int c = 0, int m = 1);
    Molecule(const std::string &coord_file, int c = 0, int m = 1);
    Molecule(const std::vector<std::string> &coord_str, int c = 0, int m = 1);
    ~Molecule();

    int getCharge() const { return this->charge; }
    int getMultiplicity() const { return this->multiplicity; }
    int getNNuclei() const { return this->nuclei.size(); }
    int getNElectrons() const;

    const double *getCenterOfMass() const { return this->COM; }
    const double *getCenterOfCharge() const { return this->COC; }

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }
    Nucleus &getNucleus(int i) { return this->nuclei[i]; }
    const Nucleus &getNucleus(int i) const { return this->nuclei[i]; }

    void printGeometry() const;
    void printProperties() const;

    void initSCFEnergy();
    void initDipoleMoment();
    void initGeometryDerivatives();
    void initMagnetizability();
    void initQuadrupoleMoment();
    void initNMRShielding(int k);
    void initHyperFineCoupling(int k);
    void initSpinSpinCoupling(int k, int l);
    void initPolarizability(double omega);
    void initOpticalRotation(double omega);

    SCFEnergy& getSCFEnergy();
    DipoleMoment& getDipoleMoment();
    QuadrupoleMoment& getQuadrupoleMoment();
    GeometryDerivatives& getGeometryDerivatives();
    Magnetizability& getMagnetizability();
    NMRShielding& getNMRShielding(int k);
    HyperFineCoupling& getHyperFineCoupling(int k);
    SpinSpinCoupling& getSpinSpinCoupling(int k, int l);
    Polarizability& getPolarizability(double omega);
    OpticalRotation& getOpticalRotation(double omega);

protected:
    int charge;
    int multiplicity;
    Nuclei nuclei;

    // Properties
    double COM[3];
    double COC[3];
    SCFEnergy *energy;
    DipoleMoment *dipole;
    QuadrupoleMoment *quadrupole;
    GeometryDerivatives *geomderiv;
    Magnetizability *magnetizability;
    NMRShielding **nmr;
    HyperFineCoupling **hfcc;
    SpinSpinCoupling ***sscc;
    std::vector<Polarizability *> polarizability;
    std::vector<OpticalRotation *> optical_rotation;

    void calcCenterOfMass();
    void calcCenterOfCharge();

    void allocNuclearProperties();
    void freeNuclearProperties();

    void clearSCFEnergy();
    void clearDipoleMoment();
    void clearQuadrupoleMoment();
    void clearGeometryDerivatives();
    void clearMagnetizability();
    void clearNMRShielding(int k);
    void clearHyperFineCoupling(int k);
    void clearSpinSpinCoupling(int k, int l);
    void clearPolarizability();
    void clearOpticalRotation();

    void readCoordinateFile(const std::string &file);
    void readCoordinateString(const std::vector<std::string> &coord_str);
};

} //namespace mrchem
