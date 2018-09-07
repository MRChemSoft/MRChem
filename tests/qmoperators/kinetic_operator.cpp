#include "catch.hpp"

#include "MRCPP/MWOperators"

#include "mrchem.h"

#include "HarmonicOscillatorFunction.h"
#include "KineticOperator.h"
#include "Orbital.h"
#include "qmfunctions/orbital_utils.h"

using namespace mrchem;
using namespace orbital;

namespace kinetic_operator {

TEST_CASE("KineticOperator", "[kinetic_operator]") {
    const double prec = 1.0e-3;
    const double thrs = prec*prec;

    int nFuncs = 3;
    OrbitalVector Phi;
    for (int n = 0; n < nFuncs; n++) {
        int nu[3] = {n,0,0};
        HarmonicOscillatorFunction f(nu);

        Orbital phi(SPIN::Paired);
        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), f);
        Phi.push_back(phi);
    }


    // reference values for harmonic oscillator eigenfunctions
    DoubleVector E_K(nFuncs);
    for (int i = 0; i < nFuncs; i++) {
        //energy = (nu + 1/2)
        double E_x = (i + 0.5);
        double E_y = (0 + 0.5);
        double E_z = (0 + 0.5);

        //virial theorem: <E_K> = <E_P> = E/2
        E_K(i) = 0.5*(E_x + E_y + E_z);
    }

    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    KineticOperator T(D);
    T.setup(prec);
    SECTION("apply") {
        Orbital Tphi_0 = T(Phi[0]);
        ComplexDouble T_00 = orbital::dot(Phi[0], Tphi_0);
        REQUIRE( T_00.real() == Approx(E_K(0)) );
        REQUIRE( T_00.imag() < thrs );
        Tphi_0.free();
    }
    SECTION("vector apply") {
        OrbitalVector TPhi = T(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble T_ii = orbital::dot(Phi[i], TPhi[i]);
            REQUIRE( T_ii.real() == Approx(E_K(i)) );
            REQUIRE( T_ii.imag() < thrs );
        }
        free(TPhi);
    }
    SECTION("expectation value") {
        ComplexDouble T_00 = T(Phi[0], Phi[0]);
        REQUIRE( T_00.real() == Approx(E_K(0)) );
        REQUIRE( T_00.imag() < thrs );
    }
    SECTION("expectation matrix ") {
        ComplexMatrix t = T(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            REQUIRE( t(i,i).real() == Approx(E_K(i)) );
            REQUIRE( t(i,i).imag() < thrs );
        }
    }
    T.clear();
    free(Phi);
}

} //namespace kinetic_operator
