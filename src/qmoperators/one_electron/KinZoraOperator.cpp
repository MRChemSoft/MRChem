#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "KinZoraOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Expectation value matrix
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix KinZoraOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T_x = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_y = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_z = ComplexMatrix::Zero(Ni, Nj);
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        OrbitalVector dBra = p_x(bra);
        OrbitalVector dKet = p_x(ket);
        OrbitalVector vzdKet = vz(dKet);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(vzdKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(vzdKet);
        T_x = orbital::calc_overlap_matrix(dBra, vzdKet);
        mrcpp::print::tree(2, "<i|p[x]*kappa*p[x]|j>", nNodes, sNodes, timer.elapsed());
    }
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        OrbitalVector dBra = p_y(bra);
        OrbitalVector dKet = p_y(ket);
        OrbitalVector vzdKet = vz(dKet);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(vzdKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(vzdKet);
        T_y = orbital::calc_overlap_matrix(dBra, vzdKet);
        mrcpp::print::tree(2, "<i|p[y]*kappa*p[y]|j>", nNodes, sNodes, timer.elapsed());
    }
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        OrbitalVector dBra = p_z(bra);
        OrbitalVector dKet = p_z(ket);
        OrbitalVector vzdKet = vz(dKet);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(vzdKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(vzdKet);
        T_z = orbital::calc_overlap_matrix(dBra, vzdKet);
        mrcpp::print::tree(2, "<i|p[z]*kappa[z]|j>", nNodes, sNodes, timer.elapsed());
    }
    return 0.5 * (T_x + T_y + T_z);
}

/** @brief Expectation value (dagger version)
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * NOT IMPLEMENTED
 */
ComplexMatrix KinZoraOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
