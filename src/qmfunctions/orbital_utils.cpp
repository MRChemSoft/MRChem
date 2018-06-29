#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "utils/RRMaximizer.h"

#include "qmfunction_utils.h"
#include "Orbital.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {

/****************************************
 * Orbital related standalone functions *
 ****************************************/

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Complicated by the fact that both bra and ket can be interpreted
 *  as complex conjugate versions of themselves. Notice that the <bra|
 *  position is already complex conjugated.
 *
 */
ComplexDouble dot(Orbital bra, Orbital ket) {
    if ((bra.spin() == SPIN::Alpha) and (ket.spin() == SPIN::Beta)) return 0.0;
    if ((bra.spin() == SPIN::Beta) and (ket.spin() == SPIN::Alpha)) return 0.0;
    
    double bra_conj(1.0), ket_conj(1.0);
    if (bra.conjugate()) bra_conj = -1.0;
    if (ket.conjugate()) ket_conj = -1.0;

    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), ket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), ket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), ket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), ket.imag());

    double real_part = rr + bra_conj*ket_conj*ii;
    double imag_part = ket_conj*ri - bra_conj*ir;
    return ComplexDouble(real_part, imag_part);
}

/** @brief Compare spin and occupancy of two orbitals
 *
 *  Returns true if orbital parameters are the same.
 *
 */
bool compare(const Orbital &orb_a, const Orbital &orb_b) {
    bool comp = true;
    if (compare_occ(orb_a, orb_b) < 0) {
        MSG_WARN("Different occupancy");
        comp = false;
    }
    if (compare_spin(orb_a, orb_b) < 0) {
        MSG_WARN("Different spin");
        comp = false;
    }
    return comp;
}

/** @brief Compare occupancy of two orbitals
 *
 *  Returns the common occupancy if they match, -1 if they differ.
 *
 */
int compare_occ(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.occ() == orb_b.occ()) comp = orb_a.occ();
    return comp;
}

/** @brief Compare spin of two orbitals
 *
 *  Returns the common spin if they match, -1 if they differ.
 *
 */
int compare_spin(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.spin() == orb_b.spin()) comp = orb_a.spin();
    return comp;
}

/** @brief out = a*inp_a + b*inp_b
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves. */
Orbital add(ComplexDouble a, Orbital inp_a, ComplexDouble b, Orbital inp_b, double prec) {
    ComplexVector coefs(2);
    coefs(0) = a;
    coefs(1) = b;

    OrbitalVector orbs;
    orbs.push_back(inp_a);
    orbs.push_back(inp_b);

    return linear_combination(coefs, orbs, prec);
}

/** @brief out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector pairwise_add(ComplexDouble a, OrbitalVector &inp_a,
                  ComplexDouble b, OrbitalVector &inp_b,
                  double prec) {
    if (inp_a.size() != inp_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) {
        Orbital out_i = add(a, inp_a[i], b, inp_b[i]);
        out.push_back(out_i);
    }
    return out;
}

/** @brief out = inp_a * inp_b
  *
  * Complicated by the fact that both inputs can be interpreted
  * as complex conjugate versions of themselves.
  *
  */
Orbital multiply(Orbital inp_a, Orbital inp_b, double prec) {
    int occ = compare_occ(inp_a, inp_b);
    int spin = compare_spin(inp_a, inp_b);
    Orbital out(spin, occ);
    qmfunction::multiply(prec, out, 1.0, inp_a, inp_a.conjugate(), inp_b, inp_b.conjugate());
    return out;
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ...
  *
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves.
  *
  */
Orbital linear_combination(const ComplexVector &c, OrbitalVector &inp, double prec) {
    if (c.size() != inp.size()) MSG_ERROR("Size mismatch");
    double thrs = mrcpp::MachineZero;

    Orbital out;
    QMFunctionVector inp_qmf;
    std::vector<bool> conj;
    // set output spin from first contributing input
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue;
        out = inp[i].paramCopy();
        break;
    }
    // check that contributing input spins must be equal and fill bool vector
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue; //
        if (out.spin() != inp[i].spin()) {
            // empty orbitals with wrong spin can occur
            if (inp[i].hasReal() or inp[i].hasImag()) MSG_FATAL("Mixing spins");
        }
        conj[i] = inp[i].conjugate();
        inp_qmf.push_back(inp[i]); //LUCA: the best would be to have a conversion from OrbitalVector to QMFunctionVector
    }
    qmfunction::linear_combination(c, inp_qmf, conj, out, prec);
    return out;
}

/** @brief Orbital transformation out_vec = U*inp_vec
 *
 * The transformation matrix is not necessarily square.
 *
 */
OrbitalVector linear_combination(const ComplexMatrix &U, OrbitalVector &inp, double prec) {
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    OrbitalVector out;
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i);
        Orbital out_i = linear_combination(c, inp, prec);
        out.push_back(out_i);
    }
    return out;
}

/** @brief Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector deep_copy(OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].deepCopy();
        out.push_back(out_i);
    }
    return out;
}

/** @brief Parameter copy
 *
 * New orbitals are constructed as parameter copies of the input set.
 *
 */
OrbitalVector param_copy(const OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].paramCopy();
        out.push_back(out_i);
    }
    return out;
}

/** @brief join two vectors
 *
 * The orbitals of the input vector are appended to
 * (*this) vector, the ownership is transferred. Leaves
 * the input vector empty.
 *
 */
OrbitalVector join(OrbitalVector &inp_a, OrbitalVector &inp_b) {
    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) out.push_back(inp_a[i]);
    for (int i = 0; i < inp_b.size(); i++) out.push_back(inp_b[i]);
    inp_a.clear();
    inp_b.clear();
    return out;
}

/** @brief Separate vector in two parts
 *
 * All orbitals of a particular spin are collected in a new vector
 * and returned. These orbitals are removed from (*this) vector,
 * and the ownership is transferred.
 *
 */
OrbitalVector split(OrbitalVector &inp, int spin) {
    OrbitalVector out;
    OrbitalVector tmp;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &orb_i = inp[i];
        if (orb_i.spin() == spin) {
            out.push_back(orb_i);
        } else {
            tmp.push_back(orb_i);
        }
    }
    inp.clear();
    inp = tmp;
    return out;
}

/** @brief Write orbitals to disk
 *
 * @param Phi: orbitals to save
 * @param file: file name prefix
 * @param n_orbs: number of orbitals to save
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Produces separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals in the
 * vector are saved.
 */
void save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs) {
    if (n_orbs < 0) n_orbs = Phi.size();
    if (n_orbs > Phi.size()) MSG_ERROR("Index out of bounds");
    for (int i = 0; i < n_orbs; i++) {
        std::stringstream orbname;
        orbname << file << "_" << i;
        Phi[i].saveOrbital(orbname.str());
    }
}

/** @brief Read orbitals from disk
 *
 * @param file: file name prefix
 * @param n_orbs: number of orbitals to read
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Reads separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals matching
 * the prefix name will be read.
 */
OrbitalVector load_orbitals(const std::string &file, int n_orbs) {
    OrbitalVector Phi;
    for (int i = 0; true; i++) {
        if (n_orbs > 0 and i >= n_orbs) break;
        std::stringstream orbname;
        orbname << file << "_" << i;
        Orbital phi_i;
        phi_i.loadOrbital(orbname.str());
        if (phi_i.hasReal() or phi_i.hasImag()) {
            Phi.push_back(phi_i);
        } else {
            break;
        }
    }
    return Phi;
}

/** @brief Frees each orbital in the vector
 *
 * Leaves an empty vector. Orbitals are freed.
 *
 */
void free(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) vec[i].free();
    vec.clear();
}

/** @brief Normalize all orbitals in the set */
void normalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].normalize();
    }
}

/** @brief Gram-Schmidt orthogonalize orbitals within the set */
void orthogonalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        Orbital &orb_i = vec[i];
        for (int j = 0; j < i; j++) {
            Orbital &orb_j = vec[j];
            orb_i.orthogonalize(orb_j);
        }
    }
}

/** @brief Orthogonalize the out orbital against all orbitals in inp */
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].orthogonalize(inp);
    }
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket) {
    return orbital::calc_overlap_matrix(braket, braket);
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &bra, OrbitalVector &ket) {
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix S(Ni, Nj);
    S.setZero();

    for (int i = 0; i < Ni; i++) {
        for (int j = 0; j < Nj; j++) {
            S(i,j) = orbital::dot(bra[i], ket[j]);
        }
    }
    return S;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi) {
    Timer timer;
    printout(1, "Calculating Löwdin orthonormalization matrix      ");

    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0/2.0);

    timer.stop();
    println(1, timer.getWallTime());
    return S_m12;
}

/** @brief Minimize the spatial extension of orbitals, by orbital rotation
 *
 * @param Phi: orbitals to localize
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 * which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix localize(double prec, OrbitalVector &Phi) {
    Printer::printHeader(0, "Localizing orbitals");
    Timer timer;

    ComplexMatrix U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rr_t;
        RRMaximizer rr(prec, Phi);
        n_it = rr.maximize();
        rr_t.stop();
        Printer::printDouble(0, "Computing Foster-Boys matrix", rr_t.getWallTime(), 5);

        if (n_it > 0) {
            println(0, " Converged after iteration   " << std::setw(30) << n_it);
            U = rr.getTotalU().transpose().cast<ComplexDouble>();
        } else {
            println(0, " Foster-Boys localization did not converge!");
        }
    } else {
        println(0, " Cannot localize less than two orbitals");
    }

    if (n_it <= 0) {
        Timer orth_t;
        U = orbital::calc_lowdin_matrix(Phi);
        orth_t.stop();
        Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);
    }

    Timer rot_t;
    OrbitalVector Psi = orbital::linear_combination(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * @param Phi: orbitals to rotate
 * @param F: Fock matrix to diagonalize
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place and Fock matrix is diagonalized in place.
 * The transformation matrix is returned.
 */
ComplexMatrix diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F) {
    Printer::printHeader(0, "Digonalizing Fock matrix");
    Timer timer;

    Timer orth_t;
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    F = S_m12.transpose()*F*S_m12;
    orth_t.stop();
    Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) math_utils::diagonalize_block(F, U, 0,       np);
    if (na > 0) math_utils::diagonalize_block(F, U, np,      na);
    if (nb > 0) math_utils::diagonalize_block(F, U, np + na, nb);
    U = U * S_m12;
    diag_t.stop();
    Printer::printDouble(0, "Diagonalizing matrix", diag_t.getWallTime(), 5);

    Timer rot_t;
    OrbitalVector Psi = orbital::linear_combination(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the Löwdin orthonormalization
 *
 * @param Phi: orbitals to orthonormalize
 *
 * Orthonormalizes the orbitals by multiplication of the Löwdin matrix S^(-1/2).
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orthonormalize(double prec, OrbitalVector &Phi) {
    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi);
    OrbitalVector Psi = orbital::linear_combination(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    return U;
}

/** @brief Returns the number of occupied orbitals */
int size_occupied(const OrbitalVector &vec) {
    int nOcc = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() > 0) nOcc++;
    }
    return nOcc;
}

/** @brief Returns the number of empty orbitals */
int size_empty(const OrbitalVector &vec) {
    int nEmpty = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 0) nEmpty++;
    }
    return nEmpty;
}

/** @brief Returns the number of singly occupied orbitals */
int size_singly(const OrbitalVector &vec) {
    int nSingly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nSingly++;
    }
    return nSingly;
}

/** @brief Returns the number of doubly occupied orbitals */
int size_doubly(const OrbitalVector &vec) {
    int nDoubly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nDoubly++;
    }
    return nDoubly;
}

/** @brief Returns the number of paired orbitals */
int size_paired(const OrbitalVector &vec) {
    int nPaired = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Paired) nPaired++;
    }
    return nPaired;
}

/** @brief Returns the number of alpha orbitals */
int size_alpha(const OrbitalVector &vec) {
    int nAlpha = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Alpha) nAlpha++;
    }
    return nAlpha;
}

/** @brief Returns the number of beta orbitals */
int size_beta(const OrbitalVector &vec) {
    int nBeta = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Beta) nBeta++;
    }
    return nBeta;
}

/** @brief Returns the spin multiplicity of the vector */
int get_multiplicity(const OrbitalVector &vec) {
    int nAlpha = get_electron_number(vec, SPIN::Alpha);
    int nBeta = get_electron_number(vec, SPIN::Beta);
    int S = std::abs(nAlpha - nBeta);
    return S + 1;
}

/** @brief Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
int get_electron_number(const OrbitalVector &vec, int spin) {
    int nElectrons = 0;
    for (int i = 0; i < vec.size(); i++) {
        const Orbital &orb = vec[i];
        if (spin == SPIN::Paired) {
            nElectrons += orb.occ();
        } else if (spin == SPIN::Alpha) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Alpha) {
                nElectrons += 1;
            }
        } else if (spin == SPIN::Beta) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Beta) {
                nElectrons += 1;
            }
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** @brief Returns a vector containing the orbital errors */
DoubleVector get_errors(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector errors = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        errors(i) = vec[i].error();
    }
    return errors;
}

/** @brief Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_errors(OrbitalVector &vec, const DoubleVector &errors) {
    if (vec.size() != errors.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setError(errors(i));
    }
}

/** @brief Returns a vector containing the orbital spins */
IntVector get_spins(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector spins = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        spins(i) = vec[i].spin();
    }
    return spins;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_spins(OrbitalVector &vec, const IntVector &spins) {
    if (vec.size() != spins.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setSpin(spins(i));
    }
}

/** @brief Returns a vector containing the orbital occupancies */
IntVector get_occupancies(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        occ(i) = vec[i].occ();
    }
    return occ;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_occupancies(OrbitalVector &vec, const IntVector &occ) {
    if (vec.size() != occ.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setOcc(occ(i));
    }
}

/** @brief Returns a vector containing the orbital square norms */
DoubleVector get_squared_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].squaredNorm();
    }
    return norms;
}

/** @brief Returns a vector containing the orbital norms */
DoubleVector get_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].norm();
    }
    return norms;
}

void print(const OrbitalVector &vec) {
    Printer::setScientific();
    printout(0, "============================================================\n");
    printout(0, " OrbitalVector:");
    printout(0, std::setw(4) << vec.size()          << " orbitals  ");
    printout(0, std::setw(4) << size_occupied(vec)  << " occupied  ");
    printout(0, std::setw(4) << get_electron_number(vec) << " electrons\n");
    printout(0, "------------------------------------------------------------\n");
    printout(0, "   n  RankID           Norm          Spin Occ      Error    \n");
    printout(0, "------------------------------------------------------------\n");
    for (int i = 0; i < vec.size(); i++) {
        println(0, std::setw(4) << i << vec[i]);
    }
    printout(0, "============================================================\n\n\n");
}

} //namespace orbital

} //namespace mrchem
