#pragma once
#include "chemistry/Cavity.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "chemistry/chemistry_fwd.h"
#include "chemistry/Nucleus.h"

using namespace mrcpp;

namespace mrchem {


//start with a cavity initialized with a geometry and a standard gaussian rho with A*exp(-B*r^2) with A = (B/pi)^(3/2) include the orbitals later

class ReactionPotential final : public QMPotential{
public:

  ReactionPotential(mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nucs, OrbitalVector *Phi);
  ~ReactionPotential(){}
  QMFunction &getPotential() {return this->V_eff_func;}
  QMFunction &getgamma() {return this->gamma_func;}
  void do_setup(double prec) {this->setup(prec);}

  friend class ReactionOperator;

protected:
  void clear(){clearApplyPrec();}


private:

  Cavity *cavity;
  Nuclei nuclei;
  OrbitalVector *orbitals;
  mrcpp::PoissonOperator *poisson;
  mrcpp::DerivativeOperator<3> *derivative;

  QMFunction rho_func;
  QMFunction cavity_func;
  QMFunction inv_eps_func;
  QMFunction rho_eff_func;
  QMFunction gamma_func;
  QMFunction V_n_func;
  mrcpp::FunctionTreeVector<3> d_cavity;
  QMFunction V_eff_func;

  void setup(double prec);

  void setup_eps();

  void calc_rho_eff();
  void calc_gamma();
};



} //namespace mrchem