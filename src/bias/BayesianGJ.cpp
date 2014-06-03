/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS BAYESIANGJ
/*
Calculate a Bayesian Score to use with the Chemical Shifts CV \ref CS2BACKBONE.

The functional form of this bias is
\f[
C=k_BT \sum_{i=1}^{N_{arg}} \log{ \left[ \frac{\pi}{\sqrt{2} \sigma} (x_i+2\sigma^2) \right]} + k_BT\log{\sigma}
\f]

where sigma is an uncertainty parameter,
sampled by a MC algorithm in the bounded interval defined by SIGMA_MIN and SIGMA_MAX.
The initial value is set at SIGMA0. The MC move is a random displacement
of maximum value equal to DSIGMA.



\par Examples
The following input tells plumed to use all the HA chemical shifts with the Bayesian Score and
to print the values of the uncertainty parameter, MC acceptance, and Bayesian score.
\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs:  CS2BACKBONE ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=0.0 NRES=13 ENSEMBLE COMPONENTS
csb: BAYESIANGJ ARG=cs.ha SIGMA0=1.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.1 KBT=2.494 NDATA=13 MC_STEPS=1 MC_STRIDE=1 MC_SEED=1234
PRINT ARG=csb.sigma,csb.accept,csb.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref PRINT).


*/
//+ENDPLUMEDOC


class BayesianGJ : public Bias
{
  const double sqrt2pi = 2.506628274631001;
  double sigma_;
  double sigma0_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  double temp_;
  int ndata_;
  int MCsteps_;
  int MCstride_;
  unsigned int MCseed_;
  unsigned int MCaccept_;
  Value* valueBias;
  Value* valueSigma;
  Value* valueAccept;
  Value* valueKappa;
 
  void doMonteCarlo();
  double getEnergy(double sigma);
  
public:
  BayesianGJ(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianGJ,"BAYESIANGJ")

void BayesianGJ::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0",   "1.0",    "initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","0.00001","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","10.0",   "maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA",   "0.1",    "maximum MC move of the uncertainty parameter");
  keys.add("compulsory","KBT",      "2.494",  "temperature of the system in energy units");
  keys.add("compulsory","NDATA",              "number of experimental data points");
  keys.add("compulsory","MC_STEPS", "1",      "number of MC steps");
  keys.add("compulsory","MC_STRIDE","1",      "MC stride");
  keys.add("compulsory","MC_SEED",  "1234",   "MC seed");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("kappa", "default","intensity of the harmonic restraint");
  keys.addOutputComponent("accept","default","MC acceptance");
}

BayesianGJ::BayesianGJ(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
sigma0_(1.0), sigma_min_(0.00001), sigma_max_(10.0), Dsigma_(0.1), temp_(2.494),
MCsteps_(1), MCstride_(1), MCseed_(1234), MCaccept_(0)
{
  parse("SIGMA0",   sigma0_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",   Dsigma_);
  parse("KBT",      temp_);
  parse("NDATA",    ndata_);
  parse("MC_STEPS", MCsteps_);
  parse("MC_STRIDE",MCstride_);
  parse("MC_SEED",  MCseed_);
  checkRead();

  log.printf("  initial value of uncertainty %f\n",sigma0_);
  log.printf("  minimum value of uncertainty %f\n",sigma_min_);
  log.printf("  maximum value of uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of the uncertainty parameter %f\n",Dsigma_);
  log.printf("  temperature of the system %f\n",temp_);
  log.printf("  number of experimental data points %d\n",ndata_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);
  log.printf("  MC seed %d\n", MCseed_);    

  addComponent("bias");   componentIsNotPeriodic("bias");
  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  addComponent("kappa");  componentIsNotPeriodic("kappa");
  valueBias=getPntrToComponent("bias");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");
  valueKappa=getPntrToComponent("kappa");


  // initialize sigma_
  sigma_ = sigma0_;
  // initialize random seed
  srand (MCseed_);
}

double BayesianGJ::getEnergy(double sigma){
 double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    ene += getArgument(i);
  }
  ene = ene/2.0/sigma/sigma + std::log(sigma) + static_cast<double>(ndata_)*std::log(sigma*sqrt2pi);
  return temp_ * ene;
}

void BayesianGJ::doMonteCarlo(){
 
 // store old energy
 double old_energy = getEnergy(sigma_);
 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move
  double r = static_cast<double>(rand()) / RAND_MAX;
  double ds = -Dsigma_ + r * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}
  // calculate new energy
  double new_energy = getEnergy(new_sigma);
  // accept or reject
  double s = static_cast<double>(rand()) / RAND_MAX;
  double delta = exp(-(new_energy-old_energy)/temp_);
  if(s<delta){
   old_energy = new_energy;
   sigma_ = new_sigma;
   MCaccept_++;
  }
 }
}

void BayesianGJ::calculate(){
  
  // do MC stuff at the right time step
  long int step = getStep();
  if((step+1)%MCstride_==0){doMonteCarlo();}

  // cycle on arguments 
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // increment energy
    ene += getArgument(i); 
    // set derivatives
    setOutputForce(i, -temp_/2.0/sigma_/sigma_);
  };
  ene = ene/2.0/sigma_/sigma_ + std::log(sigma_) + static_cast<double>(ndata_)*std::log(sigma_*sqrt2pi);

  // set value of the bias
  valueBias->set(temp_*ene);
  // set value of uncertainty
  valueSigma->set(sigma_);
  // and of the harmonic restraint on coordinates
  valueKappa->set(temp_/2.0/sigma_/sigma_);
  // calculate acceptance
  double MCtrials = std::floor(static_cast<double>(step+1) / static_cast<double>(MCstride_));
  double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCsteps_) / MCtrials; 
  // set value of acceptance 
  valueAccept->set(accept);
}


}
}


