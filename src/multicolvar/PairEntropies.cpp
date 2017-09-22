/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR PAIRENTROPIES
/*
Calculate the pair entropy of atom i using the expression:

\f[
s_i=-2\pi\rho k_B \int\limits_0^{r_{\mathrm{max}}} \left [ g(r) \ln g(r) - g(r) + 1 \right ] r^2 dr .
\f]

where \f$ g(r) $\f is the pair distribution function and \f$ r_{\mathrm{max}} $\f is a cutoff in the integration (MAXR).
For the integration the interval from 0 to  \f$ r_{\mathrm{max}} $\f is partitioned in NHIST equal intervals. 
To make the calculation of \f$ g(r) $\f differentiable, the following function is used:
\f[
g(r) = \frac{1}{4 \pi \rho r^2} \sum\limits_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(r-r_{ij})^2/(2\sigma^2)} ,
\f]
where \f$ \rho $\f is the density and \f$ sigma $\f is a broadening parameter (SIGMA).  

\par Example)

The following input tells plumed to calculate the per atom per entropy of atoms 1-250 with themselves.
The mean pair entropy is then calculated.
\verbatim
PAIRENTROPIES ...
 LABEL=s2
 SPECIES=1-250
 MAXR=0.65
 SIGMA=0.025
 NHIST=60
 MEAN
... PAIRENTROPIES
\endverbatim

*/
//+ENDPLUMEDOC


class PairEntropies : public MultiColvarBase {
private:
  double rcut2, rcut;
  double invSqrt2piSigma, sigmaSqr2, sigmaSqr;
  double maxr, sigma;
  unsigned nhist;
  double density_given;
  bool local_density;
  double deltar;
  unsigned deltaBin;
  // Integration routine
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double&der)const;
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropies(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(PairEntropies,"PAIRENTROPIES")

void PairEntropies::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","300","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.addFlag("LOCAL_DENSITY",false,"Use the local density to normalize g(r). If not specified, N/V is used");

  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
  keys.add("optional","INTEGRAND_FILE","the file on which to write the integrand");
}

PairEntropies::PairEntropies(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
  parse("MAXR",maxr);
  log.printf("Integration in the interval from 0. to %f . \n", maxr );
  parse("NHIST",nhist);
  log.printf("The interval is partitioned in %u equal parts and the integration is perfromed with the trapezoid rule. \n", nhist );
  parse("SIGMA",sigma);
  log.printf("The pair distribution function is calculated with a Gaussian kernel with deviation %f . \n", sigma);
  density_given = -1;
  parse("DENSITY",density_given);
  parseFlag("LOCAL_DENSITY",local_density);
  if (density_given>0) log.printf("The g(r) will be normalized with a density %f . \n", density_given);
  else if (local_density) log.printf("The g(r) will be normalized with the local density. Derivatives might be wrong. \n");
  else log.printf("The g(r) will be normalized with a density N/V . \n");

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();

  // Define heavily used constants
  double sqrt2piSigma = std::sqrt(2*pi)*sigma;
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/nhist;
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(3*sigma/deltar); //3*sigma is 99.7 %

  // Set the link cell cutoff
  setLinkCellCutoff( maxr + 3*sigma );
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);
  rcut = std::sqrt(rcut2);
  log.printf("Setting cut off to %f \n ", maxr + 3*sigma );
}

double PairEntropies::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   double dfunc, d2;
   Vector value;
   vector<double> gofr(nhist);
   vector<double> logGofr(nhist);
   Matrix<Vector> gofrPrime(nhist,getNumberOfAtoms());
   vector<Vector> deriv(getNumberOfAtoms());
   vector<Tensor> gofrVirial(nhist);
   Tensor virial;
   // Construct g(r)
   int countNeigh=0;
   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      Vector& distance=myatoms.getPosition(i);  
      if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor(distanceModulo/deltar);
           int minBin, maxBin;
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int j=minBin;j<maxBin+1;j+=1) {   
             double x=deltar*(j+0.5);
             gofr[j] += kernel(x-distanceModulo, dfunc);
             if (!doNotCalculateDerivatives()) {
               value = dfunc * distance_versor;
               gofrPrime[j][0] += value;
               gofrPrime[j][i] -= value;
               Tensor vv(value, distance);
               gofrVirial[j] += vv;
             }
	   } 
           countNeigh += 1;
      }
   }
   // Normalize g(r)
   double volume=getBox().determinant(); 
   double density;
   if (density_given>0) density=density_given;
   else if (local_density) density= (double) countNeigh / ( (4./3.)*pi*rcut*rcut*rcut);
   else density=getNumberOfAtoms()/volume;
   //log.printf("rcut %f \n", rcut);
   //log.printf("countNeigh %d \n", countNeigh);
   //log.printf("density %f \n", density);
   for(unsigned i=0;i<nhist;++i){
     double x=deltar*(i+0.5);
     double normConstant = 4*pi*density*x*x;
     gofr[i] /= normConstant;
     if (!doNotCalculateDerivatives()) {
       gofrVirial[i] /= normConstant;
       for(unsigned j=0;j<myatoms.getNumberOfAtoms();++j){
         gofrPrime[i][j] /= normConstant;
       }
     }
   }
   // Construct integrand
   vector<double> integrand(nhist);
   for(unsigned i=0;i<nhist;++i){
     double x=deltar*(i+0.5);
     logGofr[i] = std::log(gofr[i]);
     if (gofr[i]<1.e-10) {
       integrand[i] = x*x;
     } else {
       integrand[i] = (gofr[i]*logGofr[i]-gofr[i]+1)*x*x;
     }
   }
   // Integrate to obtain pair entropy;
   double entropy = -2*pi*density*integrate(integrand,deltar); 
   if (!doNotCalculateDerivatives()) {
     // Construct integrand and integrate derivatives
     for(unsigned i=0;i<myatoms.getNumberOfAtoms();++i) {
       vector<Vector> integrandDerivatives(nhist);
       for(unsigned j=0;j<nhist;++j){
         double x=deltar*(j+0.5);
         if (gofr[j]>1.e-10) {
           integrandDerivatives[j] = gofrPrime[j][i]*logGofr[j]*x*x;
         }
       }
       // Integrate
       deriv[i] = -2*pi*density*integrate(integrandDerivatives,deltar);
     }
     // Virial of positions
     // Construct virial integrand
     vector<Tensor> integrandVirial(nhist);
     for(unsigned i=0;i<nhist;++i){
       double x=deltar*(i+0.5);
       if (gofr[i]>1.e-10) {
         integrandVirial[i] = gofrVirial[i]*logGofr[i]*x*x;
       }
     }
     // Integrate virial
     virial = -2*pi*density*integrate(integrandVirial,deltar);
     // Virial of volume
     // Construct virial integrand
     vector<double> integrandVirialVolume(nhist);
     for(unsigned i=0;i<nhist;i+=1) {   
       double x=deltar*(i+0.5);
       integrandVirialVolume[i] = (-gofr[i]+1)*x*x;
     }
     // Integrate virial
     virial += -2*pi*density*integrate(integrandVirialVolume,deltar)*Tensor::identity();
   }
   // Assign derivatives
   for(unsigned i=0;i<myatoms.getNumberOfAtoms();++i) addAtomDerivatives( 1, i, deriv[i], myatoms );
   // Assign virial
   myatoms.addBoxDerivatives( 1, virial );
   return entropy;
}

double PairEntropies::kernel(double distance,double&der)const{
  // Gaussian function and derivative
  double result = invSqrt2piSigma*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

double PairEntropies::integrate(vector<double> integrand, double delta)const{
  // Trapezoid rule
  double result = 0.;
  for(unsigned i=1;i<(integrand.size()-1);++i){
    result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Vector PairEntropies::integrate(vector<Vector> integrand, double delta)const{
  // Trapezoid rule
  Vector result;
  for(unsigned i=1;i<(integrand.size()-1);++i){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Tensor PairEntropies::integrate(vector<Tensor> integrand, double delta)const{
  // Trapezoid rule
  Tensor result;
  for(unsigned i=1;i<(integrand.size()-1);++i){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}


}
}

