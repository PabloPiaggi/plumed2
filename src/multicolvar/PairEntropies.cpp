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
#include "tools/Communicator.h"

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
  bool local_density, one_body, no_two_body, doAverageGofr;
  unsigned averageGofrTau;
  double deltar;
  unsigned deltaBin;
  double temperature, mass, deBroglie3;
  std::vector<double> vectorX, vectorX2;
  // Integration routine
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double&der)const;
  // Output of g(r)
  bool doOutputGofr;
  void outputGofr(int index, int nat, int step, int stride, vector<double> gofr)const;
  unsigned outputStride;
  // Average g(r)
  mutable Matrix<double> avgGofr;
  mutable double avgDensity;
  mutable vector<unsigned> iteration;
  // Cut at first peak
  bool doCutAtFirstPeak;
  double cutAtFirstPeakValue;
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropies(const ActionOptions&);
  ~PairEntropies();
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
  mutable OFile gofrOfile;
};

PLUMED_REGISTER_ACTION(PairEntropies,"PAIRENTROPIES")

void PairEntropies::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","300","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("optional","TEMPERATURE","Temperature in Kelvin. It is compulsory when keyword ONE_BODY is used");
  keys.add("optional","MASS","Mass in g/mol. It is compulsory when keyword ONE_BODY is used");
  keys.addFlag("LOCAL_DENSITY",false,"Use the local density to normalize g(r). If not specified, N/V is used");
  keys.addFlag("ONE_BODY",false,"Add the one body term (S1 = 5/2 - ln(dens*deBroglie^3) ) to the entropy");
  keys.addFlag("NO_TWO_BODY",false,"Remove the two-body term. Only the one-body term is kept. This flag can only be used along with the ONE_BODY flag.");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r) to file.");
  keys.addFlag("CUT_AT_FIRST_PEAK",false,"Let the g(r) go to one after the first peak. If CUT_AT_FIRST_PEAK_VALUE is given, then the first peak is considered to be located at the first point in which the g(r) attains the given value. If CUT_AT_FIRST_PEAK_VALUE is not given, the position of the first peak is determined as the first zero of the second derivative of g(r).");
  keys.add("optional","CUT_AT_FIRST_PEAK_VALUE","Value of g(r) after which g(r) will be 1. It is optional when the keyword CUT_AT_FIRST_PEAK is used");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.addFlag("AVERAGE_GOFR",false,"Average g(r) over time.");
  keys.add("optional","AVERAGE_GOFR_TAU","Characteristic length of a window in which to average the g(r). It is in units of iterations and should be an integer. Zero corresponds to an normal average (infinite window).");
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
  log.printf("The interval is partitioned in %d equal parts and the integration is perfromed with the trapezoid rule. \n", nhist );
  parse("SIGMA",sigma);
  log.printf("The pair distribution function is calculated with a Gaussian kernel with deviation %f . \n", sigma);
  density_given = -1;
  parse("DENSITY",density_given);
  parseFlag("LOCAL_DENSITY",local_density);
  if (density_given>0) log.printf("The g(r) will be normalized with a density %f . \n", density_given);
  else if (local_density) log.printf("The g(r) will be normalized with the local density. Derivatives might be wrong. \n");
  else log.printf("The g(r) will be normalized with a density N/V . \n");
  parseFlag("ONE_BODY",one_body);
  temperature = -1.;
  mass = -1.;
  parse("TEMPERATURE",temperature);
  parse("MASS",mass);
  if (one_body) {
     if (temperature>0 && mass>0 && local_density ) log.printf("The one-body entropy will be added to the pair entropy. \n");
     if (temperature<0) error("ONE_BODY keyword used but TEMPERATURE not given. Specify a temperature greater than 0 in Kelvin using the TEMPERATURE keyword. ");
     if (mass<0) error("ONE_BODY keyword used but MASS not given. Specify a mass greater than 0 in g/mol using the MASS keyword. ");
     if (!local_density) error("ONE_BODY keyword used but LOCAL_DENSITY not given. LOCAL_DENSITY flag is compulsory with ONE_BODY.");
     double planck = 6.62607004e-16; // nm2 kg / s 
     double boltzmann = 1.38064852e-5; // nm2 kg s-2 K-1
     double avogadro= 6.0221409e23 ;
     double deBroglie = planck/std::sqrt(2*pi*(mass*1.e-3/avogadro)*boltzmann*temperature);
     deBroglie3 = deBroglie*deBroglie*deBroglie;
     log.printf("The thermal deBroglie wavelength is %f nm. Be sure to use nm as units of distance. \n", deBroglie);
  }
  parseFlag("NO_TWO_BODY",no_two_body);
  if (no_two_body) {
     if (one_body) log.printf("The two-body entropy will be removed from the pair entropy. Only the one-body term is kept. \n");
     else error("NO_TWO_BODY keyword used but ONE_BODY not specified. ONE_BODY flag is compulsory with NO_TWO_BODY.");
  }
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  if (doOutputGofr) log.printf("The g(r) of each atom will be written to a file. \n");
  if (doOutputGofr && comm.Get_size()!=1) error("OUTPUT_GOFR cannot be used with more than one MPI thread");
  outputStride=1;
  parse("OUTPUT_STRIDE",outputStride);
  if (!doOutputGofr && outputStride!=1) error("OUTPUT_STRIDE specified but OUTPUT_GOFR not given. Specify OUTPUT_GOFR or remove OUTPUT_STRIDE"); 
  if (doOutputGofr) log.printf("The output stride to write g(r) or the integrand is %d \n", outputStride);

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); 

  // Allocate space for the average g(r)
  doAverageGofr=false;
  parseFlag("AVERAGE_GOFR",doAverageGofr);
  if (doAverageGofr) {
    avgGofr.resize(getFullNumberOfTasks(),nhist);
    iteration.resize(getFullNumberOfTasks());
    avgDensity=0.;
  }
  averageGofrTau=0;
  parse("AVERAGE_GOFR_TAU",averageGofrTau);
  if (averageGofrTau!=0 && !doAverageGofr) error("AVERAGE_GOFR_TAU specified but AVERAGE_GOFR not given. Specify AVERAGE_GOFR or remove AVERAGE_GOFR_TAU");
  if (doAverageGofr && averageGofrTau==0) log.printf("The g(r) will be averaged over all frames \n");
  if (doAverageGofr && averageGofrTau!=0) log.printf("The g(r) will be averaged with a window of %d steps \n", averageGofrTau);

  doCutAtFirstPeak=false;
  cutAtFirstPeakValue=-1.;
  parseFlag("CUT_AT_FIRST_PEAK",doCutAtFirstPeak);
  parse("CUT_AT_FIRST_PEAK_VALUE",cutAtFirstPeakValue);
  if (!doCutAtFirstPeak && cutAtFirstPeakValue>0.) error("CUT_AT_FIRST_PEAK_VALUE cannot be used if CUT_AT_FIRST_PEAK has not been specified");
  if (doCutAtFirstPeak) {
    log.printf("   The g(r) will be set to 1 after the first peak \n");
    if (cutAtFirstPeakValue<0.) log.printf("   The first peak is determined as the first root of the second derivative \n");
    else log.printf("   The first peak is determined as the position where g(r) is equal to %f \n", cutAtFirstPeakValue);
  }

  checkRead();
  // Define heavily used constants
  double sqrt2piSigma = std::sqrt(2*pi)*sigma;
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/(nhist-1);
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(4*sigma/deltar); //4*sigma
  vectorX.resize(nhist);
  vectorX2.resize(nhist);
  for(unsigned i=0;i<nhist;++i){
    vectorX[i]=deltar*i;
    vectorX2[i]=vectorX[i]*vectorX[i];
  }

  // Set the link cell cutoff
  setLinkCellCutoff( maxr + 3*sigma );
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);
  rcut = std::sqrt(rcut2);
  log.printf("Setting cut off to %f \n ", maxr + 3*sigma );

  if (doOutputGofr) gofrOfile.open("gofr-multi.txt");
}

PairEntropies::~PairEntropies() {
  if (doOutputGofr) gofrOfile.close();
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
             gofr[j] += kernel(vectorX[j]-distanceModulo, dfunc);
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
   double TwoPiDensity = 2*pi*density;
   double normConstantBase = 2*TwoPiDensity;
   for(unsigned i=1;i<nhist;++i){
     double normConstant = normConstantBase*vectorX2[i];
     if (density<1.e-20) gofr[i] = 1.;
     else gofr[i] /= normConstant;
     if (!doNotCalculateDerivatives()) {
       gofrVirial[i] /= normConstant;
       for(unsigned j=0;j<myatoms.getNumberOfAtoms();++j){
         gofrPrime[i][j] /= normConstant;
       }
     }
   }
   if (doCutAtFirstPeak) {
      if (!doNotCalculateDerivatives()) error("Cannot calculate derivatives or bias using the CUT_AT_FIRST_PEAK option");
      int posFirstPeak=nhist;
      if (cutAtFirstPeakValue<0.) {
         // Check when second derivative changes sign 
         double secondDerPrev=(gofr[2]-2*gofr[1]+gofr[0])/(deltar*deltar);
         for(unsigned i=2;i<(nhist-1);++i){
            double secondDer=(gofr[i+1]-2*gofr[i]+gofr[i-1])/(deltar*deltar);
            if (secondDer*secondDerPrev<0.) {
               posFirstPeak=i;
               break;
            }
            secondDerPrev=secondDer;
         }
      } else {
         for(unsigned i=0;i<nhist;++i){
            if(gofr[i]>cutAtFirstPeakValue) {
               posFirstPeak=i;
               break;
            }
         }
      }
      for(unsigned i=posFirstPeak;i<nhist;++i){
         gofr[i] = 1.;
      }
   }
   if (doAverageGofr) {
      if (!doNotCalculateDerivatives()) error("Cannot calculate derivatives or bias using the AVERAGE_GOFR option");
      double factor;
      if (averageGofrTau==0 || (iteration[myatoms.getIndex(0)] < averageGofrTau) ) {
         iteration[myatoms.getIndex(0)] += 1;
         factor = 1./( (double) iteration[myatoms.getIndex(0)] );
      } else factor = 2./((double) averageGofrTau + 1.);
      for(unsigned i=0;i<nhist;++i){
         avgGofr[myatoms.getIndex(0)][i] += (gofr[i]-avgGofr[myatoms.getIndex(0)][i])*factor;
         gofr[i] = avgGofr[myatoms.getIndex(0)][i];
         avgDensity += (density-avgDensity)*factor;
         density = avgDensity;
      }
   }
   if (doOutputGofr && (getStep()%outputStride==0)) outputGofr(myatoms.getIndex(0)+1,getFullNumberOfTasks(),getStep(),outputStride,gofr);
   // Construct integrand
   vector<double> integrand(nhist);
   for(unsigned i=0;i<nhist;++i){
     logGofr[i] = std::log(gofr[i]);
     if (gofr[i]<1.e-10) {
       integrand[i] = vectorX2[i];
     } else {
       integrand[i] = (gofr[i]*logGofr[i]-gofr[i]+1)*vectorX2[i];
     }
   }
   // Integrate to obtain pair entropy;
   double entropy=0.;
   if (!no_two_body) {
      entropy += -TwoPiDensity*integrate(integrand,deltar); 
   }
   if (one_body) {
      entropy += 5./2. - std::log(density*deBroglie3);
   }
   if (!doNotCalculateDerivatives()) {
     // Construct integrand and integrate derivatives
     for(unsigned i=0;i<myatoms.getNumberOfAtoms();++i) {
       vector<Vector> integrandDerivatives(nhist);
       for(unsigned j=0;j<nhist;++j){
         if (gofr[j]>1.e-10) {
           integrandDerivatives[j] = gofrPrime[j][i]*logGofr[j]*vectorX2[j];
         }
       }
       // Integrate
       deriv[i] = -TwoPiDensity*integrate(integrandDerivatives,deltar);
     }
     // Virial of positions
     // Construct virial integrand
     vector<Tensor> integrandVirial(nhist);
     for(unsigned i=0;i<nhist;++i){
       if (gofr[i]>1.e-10) {
         integrandVirial[i] = gofrVirial[i]*logGofr[i]*vectorX2[i];
       }
     }
     // Integrate virial
     virial = -TwoPiDensity*integrate(integrandVirial,deltar);
     // Virial of volume
     // Construct virial integrand
     vector<double> integrandVirialVolume(nhist);
     for(unsigned i=0;i<nhist;i+=1) {   
       integrandVirialVolume[i] = (-gofr[i]+1)*vectorX2[i];
     }
     // Integrate virial
     virial += -TwoPiDensity*integrate(integrandVirialVolume,deltar)*Tensor::identity();
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

void PairEntropies::outputGofr(int index, int nat, int step, int stride, vector<double> gofr)const{
  gofrOfile.printf("# Atom index %d Step number %d Function number %d \n",index,step,((step/stride)*nat)+index-1);
  for(unsigned i=0;i<gofr.size();++i){
     gofrOfile.printField("r",vectorX[i]).printField("gofr",gofr[i]).printField();
  }
  gofrOfile.printf(" \n");
  gofrOfile.printf(" \n");
}



}
}

