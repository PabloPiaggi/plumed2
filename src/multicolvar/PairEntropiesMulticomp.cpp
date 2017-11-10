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
#include "tools/Tools.h"

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


class PairEntropiesMulticomp : public MultiColvarBase {
private:
  double rcut, rcut2, rcut3 ;
  double invSqrt2piSigma, sigmaSqr2, sigmaSqr;
  double maxr, sigma;
  unsigned nhist;
  double density_given, density_givenA, density_givenB;
  bool local_density, one_body, no_two_body;
  double deltar;
  unsigned deltaBin;
  double temperature, mass, deBroglie3;
  std::vector<double> vectorX, vectorX2;
  vector<int> atomType;
  unsigned numberOfAatoms, numberOfBatoms;
  // Integration routine
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double&der)const;
  // Output stuff
  bool doOutputGofr;
  unsigned outputStride;
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropiesMulticomp(const ActionOptions&);
  ~PairEntropiesMulticomp();
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
  void outputGofr(int index, int nat, int step, int stride, vector<double> gofrAA, vector<double> gofrAB, vector<double> gofrBB)const;
  mutable OFile gofrOfile;

};

PLUMED_REGISTER_ACTION(PairEntropiesMulticomp,"PAIRENTROPIES_MULTICOMP")

void PairEntropiesMulticomp::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("atoms","GROUPA","Atoms of type A");
  keys.add("atoms","GROUPB","Atoms of type B");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","300","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("optional","DENSITYA","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("optional","DENSITYB","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("optional","TEMPERATURE","Temperature in Kelvin. It is compulsory when keyword ONE_BODY is used");
  keys.add("optional","MASS","Mass in g/mol. It is compulsory when keyword ONE_BODY is used");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r) of AA, AB and BB pairs");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.addFlag("LOCAL_DENSITY",false,"Use the local density to normalize g(r). If not specified, N/V is used");
  keys.addFlag("ONE_BODY",false,"Add the one body term (S1 = 5/2 - ln(dens*deBroglie^3) ) to the entropy");
  keys.addFlag("NO_TWO_BODY",false,"Remove the two-body term. Only the one-body term is kept. This flag can only be used along with the ONE_BODY flag.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");

  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
  keys.add("optional","INTEGRAND_FILE","the file on which to write the integrand");
}

PairEntropiesMulticomp::PairEntropiesMulticomp(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  numberOfAatoms=ga_lista.size();
  numberOfBatoms=gb_lista.size();

  parse("MAXR",maxr);
  log.printf("Integration in the interval from 0. to %f . \n", maxr );
  parse("NHIST",nhist);
  log.printf("The interval is partitioned in %u equal parts and the integration is perfromed with the trapezoid rule. \n", nhist );
  parse("SIGMA",sigma);
  log.printf("The pair distribution function is calculated with a Gaussian kernel with deviation %f . \n", sigma);

  density_given = -1;
  density_givenA = -1;
  density_givenB = -1;
  parse("DENSITY",density_given);
  parse("DENSITYA",density_givenA);
  parse("DENSITYB",density_givenB);
  if (density_given>0 && (density_givenA<0 or density_givenB<0)) error("DENSITY has been specified but DENSITYA or DENSITYB have not. ");
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

  if (getNumberOfAtoms()!=(ga_lista.size() + gb_lista.size() )) error("Number of atoms in SPECIES is different from the sum of the number of atoms in GROUPA and GROUPB. ");
  atomType.reserve(getNumberOfAtoms());
  for(unsigned i=0;i<getNumberOfAtoms();++i){
     atomType[i]=0;
     for(unsigned j=0;j<ga_lista.size();++j){
        if (getAbsoluteIndex(i)==ga_lista[j]) atomType[i]=1;
     }
     for(unsigned j=0;j<gb_lista.size();++j){
        if (getAbsoluteIndex(i)==gb_lista[j]) atomType[i]=2;
     }
     //log.printf("Reached atom %d index %d type %d \n", i, getAbsoluteIndex(i), atomType[i]);
     if (atomType[i]==0) error("At least one atom in SPECIES doesn't have a counterpart in GROUPA or GROUPB. ");
     //log.printf("index %d, atomType %d \n", getAbsoluteIndex(i), atomType[i]);
  }

  // Define heavily used constants
  double sqrt2piSigma = std::sqrt(2*pi)*sigma;
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/(nhist-1.);
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(3*sigma/deltar); //3*sigma is 99.7 %
  vectorX.resize(nhist);
  vectorX2.resize(nhist);
  for(unsigned i=0;i<nhist;++i){
    vectorX[i]=deltar*i;
    vectorX2[i]=vectorX[i]*vectorX[i];
  }

  // Set the link cell cutoff
  setLinkCellCutoff( maxr + 3*sigma );
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);
  rcut3 = rcut2*(maxr + 3*sigma);
  rcut = std::sqrt(rcut2);
  log.printf("Setting cut off to %f \n ", maxr + 3*sigma );

  checkRead();

  if (doOutputGofr) gofrOfile.open("gofr-multi.txt");
}

PairEntropiesMulticomp::~PairEntropiesMulticomp() {
  if (doOutputGofr) gofrOfile.close();
}

double PairEntropiesMulticomp::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   double dfunc, d2;
   Vector value;
   vector<double> gofrAA(nhist), gofrAB(nhist), gofrBB(nhist);
   vector<double> logGofrAA(nhist), logGofrAB(nhist), logGofrBB(nhist);
   //Matrix<Vector> gofrPrime(nhist,getNumberOfAtoms());
   //vector<Vector> deriv(getNumberOfAtoms());
   //vector<Tensor> gofrVirial(nhist);
   Tensor virial;
   // Construct g(r)
   int countNeighA=0;
   int countNeighB=0;
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
            if  (atomType[myatoms.getIndex(0)]==1 && atomType[myatoms.getIndex(i)]==1) {
               gofrAA[j] += kernel(vectorX[j]-distanceModulo, dfunc);
            } else if  (atomType[myatoms.getIndex(0)]==2 && atomType[myatoms.getIndex(i)]==2) {
               gofrBB[j] += kernel(vectorX[j]-distanceModulo, dfunc);
            } else {
               gofrAB[j] += kernel(vectorX[j]-distanceModulo, dfunc);
            }
            /*
            if (!doNotCalculateDerivatives()) {
              value = dfunc * distance_versor;
              gofrPrime[j][0] += value;
              gofrPrime[j][i] -= value;
              Tensor vv(value, distance);
              gofrVirial[j] += vv;
            }
            */
	 } 
         if (atomType[myatoms.getIndex(0)]==1) {
            countNeighA += 1;
         } else {
            countNeighA += 1;
         }
      }
   }
   // Normalize g(r)
   double volume=getBox().determinant(); 
   double density, densityA, densityB;
   if (density_given>0) {
      density=density_given;
      densityA=density_givenA;
      densityB=density_givenB;
   } else if (local_density) {
      double volumeSphere = (4./3.)*pi*rcut3;
      density= (double) (countNeighA+countNeighB) / volumeSphere;
      densityA= (double) countNeighA / volumeSphere;
      densityB= (double) countNeighB / volumeSphere;
   } else {
      density=getNumberOfAtoms()/volume; // This is (NA+NB)/V
      densityA=numberOfAatoms/volume; // This is NA/V
      densityB=numberOfBatoms/volume; // This is NB/V
   }
   //log.printf("rcut %f \n", rcut);
   //log.printf("countNeigh %d \n", countNeigh);
   //log.printf("density %f \n", density);
   double FourPiDensityA = 4*pi*densityA;
   double FourPiDensityB = 4*pi*densityB;
   for(unsigned i=1;i<nhist;++i){
     double normConstantAA = FourPiDensityA*vectorX2[i];
     double normConstantBB = FourPiDensityB*vectorX2[i];
     double normConstantAB;
     if (atomType[myatoms.getIndex(0)]==1) {
        normConstantAB = FourPiDensityB*vectorX2[i];
     } else {
        normConstantAB = FourPiDensityA*vectorX2[i];
     }
     if (atomType[myatoms.getIndex(0)]==1) {
        gofrAA[i] /= normConstantAA;
     } else {
        gofrBB[i] /= normConstantBB;
     }
     gofrAB[i] /= normConstantAB;
     /*
     if (!doNotCalculateDerivatives()) {
       gofrVirial[i] /= normConstant;
       for(unsigned j=0;j<myatoms.getNumberOfAtoms();++j){
         gofrPrime[i][j] /= normConstant;
       }
     }
     */
   }
   if (doOutputGofr && (getStep()%outputStride==0)) outputGofr(myatoms.getIndex(0)+1,getFullNumberOfTasks(),getStep(),outputStride,gofrAA,gofrAB,gofrBB);
   // Construct integrand
   vector<double> integrandAA(nhist), integrandAB(nhist), integrandBB(nhist);
   for(unsigned i=0;i<nhist;++i){
     logGofrAB[i] = std::log(gofrAB[i]);
     if (atomType[myatoms.getIndex(0)]==1) {
        logGofrAA[i] = std::log(gofrAA[i]);
        if (gofrAA[i]<1.e-10) {
           integrandAA[i] = vectorX2[i];
        } else {
           integrandAA[i] = (gofrAA[i]*logGofrAA[i]-gofrAA[i]+1)*vectorX2[i];
        }
     } else {
        logGofrBB[i] = std::log(gofrBB[i]);
        if (gofrBB[i]<1.e-10) {
           integrandBB[i] = vectorX2[i];
        } else {
           integrandBB[i] = (gofrBB[i]*logGofrBB[i]-gofrBB[i]+1)*vectorX2[i];
        }
     }
     if (gofrAB[i]<1.e-10) {
       integrandAB[i] = vectorX2[i];
     } else {
       integrandAB[i] = (gofrAB[i]*logGofrAB[i]-gofrAB[i]+1)*vectorX2[i];
     }
   }
   // Integrate to obtain pair entropy;
   double entropy=0.;
   double prefactorAA;
   double prefactorBB;
   double prefactorAB;
   double pairAAvalue=0.;
   double pairABvalue=0.;
   double pairBBvalue=0.;
   if (atomType[myatoms.getIndex(0)]==1) {
      prefactorAA = -2*pi*(densityA*densityA/density);
      pairAAvalue =  prefactorAA*integrate(integrandAA,deltar);
   } else {
      prefactorBB = -2*pi*(densityB*densityB/density);
      pairBBvalue = prefactorBB*integrate(integrandBB,deltar);
   }
   prefactorAB = -4*pi*(densityA*densityB/density);
   pairABvalue = prefactorAB*integrate(integrandAB,deltar);
   if (!no_two_body) {
      entropy += pairAAvalue+pairBBvalue+pairABvalue;
   }
   if (one_body) {
      entropy += 5./2.;
      if (densityA>0.) entropy -= std::log(densityA*deBroglie3);
      if (densityB>0.) entropy -= std::log(densityB*deBroglie3);
   }
   /*
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
   */
   return entropy;
}

double PairEntropiesMulticomp::kernel(double distance,double&der)const{
  // Gaussian function and derivative
  double result = invSqrt2piSigma*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

double PairEntropiesMulticomp::integrate(vector<double> integrand, double delta)const{
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

Vector PairEntropiesMulticomp::integrate(vector<Vector> integrand, double delta)const{
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

Tensor PairEntropiesMulticomp::integrate(vector<Tensor> integrand, double delta)const{
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

void PairEntropiesMulticomp::outputGofr(int index, int nat, int step, int stride, vector<double> gofrAA, vector<double> gofrAB, vector<double> gofrBB)const{
  gofrOfile.printf("# Atom index %d Step number %d Function number %d \n",index,step,((step/stride)*nat)+index-1);
  for(unsigned i=0;i<gofrAA.size();++i){
     gofrOfile.printField("r",vectorX[i]).printField("gofrAA",gofrAA[i]).printField("gofrAB",gofrAB[i]).printField("gofrBB",gofrBB[i]).printField();
  }
  gofrOfile.printf(" \n");
  gofrOfile.printf(" \n");
}


}
}

