/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR PAIRENTROPY
/*
Calculate the global pair entropy using the expression:
\f[
s=-2\pi\rho k_B \int\limits_0^{r_{\mathrm{max}}} \left [ g(r) \ln g(r) - g(r) + 1 \right ] r^2 dr .
\f]
where \f$ g(r) $\f is the pair distribution function and \f$ r_{\mathrm{max}} $\f is a cutoff in the integration (MAXR).
For the integration the interval from 0 to  \f$ r_{\mathrm{max}} $\f is partitioned in NHIST equal intervals. 
To make the calculation of \f$ g(r) $\f differentiable, the following function is used:
\f[
g(r) = \frac{1}{4 \pi \rho r^2} \sum\limits_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(r-r_{ij})^2/(2\sigma^2)} ,
\f]
where \f$ \rho $\f is the density and \f$ sigma $\f is a broadening parameter (SIGMA).  
\par Example)
The following input tells plumed to calculate the pair entropy of atoms 1-250 with themselves.
\verbatim
PAIRENTROPY ...
 LABEL=s2
 GROUPA=1-250
 MAXR=0.65
 SIGMA=0.025
 NHIST=100
 NLIST
 NL_CUTOFF=0.75
 NL_STRIDE=10
... PAIRENTROPY
\endverbatim
*/
//+ENDPLUMEDOC

class PairEntropy : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  bool invalidateList;
  bool firsttime;
  bool doOutputGofr;
  bool doOutputIntegrand;
  double maxr, sigma;
  unsigned nhist;
  double rcut2;
  double invSqrt2piSigma, sigmaSqr2, sigmaSqr;
  double deltar;
  unsigned deltaBin;
  // Integration routines
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double&der)const;
  // Output gofr and integrand
  void outputGofr(vector<double> gofr);
  void outputIntegrand(vector<double> integrand);
public:
  explicit PairEntropy(const ActionOptions&);
  ~PairEntropy();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PairEntropy,"PAIRENTROPY")

void PairEntropy::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r)");
  keys.addFlag("OUTPUT_INTEGRAND",false,"Output integrand");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","1","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
}

PairEntropy::PairEntropy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  parseFlag("OUTPUT_INTEGRAND",doOutputIntegrand);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }

  addValueWithDerivatives(); setNotPeriodic();
  if(gb_lista.size()>0){
    if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
  } else {
    if(doneigh)  nl= new NeighborList(ga_lista,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,pbc,getPbc());
  }

  requestAtoms(nl->getFullAtomList());

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }


  parse("MAXR",maxr);
  log.printf("Integration in the interval from 0. to %f nm. \n", maxr );
  parse("NHIST",nhist);
  log.printf("The interval is partitioned in %u equal parts and the integration is perfromed with the trapezoid rule. \n", nhist );
  parse("SIGMA",sigma);
  log.printf("The pair distribution function is calculated with a Gaussian kernel with deviation %f nm. \n", sigma);
  double rcut = maxr + 3*sigma;
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);  // 3*sigma is hard coded
  if(doneigh){
    if(nl_cut<rcut) error("NL_CUTOFF should be larger than MAXR + 3*SIGMA");
  }

  checkRead();

  // Define heavily used expressions
  double sqrt2piSigma = std::sqrt(2*pi)*sigma;
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/nhist;
  deltaBin = std::floor(3*sigma/deltar); // 3*sigma is hard coded
}

PairEntropy::~PairEntropy(){
  delete nl;
}

void PairEntropy::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void PairEntropy::calculate()
{
  // Define output quantities
  double pairEntropy;
  vector<Vector> deriv(getNumberOfAtoms());
  Tensor virial;
  // Define intermediate quantities
  vector<double> gofr(nhist);
  vector<double> logGofr(nhist);
  Matrix<Vector> gofrPrime(nhist,getNumberOfAtoms());
  vector<Tensor> gofrVirial(nhist);
  // Setup neighbor list and parallelization
  if(nl->getStride()>0 && invalidateList){
    nl->update(getPositions());
  }
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }else{
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }
  // Loop over neighbors
  const unsigned nn=nl->size();
  for(unsigned int i=rank;i<nn;i+=stride) {
    double dfunc, d2;
    Vector distance;
    Vector distance_versor;
    unsigned i0=nl->getClosePair(i).first;
    unsigned i1=nl->getClosePair(i).second;
    if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
    if(pbc){
     distance=pbcDistance(getPosition(i0),getPosition(i1));
    } else {
     distance=delta(getPosition(i0),getPosition(i1));
    }
    if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
      double distanceModulo=std::sqrt(d2);
      Vector distance_versor = distance / distanceModulo;
      unsigned bin=std::floor(distanceModulo/deltar);
      int minBin, maxBin; // These cannot be unsigned
      // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
      minBin=bin - deltaBin;
      if (minBin < 0) minBin=0;
      if (minBin > (nhist-1)) minBin=nhist-1;
      maxBin=bin +  deltaBin;
      if (maxBin > (nhist-1)) maxBin=nhist-1;
      for(int k=minBin;k<maxBin+1;k+=1) {
        double x=deltar*(k+0.5);
        gofr[k] += kernel(x-distanceModulo, dfunc);
        Vector value = dfunc * distance_versor;
        gofrPrime[k][i0] += value;
        gofrPrime[k][i1] -= value;
        Tensor vv(value, distance);
        gofrVirial[k] += vv;
      }
    }
  }
  if(!serial){
    comm.Sum(&gofr[0],nhist);
    comm.Sum(&gofrPrime[0][0],nhist*getNumberOfAtoms());
    comm.Sum(&gofrVirial[0],nhist);
  }
  // Calculate volume and density
  double volume=getBox().determinant();
  double density=getNumberOfAtoms()/volume;
  // Normalize g(r)
  double normConstantBase = 2*pi*getNumberOfAtoms()*density;
  for(unsigned j=0;j<nhist;++j){
    double x=deltar*(j+0.5);
    double normConstant = normConstantBase*x*x;
    gofr[j] /= normConstant;
    //log.printf(" gofr after %f %f \n",x, gofr[j]);
    gofrVirial[j] /= normConstant;
    for(unsigned k=0;k<getNumberOfAtoms();++k){
      gofrPrime[j][k] /= normConstant;
    }
  }
  // Output of gofr
  if (doOutputGofr && rank==0) outputGofr(gofr);
  // Construct integrand
  vector<double> integrand(nhist);
  for(unsigned j=0;j<nhist;++j){
    double x=deltar*(j+0.5);
    logGofr[j] = std::log(gofr[j]);
    if (gofr[j]<1.e-10) {
      integrand[j] = x*x;
    } else {
      integrand[j] = (gofr[j]*logGofr[j]-gofr[j]+1)*x*x;
    }
  }
  // Output of integrands
  if (doOutputIntegrand && rank==0) outputIntegrand(integrand);
  // Integrate to obtain pair entropy;
  pairEntropy = -2*pi*density*integrate(integrand,deltar);
  // Construct integrand and integrate derivatives
  if (!doNotCalculateDerivatives() ) {
    for(unsigned int j=rank;j<getNumberOfAtoms();j+=stride) {
      vector<Vector> integrandDerivatives(nhist);
      for(unsigned k=0;k<nhist;++k){
        double x=deltar*(k+0.5);
        if (gofr[k]>1.e-10) {
          integrandDerivatives[k] = gofrPrime[k][j]*logGofr[k]*x*x;
        }
      }
      // Integrate
      deriv[j] = -2*pi*density*integrate(integrandDerivatives,deltar);
    }
    comm.Sum(&deriv[0][0],3*getNumberOfAtoms());
    // Virial of positions
    // Construct virial integrand
    vector<Tensor> integrandVirial(nhist);
    for(unsigned j=0;j<nhist;++j){
      double x=deltar*(j+0.5);
      if (gofr[j]>1.e-10) {
        integrandVirial[j] = gofrVirial[j]*logGofr[j]*x*x;
      }
    }
    // Integrate virial
    virial = -2*pi*density*integrate(integrandVirial,deltar);
    // Virial of volume
    // Construct virial integrand
    vector<double> integrandVirialVolume(nhist);
    for(unsigned j=0;j<nhist;j+=1) {
      double x=deltar*(j+0.5);
      integrandVirialVolume[j] = (-gofr[j]+1)*x*x;
    }
    // Integrate virial
    virial += -2*pi*density*integrate(integrandVirialVolume,deltar)*Tensor::identity();
  }
  // Assign output quantities
  for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (pairEntropy);
  setBoxDerivatives  (virial);
}

double PairEntropy::kernel(double distance,double&der)const{
  // Gaussian function and derivative
  double result = invSqrt2piSigma*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

double PairEntropy::integrate(vector<double> integrand, double delta)const{
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

Vector PairEntropy::integrate(vector<Vector> integrand, double delta)const{
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

Tensor PairEntropy::integrate(vector<Tensor> integrand, double delta)const{
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

void PairEntropy::outputGofr(vector<double> gofr) {
  PLMD::OFile gofrOfile;
  gofrOfile.open("gofr.txt");
  for(unsigned i=0;i<gofr.size();++i){
     double r=deltar*(i+0.5);
     gofrOfile.printField("r",r).printField("gofr",gofr[i]).printField();
  }
  gofrOfile.close();
}

void PairEntropy::outputIntegrand(vector<double> integrand) {
  PLMD::OFile gofrOfile;
  gofrOfile.open("integrand.txt");
  for(unsigned i=0;i<integrand.size();++i){
     double r=deltar*(i+0.5);
     gofrOfile.printField("r",r).printField("integrand",integrand[i]).printField();
  }
  gofrOfile.close();
}


}
}
