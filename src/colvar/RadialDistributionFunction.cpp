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
#include "tools/NeighborListParallel.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"
#include "tools/IFile.h"
#include "math.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR RADIAL_DISTRIBUTION_FUNCTION
/*

*/
//+ENDPLUMEDOC

class RadialDistributionFunction : public Colvar {
  bool pbc, serial;
  // Neighbor list stuff
  bool doneigh;
  NeighborListParallel *nl;
  vector<AtomNumber> atoms_lista;
  bool invalidateList;
  bool firsttime;
  // Others
  double minr, maxr, sigma;
  unsigned nhist;
  double rcut2;
  double sqrt2piSigma, sigmaSqr2, sigmaSqr;
  double deltar;
  unsigned deltaBin;
  double density_given;
  std::vector<double> vectorX, vectorX2;
  // Kernel to calculate g(r)
  double kernel(double distance, double invNormKernel, double&der)const;
  // Output gofr
  bool doOutputGofr;
  unsigned outputStride;
  void outputGofr(vector<double> gofr);
  mutable PLMD::OFile gofrOfile;
  // Average g(r)
  bool doAverageGofr;
  vector<double> avgGofr;
  unsigned iteration;
  // Components' names
  std::vector<string> comp_names;
public:
  explicit RadialDistributionFunction(const ActionOptions&);
  ~RadialDistributionFunction();
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(RadialDistributionFunction,"RADIAL_DISTRIBUTION_FUNCTION")

void RadialDistributionFunction::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","List of atoms");
  keys.add("compulsory","MAXR","1.","Maximum distance for the radial distribution function");
  keys.add("compulsory","SIGMA","0.01","Width of gaussians to mollify the radial distribution function");
  keys.add("optional","MINR","Minimum distance for the radial distribution function");
  keys.add("optional","NHIST","Number of bins in the radial distribution function");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list. If non specified or negative, it checks every step and rebuilds as needed.");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r)");
  keys.addFlag("AVERAGE_GOFR",false,"Average g(r) over time");
  ActionWithValue::useCustomisableComponents(keys); //unknown number of output components

}

RadialDistributionFunction::RadialDistributionFunction(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  parseAtomList("ATOMS",atoms_lista);

  log.printf("  using periodic boundary conditions\n");

// neighbor list stuff
  doneigh=false;
  bool nl_full_list=false;
  double nl_cut=0.0;
  double nl_skin;
  int nl_st=-1;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
  }
  minr=0.;
  parse("MINR",minr);
  parse("MAXR",maxr);
  log.printf("  The g(r) is defined in the interval from %f to %f distance units,\n", minr, maxr );
  parse("SIGMA",sigma);
  log.printf("  is calculated with a Gaussian kernel with deviation %f distance units,\n", sigma);

  density_given = -1;
  parse("DENSITY",density_given);
  if (density_given>0) log.printf("  and is normalized with a density %f. \n", density_given);
  else log.printf("  and is normalized with a density N/V. \n");

  double rcut = maxr + 3*sigma;
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);  // 3*sigma is hard coded
  if(doneigh){
    if(nl_cut<rcut) error("NL_CUTOFF should be larger than MAXR + 3*SIGMA");
    nl_skin=nl_cut-(maxr+2*sigma);
  }
  nhist=ceil((maxr-minr)/sigma) + 1; // Default value
  parse("NHIST",nhist);
  log.printf("  The interval is partitioned in %u equal parts.\n", nhist );
 
  doOutputGofr=false;
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  if (doOutputGofr) { 
     log.printf("  The g(r) will be written to the file gofr.txt. \n");
     gofrOfile.link(*this);
     gofrOfile.open("gofr.txt");
  }

  outputStride=1;
  parse("OUTPUT_STRIDE",outputStride);
  if (outputStride!=1 && !doOutputGofr) error("Cannot specify OUTPUT_STRIDE if OUTPUT_GOFR not used");
  if (outputStride<1) error("The output stride specified with OUTPUT_STRIDE must be greater than or equal to one.");
  if (outputStride>1) log.printf("  The output stride to write g(r)  is %d \n", outputStride);

  doAverageGofr=false;
  parseFlag("AVERAGE_GOFR",doAverageGofr);
  if (doAverageGofr) {
     iteration = 1;
     log.printf("  The g(r) will be averaged over all frames");
     avgGofr.resize(nhist);
  }

  checkRead();

  // Setup neighbor list
  if (doneigh) {
    nl= new NeighborListParallel(atoms_lista,pbc,getPbc(),comm,log,nl_cut,nl_full_list,nl_st,nl_skin);
    requestAtoms(nl->getFullAtomList());
    log.printf("  using neighbor lists with\n");
    log.printf("  cutoff %f, and skin %f\n",nl_cut,nl_skin);
    if(nl_st>=0){
      log.printf("  update every %d steps\n",nl_st);
    } else {
      log.printf("  checking every step for dangerous builds and rebuilding as needed\n");
    }
    if (nl_full_list) {
      log.printf("  using a full neighbor list\n");
    } else {
      log.printf("  using a half neighbor list\n");
    }
  } else {
    requestAtoms(atoms_lista);
  }

  // Define heavily used expressions
  sqrt2piSigma = std::sqrt(2*pi)*sigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=(maxr-minr)/(nhist-1.);
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(3*sigma/deltar); // 3*sigma is hard coded
  vectorX.resize(nhist);
  vectorX2.resize(nhist);
  for(unsigned i=0;i<nhist;i++){
    vectorX[i]=minr+deltar*i;
    vectorX2[i]=vectorX[i]*vectorX[i];
  }

  for (unsigned int i=0;i<nhist;i++){
    std::string num; Tools::convert( i, num );
    string name; name=string("bin-")+num;
    comp_names.push_back(name);
    addComponentWithDerivatives(name); componentIsNotPeriodic(name);
  }
  turnOnDerivatives();

}

RadialDistributionFunction::~RadialDistributionFunction(){
  if (doneigh) {
     nl->printStats();
     delete nl;
  }
  if (doOutputGofr) gofrOfile.close();
}

void RadialDistributionFunction::prepare(){
  if(doneigh && nl->getStride()>0){
    if(firsttime) {
      invalidateList=true;
      firsttime=false;
    } else if ( (nl->getStride()>=0) &&  (getStep()%nl->getStride()==0) ){
      invalidateList=true;
    } else if ( (nl->getStride()<0) && !(nl->isListStillGood(getPositions())) ){
      invalidateList=true;
    } else {
      invalidateList=false;
    }
  }
}

// calculator
void RadialDistributionFunction::calculate()
{
  // Define intermediate quantities
  vector<double> gofr(nhist);
  Matrix<Vector> gofrPrime(getNumberOfAtoms(),nhist);
  vector<Tensor> gofrVirial(nhist);
  // Setup parallelization
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }else{
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }
  // Normalization constant
  double volume=getBox().determinant();
  double density;
  if (density_given>0) density=density_given;
  else density=getNumberOfAtoms()/volume;
  double TwoPiDensity = 2*pi*density;
  double normConstantBase = TwoPiDensity*getNumberOfAtoms(); // Normalization of g(r)
  normConstantBase *= sqrt2piSigma; // Normalization of gaussian
  double invNormConstantBase = 1./normConstantBase; 
  // Calculation of g(r)
  if (doneigh) {
    if(invalidateList){
      nl->update(getPositions());
    }
    // Loop over all atoms
    for(unsigned int i=0;i<nl->getNumberOfLocalAtoms();i++) {
       std::vector<unsigned> neighbors;
       unsigned index=nl->getIndexOfLocalAtom(i);
       neighbors=nl->getNeighbors(index);
       Vector position_index=getPosition(index);
       // Loop over neighbors
       for(unsigned int j=0;j<neighbors.size();j++) {  
         unsigned neighbor=neighbors[j];
         Vector distance=pbcDistance(position_index,getPosition(neighbor));
         double d2;
         if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor((distanceModulo-minr)/deltar);
           int minBin, maxBin; // These cannot be unsigned
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int k=minBin;k<maxBin+1;k++) {
             double invNormKernel=invNormConstantBase/vectorX2[k];
             double dfunc;
             gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc);
             if (!doNotCalculateDerivatives()) {
                Vector value = dfunc * distance_versor;
                gofrPrime[index][k] += value;
                gofrPrime[neighbor][k] -= value;
                Tensor vv(value, distance);
                gofrVirial[k] += vv;
             }
           }
         }
       }
    }
  } else {
    for(unsigned int i=rank;i<(getNumberOfAtoms()-1);i+=stride) {
      for(unsigned int j=i+1;j<getNumberOfAtoms();j++) {
         double d2;
         Vector distance=pbcDistance(getPosition(i),getPosition(j));
         if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           Vector distance_versor = distance / distanceModulo;
           unsigned bin=std::floor((distanceModulo-minr)/deltar);
           int minBin, maxBin; // These cannot be unsigned
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(int k=minBin;k<maxBin+1;k++) {
             double invNormKernel=invNormConstantBase/vectorX2[k];
             double dfunc;
             gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc);
             if (!doNotCalculateDerivatives()) {
                Vector value = dfunc * distance_versor;
                gofrPrime[i][k] += value;
                gofrPrime[j][k] -= value;
                Tensor vv(value, distance);
                gofrVirial[k] += vv;
             }
           }
         }
      }
    }
  }
  if(!serial){
    comm.Sum(gofr);
    if (!doNotCalculateDerivatives()) {
       comm.Sum(gofrPrime);
       comm.Sum(gofrVirial);
    }
  }  
  // Add virial of volume
  if (density_given<0) { // Only if using the density N/V
    for (unsigned int i=0;i<nhist;i++){
      gofrVirial[i] -= gofr[i]*Tensor::identity();
    }
  }
  // Average g(r)
  if (doAverageGofr) {
     //if (!doNotCalculateDerivatives()) error("Cannot use the AVERAGE_GOFR keyword when biasing");
     for(unsigned i=0;i<nhist;i++){
        avgGofr[i] += (gofr[i]-avgGofr[i])/( (double) iteration);
        gofr[i] = avgGofr[i];
     }
     iteration++;
  }
  // Output
  if (doOutputGofr && (getStep()%outputStride==0)) outputGofr(gofr);
  for (unsigned int i=0;i<nhist;i++){
    Value* value=getPntrToComponent(comp_names[i].c_str());
    value->set(gofr[i]);
    for(unsigned j=0;j<getNumberOfAtoms();++j) setAtomsDerivatives(value,j,gofrPrime[j][i]);
    setBoxDerivatives(value,gofrVirial[i]);
  }
}

double RadialDistributionFunction::kernel(double distance,double invNormKernel, double&der)const{
  // Gaussian function and derivative
  double result = invNormKernel*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

void RadialDistributionFunction::outputGofr(vector<double> gofr) {
  for(unsigned i=0;i<gofr.size();i++){
     gofrOfile.printField("r",vectorX[i]).printField("gofr",gofr[i]).printField();
  }
  gofrOfile.printf("\n");
  gofrOfile.printf("\n");
}

}
}
