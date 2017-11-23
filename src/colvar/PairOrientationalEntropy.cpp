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

#include <string>
#include <math.h>

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

class PairOrientationalEntropy : public Colvar {
  bool pbc, serial, invalidateList, firsttime, doneigh;
  NeighborListParallel *nl;
  vector<AtomNumber> center_lista,start_lista,end_lista,origin_lista;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  double maxr;
  vector<int> nhist_;
  vector<double> sigma_;
  double rcut2;
  double invSqrt2piSigma, sigmaSqr2, sigmaSqr;
  double deltar, deltaAngle, deltaCosAngle;
  unsigned deltaBin, deltaBinAngle;
  // Integration routines
  double integrate(Matrix<double> integrand, vector<double> delta)const;
  //double integrate(Matrix<double> integrand, vector<double> x1, vector<double> x2)const;
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  vector<double> integrateMarginal(Matrix<double> integrand, double delta, unsigned dim)const;
  vector<double> x1;
  vector<double> x2;
  // Kernel to calculate g(r)
  double kernel(vector<double> distance, double&der)const;
  // Output gofr and integrand
  void outputGofr(Matrix<double> gofr, const char* fileName);
  void outputIntegrand(vector<double> integrand);
  void output1Dfunction(vector<double> y, vector<double> x, const char* fileName);
  int outputStride;
  bool doOutputGofr, doOutputIntegrand;
  // Average gofr
  Matrix<double> avgGofr;
  unsigned iteration;
  bool doAverageGofr;
public:
  explicit PairOrientationalEntropy(const ActionOptions&);
  ~PairOrientationalEntropy();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PairOrientationalEntropy,"PAIR_ORIENTATIONAL_ENTROPY")

void PairOrientationalEntropy::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r)");
  keys.addFlag("AVERAGE_GOFR",false,"Average g(r) over time");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.addFlag("OUTPUT_INTEGRAND",false,"Output integrand");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","ORIGIN","Define an atom that represents the origin from which to calculate the g(r,theta)");
  keys.add("atoms","CENTER","Center atoms");
  keys.add("atoms","START","Start point of vector defining orientation");
  keys.add("atoms","END","End point of vector defining orientation");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","1","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
}

PairOrientationalEntropy::PairOrientationalEntropy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  parseAtomList("CENTER",center_lista);
  parseAtomList("START",start_lista);
  parseAtomList("END",end_lista);
  if(center_lista.size()!=start_lista.size()) error("Number of atoms in START must be equal to the number of atoms in CENTER");
  if(center_lista.size()!=end_lista.size()) error("Number of atoms in START must be equal to the number of atoms in CENTER");

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// neighbor list stuff
  doneigh=false;
  double nl_cut=0.0;
  double nl_skin;
  int nl_st=-1;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   //if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  parse("MAXR",maxr);
  log.printf("Integration in the interval from 0. to %f \n", maxr );

  parseVector("SIGMA",sigma_);
  if(sigma_.size() != 2) error("SIGMA keyword takes two input values");
  log.printf("The pair distribution function is calculated with a Gaussian kernel with deviations %f and %f \n", sigma_[0], sigma_[1]);
  double rcut = maxr + 3*sigma_[0];
  rcut2 = (maxr + 3*sigma_[0])*(maxr + 3*sigma_[0]);  // 3*sigma is hard coded
  if(doneigh){
    if(nl_cut<rcut) error("NL_CUTOFF should be larger than MAXR + 3*SIGMA");
    nl_skin=nl_cut-maxr;
  }

  nhist_.resize(2);
  nhist_[0]=ceil(maxr/sigma_[0]) + 1; // Default value
  nhist_[1]=ceil(2./sigma_[1]) + 1; // Default value
  parseVector("NHIST",nhist_);
  if(nhist_.size() != 2) error("NHIST keyword takes two input values");
  log.printf("The interval is partitioned in %u equal parts in r and %u equal parts in theta. The integration is performed with the trapezoid rule. \n", nhist_[0], nhist_[1] );

  doOutputGofr=false;
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  if (doOutputGofr) { 
     log.printf("  The g(r) will be written to a file \n.");
  }
  doOutputIntegrand=false;
  parseFlag("OUTPUT_INTEGRAND",doOutputIntegrand);
  if (doOutputIntegrand) {
     log.printf("  The integrand will be written to a file \n.");
  }
  outputStride=1;
  parse("OUTPUT_STRIDE",outputStride);
  if (outputStride!=1 && !doOutputGofr && !doOutputIntegrand) error("Cannot specify OUTPUT_STRIDE if OUTPUT_GOFR or OUTPUT_INTEGRAND not used");
  if (outputStride<1) error("The output stride specified with OUTPUT_STRIDE must be greater than or equal to one.");
  if (outputStride>1) log.printf("  The output stride to write g(r) or the integrand is %d \n", outputStride);

  doAverageGofr=false;
  parseFlag("AVERAGE_GOFR",doAverageGofr);
  if (doAverageGofr) {
     iteration = 1;
     log.printf("  The g(r) will be averaged over all frames");
     avgGofr.resize(nhist_[0],nhist_[1]);
  }

  checkRead();

  // Neighbor lists
  if (doneigh) {
    nl= new NeighborListParallel(center_lista,pbc,getPbc(),comm,log,nl_cut,nl_st,nl_skin);
    log.printf("  using neighbor lists with\n");
    log.printf("  cutoff %f, and skin %f\n",nl_cut,nl_skin);
    if(nl_st>=0){
      log.printf("  update every %d steps\n",nl_st);
    } else {
      log.printf("  checking every step for dangerous builds and rebuilding as needed\n");
    }
  }
  atomsToRequest.reserve ( center_lista.size() + start_lista.size() + end_lista.size() );
  atomsToRequest.insert (atomsToRequest.end(), center_lista.begin(), center_lista.end() );
  atomsToRequest.insert (atomsToRequest.end(), start_lista.begin(), start_lista.end() );
  atomsToRequest.insert (atomsToRequest.end(), end_lista.begin(), end_lista.end() );
  requestAtoms(atomsToRequest);

  // Define heavily used expressions
  double sqrt2piSigma = std::sqrt(2*pi)*sigma_[0];
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma_[0]*sigma_[0];
  sigmaSqr = sigma_[0]*sigma_[0];
  deltar=maxr/(nhist_[0]-1);
  deltaCosAngle=2./(nhist_[1]-1);
  deltaBin = std::floor(3*sigma_[0]/deltar); // 3*sigma is hard coded
  deltaBinAngle = std::floor(3*sigma_[1]/deltaCosAngle); // 3*sigma is hard coded

  x1.resize(nhist_[0]);
  x2.resize(nhist_[1]);
  for(unsigned j=0;j<nhist_[0];++j){
     x1[j]=deltar*j;
     //x1[j]=deltar*(j+0.5);
  }
  for(unsigned j=0;j<nhist_[1];++j){
     //double theta=(pi/nhist_[1])*(j+0.5);
     x2[j]=-1+deltaCosAngle*j;
  }
}

PairOrientationalEntropy::~PairOrientationalEntropy(){
  if (doneigh) {
     nl->printStats();
     delete nl;
  }
}

void PairOrientationalEntropy::prepare(){
  if(doneigh && nl->getStride()>0){
    requestAtoms(nl->getFullAtomList());
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
void PairOrientationalEntropy::calculate()
{
  // Define output quantities
  double pairEntropy;
  vector<Vector> deriv(getNumberOfAtoms());
  Tensor virial;
  // Define intermediate quantities
  Matrix<double> gofr(nhist_[0],nhist_[1]);
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
  if (doneigh) {
    if(invalidateList){
      vector<Vector> centerPositions(getPositions().begin(),getPositions().begin() + center_lista.size());;
      nl->update(centerPositions);
    }
    // Loop over neighbors
    const unsigned nn=nl->size();
    //log.printf("nn %d \n",nn);
    for(unsigned int i=rank;i<nn;i+=stride) {
      double dfunc, d2;
      Vector distance;
      Vector distance_versor;
      unsigned i0=nl->getClosePair(i).first;
      unsigned i1=nl->getClosePair(i).second;
      if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
      //log.printf("Center1 %d Center2 %d \n", getAbsoluteIndex(i0), getAbsoluteIndex(i1) );
      if(pbc){
       distance=pbcDistance(getPosition(i0),getPosition(i1));
      } else {
       distance=delta(getPosition(i0),getPosition(i1));
      }
      if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
        double distanceModulo=std::sqrt(d2);
        Vector distance_versor = distance / distanceModulo;
        unsigned bin=std::floor(distanceModulo/deltar);
        unsigned atom1_mol1=i0+center_lista.size();
        unsigned atom2_mol1=i0+center_lista.size()+start_lista.size();
        unsigned atom1_mol2=i1+center_lista.size();
        unsigned atom2_mol2=i1+center_lista.size()+start_lista.size();
        //log.printf("Center1 %d Center2 %d atom1_mol2 %d atom1_mol2 %d \n", getAbsoluteIndex(i0), getAbsoluteIndex(i1), getAbsoluteIndex(atom1_mol2), getAbsoluteIndex(atom2_mol2) );
        Vector mol_vector1=pbcDistance(getPosition(atom1_mol1),getPosition(atom2_mol1)); 
        //Vector mol_vector2=pbcDistance(getPosition(atom1_mol2),getPosition(atom2_mol2));
        //Vector mol_vector1=pbcDistance(getPosition(i0),getPosition(i1));
        Vector mol_vector2=pbcDistance(getPosition(atom1_mol2),getPosition(atom2_mol2));
        double norm_v1 = std::sqrt(mol_vector1[0]*mol_vector1[0]+mol_vector1[1]*mol_vector1[1]+mol_vector1[2]*mol_vector1[2]);
        double norm_v2 = std::sqrt(mol_vector2[0]*mol_vector2[0]+mol_vector2[1]*mol_vector2[1]+mol_vector2[2]*mol_vector2[2]);
        double inv_v1=1./norm_v1;
        double inv_v2=1./norm_v2;
        double cosAngle=dotProduct(mol_vector1,mol_vector2)*inv_v1*inv_v2;
        double angle=acos(cosAngle);
        //log.printf("Angle %f radians %f degrees \n", angle, angle*180./pi);
        unsigned binAngle=std::floor((cosAngle+1.)/deltaCosAngle);
        int minBin, maxBin; // These cannot be unsigned
        // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
        minBin=bin - deltaBin;
        if (minBin < 0) minBin=0;
        if (minBin > (nhist_[0]-1)) minBin=nhist_[0]-1;
        maxBin=bin +  deltaBin;
        if (maxBin > (nhist_[0]-1)) maxBin=nhist_[0]-1;
        int minBinAngle, maxBinAngle; // These cannot be unsigned
        // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual angle
        minBinAngle=binAngle - deltaBinAngle;
        maxBinAngle=binAngle +  deltaBinAngle;
        //log.printf("minBinAngle %d maxBinAngle %d binAngle %d deltaBinAngle %d \n",minBinAngle,maxBinAngle,binAngle,deltaBinAngle);
        for(int k=minBin;k<maxBin+1;k+=1) {
          double x=deltar*k;
          for(int l=minBinAngle;l<maxBinAngle+1;l+=1) {
             // Include periodic effects
             int h;
             if (l<0) {
                h=-l;
             } else if (l>(nhist_[1]-1)) {
                h=2*nhist_[1]-l-2;
             } else {
                h=l;
             }
             double theta=-1+deltaCosAngle*l;
             //log.printf("l %d h %d theta %f cosAngle %f \n",l,h,theta,cosAngle);
             vector<double> pos(2);
             pos[0]=x-distanceModulo;
             pos[1]=theta-cosAngle;
             if (l==(nhist_[1]-1)) gofr[k][h] += 2*kernel(pos, dfunc);
             else if (l==0) gofr[k][h] += 2*kernel(pos, dfunc);
             else gofr[k][h] += kernel(pos, dfunc);
          }
        }
      }
    }
  } else {
    for(unsigned int i=rank;i<(center_lista.size()-1);i+=stride) {
      for(unsigned int j=i+1;j<center_lista.size();j+=1) {
        double dfunc, d2;
        Vector distance;
        Vector distance_versor;
        if(getAbsoluteIndex(i)==getAbsoluteIndex(j)) continue;
        if(pbc){
         distance=pbcDistance(getPosition(i),getPosition(j));
        } else {
         distance=delta(getPosition(i),getPosition(j));
        }
        if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
          double distanceModulo=std::sqrt(d2);
          Vector distance_versor = distance / distanceModulo;
          unsigned bin=std::floor(distanceModulo/deltar);
          unsigned atom1_mol1=i+center_lista.size();
          unsigned atom2_mol1=i+center_lista.size()+start_lista.size();
          unsigned atom1_mol2=j+center_lista.size();
          unsigned atom2_mol2=j+center_lista.size()+start_lista.size();
          //log.printf("Center1 %d Center2 %d atom1_mol2 %d atom1_mol2 %d \n", getAbsoluteIndex(i0), getAbsoluteIndex(i1), getAbsoluteIndex(atom1_mol2), getAbsoluteIndex(atom2_mol2) );
          Vector mol_vector1=pbcDistance(getPosition(atom1_mol1),getPosition(atom2_mol1)); 
          //Vector mol_vector2=pbcDistance(getPosition(atom1_mol2),getPosition(atom2_mol2));
          //Vector mol_vector1=pbcDistance(getPosition(i0),getPosition(i1));
          Vector mol_vector2=pbcDistance(getPosition(atom1_mol2),getPosition(atom2_mol2));
          double norm_v1 = std::sqrt(mol_vector1[0]*mol_vector1[0]+mol_vector1[1]*mol_vector1[1]+mol_vector1[2]*mol_vector1[2]);
          double norm_v2 = std::sqrt(mol_vector2[0]*mol_vector2[0]+mol_vector2[1]*mol_vector2[1]+mol_vector2[2]*mol_vector2[2]);
          double inv_v1=1./norm_v1;
          double inv_v2=1./norm_v2;
          double cosAngle=dotProduct(mol_vector1,mol_vector2)*inv_v1*inv_v2;
          double angle=acos(cosAngle);
          //log.printf("Angle %f radians %f degrees \n", angle, angle*180./pi);
          unsigned binAngle=std::floor((cosAngle+1.)/deltaCosAngle);
          int minBin, maxBin; // These cannot be unsigned
          // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
          minBin=bin - deltaBin;
          if (minBin < 0) minBin=0;
          if (minBin > (nhist_[0]-1)) minBin=nhist_[0]-1;
          maxBin=bin +  deltaBin;
          if (maxBin > (nhist_[0]-1)) maxBin=nhist_[0]-1;
          int minBinAngle, maxBinAngle; // These cannot be unsigned
          // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual angle
          minBinAngle=binAngle - deltaBinAngle;
          maxBinAngle=binAngle +  deltaBinAngle;
          //log.printf("minBinAngle %d maxBinAngle %d binAngle %d deltaBinAngle %d \n",minBinAngle,maxBinAngle,binAngle,deltaBinAngle);
          for(int k=minBin;k<maxBin+1;k+=1) {
            double x=deltar*k;
            for(int l=minBinAngle;l<maxBinAngle+1;l+=1) {
               // Include periodic effects
               int h;
               if (l<0) {
                  h=-l;
               } else if (l>(nhist_[1]-1)) {
                  h=2*nhist_[1]-l-2;
               } else {
                  h=l;
               }
               double theta=-1+deltaCosAngle*l;
               //log.printf("l %d h %d theta %f cosAngle %f \n",l,h,theta,cosAngle);
               vector<double> pos(2);
               pos[0]=x-distanceModulo;
               pos[1]=theta-cosAngle;
               if (l==(nhist_[1]-1)) gofr[k][h] += 2*kernel(pos, dfunc);
               else if (l==0) gofr[k][h] += 2*kernel(pos, dfunc);
               else gofr[k][h] += kernel(pos, dfunc);
            }
          }
        }
      }
    }
  }
  if(!serial){
    comm.Sum(&gofr[0][0],nhist_[0]*nhist_[1]);
  }
  vector<double> delta(2);
  delta[0]=deltar;
  delta[1]=deltaCosAngle;
  //log.printf("Integral of gaussians is %f \n", integrate(gofr,delta));
  // Calculate volume and density
  double volume=getBox().determinant();
  double density=center_lista.size()/volume;
  // Normalize g(r)
  double normConstantBase;
  if (origin_lista.size()>0) {
     normConstantBase = 4*pi*density;
  } else {
     normConstantBase = 2*pi*center_lista.size()*density;
  }
  // Take into account "volume" of angles
  double volumeOfAngles = 2.;
  normConstantBase /= volumeOfAngles;
  for(unsigned j=1;j<nhist_[0];++j){ // Starts from 1, x=0 produces a divergence
    double x=deltar*j;
    for(unsigned k=0;k<nhist_[1];++k){
       double normConstant = normConstantBase*x*x;
       gofr[j][k] /= normConstant;
    }
  }
  if (doAverageGofr) {
     for(unsigned i=0;i<nhist_[0];++i){
        for(unsigned j=0;j<nhist_[1];++j){
           avgGofr[i][j] += (gofr[i][j]-avgGofr[i][j])/( (double) iteration);
           gofr[i][j] = avgGofr[i][j];
        }
     }
     iteration += 1;
  }
  // Output of gofr
  if (doOutputGofr && (getStep()%outputStride==0) && rank==0) outputGofr(gofr,"gofr.txt");
  /*
  vector<double> gofrMarginalR(nhist_[0]);
  gofrMarginalR=integrateMarginal(gofr,deltaCosAngle,1);
  for(unsigned i=0;i<nhist_[0];++i){
     gofrMarginalR[i] /= volumeOfAngles;
  }
  if (doOutputGofr && (getStep()%outputStride==0) && rank==0) {
     output1Dfunction(gofrMarginalR, x1, "marginal.txt");
  }
  Matrix<double> ConditionalDistr(nhist_[0], nhist_[1]);
  Matrix<double> integrand(nhist_[0], nhist_[1]);
  for(unsigned i=0;i<nhist_[0];++i){
     for(unsigned j=0;j<nhist_[1];++j){
        // The 1.e-10 might seem low but leads to a smooth behavior
        if (gofrMarginalR[i]>1.e-10) ConditionalDistr[i][j] = gofr[i][j]/gofrMarginalR[i];
        if (ConditionalDistr[i][j]> 1.e-10) integrand[i][j] = ConditionalDistr[i][j]*std::log(ConditionalDistr[i][j]);
        else integrand[i][j] = -ConditionalDistr[i][j];
     }
  }
  if (doOutputGofr && (getStep()%outputStride==0)  && rank==0) outputGofr(ConditionalDistr,"conditional.txt");
  vector<double> OLE(nhist_[0]);
  OLE = integrateMarginal(integrand,deltaCosAngle,1);
  for(unsigned i=0;i<nhist_[0];++i){
     OLE[i] *= -2*pi/volumeOfAngles;
  }
  if (doOutputGofr && (getStep()%outputStride==0)  && rank==0) {
     output1Dfunction(OLE, x1, "OLE.txt");
  }
  // Build integrands for orientational and translational part of the pair entropy
  vector<double> integrandOrientation(nhist_[0]);
  vector<double> integrandTranslation(nhist_[0]);
  for(unsigned j=0;j<nhist_[0];++j){
    double x=deltar*(j+0.5);
    double logGofrMarginalR = std::log(gofrMarginalR[j]);
    if (gofrMarginalR[j]<1.e-10) {
      integrandTranslation[j] = x*x;
    } else {
      integrandTranslation[j] = (gofrMarginalR[j]*logGofrMarginalR-gofrMarginalR[j]+1)*x*x;
    }
    integrandOrientation[j] = gofrMarginalR[j]*OLE[j]*x*x;
  }
  */
  setValue           (pairEntropy);
  /*
  if (doOutputGofr && (getStep()%outputStride==0)  && rank==0) {
     output1Dfunction(integrandTranslation, x1, "integrandTranslation.txt");
  }
  if (doOutputGofr && (getStep()%outputStride==0)  && rank==0) {
     output1Dfunction(integrandOrientation, x1, "integrandOrientation.txt");
  }
  */
  // Output of integrands
  //Value* translational=getPntrToComponent("translational");
  //Value* orientational=getPntrToComponent("orientational");
  //translational->set(-2*pi*density*integrate(integrandTranslation,deltar));
  //orientational->set(density*integrate(integrandOrientation,deltar));
}

double PairOrientationalEntropy::kernel(vector<double> distance,double&der)const{
  // Gaussian function and derivative
  double result = (1./(2.*pi*sigma_[0]*sigma_[1]))*std::exp(-distance[0]*distance[0]/(2*sigma_[0]*sigma_[0])-distance[1]*distance[1]/(2*sigma_[1]*sigma_[1])) ;
  der = 0.;
  return result;
}

double PairOrientationalEntropy::integrate(Matrix<double> integrand, vector<double> delta)const{
  // Trapezoid rule
  double result = 0.;
  for(unsigned i=1;i<(nhist_[0]-1);++i){
     for(unsigned j=1;j<(nhist_[1]-1);++j){
        result += integrand[i][j];
     }
  }
  // Edges
  for(unsigned i=1;i<(nhist_[0]-1);++i){
     result += 0.5*integrand[i][0];
     result += 0.5*integrand[i][nhist_[1]-1];
  }
  for(unsigned j=1;j<(nhist_[1]-1);++j){
     result += 0.5*integrand[0][j];
     result += 0.5*integrand[nhist_[0]-1][j];
  }
  // Corners
  result += 0.25*integrand[0][0];
  result += 0.25*integrand[nhist_[0]-1][0];
  result += 0.25*integrand[0][nhist_[1]-1];
  result += 0.25*integrand[nhist_[0]-1][nhist_[1]-1];
  // Spacing
  result *= delta[0]*delta[1];
  return result;
}

vector<double> PairOrientationalEntropy::integrateMarginal(Matrix<double> integrand, double delta, unsigned dim)const{
  // First argument: Matrix defining the integrand. Only 2 dimensions are allowed.
  // Second argument: The spacing between bins in dimension to be integrated
  // Third argument: Dimension to be integrated. 0 stands for the first dimension, 1 for the second.
  // Trapezoid rule
  if (dim==0) {
     vector<double> result(nhist_[1]);
     for(unsigned i=0;i<nhist_[1];++i){
        for(unsigned j=1;j<(nhist_[0]-1);++j){
           result[i] += integrand[j][i];
        }
        result[i] += 0.5*integrand[0][i];
        result[i] += 0.5*integrand[nhist_[0]-1][i];
        result[i] *= delta;
     }
     return result;
  } else if (dim==1) {
     vector<double> result(nhist_[0]);
     for(unsigned i=0;i<nhist_[0];++i){
        for(unsigned j=1;j<(nhist_[1]-1);++j){
           result[i] += integrand[i][j];
        }
        result[i] += 0.5*integrand[i][0];
        result[i] += 0.5*integrand[i][nhist_[1]-1];
        result[i] *= delta;
     }
     return result;
  }
}

double PairOrientationalEntropy::integrate(vector<double> integrand, double delta)const{
  // Trapezoid rule
  double result=0.;
  for(unsigned i=1;i<(integrand.size()-1);++i){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}



Vector PairOrientationalEntropy::integrate(vector<Vector> integrand, double delta)const{
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

Tensor PairOrientationalEntropy::integrate(vector<Tensor> integrand, double delta)const{
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

void PairOrientationalEntropy::outputGofr(Matrix<double> gofr, const char* fileName) {
  PLMD::OFile gofrOfile;
  gofrOfile.open(fileName);
  for(unsigned i=0;i<nhist_[0];++i){
     double r=deltar*i;
     for(unsigned j=0;j<nhist_[1];++j){
        double theta=-1+deltaCosAngle*j;
        gofrOfile.printField("r",r).printField("theta",theta).printField("gofr",gofr[i][j]).printField();
     }
     gofrOfile.printf("\n");
  }
  gofrOfile.close();
}

void PairOrientationalEntropy::output1Dfunction(vector<double> y, vector<double> x, const char* fileName) {
  PLMD::OFile Outputfile;
  Outputfile.open(fileName);
  if (x.size() != y.size() ) error("Sizes of x and y arrays to be written with output1Dfunction are different.");
  for(unsigned i=0;i<y.size();++i){
     Outputfile.printField("x",x[i]).printField("y",y[i]).printField();
  }
  Outputfile.close();
}


void PairOrientationalEntropy::outputIntegrand(vector<double> integrand) {
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
