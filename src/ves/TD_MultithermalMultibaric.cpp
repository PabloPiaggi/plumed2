/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"
#include "GridIntegrationWeights.h"
#include "core/ActionRegister.h"
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cfloat>


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_MULTITHERMAL_MULTIBARIC
/*
Multithermal-multibaric target distribution (dynamic).


*/
//+ENDPLUMEDOC

class TD_MultithermalMultibaric: public TargetDistribution {
private:
  double threshold_, min_temp_, max_temp_;
  double min_press_, max_press_, press_;
  std::vector<double> sigma_;
  unsigned steps_temp_, steps_pressure_;
  bool one_region_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_MultithermalMultibaric(const ActionOptions& ao);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_MultithermalMultibaric() {}
  double GaussianSwitchingFunc(const double, const double, const double) const;
};


PLUMED_REGISTER_ACTION(TD_MultithermalMultibaric,"TD_MULTITHERMAL_MULTIBARIC")


void TD_MultithermalMultibaric::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","THRESHOLD","Maximum exploration energy in KT.");
  keys.add("compulsory","MIN_TEMP","Minimum energy.");
  keys.add("compulsory","MAX_TEMP","Maximum energy.");
  keys.add("compulsory","MIN_PRESSURE","Minimum pressure.");
  keys.add("compulsory","MAX_PRESSURE","Maximum pressure.");
  keys.add("compulsory","PRESSURE","Maximum pressure.");
  keys.add("compulsory","STEPS_PRESSURE","Minimum pressure.");
  keys.add("compulsory","STEPS_TEMP","Maximum pressure.");
  keys.add("optional","SIGMA","The standard deviation parameters of the Gaussian kernels. You should give one value for each argument. Value of 0.0 means that switch is done without a smooth switching function, this is the default behavior.");
  keys.addFlag("ONE_REGION",false,"Make one region");
}


TD_MultithermalMultibaric::TD_MultithermalMultibaric(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  threshold_(10.0),
  min_temp_(0.0),
  max_temp_(1000.0),
  sigma_(0.0)
{
  log.printf("  Multithermic-multibaric target distribution");
  log.printf("\n");
  parse("THRESHOLD",threshold_);
  if(threshold_<=0.0) {
    plumed_merror("TD_MULTICANONICAL target distribution: the value of the threshold should be positive.");
  }
  parse("MIN_TEMP",min_temp_);
  parse("MAX_TEMP",max_temp_);
  parse("MIN_PRESSURE",min_press_);
  parse("MAX_PRESSURE",max_press_);
  parse("PRESSURE",press_);
  parseVector("SIGMA",sigma_);
  if(sigma_.size()!=2) {plumed_merror(getName()+": SIGMA takes 2 values as input, e.g. -1000,100");}
  parse("STEPS_TEMP",steps_temp_);
  parse("STEPS_PRESSURE",steps_pressure_);
  steps_temp_ += 1;
  steps_pressure_ += 1;

  one_region_=false;
  parseFlag("ONE_REGION",one_region_);

  setDimension(2);
  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_MultithermalMultibaric::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_MultithermalMultibaric");
  return 0.0;
}


void TD_MultithermalMultibaric::updateGrid() {
  if (getStep() == 0) {
    // Use uniform TD
    std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = 1.0;
      norm += integration_weights[l]*value;
      targetDistGrid().setValue(l,value);
    }
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logTargetDistGrid().setMinToZero();
  } else {
    double beta = getBeta();
    double beta_prime_min = 1./(plumed.getAtoms().getKBoltzmann()*min_temp_);
    double beta_prime_max = 1./(plumed.getAtoms().getKBoltzmann()*max_temp_);
    plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_MultithermalMultibaric!");
    // Set all to zero
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = 0.0;
      targetDistGrid().setValue(l,value);
    }
    // Loop over pressures and temperatures   
    for(unsigned i=0; i<steps_temp_; i++) {
      double beta_prime=beta_prime_min + (beta_prime_max-beta_prime_min)*i/(steps_temp_-1);
      for(unsigned j=0; j<steps_pressure_; j++) {
        double pressure_prime=min_press_ + (max_press_-min_press_)*j/(steps_pressure_-1);
        // Find minimum for this pressure and temperature
        double minval=DBL_MAX;
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double volume = targetDistGrid().getPoint(l)[1];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy + (beta_prime*pressure_prime-beta*press_)*volume;
          if(value<minval) {
            minval=value;
          }
        }
        // Now check which energies and volumes are below X kt
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double volume = targetDistGrid().getPoint(l)[1];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy + (beta_prime*pressure_prime-beta*press_)*volume - minval;
          if (value<threshold_) {
            double value = 1.0;
            targetDistGrid().setValue(l,value);
          }
        }
      }
    }
    std::vector<unsigned> nbin=targetDistGrid().getNbin();
    std::vector<double> dx=targetDistGrid().getDx();
    if (one_region_) {
      OFile grid1; grid1.link(*this);
      std::string grid1string="grid1.txt";
      grid1.open(grid1string);
      targetDistGrid().writeToFile(grid1);
      // Make one region
      for(unsigned i=0; i<nbin[0]; i++) {
        unsigned jmin=nbin[1]-1;
        unsigned jmax=0;
        // Find minimum and maximum j for which this i is 1 
        for(unsigned j=0; j<nbin[1]; j++) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          //Grid::index_t index = targetDistGrid().getIndex(indices);
          double value = targetDistGrid().getValue(indices);
          if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
            jmin=j;
            break;
          }
        }
        for(unsigned j=nbin[1]-1; j>0; j--) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          double value = targetDistGrid().getValue(indices);
          if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
            jmax=j;
            break;
          }
        }
        // Now make 1 all values in between
        for(unsigned j=jmin; j<jmax; j++) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          double value = 1.0;
          targetDistGrid().setValue(indices,value);
        }
      }
      for(unsigned j=0; j<nbin[1]; j++) {
        unsigned imin=nbin[0]-1;
        unsigned imax=0;
        // Find minimum and maximum j for which this i is 1 
        for(unsigned i=0; i<nbin[0]; i++) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          double value = targetDistGrid().getValue(indices);
          if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
            imin=i;
            break;
          }
        }
        for(unsigned i=nbin[0]-1; i>0; i--) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          double value = targetDistGrid().getValue(indices);
          if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
            imax=i;
            break;
          }
        }
        // Now make 1 all values in between
        for(unsigned i=imin; i<imax; i++) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          double value = 1.0;
          targetDistGrid().setValue(indices,value);
        }
      }
    }
    // Smoothening
    for(unsigned i=0; i<nbin[0]; i++) {
      for(unsigned j=0; j<nbin[1]; j++) {
        std::vector<unsigned> indices(2);
        indices[0]=i;
        indices[1]=j;
        Grid::index_t index = targetDistGrid().getIndex(indices);
        double energy = targetDistGrid().getPoint(index)[0];
        double volume = targetDistGrid().getPoint(index)[1];
        double value = targetDistGrid().getValue(index);
        if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
          // Apply gaussians around
          std::vector<int> minBin(2), maxBin(2), deltaBin(2); // These cannot be unsigned
          // Only consider contributions less than n*sigma bins apart from the actual distance
          deltaBin[0]=std::floor(5*sigma_[0]/dx[0]);;
          deltaBin[1]=std::floor(5*sigma_[1]/dx[1]);;
          // For energy
          minBin[0]=i - deltaBin[0];
          if (minBin[0] < 0) minBin[0]=0;
          if (minBin[0] > (nbin[0]-1)) minBin[0]=nbin[0]-1;
          maxBin[0]=i +  deltaBin[0];
          if (maxBin[0] > (nbin[0]-1)) maxBin[0]=nbin[0]-1;
          // For volume
          minBin[1]=j - deltaBin[1];
          if (minBin[1] < 0) minBin[1]=0;
          if (minBin[1] > (nbin[1]-1)) minBin[1]=nbin[1]-1;
          maxBin[1]=j +  deltaBin[1];
          if (maxBin[1] > (nbin[1]-1)) maxBin[1]=nbin[1]-1;
          for(unsigned l=minBin[0]; l<maxBin[0]+1; l++) {
            for(unsigned m=minBin[1]; m<maxBin[1]+1; m++) {
              std::vector<unsigned> indices_prime(2);
              indices_prime[0]=l;
              indices_prime[1]=m;
              Grid::index_t index_prime = targetDistGrid().getIndex(indices_prime);
              double energy_prime = targetDistGrid().getPoint(index_prime)[0];
              double volume_prime = targetDistGrid().getPoint(index_prime)[1];
              double value_prime = targetDistGrid().getValue(index_prime);
              // Apply gaussian
              double gaussian_value = GaussianSwitchingFunc(energy_prime,energy,sigma_[0])*GaussianSwitchingFunc(volume_prime,volume,sigma_[1]);
              if (value_prime<gaussian_value) {
                targetDistGrid().setValue(index_prime,gaussian_value);
              }
            }
          }
        }
      }
    }
    // Normalize
    std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = targetDistGrid().getValue(l);
      norm += integration_weights[l]*value;
    }
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logTargetDistGrid().setMinToZero();
  }
}

inline
double TD_MultithermalMultibaric::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
  if(sigma>0.0) {
    double arg=(argument-center)/sigma;
    return exp(-0.5*arg*arg);
  }
  else {
    return 0.0;
  }
}


}
}
