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

//+PLUMEDOC VES_TARGETDIST TD_MULTICANONICAL
/*
Multicanonical target distribution (dynamic).


*/
//+ENDPLUMEDOC

class TD_Multicanonical: public TargetDistribution {
private:
  double threshold_, min_temp_, max_temp_;
  double sigma_min_, sigma_max_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Multicanonical(const ActionOptions& ao);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_Multicanonical() {}
  double GaussianSwitchingFunc(const double, const double, const double) const;
};


PLUMED_REGISTER_ACTION(TD_Multicanonical,"TD_MULTICANONICAL")


void TD_Multicanonical::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","THRESHOLD","Maximum exploration energy in KT.");
  keys.add("compulsory","MIN_TEMP","Minimum energy of the multicanonical ensemble.");
  keys.add("compulsory","MAX_TEMP","Maximum energy of the multicanonical ensemble.");
  keys.add("optional","SIGMA_MINIMA","The standard deviation parameters of the Gaussian switching functions for the minima of the intervals. You should give one value for each argument. Value of 0.0 means that switch is done without a smooth switching function, this is the default behavior.");
  keys.add("optional","SIGMA_MAXIMA","The standard deviation parameters of the Gaussian switching functions for the maxima of the intervals. You should give one value for each argument. Value of 0.0 means that switch is done without a smooth switching function, this is the default behavior.");

}


TD_Multicanonical::TD_Multicanonical(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  threshold_(10.0),
  min_temp_(0.0),
  max_temp_(1000.0)
{
  log.printf("  Multicanonical target distribution");
  log.printf("\n");
  parse("THRESHOLD",threshold_);
  if(threshold_<=0.0) {
    plumed_merror("TD_MULTICANONICAL target distribution: the value of the threshold should be positive.");
  }
  parse("MIN_TEMP",min_temp_);
  parse("MAX_TEMP",max_temp_);
  parse("SIGMA_MINIMA",sigma_min_);
  parse("SIGMA_MAXIMA",sigma_max_);
  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_Multicanonical::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_Multicanonical");
  return 0.0;
}


void TD_Multicanonical::updateGrid() {
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
    plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_Multicanonical!");
 
    // Find minimum of F(U) at temperature min
    double minval=DBL_MAX;
    Grid::index_t minindex = (targetDistGrid().getSize())/2;
    double minpos = targetDistGrid().getPoint(minindex)[0];
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = getFesGridPntr()->getValue(l);
      double argument = targetDistGrid().getPoint(l)[0];
      value = beta*value + (beta_prime_min-beta)*argument;
      //log.printf(" l arg val %d %f %f \n",l,argument,value);
      if(value<minval) {
        minval=value;
        minpos=argument;
        minindex=l;
      }
    }
 
    //log.printf("  minval %f minpos %f minindex %d \n",minval,minpos,minindex);

    // Find minimum energy at low temperature
    double minimum_low = minpos;
    for(Grid::index_t l=minindex; l>1; l-=1) {
      double argument = targetDistGrid().getPoint(l)[0];
      double argument_next = targetDistGrid().getPoint(l-1)[0];
      double value = getFesGridPntr()->getValue(l);
      double value_next = getFesGridPntr()->getValue(l-1);
      value = beta*value + (beta_prime_min-beta)*argument - minval;
      value_next = beta*value_next + (beta_prime_min-beta)*argument_next - minval;
      if (value<threshold_ && value_next>threshold_) {
        minimum_low = argument_next;
        break;
      }
    }
  
    // Find maximum energy at low temperature
    double maximum_low = minpos;
    for(Grid::index_t l=minindex; l<(targetDistGrid().getSize()-1); l++) {
      double argument = targetDistGrid().getPoint(l)[0];
      double argument_next = targetDistGrid().getPoint(l+1)[0];
      double value = getFesGridPntr()->getValue(l);
      double value_next = getFesGridPntr()->getValue(l+1);
      value = beta*value + (beta_prime_min-beta)*argument - minval;
      value_next = beta*value_next + (beta_prime_min-beta)*argument_next - minval;
      if (value<threshold_ && value_next>threshold_) {
        maximum_low = argument_next;
        break;
      }
    }
 
    // Find minimum of F(U) at temperature max
    //minindex = (targetDistGrid().getSize())/2;
    //minpos = targetDistGrid().getPoint(minindex)[0];
    minval=DBL_MAX;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = getFesGridPntr()->getValue(l);
      double argument = targetDistGrid().getPoint(l)[0];
      value = beta*value + (beta_prime_max-beta)*argument;
      //log.printf(" l arg val %d %f %f \n",l,argument,value);
      if(value<minval) {
        minval=value;
        minpos=argument;
        minindex=l;
      }
    }

    //log.printf("  minval %f minpos %f minindex %d \n",minval,minpos,minindex);
 
    // Find minimum energy at high temperature
    double minimum_high = minpos;
    for(Grid::index_t l=minindex; l>1; l-=1) {
      double argument = targetDistGrid().getPoint(l)[0];
      double argument_next = targetDistGrid().getPoint(l-1)[0];
      double value = getFesGridPntr()->getValue(l);
      double value_next = getFesGridPntr()->getValue(l-1);
      value = beta*value + (beta_prime_max-beta)*argument - minval;
      value_next = beta*value_next + (beta_prime_max-beta)*argument_next - minval;
      if (value<threshold_ && value_next>threshold_) {
        minimum_high = argument_next;
        break;
      }
    }
  
    // Find maximum energy at high temperature
    double maximum_high = minpos;
    for(Grid::index_t l=minindex; l<(targetDistGrid().getSize()-1); l++) {
      double argument = targetDistGrid().getPoint(l)[0];
      double argument_next = targetDistGrid().getPoint(l+1)[0];
      double value = getFesGridPntr()->getValue(l);
      double value_next = getFesGridPntr()->getValue(l+1);
      value = beta*value + (beta_prime_max-beta)*argument - minval;
      value_next = beta*value_next + (beta_prime_max-beta)*argument_next - minval;
      if (value<threshold_ && value_next>threshold_) {
        maximum_high = argument_next;
        break;
      }
    }

    double minimum = minimum_low;
    if (minimum_high<minimum_low) minimum=minimum_high;

    double maximum = maximum_low;
    if (maximum_high>maximum_low) maximum=maximum_high;
  
    //log.printf("  Minimum H %f L %f final %f \n", minimum_high, minimum_low, minimum);
    //log.printf("  Maximum H %f L %f final %f \n", maximum_high, maximum_low, maximum);

    // Construct uniform TD in the interval between minimum and maximum
    std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double argument = targetDistGrid().getPoint(l)[0];
      double value = 1.0;
      double tmp;
      if(argument < minimum) {
        tmp = GaussianSwitchingFunc(argument,minimum,sigma_min_);
      }
      else if(argument > maximum) {
        tmp = GaussianSwitchingFunc(argument,maximum,sigma_max_);
      }
      else {
        tmp = 1.0;
      }
      value *= tmp;
      //log.printf("  Argument %f \n", argument);
      norm += integration_weights[l]*value;
      targetDistGrid().setValue(l,value);
    }
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logTargetDistGrid().setMinToZero();
  }
}

inline
double TD_Multicanonical::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
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
