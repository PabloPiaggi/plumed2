/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_TEMP
/*
Calculate weights for ensemble averages allow for the computing of ensemble averages at temperatures lower/higher than that used in your original simulation.

We can use our knowledge of the Boltzmann distribution in the cannonical ensemble to reweight the data
contained in trajectories.  Using this procedure we can take trajectory at temperature \f$T_1\f$ and use it to
extract probabilities at a different temperature, \f$T_2\f$, using:

\f[
P(s',t) = \frac{ \sum_{t'}^t \delta( s(x) - s' ) \exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }{ \sum_{t'}^t \exp\left( +\left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }
\f]

The weights calculated by this action are equal to \f$\exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right)\f$ and take
the effect the bias has on the system into account.  These weights can be used in any action
that computes ensemble averages.  For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE.

\par Examples

The following input can be used to postprocess a molecular dynamics trajectory calculated at a temperature of 500 K.
The \ref HISTOGRAM as a function of the distance between atoms 1 and 2 that would have been obtained if the simulation
had been run at the lower temperature of 300 K is estimated using the data from the higher temperature trajectory and output
to a file.

\plumedfile
x: DISTANCE ATOMS=1,2
aa: REWEIGHT_TEMP TEMP=500 REWEIGHT_TEMP=300
hB: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 LOGWEIGHTS=aa
DUMPGRID GRID=hB FILE=histoB
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightPressure : public ReweightBase {
private:
///
  double rpress_, press_;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightPressure(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightPressure,"REWEIGHT_PRESSURE")

void ReweightPressure::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.add("compulsory","REWEIGHT_PRESSURE","reweight data from a trajectory at one pressure and output the probability "
           "distribution at a second temperature. This is not possible during postprocessing.");
  keys.add("compulsory","PRESSURE","Pressure of the trajectory");
  keys.remove("ARG"); keys.add("compulsory","ARG","Volume");
}

ReweightPressure::ReweightPressure(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  parse("REWEIGHT_PRESSURE",rpress_);
  parse("PRESSURE",press_);
  //rtemp*=plumed.getAtoms().getKBoltzmann();
  if(getNumberOfArguments()>1) error(" Number of arguments should be 1. Only the volume should be given.");
}

double ReweightPressure::getLogWeight() {
  double volume=getArgument(0);
  return - (1./simtemp) * ( rpress_ - press_ )*volume;
}

}
}
