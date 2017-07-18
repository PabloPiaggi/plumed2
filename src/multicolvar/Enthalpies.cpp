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

//+PLUMEDOC MCOLVAR ENERGIES
/*
*/
//+ENDPLUMEDOC


class Enthalpies : 
   public MultiColvarBase
   {
private:
  double rcut2, pressure;
public:
  static void registerKeywords( Keywords& keys );
  explicit Enthalpies(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Enthalpies,"ENTHALPIES")

void Enthalpies::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","PRESSURE","0.","Pressure");
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
}

Enthalpies::Enthalpies(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
  // Parse the pressure
  parse("PRESSURE",pressure);
  log.printf("The pressure is %f . \n", pressure );

  // And setup the ActionWithVessel
  setLinkCellCutoff( 10.0 );
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
}

double Enthalpies::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   double volume=getBox().determinant();
   return getEnergyAtom(myatoms.getIndex(0)) + pressure*volume/getNumberOfAtoms();
}

}
}

