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

//+PLUMEDOC MCOLVAR TETRAHEDRAL ORDER PARAMETER
/*
Calculate the degree of tethaedrality around a central atom, by computing the order parameter proposed by: 
'A new order parameter for tetrahedral configurations' (Chau, Hardwick, Molecular Physics 1998)
rescaled as done in 'Relationship between structural order and the anomalies of liquid water' (Errington, Debenedetti, Nature 2001).

\f[
s = 1 - \frac{3}{32} \sum_{i=0} ^{N-1} \sum_{j=i+1} ^{N} (\cos(\theta_{ikj}+1/3))^2
\f]

Then you can calculate functions of the distribution of this local tethraedral order parameter such as
the mean, minimum, the number less than a certain quantity and so on.

\par Examples

The following input tells plumed to calculate the tethraedral paramter of atoms 1-100 with themselves.
The mean of these parameter is then calculated:
\verbatim
tb: TETRABRICK SPECIES=1-216 D_0=0.25 R_0=0.05 MEAN
\endverbatim

*/
//+ENDPLUMEDOC


class Tetrabrick : public MultiColvarBase {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction; 
  bool entropyFlag, normalizeFlag;
public:
  static void registerKeywords( Keywords& keys );
  explicit Tetrabrick(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Tetrabrick,"TETRABRICK")

void Tetrabrick::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("ENTROPY",false,"Calculate the tetrahedral entropy instead of the entropy order parameter.");
  keys.addFlag("NORMALIZE_NUM_NEIGH",false,"Normalize the order parameter with the number of neighbors.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

Tetrabrick::Tetrabrick(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();

  parseFlag("ENTROPY",entropyFlag);
  if (entropyFlag) log.printf("The tetrahedral entropy will be calculated instead of the tetrahedral order paramater. Derivatives are not available. \n");

  parseFlag("NORMALIZE_NUM_NEIGH",normalizeFlag);
  if (normalizeFlag) log.printf("The tetrahedral parameter will be normalized by the number of neighbors. \n");

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
}

double Tetrabrick::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   // --- Calculate the tetrahedral order parameter ---

   // Define output quantities
   double tetra=0;
   vector<Vector> deriv(getNumberOfAtoms());
   Tensor virial;
   // Define temp quantities
   double d2i, d2j;								//square distances
   double di_mod, dj_mod;						//modulo of distances
   double di_inv, dj_inv;						//inverse of distances
   double cos, cos2, versors; 							//cosine term, and squared version, and scalar product between versors
   double sw_i, sw_j, df_i, df_j;				//switch functions and derivatives
   double prefactor;							//prefactor of the derivative of the cosine part 
   Vector der_sw_i, der_sw_j;					//switch derivative parts
   Vector der_i, der_j;							//derivatives i and j

   // Loop on nearest neighbors and load distances di and dj
	double numNeigh;

   for(unsigned i=1;i<(myatoms.getNumberOfAtoms()-1);++i){
      Vector& di=myatoms.getPosition(i); 		//relative position of atom i (with respect to k)
      if ( (d2i=di[0]*di[0])<rcut2 &&
           (d2i+=di[1]*di[1])<rcut2 &&
           (d2i+=di[2]*di[2])<rcut2) {
         for(unsigned j=i+1;j<myatoms.getNumberOfAtoms();++j){
            Vector& dj=myatoms.getPosition(j); 	//relative position of atom j (with respect to k)
            if ( (d2j=dj[0]*dj[0])<rcut2 &&
                 (d2j+=dj[1]*dj[1])<rcut2 &&
                 (d2j+=dj[2]*dj[2])<rcut2) {
					 
			// -- (1) COMPUTE --
			// useful terms
			di_mod=sqrt(d2i);
			dj_mod=sqrt(d2j);
			di_inv=1./di_mod;
			dj_inv=1./dj_mod;
			// cosine term
			versors = dotProduct( di,dj )/( di_mod * dj_mod );
			cos = versors + 1./3.;
			cos2 = cos*cos;
			// compute switching functions (derivatives return der/dist)
			sw_i = switchingFunction.calculateSqr( d2i, df_i );
			sw_j = switchingFunction.calculateSqr( d2j, df_j );
			// order parameter (relative to ikj triplet)
            tetra += cos2*sw_i*sw_j;
			numNeigh += sw_i*sw_j;

			// -- (2) DERIVATE --
			// useful terms
			der_sw_i = di * df_i * sw_j * cos2;		// sw part relative to i derivative
			der_sw_j = sw_i * dj * df_j * cos2;		// sw part relative to j derivative
			prefactor = 2*sw_i*sw_j*cos;			// der prefactor
      
			//deriv[0] -= der_sw_i + prefactor * ( di_inv*dj_inv - versors/d2i) * di + der_sw_j + prefactor * ( di_inv*dj_inv - versors/d2j) * dj;
			der_i = der_sw_i + prefactor * ( di_inv*dj_inv * dj - (versors/d2i)* di);
			der_j = der_sw_j + prefactor * ( di_inv*dj_inv * di - (versors/d2j)* dj);
			
			deriv[i] += der_i;	
			deriv[j] += der_j;	
			
			// -- (3) VIRIAL --
			 myatoms.addBoxDerivatives( 1, -(Tensor(di,-3./8.*der_i)+Tensor(dj,-3./8.*der_j)) );		
			}
         }
      }
   }

/* OLD VERSION!
   for(unsigned i=1;i<(myatoms.getNumberOfAtoms()-1);++i){
      Vector& di=myatoms.getPosition(i); 		//relative position of atom i (with respect to k)
      if ( (d2i=di[0]*di[0])<rcut2 &&
           (d2i+=di[1]*di[1])<rcut2 &&
           (d2i+=di[2]*di[2])<rcut2) {
         for(unsigned j=i+1;j<myatoms.getNumberOfAtoms();++j){
            Vector& dj=myatoms.getPosition(j); 	//relative position of atom j (with respect to k)
            if ( (d2j=dj[0]*dj[0])<rcut2 &&
                 (d2j+=dj[1]*dj[1])<rcut2 &&
                 (d2j+=dj[2]*dj[2])<rcut2) {

			// -- (1) COMPUTE --
			// useful terms
			di_mod=sqrt(d2i);
			dj_mod=sqrt(d2j);
			di_inv=1./di_mod;
			dj_inv=1./dj_mod;
			// cosine term
			cos = dotProduct( di,dj )/( di_mod * dj_mod ) +1./3.;
			cos2=cos*cos;
			// order parameter (relative to ikj triplet)
            tetra += cos2;

			// -- (2) DERIVATE --
			// useful terms
			prefactor = 2*cos;			// der prefactor
      
			for(unsigned k=0;k<3;++k){
				c_i[k] = ( di[k]*di[k] * di_inv*di_inv -1)*di_inv *dj[k]*dj_inv; //!!!
				c_j[k] = ( dj[k]*dj[k] * dj_inv*dj_inv -1)*dj_inv *di[k]*di_inv;
			}
			// central atom
			deriv[0] += prefactor *(c_i + c_j); 		//derivative with respect to central atom has a minus in both the sw terms
			deriv[i] -= prefactor * c_i;							//the derivative part of the cosine has a minus with respect to the central one
			deriv[j] -= prefactor * c_j;				
			
			// -- (3) VIRIAL --
			//Tensor vv(value, distance);	
			//myatoms.addBoxDerivatives( 1, virial );		
			}
         }
      }
   }
*/
/*   DOUBLE SUM WITH SWs
 * 	for(unsigned i=1;i<(myatoms.getNumberOfAtoms()-1);++i){
      Vector& di=myatoms.getPosition(i); 		//relative position of atom i (with respect to k)
      if ( (d2i=di[0]*di[0])<rcut2 &&
           (d2i+=di[1]*di[1])<rcut2 &&
           (d2i+=di[2]*di[2])<rcut2) {
         for(unsigned j=i+1;j<myatoms.getNumberOfAtoms();++j){
            Vector& dj=myatoms.getPosition(j); 	//relative position of atom j (with respect to k)
            if ( (d2j=dj[0]*dj[0])<rcut2 &&
                 (d2j+=dj[1]*dj[1])<rcut2 &&
                 (d2j+=dj[2]*dj[2])<rcut2) {

			// -- (1) COMPUTE --
			// useful terms
			di_mod=sqrt(d2i);
			dj_mod=sqrt(d2j);
			di_inv=1./di_mod;
			dj_inv=1./dj_mod;
			// cosine term
			cos = dotProduct( di,dj )/( di_mod * dj_mod ) +1./3.;
			cos2=cos*cos;
			// compute switching functions (derivatives return der/dist)
			sw_i = switchingFunction.calculateSqr( d2i, df_i );
			sw_j = switchingFunction.calculateSqr( d2j, df_j );
			// order parameter (relative to ikj triplet)
            tetra += sw_i*sw_j;

			// -- (2) DERIVATE --
			// useful terms
			der_sw_i = di * df_i * sw_j ;		// sw part relative to i derivative
			der_sw_j = sw_i * dj * df_j ;		// sw part relative to j derivative
			prefactor = 2*sw_i*sw_j*cos;		// der prefactor
      
			for(unsigned k=0;k<3;++k){
				c_i[k] = ( di[k]*di[k] * di_inv*di_inv -1)*di_inv *dj[k];
				c_j[k] = ( dj[k]*dj[k] * dj_inv*dj_inv -1)*dj_inv *di[k];
			}
			// central atom
			deriv[0] -= der_sw_i + der_sw_j; 		//derivative with respect to central atom has a minus in both the sw terms
			deriv[i] += der_sw_i;							//the derivative part of the cosine has a minus with respect to the central one
			deriv[j] += der_sw_j;				
			
			// -- (3) VIRIAL --
			//Tensor vv(value, distance);	
			//myatoms.addBoxDerivatives( 1, virial );		
			}
         }
      }
   }
*/
/*	COORDINATION NUMBER
 * for(unsigned i=1;i<(myatoms.getNumberOfAtoms()-1);++i){
      Vector& di=myatoms.getPosition(i); 		//relative position of atom i (with respect to k)
      if ( (d2i=di[0]*di[0])<rcut2 &&
           (d2i+=di[1]*di[1])<rcut2 &&
           (d2i+=di[2]*di[2])<rcut2) {
        
			// -- (1) COMPUTE --

			// compute switching functions (derivatives return der/dist)
			sw_i = switchingFunction.calculateSqr( d2i, df_i );
			// order parameter (relative to ikj triplet)
            tetra += sw_i;

			// -- (2) DERIVATE --
			// useful terms
			der_sw_i = di * df_i ;		// sw part relative to i derivative
      
			// central atom
			deriv[0] -= der_sw_i; 		//derivative with respect to central atom has a minus in both the sw terms
			deriv[i] += der_sw_i;		//the derivative part of the cosine has a minus with respect to the central one	

      }
   }
*/

   // output quantities
   if (tetra<1.e-10) tetra=1.e-10;
	tetra /= numNeigh;
   tetra = 1 - 3./8.*tetra;
   if (entropyFlag) tetra = (3./2.)*std::log(1-tetra);

   // Assign derivatives
   for(unsigned i=0;i<myatoms.getNumberOfAtoms();++i){
	  deriv[0]+=3./8.*deriv[i];
	  addAtomDerivatives(1,i,-3./8.*deriv[i],myatoms);
   }
   addAtomDerivatives(1,0,deriv[0],myatoms);
   
   // Assign virial	
   //myatoms.addBoxDerivatives( 1, virial );
   
   return tetra;
}

}
}
