/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "ActionWithArguments.h"
#include "ActionWithValue.h"
#include "tools/PDB.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include <iostream>

using namespace std;
namespace PLMD{

void ActionWithArguments::registerKeywords(Keywords& keys){
  keys.reserve("compulsory","ARG","the input for this action is the scalar output from one or more other actions. The particular scalars that you will use "
                                  "are referenced using the label of the action. If the label appears on its own then it is assumed that the Action calculates "
                                  "a single scalar value.  The value of this scalar is thus used as the input to this new action.  If * or *.* appears the "
                                  "scalars calculated by all the proceding actions in the input file are taken.  Some actions have multi-component outputs and "
                                  "each component of the output has a specific label.  For example a \\ref DISTANCE action labelled dist may have three componets "
                                  "x, y and z.  To take just the x component you should use dist.x, if you wish to take all three components then use dist.*."
                                  "More information on the referencing of Actions can be found in the section of the manual on the PLUMED \\ref Syntax.  "
                                  "Scalar values can also be "
                                  "referenced using POSIX regular expressions as detailed in the section on \\ref Regex. To use this feature you you must compile "
                                  "PLUMED with the appropriate flag.");
}

void ActionWithArguments::parseArgumentList(const std::string&key,std::vector<Value*>&arg){
  vector<string> c; arg.clear(); parseVector(key,c); interpretArgumentList(c,arg);
}

void ActionWithArguments::interpretArgumentList(const std::vector<std::string>& c, std::vector<Value*>&arg){
  plumed.getActionSet().interpretArgumentList( c, arg, dynamic_cast<Action*>(this) );
}

void ActionWithArguments::expandArgKeywordInPDB( PDB& pdb ){
  std::vector<std::string> pdb_remark=pdb.getRemark();
  std::vector<std::string> arg_names;
  bool found=Tools::parseVector(pdb_remark,"ARG",arg_names);
  if( found ){ 
      std::vector<Value*> arg_vals;
      interpretArgumentList( arg_names, arg_vals );
      std::string new_args="ARG=" + arg_vals[0]->getName();
      for(unsigned i=1;i<arg_vals.size();++i) new_args = new_args + "," + arg_vals[i]->getName();
      pdb.setArgKeyword( new_args );
  }
}

void ActionWithArguments::requestArguments(const vector<Value*> &arg){
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  arguments=arg;
  clearDependencies();
  std::string fullname,name;
  for(unsigned i=0;i<arguments.size();i++){
     fullname=arguments[i]->getName();
     if(fullname.find(".")!=string::npos){
       std::size_t dot=fullname.find_first_of('.');
       name=fullname.substr(0,dot);
     } else {
       name=fullname;
     }
     ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(name);
     plumed_massert(action,"cannot find action named (in requestArguments - this is weird)" + name);
     addDependency(action);
  }
}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false)
{
  if( keywords.exists("ARG") ){
     vector<Value*> arg;
     parseArgumentList("ARG",arg);

     if(!arg.empty()){
       log.printf("  with arguments");
       for(unsigned i=0;i<arg.size();i++) log.printf(" %s",arg[i]->getName().c_str());
       log.printf("\n");
     }
     requestArguments(arg);
  }
}

void ActionWithArguments::calculateNumericalDerivatives( ActionWithValue* a ){
  if(!a){
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }

  const int nval=a->getNumberOfComponents();
  const int npar=arguments.size();
  std::vector<double> value (nval*npar);
  for(int i=0;i<npar;i++){
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+sqrt(epsilon));
    a->calculate();
    arguments[i]->set(arg0);
    for(int j=0;j<nval;j++){
      value[i*nval+j]=a->getOutputQuantity(j);
    }
  }
  a->calculate();
  a->clearDerivatives();
  for(int j=0;j<nval;j++){
    Value* v=a->copyOutput(j);
    if( v->getNumberOfDerivatives()>0 ) for(int i=0;i<npar;i++) v->addDerivative(i,(value[i*nval+j]-a->getOutputQuantity(j))/sqrt(epsilon));
  }
}

double ActionWithArguments::getProjection(unsigned i,unsigned j)const{
  plumed_massert(i<arguments.size()," making projections with an index which  is too large");
  plumed_massert(j<arguments.size()," making projections with an index which  is too large");
  const Value* v1=arguments[i];
  const Value* v2=arguments[j];
  return Value::projection(*v1,*v2);
}

void ActionWithArguments::addForcesOnArguments( const std::vector<double>& forces ){
  for(unsigned i=0;i<arguments.size();++i) arguments[i]->addForce( forces[i] );
}

}
