/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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

#include "Function.h"
#include "ActionRegister.h"
#include "tools/IFile.h"
#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION KL_DIVERGENCE
/*

*/
//+ENDPLUMEDOC

class KLdivergence :
  public Function
{
  std::vector<double> reference_;
  std::vector<double> weight_;
  double eps=1.e-10;
public:
  explicit KLdivergence(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(KLdivergence,"KL_DIVERGENCE")

void KLdivergence::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","REF","the reference values for caluclating the divergence");
  keys.add("optional","DOMAIN","the domain of integration to be used squared as weight");
  keys.add("optional","REF_FILE","the name of the file with the reference");
}

KLdivergence::KLdivergence(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  const unsigned narg=getNumberOfArguments();

  std::string referenceFileName;
  parse("REF_FILE",referenceFileName); 
  parseVector("REF",reference_);
  plumed_massert(!(!referenceFileName.empty() && reference_.size()>0),"Cannot use REF and REF_FILE keywords simultaneously");

  if (!referenceFileName.empty() ) {
    // Use file as reference
    log.printf("  Reading a reference from the file %s . \n", referenceFileName.c_str() );
    IFile ifile; 
    ifile.link(*this);
    ifile.open(referenceFileName);
    reference_.resize(narg);
    weight_.resize(narg);
    for(unsigned int i=0;i<narg;i++) {
      ifile.scanField("r",weight_[i]).scanField("ref",reference_[i]).scanField();
    }
  } else if (reference_.size()>0) {
    // Use values in REF as reference
    plumed_massert(reference_.size()==narg,"Size of REF array should be the same as number for arguments");
    parseVector("DOMAIN",weight_);
    if (weight_.size()==0)
      weight_.resize(narg,1);
    plumed_massert(weight_.size()==narg,"Size of DOMAIN array should be the same as number for arguments");
  }
  for (unsigned i=0; i<narg; i++)
    if (!(reference_[i]>0)) reference_[i]=eps;
  for (unsigned i=0; i<narg; i++)
    weight_[i]*=weight_[i];

  addValueWithDerivatives();
  checkRead();
  setNotPeriodic();
  log.printf("  with reference:");
  for(unsigned i=0; i<narg; i++) 
    log.printf("%f,",reference_[i]);
  log.printf("  with weights:");
  for(unsigned i=0; i<narg; i++) 
    log.printf("%f,",weight_[i]);
  log.printf("\n");
}

void KLdivergence::calculate() {
  double KL=0.0;
  for(unsigned i=0; i<reference_.size(); i++) 
  {
    const double arg_i=getArgument(i);
    const double log_ratio=std::log(arg_i/reference_[i]);
    KL+=weight_[i]*(arg_i*log_ratio-arg_i+reference_[i]);
    setDerivative(i,weight_[i]*log_ratio);
  }
  setValue(KL);
}

}
}
