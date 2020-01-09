/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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

#include <string>
#include <cmath>
#include <stdio.h>

#define THREADS_PER_BLOCK 32

using namespace std;

namespace PLMD {
namespace gpucolvar {

//+PLUMEDOC COLVAR COORDINATIONNUMBERGPU
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

\plumedfile
# This should be a sample input.
t: COORDINATIONNUMBERGPU ATOMS=1,2
PRINT ARG=t STRIDE=100 FILE=COLVAR
\endplumedfile
<!---You should reference here the other actions used in this example--->
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class CoordinationNumberGPU : public Colvar {
  bool pbc;

public:
  explicit CoordinationNumberGPU(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(CoordinationNumberGPU,"COORDINATIONNUMBERGPU")

void CoordinationNumberGPU::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.addFlag("COORDINATIONNUMBERGPU_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
  keys.addFlag("COORDINATIONNUMBERGPU_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  //keys.add("compulsory","COORDINATIONNUMBERGPU_COMPULSORY","all compulsory keywords should be added like this with a description here");
  keys.add("optional","COORDINATIONNUMBERGPU_OPTIONAL","all optional keywords that have input should be added like a description here");
  keys.add("atoms","ATOMS","the keyword with which you specify what atoms to use should be added like this");
}

CoordinationNumberGPU::CoordinationNumberGPU(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  requestAtoms(atoms);
}


__device__ void pbcDistanceGPU(float *position1, float *position2, float *box, float *pbcDistance) {
  //float *pbcDistance = new float[3];
  pbcDistance[0] = position2[0]-position1[0];
  pbcDistance[1] = position2[1]-position1[1];
  pbcDistance[2] = position2[2]-position1[2];
  // Only for orthorombic boxes
  /*
  float *diag = new float[3]; float *hdiag = new float[3]; float *mdiag = new float[3];
  for(unsigned i=0; i<3; i++) {
    diag[i]=box[4*i];
    hdiag[i]=0.5*box[4*i];
    mdiag[i]=-0.5*box[4*i];
    printf(" i is %d; diag is %f; hdiag is %f; mdiag is %f; box is %f \n",i,diag[i],hdiag[i],mdiag[i],box[4*i]);
  }
  */
  for(unsigned i=0; i<3; i++) {
    while(pbcDistance[i]>0.5*box[4*i]) pbcDistance[i]-=box[4*i];
    while(pbcDistance[i]<=-0.5*box[4*i]) pbcDistance[i]+=box[4*i];
  }
  //return pbcDistance;
  //delete [] pbcDistance;
  //delete [] diag; delete [] hdiag; delete [] mdiag;
}

__global__ void GPUFunction(float *value, float *positions,int n, float *box) {
  *value=0.0;
  __shared__ float coordinationNumbers[THREADS_PER_BLOCK];
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    for (int j = 0; j < n; j++) {
      float *distance = new float[3];
      float *position1 = new float[3];
      float *position2 = new float[3];
      for (int k = 0; k < 3; k++) position1[k]=positions[3*i+k]; 
      for (int k = 0; k < 3; k++) position2[k]=positions[3*j+k]; 
      pbcDistanceGPU(position1,position2,box,distance);
      float distance2; 
      distance2 = distance[0]*distance[0];
      distance2 += distance[1]*distance[1];
      distance2 += distance[2]*distance[2];
      if (distance2 < 0.5) coordinationNumbers[threadIdx.x] += 1;
      delete [] distance; 
      delete [] position1; 
      delete [] position2; 
    }
  }

  __syncthreads(); 

  if( 0 == threadIdx.x ) {
    float sum = 0.0;
    for( int i = 0; i < THREADS_PER_BLOCK; i++ )
      sum += coordinationNumbers[i];
    atomicAdd( value , sum );
  }
  
  // Reset to 0 coordinationNumbers
  for (int i = 0; i < THREADS_PER_BLOCK; i++) {
    coordinationNumbers[threadIdx.x] = 0.0;
  }
}

// calculator
void CoordinationNumberGPU::calculate() {
  float *value;
  float *positions;
  float *box;
  cudaMallocManaged(&positions, getNumberOfAtoms()*3*sizeof(float));
  cudaMallocManaged(&box, 9*sizeof(float));
  cudaMallocManaged(&value, sizeof(float));

  for (unsigned i=0;i<getNumberOfAtoms();i++) {
    for (unsigned j=0;j<3;j++) {
      positions[3*i+j] = getPosition(i)[j];
    }
  }
  for (unsigned i=0;i<3;i++) {
    for (unsigned j=0;j<3;j++) {
      box[3*i+j] = getBox()[i][j];
    }
  }

  int N = getNumberOfAtoms();
  int numBlocks = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  // function to run on the gpu
  GPUFunction<<<numBlocks, THREADS_PER_BLOCK>>>(value, positions, N, box);
  //GPUFunction<<<1, 1>>>(positions,getNumberOfAtoms()*3);
  // kernel execution is asynchronous so sync on its completion
  cudaDeviceSynchronize();

  cudaFree(positions);
  cudaFree(box);

  /*
  Vector distance;
  if(pbc) {
    distance=pbcDistance(getPosition(0),getPosition(1));
  } else {
    distance=delta(getPosition(0),getPosition(1));
  }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  setAtomsDerivatives(0,-invvalue*distance);
  setAtomsDerivatives(1,invvalue*distance);
  setBoxDerivatives  (-invvalue*Tensor(distance,distance));
  */
  setValue           (*value);
  cudaFree(value);
}

}
}



