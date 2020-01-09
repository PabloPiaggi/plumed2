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
#include "GpuColvar.h"
#include "core/Colvar.h"

#include <string>
#include <cmath>
#include <stdio.h>

#define THREADS_PER_BLOCK 32

using namespace std;

namespace PLMD {
namespace gpucolvar {

GpuColvar::GpuColvar(const ActionOptions&ao) :
  Action(ao),
  Colvar(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("SPECIES",atoms);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  requestAtoms(atoms);
}

void GpuColvar::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","SPECIES","the keyword with which you specify what atoms to use should be added like this");
}


__host__ __device__ void GpuColvar::pbcDistanceGPU(float *position1, float *position2, float *box, float *pbcDistance) {
  pbcDistance[0] = position2[0]-position1[0];
  pbcDistance[1] = position2[1]-position1[1];
  pbcDistance[2] = position2[2]-position1[2];
  // Only for orthorombic boxes
  for(unsigned i=0; i<3; i++) {
    while(pbcDistance[i]>0.5*box[4*i]) pbcDistance[i]-=box[4*i];
    while(pbcDistance[i]<=-0.5*box[4*i]) pbcDistance[i]+=box[4*i];
  }
}

// See the use of many files here:
// https://devblogs.nvidia.com/separate-compilation-linking-cuda-device-code/
// https://devtalk.nvidia.com/default/topic/1014106/separate-compilation-of-cuda-code-into-library-for-use-with-existing-code-base/
// Maybe pass the class GpuColvar as argument
/*
__global__ void GPUFunction(float *value, float *positions,int n, float *box) {
  *value=0.0;
  //__shared__ float coordinationNumbers[THREADS_PER_BLOCK];
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
      printf("%f %f %f \n",distance[0],distance[1],distance[2]);
    }
  }
}
*/
__global__ void GPUFunction(Bundle *mysystem) {
  // This function should call inside the loop a virtual function
  // that is defined in the derived class.
  //mygpucolvar.*value=1.;
  //if( 0 == threadIdx.x ) {
    //if( 0 == blockIdx.x ) {
      float* value = new float(1);
      value[0]=1.0;
      mysystem->setValue(value);
      delete [] value;
    //}
  //}
  /*
  //__shared__ float coordinationNumbers[THREADS_PER_BLOCK];
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
      printf("%f %f %f \n",distance[0],distance[1],distance[2]);
    }
  }
  */
}

void GpuColvar::calculate() {
  cudaMallocManaged(&positions, getNumberOfAtoms()*3*sizeof(float));
  cudaMallocManaged(&box, 9*sizeof(float));
 // cudaMallocManaged(&value, sizeof(float));
  Bundle *devsystem;
  cudaMallocManaged(&devsystem, sizeof(Bundle));
  devsystem->initValue();

  for (unsigned i=0;i<getNumberOfAtoms();i++) {
    for (unsigned j=0;j<3;j++) {
      positions[3*i+j] = Colvar::getPosition(i)[j];
    }
  }
  for (unsigned i=0;i<3;i++) {
    for (unsigned j=0;j<3;j++) {
      box[3*i+j] = Colvar::getBox()[i][j];
    }
  }

  int N = getNumberOfAtoms();
  int numBlocks = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  // function to run on the gpu
  GPUFunction<<<1, 1>>>(devsystem);
  //GPUFunction<<<numBlocks, THREADS_PER_BLOCK>>>(this);
  //GPUFunction<<<numBlocks, THREADS_PER_BLOCK>>>(value, positions, N, box);
  cudaDeviceSynchronize();

  cudaFree(positions);
  cudaFree(box);

  Colvar::setValue           (devsystem->getValue());
  cudaFree(devsystem);
}

/*
GpuColvar::distances(thrust::device_vector<float> *positions1, thrust::device_vector<float> *positions2) {
 int N=positions1.size()*positions2.size()/3 ;
 thrust::device_vector<float> d(N);
 
 return d;
}
*/
}
}
