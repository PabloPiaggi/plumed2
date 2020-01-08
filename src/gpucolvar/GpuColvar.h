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

#ifndef __PLUMED_gpucolvar_GpuColvar_h
#define __PLUMED_gpucolvar_GpuColvar_h

#include "core/Colvar.h"

#include <string>
#include <cmath>
#include <stdio.h>
#include <thrust/device_vector.h>

#define PLUMED_GPUCOLVAR_GPUCOLVAR_INIT(ao) GpuColvar(ao)

namespace PLMD {
namespace gpucolvar {

class GpuColvar : public Colvar {
  bool pbc;
  float *positions;
  float *box;
  float *value;
public:
  explicit GpuColvar(const ActionOptions&);
// active methods:
  float* getPositions();
  __host__ __device__ void setValue(float *myvalue);
  void calculate() override;
  //thrust::device_vector<float> getPosition(int i);
  //thrust::device_vector<float> getBox();
  __host__ __device__ void pbcDistanceGPU(float *position1, float *position2, float *box, float *pbcDistance);
};

inline
float* GpuColvar::getPositions() {
  float* pos = new float(3*getNumberOfAtoms());
  for (unsigned i=0;i<getNumberOfAtoms();i++) {
    for (unsigned j=0;j<3;j++) {
      pos[3*i+j] = Colvar::getPosition(i)[j];
    }
  }
  return pos;
}

inline
__host__ __device__ void GpuColvar::setValue(float *myvalue) {
  value = myvalue;
}


/*
inline
thrust::device_vector<float> GpuColvar::getPosition(int i) {
  thrust::device_vector<float> pos(3);
  for (unsigned j=0;j<3;j++) {
    pos[j] = Colvar::getPosition(i)[j];
  }
  return pos;
}

inline
thrust::device_vector<float> GpuColvar::getBox() {
  thrust::device_vector<float> box(9);
  for (unsigned i=0;i<3;i++) {
    for (unsigned j=0;j<3;j++) {
      box[3*i+j] = Colvar::getBox()[i][j];
    }
  }
  return box;
}
*/


}
}
#endif
