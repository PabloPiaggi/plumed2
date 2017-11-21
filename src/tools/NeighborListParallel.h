/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_tools_NeighborListParallel_h
#define __PLUMED_tools_NeighborListParallel_h

#include "Vector.h"
#include "AtomNumber.h"
#include "Communicator.h"
#include "Log.h"

#include <vector>

namespace PLMD {

class Pbc;

/// \ingroup TOOLBOX
/// A class that implements neighbor lists from two lists or a single list of atoms
class NeighborListParallel
{
  bool reduced;
  bool do_pair_,do_pbc_,twolists_,firsttime_;
  const PLMD::Pbc* pbc_;
  std::vector<PLMD::AtomNumber> fullatomlist_,requestlist_;
  std::vector<PLMD::Vector> positions_old_;
  std::vector<std::pair<unsigned,unsigned> > neighbors_;
  double distance_, skin_;
  int stride_;
  unsigned nlist0_,nlist1_,nallpairs_,lastupdate_;
  unsigned dangerousBuilds_, numberOfBuilds_;
  double maxLoadImbalance_,avgLoadImbalance_, avgTotalNeighbors_;
/// Return the pair of indexes in the positions array
/// of the two atoms forming the i-th pair among all possible pairs
  //std::pair<unsigned,unsigned> getIndexPair(unsigned i);
/// Communicator
  Communicator& mycomm;
/// Log
  Log& mylog;
public:
  NeighborListParallel(const std::vector<PLMD::AtomNumber>& list0,
               const std::vector<PLMD::AtomNumber>& list1,
               const bool& do_pair, const bool& do_pbc, const PLMD::Pbc& pbc, Communicator& cc,
               Log& log, const double& distance=1.0e+30, const int& stride=0, const double& skin=0.1);
  NeighborListParallel(const std::vector<PLMD::AtomNumber>& list0, const bool& do_pbc,
               const PLMD::Pbc& pbc, Communicator& cc, Log& log, const double& distance=1.0e+30,
               const int& stride=0, const double& skin=0.1);
/// Return the list of all atoms. These are needed to rebuild the neighbor list.
  std::vector<PLMD::AtomNumber>& getFullAtomList();
/// Check if the nieghbor list must be rebuilt
  bool isListStillGood(const std::vector<Vector>& positions);
/// Update the indexes in the neighbor list to match the
/// ordering in the new positions array
/// and return the new list of atoms that must be requested to the main code
//  std::vector<PLMD::AtomNumber>& getReducedAtomList();
/// Update the neighbor list and prepare the new
/// list of atoms that will be requested to the main code
  void update(const std::vector<PLMD::Vector>& positions);
/// Get the update stride of the neighbor list
  int getStride() const;
/// Get the last step in which the neighbor list was updated
  unsigned getLastUpdate() const;
/// Set the step of the last update
  void setLastUpdate(unsigned step);
/// Get the size of the neighbor list
  unsigned size() const;
/// Get the i-th pair of the neighbor list
  std::pair<unsigned,unsigned> getClosePair(unsigned i) const;
/// Get the list of neighbors of the i-th atom
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborListParallel() {}
/// Print statistics of neighbor list
  void printStats();
/// Gather statistics of neighbor list
  void gatherStats(const std::vector<PLMD::Vector>& positions);

};

}

#endif
