#!/usr/bin/env python
import os
import sys
import time

import numpy as np
from plumed import *
from tools import *

# Read XYZ
filename= os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    os.pardir,
    os.pardir,
    "regtest", "crystallization", "rt-q6", "64.xyz")
atom_type, pos = read_xyz(filename)

# Define a few things
step=0
box=np.diag(12.41642*np.ones(3,dtype=float))
virial=np.zeros((3,3),dtype=float)

# Repeat the call multiple times to see if 
# there are problems with multiple start/stop of the environment, and
# for more accurate timing. Remember to delete the bck* files before running!
num_loops = 1

# Create the class only once
plumed = Plumed()

# Start time
t1 = time.time()

for idx in range(num_loops):
    # Create appropriate arrays and values on the fly;
    # needed if the structure change at each loop (that is
    # the correct behavior)
    num_atoms=pos.shape[0]
    # Unfortunately, these need to be defined otherwise plumed complains
    masses=np.ones(num_atoms,dtype=float)
    forces=np.zeros((num_atoms,3),dtype=float)
    charges=np.zeros(num_atoms,dtype=float)
    
    # Start the environment
    plumed.start_plumed()

    # Pre-init settings 
    ## This one should be probably hiddend in the Plumed class
    plumed.cmd("setRealPrecision", 8) # float64

    plumed.cmd("setMDEngine","none")
    plumed.cmd("setTimestep", 1.) # Not used but must be defined
    plumed.cmd("setKbT", 1.) # Not used but must be defined
    plumed.cmd("setNatoms",num_atoms)
    plumed.cmd("setPlumedDat","") # Empty, will use the 'action' command
    # TODO: write to memory, or disable completely logging
    plumed.cmd("setLogFile","test.log") # To avoid printing on screen

    # Init
    plumed.cmd("init")
    # New command by Pablo Piaggi to send direclty commands that would be written
    # in the input file.
    # Check if atoms were correctly received by Plumed.
    # plumed.cmd("action", "DUMPATOMS ATOMS=1-64 FILE=testout.xyz")
    # Calculate Q6
    plumed.cmd("createAction","Q6 SPECIES=1-64 D_0=3.0 R_0=1.5 MEAN LABEL=q6")
    # plumed.cmd("action","PRINT ARG=q6.* FILE=colv")
    # Post-init settings
    plumed.cmd("setStep",step)
    plumed.cmd("setBox",box)
    plumed.cmd("setPositions",pos)
    plumed.cmd("setMasses", masses)
    plumed.cmd("setCharges", charges)
    plumed.cmd("setForces", forces)
    plumed.cmd("setVirial", virial)
    # Calculate
    # plumed.cmd("calc")

    # TODO: here get values instead of reading files
    positions = plumed.grab("q6.mean")
    print 'outshape:', positions.shape
    print 'outval:', positions
    q6data = plumed.grab("q6")
    print 'outq6s:', q6data.shape
    print 'outq6:', q6data 
   
    # Stop the environment, delete 
    plumed.stop_plumed()

# End time
t2 = time.time()
dt = (t2-t1)*1000. # ms

print >>sys.stderr, "Total time: {:4.1f} ms".format(dt)
print >>sys.stderr, "Num loops:  {:4d}".format(num_loops)
print >>sys.stderr, "Time/loop:  {:4.1f} ms".format(dt/float(num_loops))

# TO DO: Write wrapper for actions. For instance "DISTANCE LABEL=l6 ATOMS=3,4" = Distance(label="l6",atoms)
# TO DO: Wrappers to get values back from actions. Very important!
# TO DO: Perhaps have two options. plumed.setMasses(masses) = plumed.cmd("setMasses",masses) . The first would be more pythonic...
# IDEAS AND TO DO: Retrieve results of calculation directly in plumed. Check how forces, etc. are passed back to the md code (I think that they are pointers). Remember that the bias is passed back to gromacs for replica exchange so we can also passed it to python (perhaps do the same for CVs?). Pass things back according to a label. Check how to delete actions.
