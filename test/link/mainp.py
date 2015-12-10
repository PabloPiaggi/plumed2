#!/usr/bin/env python
import os
import sys
import time

import numpy as np
import plumed
from plumed import tools

# Read XYZ
filename= os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    os.pardir,
    os.pardir,
    "regtest", "crystallization", "rt-q6", "64.xyz")
atom_type, pos = tools.read_xyz(filename)

# Define a few things
step=0
box=np.diag(12.41642*np.ones(3,dtype=float))
virial=np.zeros((3,3),dtype=float)

# Repeat the call multiple times to see if 
# there are problems with multiple start/stop of the environment, and
# for more accurate timing. Remember to delete the bck* files before running!
num_loops = 1#50

# Create the class only once
plumed = plumed.Plumed()

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
    # TODO: write to memory
    # Avoid printing on screen
    #plumed.cmd("setLogFile","test.log") 
    # Disable plumed's output completely
    plumed.cmd("setLogFile","/dev/null")

    # Init
    plumed.cmd("init")

    # "action": new command by Pablo Piaggi to send directly commands 
    # that would be normally written in the input file.

    # Check if atoms were correctly received by Plumed.
    # plumed.cmd("action","DUMPATOMS ATOMS=1-64 FILE=testout.xyz")

    # Calculate Q6
    plumed.cmd("action",
               "Q6 SPECIES=1-32 D_0=3.0 R_0=1.5 MEAN MIN={BETA=0.1} LABEL=q6") 
    # Write on file
    plumed.cmd("action","PRINT ARG=q6.* FILE=colv")
    plumed.cmd("action","DUMPMULTICOLVAR DATA=q6 FILE=MULTICOLVAR.xyz")

    # Post-init settings
    plumed.cmd("setStep",step)
    plumed.cmd("setBox",box)
    plumed.cmd("setPositions",pos)
    plumed.cmd("setMasses", masses)
    plumed.cmd("setCharges", charges)
    plumed.cmd("setForces", forces)
    plumed.cmd("setVirial", virial)

    # Calculate
    plumed.cmd("calc")

    # Get values instead of reading files

    # Test with positions 
    #positions = plumed.grab("positions")
    #print 'outshape:', positions.shape
    #print 'outval:', positions
   
    # Grab single-valued actions and actions with components 
    outGrab = plumed.grab("q6.mean,q6.min,q6.*")
    print 'outGrab:', outGrab

    # Grab actions with vessels
    outVessels = plumed.grab("vessels q6")
    print 'outVessels:', outVessels

    # Stop the environment, delete 
    plumed.stop_plumed()

# End time
t2 = time.time()
dt = (t2-t1)*1000. # ms

print >>sys.stderr, "Total time: {:4.1f} ms".format(dt)
print >>sys.stderr, "Num loops:  {:4d}".format(num_loops)
print >>sys.stderr, "Time/loop:  {:4.1f} ms".format(dt/float(num_loops))

# TO DO: Write wrapper for actions. For instance "DISTANCE LABEL=l6 ATOMS=3,4" = Distance(label="l6",atoms)
# TO DO: Perhaps have two options. plumed.setMasses(masses) = plumed.cmd("setMasses",masses) . The first would be more pythonic...
