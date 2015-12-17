import numpy as np

# Read xyz

def read_xyz(filename):
   # Based on http://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/03/IntroductiontoPython/p03_instructor.html and change to numpy style.
   # It only reads one frame.
   counter = 0
   xyz = open(filename)
   n_atoms = int(xyz.readline())
   atom_type = np.zeros(n_atoms).astype(str)
   coordinates = np.zeros([n_atoms,3])
   title = xyz.readline()
   for line in xyz:
       if ( (counter+1) > n_atoms ):
		raise ValueError("File says %d atoms but the number of atom lines is greater." % (n_atoms))
       atom,x,y,z = line.split()
       atom_type[counter]=atom
       coordinates[counter,:]=np.array([float(x),float(y),float(z)])
       counter += 1
   xyz.close()
   return atom_type, coordinates
