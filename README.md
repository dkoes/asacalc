This is a simple little program that uses Shake-Rupley and OpenBabel
to calculate the solvent accessible surface area of a molecule.

It only works for PDBs and takes a file as its only argument or reads 
from stdin.  The SASA is output in the occupancy field. 
The total is provided as a comment at the end.

If you give it something that isn't a PDB, bad things will probably happen.
