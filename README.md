##### This is small python script to calculate distance matrix from the lmpdat file (LAMMPS).

The idea is to create matrix which contains all paired distances between chain monomers. The numbers of chain beads should be a row.

Input parameters:  matrix size (beads), path to lmp file, path to output matrix and path to output mol2 file. Matrix size can be smaller than the length of the chain. In this case, the matrix will be averaged on set of segments. For example, if length of the polymer is 10 000 beads and you set matrix size as 1 000. So the polymer chain will be splited in 10 segments of 1 000 each and the final result will be averaged on 10 matrices.

In the case of using this code, please, cite the repository.
