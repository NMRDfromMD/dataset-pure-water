"""This Python script creates a box water molecules"""

from utilities import place_molecules, write_lammps, prepare_lammps

# fix the number of water
Number_water = 800

atoms, atoName, resName, Lx, Ly, Lz = place_molecules(Number_water)

# Convert to LAMMPS
atoms, bonds, angles = prepare_lammps(atoms)
print('The total number of atoms is', str(len(atoms)))
write_lammps(atoms, bonds, angles, Lx, Ly, Lz)
