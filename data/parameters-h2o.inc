# LAMMPS Input Script
# LAMMPS version: 30Mar2026
# Author: Simon Gravelle
# License: GNU General Public License v3.0
# Project: NMRDfromMD – https://nmrdfrommd.github.io
# Description: Parameter file for simulating water

# Mass
mass 		1 15.9994 # H2O O
mass 		2 1.008 # H2O H

# Pair Coeff
pair_coeff 	1 1 0.185199 3.1589
pair_coeff 	2 2 0.0 0.0

# Bond coeff
bond_coeff 	1 0 0.9572

# Angle coeff
angle_coeff     1 0 104.52
