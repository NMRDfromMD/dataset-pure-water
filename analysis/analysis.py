
import os, time
import numpy as np
import MDAnalysis as mda
from nmrdfrommd import NMRD

from utilities import save_result, get_git_repo_path

git_path = get_git_repo_path()
data_dir = os.path.join(git_path, "data")
topology_file = os.path.join(data_dir, "equilibrate.data")
trajectory_file = os.path.join(data_dir, "production.xtc")

u = mda.Universe(topology_file, trajectory_file)

n_TOT = u.atoms.n_residues
n_H2O = u.select_atoms("type 1 2").n_residues

print(f"The total number of H2O molecules is {n_H2O}")

timestep = np.int32(u.trajectory.dt)
print(f"The timestep is {timestep} ps")
total_time = np.int32(u.trajectory.totaltime)
print(f"The total simulation time is {total_time//1000} ns")

H_H2O = u.select_atoms("type 2")

for n in [5, 20, 80, 320, 1280]:

    ti = time.time()

    nmr_intra_H2O = NMRD(
        u=u,
        atom_group=H_H2O,
        type_analysis="full",
        number_i=n)
    nmr_intra_H2O.run_analysis()
    save_result(nmr_intra_H2O, n, f"nmr_full_H2O")

    tf = time.time()
    print("time", np.round(tf-ti,1))