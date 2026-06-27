
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

for n, repet in zip([5, 25, 50], [10, 2, 1]): # , 320, 1280]:

    for iteration in range(repet):

        ti = time.time()

        nmr_intra_H2O = NMRD(
            u=u,
            atom_group=H_H2O,
            type_analysis="intra_molecular",
            number_i=n)
        results_nmr_intra_H2O = nmr_intra_H2O.run_analysis()
        save_result(results_nmr_intra_H2O, n, iteration, f"nmr_intra")

        nmr_inter_H2O = NMRD(
            u=u,
            atom_group=H_H2O,
            type_analysis="inter_molecular",
            number_i=n)
        results_nmr_inter_H2O = nmr_inter_H2O.run_analysis()
        save_result(results_nmr_inter_H2O, n, iteration, f"nmr_inter")

        nmr_intra_H2O = NMRD(
            u=u,
            atom_group=H_H2O,
            type_analysis="full",
            isotropic=False,
            number_i=n)
        results_nmr_intra_H2O = nmr_intra_H2O.run_analysis()
        save_result(results_nmr_intra_H2O, n, iteration, f"nmr_full")

        tf = time.time()
        print("time", np.round(tf-ti,1))
