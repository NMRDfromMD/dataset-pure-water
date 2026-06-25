import os, git
import numpy as np


def get_git_repo_path():
    current_path = os.getcwd()
    repo = git.Repo(current_path, search_parent_directories=True)
    return repo.git.rev_parse("--show-toplevel")

def save_result(data, n, name):
    """Save or update correlation function results into a .npy file."""
    
    # Create directory if it doesn't exist
    output_dir = f"{name}"
    os.makedirs(output_dir, exist_ok=True)
    saving_file = os.path.join(output_dir, f"result{n}.npy")
    
    # Extract relevant data
    t = data.results["t"]
    f = data.results["f"]
    C = data.results["gij"]
    R1 = data.results["R1"]
    R2 = data.results["R2"]
    N = data.group_j.atoms.n_atoms
    
    # If file exists, load and combine with previous data
    if os.path.exists(saving_file):
        try:
            previous_data = np.load(saving_file, allow_pickle=True).item()
            t_prev = np.real(previous_data["t"])
            
            if len(t_prev) == len(t):
                N_prev = np.real(previous_data["N"])
                C_prev = np.real(previous_data["C"])
                R1_prev = np.real(previous_data["R1"])
                R2_prev = np.real(previous_data["R2"])

                # Weighted averaging
                total_N = N + N_prev
                C = (C * N + C_prev * N_prev) / total_N
                R1 = (R1 * N + R1_prev * N_prev) / total_N
                R2 = (R2 * N + R2_prev * N_prev) / total_N
                N = total_N
        except Exception as e:
            print(f"Warning: Could not load previous data due to: {e}")
    
    # Save updated data
    result = {
        "t": t,
        "f": f,
        "C": C,
        "N": N,
        "R1": R1,
        "R2": R2
    }
    np.save(saving_file, result)


def log_bin(x, y, num_bins=50):
    x = np.array(x)
    y = np.array(y)

    assert np.all(x > 0)

    log_min = np.log10(x.min())
    log_max = np.log10(x.max())

    log_bins = np.logspace(log_min, log_max, num_bins + 1)
    bin_centers = np.sqrt(log_bins[:-1] * log_bins[1:])  # geometric mean

    indices = np.digitize(x, log_bins)

    binned_y = [
        y[indices == i].mean() if np.any(indices == i) else np.nan
        for i in range(1, len(log_bins))
    ]

    bin_centers = bin_centers[~np.isnan(binned_y)]
    binned_y = np.array(binned_y)[~np.isnan(binned_y)]

    return bin_centers, binned_y