#!/usr/bin/env python
"""
Gen_Sim_CryoEM_vol_multi.py

Description:
    This script generates synthetic cryo-EM volumes from a set of PDB files. Each PDB file is read to 
    generate a 3D density volume by depositing a precomputed Gaussian kernel at the location of every atom.
    The script can optionally filter out non-protein atoms (e.g., include only standard amino acids), and 
    simulate experimental imperfections such as alignment errors (random translations and rotations) and 
    additive Gaussian noise. Finally, the script sums the individual volumes and computes an average volume,
    which is then saved in MRC format for further visualization (e.g., in ChimeraX).

Usage:
    Adjust the parameters below (grid dimensions, voxel_size, sigma, directory paths, and simulation options)
    as needed, then run the script from the command line:
    
        python Gen_Sim_CryoEM_vol_multi.py

Requirements:
    - Python 3.x
    - NumPy
    - mrcfile
    - SciPy

Author:
    Zachary Berndsen
    ChatGPT o3-mini-high
    
Date:
    03-04-2025

License:
    MIT License
"""

import os
import numpy as np
import mrcfile
from scipy.ndimage import shift, rotate

# ------------------------------------------------
# PARAMETERS (adjust these as needed)
# ------------------------------------------------
voxel_size = 2.0             # Angstrom per voxel
nx, ny, nz = 128, 128, 128   # Dimensions of the 3D grid
sigma = 1.5                  # Gaussian width in Angstrom
grid_center = np.array([nx / 2, ny / 2, nz / 2]) # define grid center

# input PDB file directory
pdb_directory = "/path/to/pdb/files/"  # Update this path

# output file name and directory
output_filename = "/path/to/output/volume.mrc"

# Option flags:
include_only_amino_acids = False   # If True, only include standard amino acids
simulate_alignment_error = False   # If True, add random translation and rotation
simulate_noise = True             # If True, add Gaussian noise

# Alignment error parameters:
max_translation = 2         # Maximum translation error (in voxels)
max_rotation_deg = 5        # Maximum rotation error (in degrees)

# Noise parameters:
noise_std = 0.05            # Standard deviation of the noise

# List of standard 3-letter amino acid codes.
amino_acids = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

# ------------------------------------------------
# PRECOMPUTE THE GAUSSIAN KERNEL (for all atoms)
# ------------------------------------------------
def precompute_gaussian_kernel(sigma, voxel_size):
    """
    Precompute a 3D Gaussian kernel for a given sigma (in Angstrom) and voxel size.
    The kernel will extend to ±3σ (in voxel units).
    
    Returns:
      kernel: a 3D numpy array of Gaussian values.
      r: radius of the kernel (in voxels).
    """
    sigma_voxel = sigma / voxel_size
    r = int(np.ceil(3 * sigma_voxel))
    coords = np.arange(-r, r + 1)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
    kernel = np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma_voxel**2))
    return kernel, r

kernel, r = precompute_gaussian_kernel(sigma, voxel_size)
print(f"Precomputed Gaussian kernel with shape: {kernel.shape}")

# ------------------------------------------------
# FUNCTION: Generate a volume from a single PDB file
# ------------------------------------------------
def generate_volume_for_file(pdb_file, include_only_amino_acids=False):
    """
    Generates a density volume from a single PDB file.
      - Reads the file line by line.
      - If include_only_amino_acids is True, only processes lines
        whose residue name (columns 18–20) is in the set of standard amino acids.
      - Computes the geometric center so that the molecule is centered.
      - Deposits a precomputed Gaussian kernel at each atom's (centered) location.
    
    Returns:
      volume: a 3D numpy array with shape (nx, ny, nz).
    """
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            # Process only ATOM and HETATM lines.
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            # If filtering for amino acids, check the residue name (columns 18-20).
            if include_only_amino_acids:
                resname = line[17:20].strip().upper()
                if resname not in amino_acids:
                    continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append([x, y, z])
    coords = np.array(coords)

    if coords.shape[0] == 0:
        # No atoms found; return an empty volume.
        return np.zeros((nx, ny, nz), dtype=np.float32)

    # Compute the geometric center so that the molecule is centered.
    center = coords.mean(axis=0)

    # Create an empty volume.
    volume = np.zeros((nx, ny, nz), dtype=np.float32)

    # Deposit the Gaussian kernel at each atom.
    for coord in coords:
        # Shift coordinate so that the molecule is centered.
        shifted_coord = coord - center
        # Convert from Angstroms to voxel coordinates.
        voxel_coord = shifted_coord / voxel_size + grid_center
        # Round to the nearest voxel index.
        ix, iy, iz = np.round(voxel_coord).astype(int)

        # Determine the sub-region of the volume to update.
        x_min = max(0, ix - r)
        x_max = min(nx, ix + r + 1)
        y_min = max(0, iy - r)
        y_max = min(ny, iy + r + 1)
        z_min = max(0, iz - r)
        z_max = min(nz, iz + r + 1)

        # Determine the corresponding indices in the kernel.
        kx_min = 0 if ix - r >= 0 else r - ix
        ky_min = 0 if iy - r >= 0 else r - iy
        kz_min = 0 if iz - r >= 0 else r - iz
        kx_max = kernel.shape[0] - ((ix + r + 1) - x_max)
        ky_max = kernel.shape[1] - ((iy + r + 1) - y_max)
        kz_max = kernel.shape[2] - ((iz + r + 1) - z_max)

        # Deposit the kernel into the volume.
        volume[x_min:x_max, y_min:y_max, z_min:z_max] += \
            kernel[kx_min:kx_max, ky_min:ky_max, kz_min:kz_max]

    return volume

# ------------------------------------------------
# LOOP OVER ALL PDB FILES, SUM THE INDIVIDUAL VOLUMES
# ------------------------------------------------
pdb_files = sorted([os.path.join(pdb_directory, f)
                    for f in os.listdir(pdb_directory) if f.endswith(".pdb")])
n_files = 0
volume_sum = np.zeros((nx, ny, nz), dtype=np.float32)

print(f"Found {len(pdb_files)} PDB files.")
for pdb in pdb_files:
    print("Processing:", pdb)
    vol = generate_volume_for_file(pdb, include_only_amino_acids=include_only_amino_acids)

    # Simulate alignment error if enabled.
    if simulate_alignment_error:
        # Apply a random translation.
        random_translation = np.random.uniform(-max_translation, max_translation, size=3)
        vol = shift(vol, shift=random_translation, mode='nearest')
        # Apply a random rotation about the z-axis.
        angle = np.random.uniform(-max_rotation_deg, max_rotation_deg)
        vol = rotate(vol, angle=angle, axes=(0, 1), reshape=False, mode='nearest')

    # Simulate noise if enabled.
    if simulate_noise:
        noise = np.random.normal(0, noise_std, vol.shape)
        vol = vol + noise

    volume_sum += vol
    n_files += 1

if n_files == 0:
    raise ValueError("No valid PDB files were processed!")

# Compute the average volume.
average_volume = volume_sum / n_files
print(f"Processed {n_files} PDB files. Calculated average volume.")

# ------------------------------------------------
# SAVE THE AVERAGED VOLUME AS AN MRC FILE
# ------------------------------------------------
with mrcfile.new(output_filename, overwrite=True) as mrc:
    mrc.set_data(average_volume)
    mrc.voxel_size = voxel_size  # Store voxel size metadata

print(f"Averaged volume saved as '{output_filename}'.")
