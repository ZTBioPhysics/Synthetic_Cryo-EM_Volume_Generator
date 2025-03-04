Figure legend | Simulated cryo-EM maps. Synthetic cryo-EM density maps generated from the 1000 atomistic models (Figure 3) of fully glycosylated HK/68 (blue) and SING/16 (green) HAs and from protein residues only (grey) shown at high and low contour levels from the top and side.

Methods

The synthetic cryo-EM maps were generated with a custom python script available at: 

This script generates synthetic cryo-EM volumes from a set of PDB files. Each PDB file is read to generate a 3D density volume by depositing a precomputed Gaussian kernel at the location of every atom. The script can optionally filter out non-protein atoms (e.g., include only standard amino acids), and simulate experimental imperfections such as alignment errors (random translations and rotations) and additive Gaussian noise. Finally, the script sums the individual volumes and computes an average volume, which is then saved in MRC format for further visualization (e.g., in ChimeraX).

# Synthetic Cryo-EM Volume Generator

This repository contains a Python script that generates synthetic cryo-EM volumes from atomic models (PDB files). The script creates a 3D density map for each PDB by depositing a precomputed Gaussian kernel at each atom's position, with optional filtering to include only standard amino acids (to exclude glycans and other non-protein residues). It can also simulate realistic experimental errors such as alignment errors (via random translations and rotations) and additive Gaussian noise. The final output is an averaged MRC volume that can be viewed in tools such as ChimeraX.

## Features

- **Volume Generation:**  
  Creates a 3D density map from a set of PDB files by placing a Gaussian kernel at every atom.

- **Residue Filtering:**  
  Optionally include only standard amino acid residues, excluding non-protein molecules.

- **Simulation of Experimental Errors:**  
  - **Alignment Errors:** Simulates random translations and rotations to mimic misalignment.
  - **Noise Simulation:** Adds Gaussian noise to the volume.

- **Adjustable Parameters:**  
  Customize grid size, voxel size (to maintain physical dimensions), and simulation parameters to suit your experiment.

- **Output:**  
  Produces an averaged MRC volume from all input PDB files.

## Requirements

- Python 3.x
- NumPy
- mrcfile
- SciPy

Install the required packages using pip:

```bash
pip install numpy mrcfile scipy
![image](https://github.com/user-attachments/assets/4ae59079-7c88-4342-8f8e-f087e9acb2a4)
