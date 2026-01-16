# Synthetic Cryo-EM Volume Generator

A Python script that generates synthetic cryo-EM volumes from atomic models (PDB files). The script creates a 3D density map for each PDB by depositing a precomputed Gaussian kernel at each atom's position, with optional filtering to include only standard amino acids (to exclude glycans and other non-protein residues). It can also simulate realistic experimental errors such as alignment errors (via random translations and rotations) and additive Gaussian noise. The final output is an averaged MRC volume that can be viewed in tools such as ChimeraX.

## Repository Structure

```
├── python/
│   └── Gen_Synth_CryoEM_Volume.py    # Main volume generation script
├── requirements.txt                   # Python dependencies
└── LICENSE
```

## Features

- **Volume Generation:** Creates a 3D density map from a set of PDB files by placing a Gaussian kernel at every atom.

- **Residue Filtering:** Optionally include only standard amino acid residues, excluding non-protein molecules.

- **Simulation of Experimental Errors:**
  - **Alignment Errors:** Simulates random translations and rotations to mimic misalignment.
  - **Noise Simulation:** Adds Gaussian noise to the volume.

- **Adjustable Parameters:** Customize grid size, voxel size (to maintain physical dimensions), and simulation parameters to suit your experiment.

- **Output:** Produces an averaged MRC volume from all input PDB files.

## Requirements

- Python 3.9+
- NumPy
- mrcfile
- SciPy

### Installation

```bash
pip install -r requirements.txt
```

## Usage

1. Update the paths in `python/Gen_Synth_CryoEM_Volume.py`:
   - `pdb_directory`: Path to folder containing input PDB files
   - `output_filename`: Path for output MRC file

2. Adjust parameters as needed:
   - `voxel_size`: Angstrom per voxel (default: 2.0)
   - `nx, ny, nz`: Grid dimensions (default: 128x128x128)
   - `sigma`: Gaussian width in Angstrom (default: 1.5)
   - `include_only_amino_acids`: Filter out non-protein atoms
   - `simulate_alignment_error`: Add random translations/rotations
   - `simulate_noise`: Add Gaussian noise

3. Run from the repository root:

```bash
python python/Gen_Synth_CryoEM_Volume.py
```

## Output

The script produces an MRC file containing the averaged 3D density volume, which can be visualized in molecular graphics programs such as UCSF ChimeraX.

## License

MIT License - see [LICENSE](LICENSE) file.
