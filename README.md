# `mdx-examples`

Workflows for diffuse scattering data reduction and lattice disorder modeling with MATLAB and [mdx-lib](https://github.com/ando-lab/mdx-lib).

The lysozyme examples reproduce the workflow described in [Meisburger, Case, & Ando, 2022].

## Version history

### Version 0.1.1

Version to accompany the publication of our manuscript in Nature Communications.

### Version 0.1.0

Current development version, accompanies our [preprint on bioRxiv](<https://doi.org/10.1101/2022.08.22.504832>).

## Repository Contents

Each example consists of a series of MATLAB scripts (starting with job##) that are run in numerical order. A "setup" job downloads source data from the SBGrid databank, and all other input files are included in the repository (such as XDS input files, structure factor files, and atomic coordinate files).

In each example, the workflow is divided into two parts organized within subfolders of this repository. In subfolders ending in `_map`, the diffraction images are processed to produce a three-dimensional map of diffuse scattering. In subfolders ending in `_model`, a lattice disorder model is fit to the diffuse scattering, and the data are further analyzed following the workflow described in [Meisburger, Case, & Ando, 2022].

### Example 1: Triclinic lysozyme

- Part 1: [lys_tri_map](lys_tri_map)
- Part 2: [lys_tri_model](lys_tri_model)

The dataset for lysozyme in the triclinic space group (P1) is described in [Meisburger, Case, & Ando, 2020] and in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 747]), and the structure has been deposited in the Protein Data Bank ([PDB 6o2h]).

### Example 2: Orthorhombic lysozyme

- Part 1: [lys_ortho_map](lys_ortho_map)
- Part 2: [lys_ortho_model](lys_ortho_model)

The dataset for lysozyme in orthorhombic space group (P2~1~2~1~2~1~) is described in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 958]), and the structure has been deposited in the Protein Data Bank ([PDB 8dz7]).

### Example 3: Tetragonal lysozyme

- Part 1: [lys_tet_map](lys_tet_map)
- Part 2: [lys_tet_model](lys_tet_model)

This dataset for lysozyme in tetragonal space group (P4~3~2~1~2) is described in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 957]), and the structure has been deposited in the Protein Data Bank ([PDB 8dyz]).

## Requirements

### Software

- Unix-like operating system (tested on Mac)
- [MATLAB](https://www.mathworks.com) and the following toolboxes installed:
  - Parallel Processing Toolbox is used, but can be disabled by setting `useParallel` to false in individual scripts
  - Statistics and Machine Learning Toolbox
- [mdx-lib](https://github.com/ando-lab/mdx-lib), **version 1.2**

The examples have been tested using MATLAB version R2021a on Mac OS 10.14

### Hardware

The `_map` examples can make use of multi-core processing. The required RAM can be very large, especially when building the fine maps. The examples were tested on a desktop computer with 64 Gb of RAM and a 10-core CPU. On the test system, each map example required several hours to run.

The `_model` examples are less resource intensive, however sufficient RAM is needed for storing the diffuse maps. On the test system, each example required about 15 minutes.

## Instructions

After downloading or cloning this repository, open MATLAB and navigate to the example directory (such as `lys_tri_map`). Then,

1. Add mdx-lib directory to the MATLAB path

2. Run the jobs in numerical order.

Currently, the `_model` jobs require that you have run the corresponding `_map` job first, since the experimental maps are needed for modeling. If you want to skip the `_map` step, the processed maps can be downloaded from a public repository (_not yet: coming soon_).

**Notes:**

- In general, the job m-files should be inspected before they are run, in case there are variables that need to be modified (such as `useParallel`).
- Additional information can be found in the m-file comments
- Jobs can be run all at once, or step-wise using cell mode in the editor.
- The workspace can be cleared between jobs. Information that needs to be carried over is saved in files, typically in the `proc` or `export` subdirectories.

## Expected Results

The final modeling job generates a report summarizing the results. Reports are included here:

- [report_lys_tri_model.md](lys_tri_model/report/report_lys_tri_model.md)
- [report_lys_ortho_model.md](lys_ortho_model/report/report_lys_ortho_model.md)
- [report_lys_tet_model.md](lys_tet_model/report/report_lys_tet_model.md)

## References

[Meisburger, Case, & Ando, 2020]: https://doi.org/10.1038/s41467-020-14933-6
Meisburger, SP, Case, DA & Ando, N. (2020). Diffuse X-ray scattering from correlated motions in a protein crystal. _Nat Commun_ **11**, 1271. <https://doi.org/10.1038/s41467-020-14933-6>

[SBGrid 747]: https://doi.org/10.15785/SBGRID/747

[PDB 6o2h]: http://doi.org/10.2210/pdb6O2H/pdb

[SBGrid 958]: https://doi.org/10.15785/SBGRID/958

[PDB 8dz7]: http://doi.org/10.2210/pdb8DZ7/pdb

[SBGrid 957]: https://doi.org/10.15785/SBGRID/957

[PDB 8dyz]: http://doi.org/10.2210/pdb8DYZ/pdb

[Meisburger, Case, & Ando, 2022]: https://doi.org/10.1101/2022.08.22.504832
Meisburger, SP, Case, DA & Ando, N. (2022). Robust total X-ray scattering workflow to study correlated motion of proteins in crystals. _bioRxiv_ 2022.08.22.504832; doi: <https://doi.org/10.1101/2022.08.22.504832>
