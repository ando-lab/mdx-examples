# `mdx-examples`

Workflows for diffuse scattering data reduction and lattice disorder modeling with MATLAB and [mdx-lib](https://github.com/ando-lab/mdx-lib).

The lysozyme examples reproduce the workflow described in [Meisburger, Case, & Ando, 2022].

## Version history

## Repository Contents

Each example consists of a series of MATLAB scripts (starting with job##) that are run in numerical order. A "setup" job downloads source data from the SBGrid databank, and all other input files are included in the repository (such as XDS input files, structure factor files, and atomic coordinate files).

In each example, the workflow is divided into two parts organized within subfolders of this repository. In subfolders ending in `_map`, the diffraction images are processed to produce a three-dimensional map of diffuse scattering. In subfolders ending in `_model`, a lattice disorder model is fit to the diffuse scattering, and the data are further analyzed following the workflow described in [Meisburger, Case, & Ando, 2022].

### Example 1: Triclinic lysozyme

- Part 1: [lys_tri_map](lys_tri_map)
- Part 2: [lys_tri_model](lys_tri_model)

The dataset for lysozyme in the triclinic space group (P1) is described in [Meisburger, Case, & Ando, 2020] and in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 747]), and the structure has been deposited in the Protein Data Bank ([PDB 6o2h]).

The final job generates a report summarizing the results. A copy is included here: [report_lys_tri_model.md](lys_tri_model/report/report_lys_tri_model.md)

### Example 2: Orthorhombic lysozyme

- Part 1: [lys_ortho_map](lys_ortho_map)
- Part 2: [lys_ortho_model](lys_ortho_model)

The dataset for lysozyme in orthorhombic space group (P2~1~2~1~2~1~) is described in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 958]), and the structure has been deposited in the Protein Data Bank ([PDB 8dz7]).

The final job generates a report summarizing the results. A copy is included here: [report_lys_ortho_model.md](lys_ortho_model/report/report_lys_ortho_model.md)

### Example 3: Tetragonal lysozyme

- Part 1: [lys_tet_map](lys_tet_map)
- Part 2: [lys_tet_model](lys_tet_model)

This dataset for lysozyme in tetragonal space group (P4~3~2~1~2) is described in [Meisburger, Case, & Ando, 2022].  Diffraction images are in the SBGrid databank ([SBGrid 957]), and the structure has been deposited in the Protein Data Bank ([PDB 8dyz]).

The final job generates a report summarizing the results. A copy is included here: [report_lys_tet_model.md](lys_ortho_model/report/report_lys_tet_model.md)

## Requirements

- MATLAB (<https://www.mathworks.com/>). Tested using version R2021a on Mac.
- mdx-lib (<https://github.com/ando-lab/mdx-lib>). Tested using version



## Instructions

## Expected Results

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
