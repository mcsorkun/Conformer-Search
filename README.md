# Conformer Search
Conformers are different configurations of a molecule (same formula and connections) by rotating the single bonds. 

This repository contains a simple workflow for minimum energy conformer search based on force field optimization methods (UFF and MMFF94) available from RDKit.

Workflow consists of the following steps:

- Import molecule from SMILES
- Generate conformers for given molecule
- Optimiza conformers by force field optimization (UFF and MMFF94)
- Select the minimum energy conformer and write it into a SDF file

![alt text](https://raw.githubusercontent.com/mcsorkun/Conformer-Search/master/conformer-search.png)

Requirements
- RDKit

Disclaimer/Warning from RDKit: Conformation generation is a difficult and subtle task. The original 2D->3D conversion provided with the RDKit was not intended to be a replacement for a “real” conformational analysis tool; it merely provides quick 3D structures for cases when they are required. We believe, however, that the newer ETKDG method[#riniker2]_ should be adequate for most purposes.

Note: It is recommended to use MMFF94 for organic-like molecules and UFF for metal contained molecules.


References:
[1] Riniker, S.; Landrum, G. A. (2015). Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation J. Chem. Inf. Comp. Sci. 55:2562-74.
[2] Ebejer, J. P., Morris, G. M., & Deane, C. M. (2012). Freely available conformer generation methods: how good are they?. Journal of chemical information and modeling, 52(5), 1146-1158.
