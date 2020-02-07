#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 16:14:17 2020

@author: Murat Cihan Sorkun

A basic search for the minimum energy conformer
Note: It is more efficient to align similar conformers first, but this example does not include it
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Creating mol object from SMILES and add H's
mol = Chem.MolFromSmiles("O=C(NC1CCCCC1)NS(=O)(c1ccc(CCNC(c2cnc(cn2)C)=O)cc1)=O")
mol_h_UFF = Chem.AddHs(mol)
mol_h_MMFF = Chem.AddHs(mol)

# Number of conformers to be generated
num_of_conformer=100
max_iter=500
# Default values for min energy conformer
min_energy_UFF=10000
min_energy_index_UFF=0
min_energy_MMFF=10000
min_energy_index_MMFF=0

# Generate conformers (stored in side the mol object)
#cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=num_of_conformer)
cids = AllChem.EmbedMultipleConfs(mol_h_UFF, numConfs=num_of_conformer,params=AllChem.ETKDG())
cids = AllChem.EmbedMultipleConfs(mol_h_MMFF, numConfs=num_of_conformer,params=AllChem.ETKDG())
ids = list(cids) #You can reach conformers by ids


results_UFF = AllChem.UFFOptimizeMoleculeConfs(mol_h_UFF,maxIters=max_iter)
results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,maxIters=max_iter)


# Search for the min energy conformer from results(tuple(is_converged,energy))
print("Searching conformers by UFF ")       
for index, result in enumerate(results_UFF):
    if(min_energy_UFF>result[1]):       
        min_energy_UFF=result[1]
        min_energy_index_UFF=index
        print(min_energy_index_UFF,":",min_energy_UFF)
		
print("\nSearching conformers by MMFF ")   
for index, result in enumerate(results_MMFF):
    if(min_energy_MMFF>result[1]):       
        min_energy_MMFF=result[1]
        min_energy_index_MMFF=index
        print(min_energy_index_MMFF,":",min_energy_MMFF)



# Write minimum energy conformers into a SDF file
w = Chem.SDWriter('minimum-energy-conformer-UFF.sdf')
w.write(Chem.Mol(mol_h_UFF,False,min_energy_index_UFF))
w.flush()  
w.close()


# Write minimum energy conformer into a SDF file
w = Chem.SDWriter('minimum-energy-conformer-MMFF.sdf')
w.write(Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF))
w.flush()  
w.close()

