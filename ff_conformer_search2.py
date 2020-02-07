#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 16:14:17 2020

@author: Murat Cihan Sorkun

An alternative method (same FFs) that generates and calculates all the conformers together
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Creating mol object from SMILES and add H's
mol = Chem.MolFromSmiles("O=C(NC1CCCCC1)NS(=O)(c1ccc(CCNC(c2cnc(cn2)C)=O)cc1)=O")
mol_h = Chem.AddHs(mol)

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
cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=num_of_conformer,params=AllChem.ETKDG())
ids = list(cids) #You can reach conformers by ids

# copy molecule for UFF and MMFF optimizations
mol_h_UFF=mol_h
mol_h_MMFF=mol_h

results_UFF = AllChem.UFFOptimizeMoleculeConfs(mol_h_UFF,maxIters=max_iter)
results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,maxIters=max_iter)


# Search for the min energy conformer from results(tuple(is_converged,energy))
print("\nSearching conformers by UFF ")       
for index, result in enumerate(results_UFF):
    if(min_energy_UFF>result[1]):       
        min_energy_UFF=result[1]
        min_energy_index_UFF=index
        print(min_energy_index_UFF,":",min_energy_UFF)
		
print("Searching conformers by MMFF ")   
for index, result in enumerate(results_MMFF):
    if(min_energy_MMFF>abs(result[1])):       
        min_energy_MMFF=abs(result[1])
        min_energy_index_MMFF=index
        print(min_energy_index_MMFF,":",min_energy_MMFF)


# Compare UFF and MMFF min energies 
if(min_energy_MMFF>min_energy_index_UFF):
    mol_to_save = Chem.Mol(mol_h_UFF,False,min_energy_index_UFF)
else:
    mol_to_save = Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF)

# Write minimum energy conformer into a SDF file
w = Chem.SDWriter('minimum-energy-conformer.sdf')
w.write(mol_to_save)
w.flush()  
w.close()