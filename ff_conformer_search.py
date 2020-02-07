#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:13:55 2020

@author: Murat Cihan Sorkun

A basic search for the minimum energy conformer
Note: It is more efficient to align similar conformers first, but this example does not include it
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


# Here we are going to create copies of molecule into a list
molecules_UFF=[]
molecules_MMFF=[]
for x in range(num_of_conformer):
    molecules_UFF.append(mol_h)
    molecules_MMFF.append(mol_h)


  
# Search for the min energy conformer by UFF (universal force field) optimization
print("Searching conformers by UFF ") 
for index, mol in enumerate(molecules_UFF):
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    ff = AllChem.UFFGetMoleculeForceField(mol)
    ff.Initialize()
    ff.Minimize()
    if(min_energy_UFF>ff.CalcEnergy()):       
        min_energy_UFF=ff.CalcEnergy()
        min_energy_index_UFF=index
        print(min_energy_index_UFF,":",min_energy_UFF)
        
 
    
# Search for the min energy conformer by MMFF94 force field optimization
print("\nSearching conformers by MMFF ")   
for index, mol in enumerate(molecules_MMFF):
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    prop = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = AllChem.MMFFGetMoleculeForceField(mol,prop)
    ff.Initialize()
    ff.Minimize()
    if(min_energy_MMFF>abs(ff.CalcEnergy())):       
        min_energy_MMFF=abs(ff.CalcEnergy())
        min_energy_index_MMFF=index
        print(min_energy_index_MMFF,":",min_energy_MMFF)

# Compare UFF and MMFF min energies 
if(min_energy_MMFF>min_energy_index_UFF):
    mol_to_save = molecules_UFF[min_energy_index_UFF]
else:
    mol_to_save = molecules_MMFF[min_energy_index_MMFF]

# Write minimum energy conformer into a SDF file
w = Chem.SDWriter('minimum-energy-conformer.sdf')
w.write(mol_to_save)
w.flush()  
w.close()


    
