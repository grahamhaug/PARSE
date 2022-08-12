import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from itertools import chain

"""
- Validate SMARTS definitition for Carboxylic Acids bound to non-aromatic carbons
- Parses a list of CA SMILES for CA substructures
- Outputs an image with matched substructures highlighted
"""

def match_n_sort_CAs(canon_smiles):
    """
    Input: 
        list of canonical SMILES for "carboxylic acids"
        -might contain all sorts of strange FGs from PubChem
    Returns:
        list of SMILES with at least one confirmed SMARTS match
    Prints:
        Count of mono-/multi-substituted CAs returned/and # excluded structures
    """

    #SMARTS definition for carboxylic acids; skips CAs bound to aromatic carbons
    carboxacid_smarts = Chem.MolFromSmarts('[CX3&!$(Cc)](=O)[OX2H1]')

    #convert canon_smiles to list of mols
    mols = [Chem.MolFromSmiles(smile) for smile in canon_smiles]

    #lists
    mono_subs=[] #collect the mono-substituted CAs
    multi_subs=[] #collect the multi-substituted CAs
    chaff=[] #collect anything else
    highlights=[] #collect the sulfide matches for images
    
    #loop over each of the mols in the mols list, check for sulfides
    for mol in mols:
        #Print the SMILE for debug:
        current_smile = Chem.MolToSmiles(mol)
        print(f'\nSMILE: {current_smile}:')

        #check how many sulfides are present in current_smile
        matches = mol.GetSubstructMatches(carboxacid_smarts)

        #count the number of matches
        num_subs = len(matches)
        print(f'Number of Carboxylic Acids: {num_subs}')

        #highlight matched atoms
        allsubs = tuple(chain.from_iterable(matches))
        highlights.append(allsubs)

        #retain the mono-substituted CAs
        if num_subs == 1:
            mono_subs.append(Chem.MolToSmiles(mol))
        elif num_subs >= 1:
            multi_subs.append(Chem.MolToSmiles(mol))
        #everything else gets stuck into "chaff"
        elif num_subs == 0:
            chaff.append(Chem.MolToSmiles(mol))

    #output an image of the matched acids
    image = Draw.MolsToGridImage(mols, highlightAtomLists=highlights, molsPerRow=5, subImgSize=(600, 300))
    image.save(f'Carboxylic_Acids_Highlighted.png')

    print(f'\nTotal # mono Carboxylic Acids: {len(mono_subs)}')
    print(f'\nTotal # multi Carboxylic Acids: {len(multi_subs)}')
    print(f'Total # non-Carboxylic Acids: {len(chaff)}')

    return mono_subs


#list of carboxylic acid SMILES
acids_raw = ['OC(C1=CC=CC(CC(C)C(O)=O)=C1)=O',
                'O=C(OC)CC(O)=O',
                'OC(CC(N(C)C)=O)=O',
                'OC(CC(O)=O)=O'
                ]

#first convert inputs to canonical smiles
canon_smiles = []
for smile in acids_raw:
    canon = Chem.CanonSmiles(smile)
    canon_smiles.append(canon)

print(f'\nCanonical SMILES:')
for smile in canon_smiles:
    print(smile)

processed_CAs = match_n_sort_CAs(canon_smiles)
print(f"\nConfirmed CAs:")
for acid in processed_CAs:
    print(acid)
