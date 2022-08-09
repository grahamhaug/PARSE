import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from itertools import chain



def match_n_sort_sulfides(canon_smiles):
    """
    Input: 
        list of canonical SMILES for "sulfides"
        -might contain thiols/disulfides, whatnots
    Returns:
        list of SMILES with at least one confirmed sulfide
    Prints:
        Count of sulfides returned/"sulfides" excluded as "chaff"
    """

    #SMARTS definition for sulfides; skips thiols, skips disulfides
    sulfide_smarts = Chem.MolFromSmarts('[#16X2H0&!$(SS)]')

    #convert canon_smiles to list of mols
    mols = [Chem.MolFromSmiles(smile) for smile in canon_smiles]

    #lists
    sulfides=[] #collect the confirmed sulfides
    chaff=[] #collect anything else (thiols, disulfides, extraneous)
    highlights=[] #collect the sulfide matches for images
    
    #loop over each of the mols in the mols list, check for sulfides
    for mol in mols:
        #Print the SMILE for debug:
        current_smile = Chem.MolToSmiles(mol)
        print(f'\nSMILE: {current_smile}:')

        #check how many sulfides are present in current_smile
        matches = mol.GetSubstructMatches(sulfide_smarts)

        #count the number of matches
        num_subs = len(matches)
        print(f'Number of Sulfides: {num_subs}')

        #highlight matched sulfide S atoms
        allsubs = tuple(chain.from_iterable(matches))
        highlights.append(allsubs)

        #retain the sulfides
        if num_subs >= 1:
            sulfides.append(Chem.MolToSmiles(mol))
        #everything else gets stuck into "chaff"
        elif num_subs == 0:
            chaff.append(Chem.MolToSmiles(mol))

    #output an image of the matched sulfides
    image = Draw.MolsToGridImage(mols, highlightAtomLists=highlights, molsPerRow=5, subImgSize=(600, 300))
    image.save(f'Sulfides_highlighted.png')

    print(f'\nTotal # sulfides: {len(sulfides)}')
    print(f'Total # non-sulfides: {len(chaff)}')

    return sulfides


def oxidize_sulfides(processed_sulfides, rxn_oxidize_sulfide):
    """
    Inputs: 
        processed_sulfides - a list of confirmed sulfide-containing canonical SMILES
        rxn_oxidize_sulfide - reaction SMARTS def for sulfide => sulfoxide oxidation
    Returns:
        list of unique sulfoxide products

    - where multiple sulfides exist, returns only the maximally "sulfoxidized" product    
    """

    sulfide_smarts = Chem.MolFromSmarts('[#16X2H0&!$(SS)]')

    #convert incoming sulfides to list of mols
    mols = [Chem.MolFromSmiles(smile) for smile in processed_sulfides]

    #retain only unique products via set() method
    products = set()

    #store these highlighted atoms for vis.
    highlights = []

    for mol in mols:
        #check how many sulfide matches in current mol
        matches = mol.GetSubstructMatches(sulfide_smarts)
        #count the number of matches
        num_subs = len(matches)
        allsubs = tuple(chain.from_iterable(matches))
        highlights.append(allsubs)

        #apply reaction SMARTS to convert sulfide ==> sulfoxide
        ps = rxn_oxidize_sulfide.RunReactants((mol,))

        #if ps has more than one solution => multiple substitutions possible
        #Iterate until possible substitutions are exhausted
        while len(ps) > 1:
            #IDK if I'm just missing something but I can't pass tuple directly back in to RunReactants()
            #unless I do this back conversion to SMILE and then convert to Mol again. Something strange
            #going on under the hood of th is RunReactants() method. For now, I have to do this iterative
            #approach with multi-substituted sulfides. 
            # convert mol to smile
            tmp_smile = Chem.MolToSmiles(ps[0][0],isomericSmiles=True)
            # convert back to mol
            tmp_mol = Chem.MolFromSmiles(tmp_smile)

            #run rxn again on the iterated product
            ps = rxn_oxidize_sulfide.RunReactants((tmp_mol,))

        #now ps should be len 1 no matter what
        #add the sulfoxide to the list of products
        products.add(Chem.MolToSmiles(ps[0][0],isomericSmiles=True))

    #output image of confirmed sulfides
    image = Draw.MolsToGridImage(mols, highlightAtomLists=highlights, molsPerRow=5, subImgSize=(600, 300))
    image.save(f'Confirmed_Sulfides.png')

    return products


#list of sulfides to convert to sulfoxides
sulfides_raw = ['CSCC',
                '[H]SC(C)Cl',
                'CC(CC(SCC)C1=CC=CC=C1)S[H]',
                'CC(C)SC(C1=CC=CC=C1)C',
                'CCSC(SC(C)C)C',
                'CSSC(C)C'
                ]

#first convert inputs to canonical smiles
canon_smiles = []
for smile in sulfides_raw:
    canon = Chem.CanonSmiles(smile)
    canon_smiles.append(canon)

print(f'\nCanonical SMILES:')
for smile in canon_smiles:
    print(smile)

processed_sulfides = match_n_sort_sulfides(canon_smiles)
print(f"\nConfirmed Sulfides:")
for sulfide in processed_sulfides:
    print(sulfide)

rxn_oxidize_sulfide = AllChem.ReactionFromSmarts('[#16X2H0&!$(SS):1]>>[#16X430&!$(SS):1]=[O]')

products = oxidize_sulfides(processed_sulfides, rxn_oxidize_sulfide)

product_mols = list(map(Chem.MolFromSmiles, products))
product_image = Draw.MolsToGridImage(product_mols, molsPerRow=5, subImgSize=(600, 300))
product_image.save(f'Oxidized_to_Sulfoxides.png')
