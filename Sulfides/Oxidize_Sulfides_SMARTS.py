import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from itertools import chain
from Functions_for_Parsing import process_in_batches
import time

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
    sulfide_smarts = Chem.MolFromSmarts('[#16X2H0&!$(s1cccc1)&!$(SS)]')

    #convert canon_smiles to list of mols
    mols = [Chem.MolFromSmiles(smile) for smile in canon_smiles]

    #lists
    sulfides=[] #collect the confirmed sulfides
    chaff=[] #collect anything else (thiols, disulfides, extraneous)
    # highlights=[] #collect the sulfide matches for images
    
    #loop over each of the mols in the mols list, check for sulfides
    for mol in mols:
        #Print the SMILE for debug:
        # current_smile = Chem.MolToSmiles(mol)
        # print(f'\nSMILE: {current_smile}:')

        #check how many sulfides are present in current_smile
        matches = mol.GetSubstructMatches(sulfide_smarts)

        #count the number of matches
        num_subs = len(matches)
        print(f'Number of Sulfides: {num_subs}')

        # #highlight matched sulfide S atoms
        # allsubs = tuple(chain.from_iterable(matches))
        # highlights.append(allsubs)

        #retain the sulfides
        if num_subs >= 1:
            sulfides.append(Chem.MolToSmiles(mol))
        #everything else gets stuck into "chaff"
        elif num_subs == 0:
            chaff.append(Chem.MolToSmiles(mol))

    # #output an image of the matched sulfides
    # image = Draw.MolsToGridImage(mols, highlightAtomLists=highlights, molsPerRow=5, subImgSize=(600, 300))
    # image.save(f'Sulfides_highlighted.png')

    print(f'\nTotal # sulfides: {len(sulfides)}')
    print(f'Total # non-sulfides: {len(chaff)}')

    return sulfides


def oxidize_sulfides(smiles_list, rxn_oxidize_sulfide,csv_name):
    """
    Inputs: 
        processed_sulfides - a list of confirmed sulfide-containing canonical SMILES
        rxn_oxidize_sulfide - reaction SMARTS def for sulfide => sulfoxide oxidation
        csv_name - name for the output file(s)
    Returns:
        list of unique sulfoxide products
        a series of csv files containing the products of reaction

    - where multiple sulfides exist, returns only the maximally "sulfoxidized" product    
    """

    sulfide_smarts = Chem.MolFromSmarts('[#16X2H0&!$(s1cccc1)&!$(SS)]')

    #for naming
    # output_chunk = 10000
    # chunk_count = 0
    # append_count = 1

    #retain only unique products via set() method
    products = []

    #slice smiles_list into smaller chunks of 10k SMILES for memory
    chunked_list = process_in_batches(smiles_list)   
    num_chunks = len(chunked_list)
    print(f'\nOxidizing sulfides to sulfoxides in {num_chunks} batches...')
    batch_num = 1

    #now working through an individual chunk
    for chunk in chunked_list:
        print(f'\tBatch {batch_num}')

        #convert SMILES in chunk to list of mols
        mols = [Chem.MolFromSmiles(smile) for smile in chunk]

        #now loop through the chunked mols
        for mol in mols:

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
            products.append(Chem.MolToSmiles(ps[0][0],isomericSmiles=True))


        #write the products to a csv
        csv_label = f"{csv_name}_{batch_num}"
        outData = pd.DataFrame()
        outData['SMILEs'] = products
        outData.to_csv(f'{csv_label}.csv', index=False, header=True)
        print(f'\t\t{len(products)} Sulfoxides written to {csv_label}.csv')

        #reset products
        products = []

        batch_num += 1      

    #output remaining sulfoxides to a csv
    outDatac = pd.DataFrame(columns=['SMILEs'])
    outDatac['SMILEs'] = products
    outDatac.to_csv(f'{csv_name}_last.csv', index=False, header=True)

#timeit
start=time.time()

### Data Import ###
import_file_name = 'sulfide_SMILES.csv'
chop_name = import_file_name.split('.',1)[0]
#using fields to not pull in the full mega frame
fields = ['SMILES']
incoming_data = pd.read_csv(f'{import_file_name}', usecols=fields)
#convert 'isosmiles' column to list, drop null values
raw_smiles = incoming_data[incoming_data['SMILES'].notnull()]
incoming_smiles = raw_smiles['SMILES'].to_list()

#count number of incoming sulfides:
incoming_length = len(incoming_smiles)
print(f'\nOxidizing sulfides contained in: {import_file_name}')
print(f'Incoming sulfide SMILES: {incoming_length}')

#Reaction SMARTS for oxidation of sulfide ==> sulfoxide
rxn_oxidize_sulfide = AllChem.ReactionFromSmarts('[#16X2H0&!$(s1cccc1)&!$(SS):1]>>[#16X430&!$(SS):1]=[O]')

#Function appends to a .csv 
oxidize_sulfides(incoming_smiles,rxn_oxidize_sulfide,'oxidized_sulfoxides')

#timeit
end = time.time()
elapsed = end-start
format_elapsed = "{:.2f}".format(elapsed)
print(f'Sec Elapsed: {format_elapsed}')


### Use below for small sets of sulfides ### 
#list of sulfides to convert to sulfoxides
# sulfides_raw = ['CSCC',
#                 '[H]SC(C)Cl',
#                 'CC(CC(SCC)C1=CC=CC=C1)S[H]',
#                 'CC(C)SC(C1=CC=CC=C1)C',
#                 'CCSC(SC(C)C)C',
#                 'CSSC(C)C',
#                 'CC1(c2ccccc2)OC(CC(SCCO)c2ccsc2)=CC1=O'
#                 ]             

# #first convert inputs to canonical smiles
# canon_smiles = []
# for smile in sulfides_raw:
#     canon = Chem.CanonSmiles(smile)
#     canon_smiles.append(canon)

# print(f'\nCanonical SMILES:')
# for smile in canon_smiles:
#     print(smile)

# processed_sulfides = match_n_sort_sulfides(canon_smiles)
# print(f"\nConfirmed Sulfides:")
# for sulfide in processed_sulfides:
#     print(sulfide)

# #rxn_oxidize_sulfide = AllChem.ReactionFromSmarts('[#16X2H0&!$(SS):1]>>[#16X430&!$(SS):1]=[O]')
# rxn_oxidize_sulfide = AllChem.ReactionFromSmarts('[#16X2H0&!$(s1cccc1)&!$(SS):1]>>[#16X430&!$(SS):1]=[O]')

# products = oxidize_sulfides(processed_sulfides, rxn_oxidize_sulfide)

# product_mols = list(map(Chem.MolFromSmiles, products))
# product_image = Draw.MolsToGridImage(product_mols, molsPerRow=5, subImgSize=(600, 300))
# product_image.save(f'Oxidized_to_Sulfoxides.png')

# #output a csv
# oxidized_sulfoxides = pd.DataFrame(columns=['SMILES'])
# oxidized_sulfoxides['SMILES'] = products
# oxidized_sulfoxides.to_csv(f'sulfides_to_sulfoxides.csv', index=False)
