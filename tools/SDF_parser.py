import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

#these are for drawing, not really necessary because img is terrible (really small)
import rdkit.Chem.Draw
from rdkit.Chem.Draw import rdMolDraw2D
try:
    import Image
except ImportError:
    from PIL import Image
from io import BytesIO

#can optionally output an image - it's kind of janky, though
def DrawMolsZoomed(mols, molsPerRow=10, subImgSize=(600, 600)):
	nRows = len(mols) // molsPerRow
	if len(mols) % molsPerRow: nRows += 1
	fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
	full_image = Image.new('RGBA', fullSize )
	for ii, mol in enumerate(mols):
		if mol.GetNumConformers() == 0:
			AllChem.Compute2DCoords(mol)
		column = ii % molsPerRow
		row = ii // molsPerRow
		offset = ( column*subImgSize[0], row * subImgSize[1] )
		d2d = rdMolDraw2D.MolDraw2DCairo(subImgSize[0], subImgSize[1])
		d2d.DrawMolecule(mol)
		d2d.FinishDrawing()
		sub = Image.open(BytesIO(d2d.GetDrawingText()))
		full_image.paste(sub,box=offset)
	return full_image

#specify an SDF file
fileName = 'sulfinate_salts.sdf'

#parse the sdf file
molecules = Chem.SDMolSupplier(fileName)

#count number of mols contained in the sdf
print(f'Number of mols in sdf: {len(molecules)}')

#for drawing if necessary
#img = DrawMolsZoomed(molecules)
#img.save('testing.png')

#make some empty lists to store data
numData = []
smiles = []
#inchis = []
molIDs = []
numAtoms = []
molWeights = []
molFPS3s = []
molcomplexities = []

#count the #mols in sdf to index the csv
molcount = len(molecules)
for value in range(1,molcount+1):
	numData.append(value)

#convert mols to smiles and inchis
for molecule in molecules:
		#
		if molecule is None: continue

		#get the name
		mol_id = molecule.GetProp('_Name')
		molIDs.append(mol_id)

		#convert to smiles
		smi = Chem.MolToSmiles(molecule)
		smiles.append(smi)

		#count numAtoms
		atomCount = molecule.GetNumAtoms()
		numAtoms.append(atomCount)

		#inch = rdkit.Chem.inchi.MolToInchi(molecule, options='-SNon', logLevel=None, treatWarningAsError=False)
		#trimInch = inch.replace("InChI=", "")
		#inchis.append(trimInch)

#calculate Bottcher Complexities using modified BottchScore.py from Forli Lab
Bscores = BottchScore.runBS(fileName)

#calculate Bottcher per Atom
BPA = [float(bottch) / int(atoms) for bottch,atoms in zip(Bscores, numAtoms)]

#calculate MW and FSP3 from extracted SMILEs
for smile in smiles:
	#calculate MW from SMILEs
	calcMW = Descriptors.MolWt(Chem.MolFromSmiles(smile))
	molWeights.append(calcMW)

	#calculate FSP3 - might need to do molfromsmiles
	calcFSP3 = rdkit.Chem.Lipinski.FractionCSP3(Chem.MolFromSmiles(smile))
	molFPS3s.append(calcFSP3)

#write data to a CSV
outData = pd.DataFrame(columns=['Num','molID', 'SMILEs', 'NumAtoms', 'MW', 'FSP3', 'Bottcher', 'BottchPA'])
outData['Num'] = numData
outData['molID'] = molIDs
outData['SMILEs'] = smiles
outData['NumAtoms'] = numAtoms
outData['MW'] = molWeights
outData['FSP3'] = molFPS3s
outData['Bottcher'] = Bscores
outData['BottchPA'] = BPA
#outData['InChIs'] = inchis

outData.to_csv('sulfinate1.csv', index=False)

