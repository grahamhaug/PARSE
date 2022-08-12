#### Python scripts for parsing Carboxylic Acids (at scale).
SMARTS substructure search for Carboxylic Acids (excluding other carbonyl compounds and CAs bound to aromatic carbons):  
![image](https://user-images.githubusercontent.com/49004818/184446901-f5c897b5-dc60-4725-a385-9a6e169b97dc.png)  

*Parse_PubChem_Sulfides.py*
-------------------------------
- Process isosmiles downloaded from Pubchem.  
- Splits multicomponents, removes bad FGs, removes SMILES with nonzero net charge, removes SMILES with MW > 1000.  
- Performs SMARTS substructure search to confirm each structure contains exactly one Carboxylic Acid (structures with multiple CAs are stored separately)    
- Batches the above operations to not crash local PCs due to memory when working with large datasets (tested up to ~5 million structures)   
- Returns a .csv file containing canonical SMILES for non-duplicate structures with a single valid Carboxylic Acid substructure   
- Has dependencies for functions in Functions_for_Parsing.py 
