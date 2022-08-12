#### Python scripts for the oxidation of sulfides to sulfoxides (at scale) using Reaction SMARTS.
SMARTS substructure search for sulfides (exclude thiols, disulfides, thiophene derivatives, thioesters):
![image](https://user-images.githubusercontent.com/49004818/184387077-8680aaa2-f22c-4950-9893-2cf3f3271874.png)

Selectively oxidize confirmed sulfides to sulfoxides using Reaction SMARTS:
![image](https://user-images.githubusercontent.com/49004818/184387147-774e622c-6269-4b6c-a056-dfd77476d885.png)

*Parse_PubChem_Sulfides.py*
-------------------------------
- Process isosmiles downloaded from Pubchem.  
- Splits multicomponents, removes bad FGs, removes SMILES with nonzero net charge, removes SMILES with MW > 1000.  
- Performs SMARTS substructure search to confirm each structure contains at least one sulfide (does not match thiols, disulfides, thioesters, nor thiophenes).  
- Batches the above operations to not crash local PCs due to memory when working with large datasets (tested up to ~5 million structures)   
- Returns a .csv file containing canonical SMILES for non-duplicate sulfides  
- Has dependencies for functions in Functions_for_Parsing.py  
![image](https://user-images.githubusercontent.com/49004818/184448982-4e1bfd21-708f-4cf8-bc8e-c3ef43370ea1.png)

*Oxidize_Sulfides_SMARTS.py*
-------------------------------
- Oxidize sulfides SMILES to sulfoxide SMILES using Reaction SMARTS  
- Processes data in batches of ~10,000 SMILES for memory  
- Outputs a series of .csv files containing the resulting sulfoxide SMILES (ex: "oxidized_sulfides.csv")  
- A sample .csv is provided in the current directory ("sulfide_SMILES.csv")  
![image](https://user-images.githubusercontent.com/49004818/184449161-a2b6f28a-4571-4572-90d3-e986cd00e8d8.png)

*Calculate_Complexity_Descriptors.py*
-------------------------------------
- Calculates MW, FSP3, and Cm chemical descriptors  
- Has dependencies for Functions_for_PropCalc.py and BottchScore.py  
- Outputs a .csv file containing structure Number, SMILES, MW, FSP3, and Cm data  
- Example output: "oxidized_sulfides_Properties.csv"  
- Some OpenBabel memory printouts are likely but inoffensive in the scheme of things  
