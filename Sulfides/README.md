#### Python scripts for the oxidation of sulfides to sulfoxides.
SMARTS substructure search for sulfides (exclude thiols, disulfides, thiophene derivatives):
![image](https://user-images.githubusercontent.com/49004818/184214576-20cecf8b-39e9-43a9-a098-bdda980b4d3c.png)

Oxidize confirmed sulfides to sulfoxides using Reaction SMARTS:
![image](https://user-images.githubusercontent.com/49004818/184214668-0928a15a-4067-48ba-9776-04873081e246.png)

*Parse_PubChem_Sulfides.py*
-------------------------------
-Process isosmiles downloaded from Pubchem.  
-Splits multicomponents, removes bad FGs, removes SMILES with nonzero net charge, removes SMILES with MW > 1000.  
-Performs SMARTS substructure search to confirm each structure contains at least one sulfide (does not match thiols, disulfides, nor thiophenes).  
-Batches the above operations to not crash local PCs due to memory  
-Returns a .csv file containing canonical SMILES for non-duplicate sulfides  
![image](https://user-images.githubusercontent.com/49004818/184212482-d9021d20-f924-40df-a2d0-d993e01ac6b9.png)  

*Oxidize_Sulfides_SMARTS.py*
-------------------------------
-Oxidize sulfides SMILES to sulfoxide SMILES using Reaction SMARTS  
-Processes data in batches of ~10,000 SMILES for memory  
-Outputs a series of .csv files containing the resulting sulfoxide SMILES
![image](https://user-images.githubusercontent.com/49004818/184225353-41fbcf7f-f6a1-49e5-abb7-ebeb1307d80b.png)

