# PARSE
Prospective Analysis of Reaction ScopE
Python tools for mapping chemical space.
Technical Supplemental Information for DOI: 

This repo contains the Python code and data for reproducing the results within the above article. Each python script contains details for use and sample output files.

Sulfoxide SMILES can be generated from SMILES of carboxylic acids and sulfinic salts by designating linkages within the substructure SMILEs. For instance, COOH groups of carboxylic acids or 'S(=O)O' within sulfinates can be converted to "joints"  and then combined to designate sulfoxides:

![image](https://user-images.githubusercontent.com/49004818/172218557-d3a3aa43-17b3-420b-8ce7-5397cb693999.png)

Once prepared, the two sets of substructures are systematically combined to generate a set of possible sulfoxides: 

![image](https://user-images.githubusercontent.com/49004818/172219395-a4199b4c-9e41-4287-af4c-60c592f2fa62.png)
![image](https://user-images.githubusercontent.com/49004818/172219417-bb1fd827-e8f8-4fc5-a644-156cc43da0b4.png)

Chemical descriptors for molecular complexity are then calculated for each sulfoxide: molecular weight (MW), fraction of SP3-hybridized carbon atoms (FSP3), and Bottcher Complexity (Cm, calculated with a modified version of the Forli group's BottchScore: https://github.com/forlilab/bottchscore): The set of generated sulfoxides can then be plotted relative to the set of known sulfoxides in various manners:

First, as a 3D scatter plot colored according to density of structures:

![image](https://user-images.githubusercontent.com/49004818/172219857-4d566497-09e4-4396-8cf1-0ac01f121a7c.png)

Alternatively, the FSP3/Cm data can be presented as 3D surfaces with amplitude determined by the value of the probability density function (PDF): 

![image](https://user-images.githubusercontent.com/49004818/172220046-08d122fa-bdb7-4753-ad22-b0264f101351.png)

Finally, the data can be presented in plots of the first two normalized Principal Moments of Inertia to determine the shape space of molecules within each set: 

![image](https://user-images.githubusercontent.com/49004818/172220219-1a57c2c2-ddbc-4a6e-95f7-8d5087132b79.png)

