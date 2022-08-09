Reviewers commented that Reaction SMARTS could be a more efficient option for operations with SMILES vs. operating using SMILES only. 

Reviewers also mentioned that it would be more meaningful to compare sulfoxide generation from Carboxylic Acids + Sulfinates vs. existing methods of sulfoxide synthesis.
Here is an implementation of using reaction SMARTS to oxidize sulfides to sulfoxides (oxidation of sulfides). 

Candidate Sulfides are parsed for sulfide substructures using reaction SMARTS, excluding thiols and disulfides:
![image](https://user-images.githubusercontent.com/49004818/183739427-1ca2b405-fba6-4b90-9f9b-4ca54d4dab9f.png)

Structures containing confirmed sulfides are oxidized to produce sulfoxides. In instances of multiple sulfides, all sulfides are oxidized.
![image](https://user-images.githubusercontent.com/49004818/183739652-40d0051e-482a-40e0-b9b0-a86b55901084.png)
