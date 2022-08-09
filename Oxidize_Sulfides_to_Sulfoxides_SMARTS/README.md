Reviewers suggested that using Reaction SMARTS for chemical transformation could be a more efficient option vs. operating using SMILES only. SMARTS definitely results in "cleaner" code but comes at the cost of having to ensure a valid reaction SMARTS definintion which can be tricky to implement. Additionally, as far as I can tell, there's no satisfying way of iterating over a reaction multiple times in cases where multiple substitutions are possible in RDkit. It's doable by converting each iteration's product to a SMILE and then back into a mol, but it would be better (faster/cleaner) to pass back the tuple directly. RDkit seems to do something under the hood to the product mol object that prevents it from being passed back in directly. It's strange. Maybe someone knows a better solution. 

Reviewers also mentioned that it would be more meaningful to compare sulfoxide generation from Carboxylic Acids + Sulfinates vs. existing methods of sulfoxide synthesis.
Here is an implementation of using reaction SMARTS to oxidize sulfides to sulfoxides (oxidation of sulfides). 

Candidate Sulfides are parsed for sulfide substructures using reaction SMARTS, excluding thiols and disulfides:
![image](https://user-images.githubusercontent.com/49004818/183739427-1ca2b405-fba6-4b90-9f9b-4ca54d4dab9f.png)

Structures containing confirmed sulfides are oxidized to produce sulfoxides. In instances of multiple sulfides, all sulfides are oxidized.
![image](https://user-images.githubusercontent.com/49004818/183739652-40d0051e-482a-40e0-b9b0-a86b55901084.png)
