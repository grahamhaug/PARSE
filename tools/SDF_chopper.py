import os
'''
Reads in a raw SDF
Chops # of molecules off the end of input SDF
Outputs new SDF with user-selected # of mols
'''

#specify the input sdf
input_sdf = 'Sulfinate_Salts/Substance_20220322_1912_Li500.sdf'

#this delimiter at end of each mol definition
target_string = '$$$$'

#determine how many molecules are in the SDF
with open(input_sdf, 'rt') as input_text:
	contents = input_text.read()
	count = contents.count(target_string)

#report how many molecules were found
print(f'\nCounted {count} molecules in {input_sdf}')

#ask user for an output name - will create example.png in local dir. 
output_sdf = input("\nEnter desired filename for output sdf: ")

#user selects how many molecules to retain
retain_prompt = "\nHow many molecules (int) to include in final sdf? "
num_to_keep = input(retain_prompt)
num_to_keep = int(num_to_keep)
print(f"{num_to_keep} molecules will be captured...")

#copy up to the specified occurrence of $$$$
with open(input_sdf, 'r') as file_to_chop, open(output_sdf, 'w') as file_to_copy:
	target_counter = 0
	for line in file_to_chop:
		# copy line to new sdf file
		file_to_copy.write(line)
		# search for target string in line 
		found = line.find(target_string)
		# if 'found' is not -1, the string was found
		if found != -1:
			target_counter = target_counter + 1
			# check if this is the nth occurrence:
		if target_counter == num_to_keep:
			# stop the loop
			break






