import os

sdf_files = [file for file in os.listdir('.') if os.path.isfile(file) and file.endswith('.sdf')]

if sdf_files:
	print(f"\nFound {len(sdf_files)} .sdf files in local directory:")
	for sdf_file in sdf_files:
		print(sdf_file)
else:
	print("No .sdf files found in local directory.")

#ask user for an output name - will create example.png in local dir. 
outputName = input("\nEnter desired filename for output sdf (no extension): ")
#sdfOutput = (f'{outputName}.sdf')

with open('%s.sdf' % outputName, 'w') as outfile:
	for sdf_file in sdf_files:
		with open(sdf_file) as infile:
			outfile.write(infile.read())



