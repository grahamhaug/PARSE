import os
import pandas as pd

'''
imports a csv file
reads the SMILEs column
finds any SMILE that starts with the flagged atom 'B'
Makes a new csv from the truncated dataset
'''

def csv_chopper():

	#pull in an existing csv
	incoming_csv = pd.read_csv('xinchi_cacids.csv')

	#read the SMILEs column; retain row if SMILEs has 
	#the substituted group a the start
	chopped_df = incoming_csv[incoming_csv['SMILEs'].str.startswith('B')]
	#print(test_df)

	#write to a new csv
	chopped_df.to_csv('chopped-SAs.csv', index=False)