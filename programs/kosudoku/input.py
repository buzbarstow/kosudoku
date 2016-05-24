# ------------------------------------------------------------------------------------------------ #
# kosudoku.py
# Created by Buz Barstow 2015-6-14
# Last updated by Buz Barstow 2016-02-03
# Utility code for analysis of Sudoku sequencing data
# ------------------------------------------------------------------------------------------------ #

# Input parsing routines

# ------------------------------------------------------------------------------------------------ #
# Function to search input data for a given key
def inputParameterSearch(inputData, inputItem):
	import re
	import pdb
	import numpy
	
	#from string import rstrip
	regex = re.compile('\s*' + inputItem + '\s*' + '='+ '\s*', re.IGNORECASE)
	
	inputFound = False
	continuation = False
	
	for line in inputData:
		if line[0] != '#':
			if continuation == False:
				if regex.match(line) != None:
					if inputFound == False:
						inputLine = line.split('=')
					
						inputFound = True
					
						if line.strip()[-1] == '\\':
	# 						pdb.set_trace()
							continuation = True
							inputParameter = inputLine[1].strip()[0:-1]
						else:
							continuation = False
							inputParameter = inputLine[1].strip()
				
					elif inputFound == True:
						print("An input parameter can only be defined once")
						inputParameter = None
			
			elif continuation == True:
				try:
					if line.strip()[-1] == '\\':
						continuation = True
						inputParameter += line.strip()[0:-1]
					else:
						continuation = False
						inputParameter += line.strip()
				except:
					pdb.set_trace()
			
	if inputFound == False:
		inputParameter = None
	

	return inputParameter
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Input parser
def get_input(argv, inputParameters):
	import re
	import pdb
	
	if len(argv) != 2:
		print("IndexSummarize takes 1 argument, the input file name")
		return
	
	file = open(argv[1],'r')
	inputData = file.readlines()
	
	# Go through all lines and search for input parameters	
	
	inputParameterValues = {}
	
	for item in inputParameters:
		result = inputParameterSearch(inputData, item)
		if result != None:
			inputParameterValues[item] = result
	
	return inputParameterValues
# ------------------------------------------------------------------------------------------------ #