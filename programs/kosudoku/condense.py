
from kosudoku.grid import SudokuWell, SudokuGenomicCoord


# ------------------------------------------------------------------------------------------------ #
def FindFeaturesWithDisruptionsInProgenitorButNotInQC(featureArray): 
# Take an updated feature array and check it for features that have representative mutants in the 
# progenitor but not in the QC collection

	i = 0

	featuresWithDisruptionsInQCAndProgenitor = []
	featuresWithDisruptionsInProgenitorOnly = []
	

	while i < len(featureArray):
		if len(featureArray[i].sudokuGridEntries) > 0:
			
			featureHasAtLeastOneRepInQCCollection = False
			featureHasAtLeastOneRepInProgenitorCollection = False
			
			sudokuGridEntries = featureArray[i].sudokuGridEntries
			j = 0
			while j < len(sudokuGridEntries):
				
				addressDictKeys = list(sudokuGridEntries[j]['addressDict'].keys())
				
				if 'colonyPurified' in addressDictKeys:
					featureHasAtLeastOneRepInQCCollection = True
				
				if 'progenitor' in addressDictKeys:
					featureHasAtLeastOneRepInProgenitorCollection = True
				
				j += 1
				
			if featureHasAtLeastOneRepInProgenitorCollection == True \
			and featureHasAtLeastOneRepInQCCollection == True:
				featuresWithDisruptionsInQCAndProgenitor.append(featureArray[i])
			elif featureHasAtLeastOneRepInProgenitorCollection == True \
			and featureHasAtLeastOneRepInQCCollection == False:
				featuresWithDisruptionsInProgenitorOnly.append(featureArray[i])

		
		i += 1

	return [featuresWithDisruptionsInQCAndProgenitor, featuresWithDisruptionsInProgenitorOnly]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignProgenitorWellsToCondensedCollection(wellsToBeAssignedInput, startPlateNumber, \
startRowNumber, startColNumber, rows, columns, fillOrder='columns', fillPattern='checkerboard', \
plateRowForCondensedWells=1, fillLastPlateWithBlanks=True):

# Note on input parameters: rows and columns are ordered listings of all of the rows and columns
# on a plate. 

	from copy import deepcopy
	import pdb
	
	wellsToBeAssigned = deepcopy(wellsToBeAssignedInput)
	wellsToBeAssigned.reverse()
	totalWellsForRearray = len(wellsToBeAssigned)
	
	currentRowNumber = startRowNumber
	currentColNumber = startColNumber
	currentPlateNumber = startPlateNumber
	
	wellsAssignedWithMutants = 0
	wellsAssignedOnCurrentPlate = 0
	condensedWellArray = []
	
	while wellsAssignedWithMutants < totalWellsForRearray:
	
		currentRow = rows[currentRowNumber - 1]
		currentCol = columns[currentColNumber - 1]
	
		newWell = SudokuCondensedWell(currentPlateNumber, plateRowForCondensedWells, \
		currentPlateNumber, currentRow, currentCol, addressSystem='condensed')
	
		# Test if this is supposed to be a blank well
		blank = TestIfCondensedCollectionPlateWellIsBlank(currentPlateNumber, currentRowNumber, \
		currentColNumber, fillPattern)
	
		if blank:
		# If the well is supposed to be blank, set readalignment coords to blank
			newWell.readAlignmentCoords.append(deepcopy(blankSudokuCoord))
			newWell.condensationData = blankSudokuGridEntrySummary
			newWell.addressDict['progenitor'] = \
			blankSudokuGridEntrySummary['addressDict']['progenitor']
		else:
		# If the well isn't supposed to be blank, copy over information from wells to be 
		# assigned
			wellToBeAssigned = wellsToBeAssigned.pop()
			addressDictData = deepcopy(wellToBeAssigned[1]['well'].addressDict['progenitor'])
			readAlignmentCoords = deepcopy(wellToBeAssigned[1]['well'].readAlignmentCoords)
			
			newWell.addressDict['progenitor'] = addressDictData
			newWell.readAlignmentCoords = readAlignmentCoords
			newWell.condensationData = wellToBeAssigned[1]['sudokuGridEntry']
		
			wellsAssignedWithMutants += 1
	
		# Add the new well to the growing list of condensed collection wells
		condensedWellArray.append(newWell)
		wellsAssignedOnCurrentPlate += 1
	
		# Increment the row, column and plate numbers
		[nextPlateNumber, nextRowNumber, nextColNumber] = \
		IncrementWell(currentPlateNumber, currentRowNumber, currentColNumber, rows, columns, \
		fillOrder)
	
		if nextPlateNumber != currentPlateNumber:
			wellsAssignedOnCurrentPlate = 0
	
		currentPlateNumber = nextPlateNumber
		currentRowNumber = nextRowNumber
		currentColNumber = nextColNumber
		
	
	# If desired, fill out the remainder of the last plate with blank entries
	if fillLastPlateWithBlanks == True:
					
		lastPlateNumberWithEntries = deepcopy(currentPlateNumber)
		
		while lastPlateNumberWithEntries == currentPlateNumber:
			currentRow = rows[currentRowNumber - 1]
			currentCol = columns[currentColNumber - 1]
	
			newWell = SudokuCondensedWell(currentPlateNumber, plateRowForCondensedWells, \
			currentPlateNumber, currentRow, currentCol, addressSystem='condensed')
			
			newWell.addressDict['progenitor'] = \
			blankSudokuGridEntrySummary['addressDict']['progenitor']
			
			newWell.condensationData = blankSudokuGridEntrySummary
			
			newWell.readAlignmentCoords.append(deepcopy(blankSudokuCoord))
			condensedWellArray.append(newWell)
			
			# Increment the row, column and plate numbers
			[nextPlateNumber, nextRowNumber, nextColNumber] = \
			IncrementWell(currentPlateNumber, currentRowNumber, currentColNumber, rows, columns, \
			fillOrder)
	
			if nextPlateNumber != currentPlateNumber:
				wellsAssignedOnCurrentPlate = 0
	
			currentPlateNumber = nextPlateNumber
			currentRowNumber = nextRowNumber
			currentColNumber = nextColNumber

	return [condensedWellArray, currentPlateNumber, currentRowNumber, currentColNumber]

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SortBestSudokuWellsIntoSinglyAndMultiplyOccupiedSets(featuresAndBestWellArray, \
rearrayPickOrder='columns'):
	
	from copy import deepcopy
	
	featuresAndBestWellArrayCopy = deepcopy(featuresAndBestWellArray)

	singlyOccupiedRepresenativeWells = []
	multiplyOccupiedRepresentativeWells = []
	totalWellsForRearrayOfSinglyOccupiedWells = 0
	totalWellsForRearrayOfMultiplyOccupiedWells = 0

	for featureAndBestWell in featuresAndBestWellArrayCopy:
		locationSummary = featureAndBestWell[1]['sudokuGridEntry']
		multipleOccupancy = locationSummary['multipleOccupancy']
		coloniesToPick = featureAndBestWell[1]['minimumColoniesToPick']
	
		if multipleOccupancy:
			j = 0
			while j < coloniesToPick:
				multiplyOccupiedRepresentativeWells.append(deepcopy(featureAndBestWell))
				multiplyOccupiedRepresentativeWells[-1][1]['sudokuGridEntry']['condensationType'] \
				= 'Colony Purification'
				j += 1
			totalWellsForRearrayOfMultiplyOccupiedWells += coloniesToPick
		else:
			singlyOccupiedRepresenativeWells.append(deepcopy(featureAndBestWell))
			singlyOccupiedRepresenativeWells[-1][1]['sudokuGridEntry']['condensationType'] = \
			'Rearray'
			totalWellsForRearrayOfSinglyOccupiedWells += coloniesToPick

	totalWellsForRearrayOfSinglyOccupiedWells = int(totalWellsForRearrayOfSinglyOccupiedWells)
	totalWellsForRearrayOfMultiplyOccupiedWells = int(totalWellsForRearrayOfMultiplyOccupiedWells)

	if rearrayPickOrder == 'columns':
		secondSortIndex = 'col'
		thirdSortIndex = 'row'
	elif rearrayPickOrder == 'rows':
		secondSortIndex = 'row'
		thirdSortIndex = 'col'
	else:
		secondSortIndex = 'col'
		thirdSortIndex = 'row'


	# Sort the singly occupied wells array
	singlyOccupiedRepresenativeWells = \
	sorted(singlyOccupiedRepresenativeWells, key = lambda x: \
	(x[1]['sudokuGridEntry']['plateName'], x[1]['sudokuGridEntry'][secondSortIndex], \
	x[1]['sudokuGridEntry'][thirdSortIndex]))

	# Sort the multiply occupied wells array
	multiplyOccupiedRepresentativeWells = \
	sorted(multiplyOccupiedRepresentativeWells, key = lambda x: \
	(x[1]['sudokuGridEntry']['plateName'], x[1]['sudokuGridEntry'][secondSortIndex], \
	x[1]['sudokuGridEntry'][thirdSortIndex]))

	
	return [singlyOccupiedRepresenativeWells, multiplyOccupiedRepresentativeWells]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #










# ------------------------------------------------------------------------------------------------ #
# Pick Best Mutant To Disrupt Feature
def PickBestSudokuWellToDisruptFeature(featureArray, sudokuGridLookupDict, pickProbability=0.95):
	
	import pdb
	from operator import itemgetter
	from math import floor, ceil
	from numpy import log
	
	featuresAndBestWellArray = []
	
	i = 0
	while i < len(featureArray):
		
		feature = featureArray[i]
		featureLocusTag = feature.tagDict['locus_tag']

		sudokuGridEntries = feature.sudokuGridEntries
		sudokuWellObjectsToScore = []
		
		j = 0
		while j < len(sudokuGridEntries) and len(sudokuGridEntries) > 0:
				
			addressDictKeys = list(sudokuGridEntries[j]['addressDict'].keys())
			
			if 'progenitor' in addressDictKeys:
				plateCol = sudokuGridEntries[j]['addressDict']['progenitor']['plateCol']
				plateRow = sudokuGridEntries[j]['addressDict']['progenitor']['plateRow']
				row = sudokuGridEntries[j]['addressDict']['progenitor']['row']
				col = sudokuGridEntries[j]['addressDict']['progenitor']['col']
				wellOccupants = sudokuGridEntries[j]['wellOccupants']
				fracDistFromTranslationStart = sudokuGridEntries[j]['fracDistFromTranslationStart']
				locatability = sudokuGridEntries[j]['locatability']
				
				sudokuWell = sudokuGridLookupDict[plateRow][plateCol].wellGrid[row][col]
				
				if locatability == 'unambiguous':
					locatabilityDiscountFactor = 1.0
				else:
					locatabilityDiscountFactor = 0.5
				
				purifiability = \
				(1.0 - fracDistFromTranslationStart)*locatabilityDiscountFactor*(1.0/wellOccupants)
				
				pA = (1.0/(float(wellOccupants)))
		
				if pA == 1.0:
					minimumColoniesToPick = 1.0
				else:
					minimumColoniesToPick = floor(log(1-pickProbability)/log(1-pA))
				
				sudokuWellObjectsToScore.append({'well':sudokuWell, 'purifiability':purifiability,\
				'sudokuGridEntry':sudokuGridEntries[j], \
				'minimumColoniesToPick':minimumColoniesToPick})
				
			j += 1
				
		
		sudokuWellObjectsToScore = sorted(sudokuWellObjectsToScore, reverse=True, \
		key=itemgetter('purifiability'))
		bestWell = sudokuWellObjectsToScore[0]

		featuresAndBestWellArray.append([feature, bestWell])
		
		i += 1

	return featuresAndBestWellArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignPlateGridCoordinatesToCondensedCollectionPlates(condensedWellArray):
# This function assigns plate row and plate column coordinates to plates in the condensed 
# collection. 

	from copy import deepcopy
	from scipy import unique
	from numpy import sqrt
	from math import ceil, floor


	wellList = deepcopy(condensedWellArray)

	# Figure out how many plates we have in the condensed collection. The names don't need to 
	# be continuous or even numbers, just sortable. 
	j = 0
	plateNames = []
	while j < len(wellList):
		plateNames.append(wellList[j].plateName)
		j += 1

	uniquePlateNames = unique(plateNames)
	lengthPlateNames = len(uniquePlateNames)
	
	# Figure out the dimension of the plate grid
	sqrtLenPlates = sqrt(lengthPlateNames)
	plateCols = ceil(sqrtLenPlates)
	plateRows = floor(sqrtLenPlates)
	
	while plateRows*plateCols < lengthPlateNames:
		plateRows += 1


	# Generate the keys for the plate rows and plate columns
	i = 1
	plateRowKeys = []
	while i <= plateRows:
		plateRowKeys.append('PR' + str(i).zfill(2))
		i += 1

	i = 1
	plateColKeys = []
	while i <= plateCols:
		plateColKeys.append('PC' + str(i).zfill(2))
		i += 1
		
	# Assign the plate names to positions in the new plate grid. 
	uniquePlateNames = sorted(uniquePlateNames, reverse=True)
	uniquePlateNamesCopy = deepcopy(uniquePlateNames)

	plateGridLookupDict = {}

	continueAdding = True
	i = 0
	while i < len(plateRowKeys):
		plateGridLookupDict[plateRowKeys[i]] = {}
		j = 0 
		while j < len(plateColKeys) and len(uniquePlateNamesCopy) > 0:
			plateGridLookupDict[plateRowKeys[i]][plateColKeys[j]] = uniquePlateNamesCopy.pop()
			j += 1
		i += 1

	# Use the plate grid lookup table to assign plate grid coordinates to each well in the 
	# condensed collection. 
	i = 0
	while i < len(wellList):
		plateName = wellList[i].plateName
	
		plateNameMatchFound = False
		j = 0
		while j < len(plateRowKeys) and plateNameMatchFound == False:
			k = 0
		
			tempPlateColKeys = sorted(list(plateGridLookupDict[plateRowKeys[j]].keys()))
		
			while k < len(tempPlateColKeys):
				testPlateName = str(plateGridLookupDict[plateRowKeys[j]][tempPlateColKeys[k]])
				if str(testPlateName) == str(plateName):
					wellList[i].plateCol = plateColKeys[k]
					wellList[i].plateRow = plateRowKeys[j]
					wellList[i].addressDict['condensed']['plateCol'] = plateColKeys[k]
					wellList[i].addressDict['condensed']['plateRow'] = plateRowKeys[j]
					plateNameMatchFound = True
					break
				k += 1
			j += 1
		if plateNameMatchFound == False:
			print('Plate name match not found for wellList entry ' + str(i))
		i += 1
	
	return wellList
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def WriteSudokuGridCondensationInstructions(condensationFileName, condensedWellArray):
	
	from kosudoku.xml import ExportListForXML
	import pdb
	
	headers = [\
	'#Condensed Collection Plate', 'Condensed Collection Row', \
	'Condensed Collection Column', 'Condensed Collection Well', \
	'Condensed Collection Plate Row', 'Condensed Collection Plate Column', \
	'Progenitor Collection Plate', 'Progenitor Collection Row', \
	'Progenitor Collection Column', 'Progenitor Collection Well', \
	'Condensation Type', \
	'Desired Transposon Coordinate', \
	'Desired Feature', 'Desired Locus Tag']
	
	headerStr = ExportListForXML(headers) + '\n'
	
	outputStr = ''
	i = 0
	while i < len(condensedWellArray):
		
		condensedWell = condensedWellArray[i]
		
		addressDict = condensedWell.addressDict
		try:
			progAddressDict = addressDict['progenitor']
		except:
			pdb.set_trace()
		
		
		try:
			condAddressDict = addressDict['condensed']
		except:
			pdb.set_trace()
		
		condCollectionPlate = condAddressDict['plateName']
		condCollectionRow = condAddressDict['row']
		condCollectionCol = condAddressDict['col']
		condCollectionWell = str(condCollectionRow) + str(condCollectionCol).zfill(2)
		condCollectionPR = condAddressDict['plateRow']
		condCollectionPC = condAddressDict['plateCol']
		
		progCollectionPlate = progAddressDict['plateName']
		progCollectionRow = progAddressDict['row']
		progCollectionCol = progAddressDict['col']
		if progCollectionRow == 'NA' and progCollectionCol == 'NA':
			progCollectionWell = 'NA'
		else:	
			progCollectionWell = str(progCollectionRow) + str(progCollectionCol).zfill(2)
		
		
		try:
			condensationType = condensedWell.condensationData['condensationType']
		except:
			pdb.set_trace()
			
		desiredTransposonCoord = condensedWell.condensationData['readAlignmentCoord']
		readAlignmentCoords = condensedWell.readAlignmentCoords
		
		desiredFeature = ''
		desiredLocusTag = ''
		for coord in readAlignmentCoords:
			if desiredTransposonCoord == coord.coord:
				desiredFeature = coord.featureName[0]
				desiredLocusTag = coord.locusTag[0]


		outputLineArray = [condCollectionPlate, condCollectionRow, condCollectionCol, \
		condCollectionWell, condCollectionPR, condCollectionPC, \
		progCollectionPlate, progCollectionRow, progCollectionCol, progCollectionWell, \
		condensationType, desiredTransposonCoord, desiredFeature, \
		desiredLocusTag]
		
		lineStr = ExportListForXML(outputLineArray) + '\n'
		outputStr += lineStr
		
		i += 1
	
	outputHandle = open(condensationFileName, 'w')
	outputHandle.write(headerStr)
	outputHandle.write(outputStr)
	outputHandle.close()
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportPredictedCondensedCollectionCatalog(filename):
# Imports data from the kosudoku-condense run for use in kosudoku-isitinthere. 
# Assumes that data can be found in the order:

# Outputs a dictionary of the condensed collection:

# condensedCatalog[condensed collection location key] = 
# {progenitor plate:, progenitor row:, progenitor col:, hoped for transposon:}

# At this stage, the progenitor contents don't need to be filled in. 
# This will be with a code called UpdateColonyPurificationDictWithProgenitorCollection


	import pdb
	from kosudoku.utils import ParseCSVLine
	from ast import literal_eval
	
	fileHandle = open(filename, 'r')
	
	data = fileHandle.readlines()
	
	sudokuCatalog = {}
	
	i = 0
	while i < len(data):
		line = data[i]
		
		if line[0] != '#':
			
			columns = ParseCSVLine(data[i])
			entryDict = {}
			
			entryDict['cpPlate'] = columns[0] 
			entryDict['cpRow'] = columns[1]
			entryDict['cpCol'] = columns[2]
			entryDict['cpWell'] = columns[3]
			entryDict['cpPlateRow'] = columns[4]
			entryDict['cpPlateCol'] = columns[5]
				
			entryDict['progenitorPlate'] = columns[6]
			entryDict['progenitorRow'] = columns[7]
			entryDict['progenitorCol'] = columns[8]
			entryDict['progenitorWell'] = columns[9]
			entryDict['condensationType'] = columns[10]
			entryDict['hopedForTransposonCoord'] = columns[11]
			entryDict['hopedForFeature'] = columns[12]
			entryDict['hopedForLocusTag'] = columns[13]
			entryDict['progenitorContents'] = []
			
			cpKey = entryDict['cpPlate'] + '_' + entryDict['cpWell']

			sudokuCatalog[cpKey] = entryDict
		
		i += 1
	
	return sudokuCatalog
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class SudokuCondensedWell(SudokuWell):

	def __init__(self, sourcePlateName, sourcePlateRow, sourcePlateCol, sourceRow, \
	sourceCol, addressSystem='condensed', OD=1.0):
				
		SudokuWell.__init__(self, sourcePlateName, sourcePlateRow, sourcePlateCol, sourceRow, \
		sourceCol, addressSystem=addressSystem, OD=1.0)

		self.condensationData = {}
	
# ------------------------------------------------------------------------------------------------ #

# ----------------------------------------------------------------------------------------------- #
def TestIfCondensedCollectionPlateWellIsBlank(currentPlate, currentRowNumber, currentColNumber, \
pattern):

	# Returns if true is well is supposed to be blank, false if it is supposed to have something in
	# it.
	
	if pattern == 'checkerboard':
		# Currently, with a checkerboard, the plate number doesn't matter at all. Basically,
		# wells that an have odd row and odd column or even row and even column numbers (e.g. A1 
		# or B2) are blank, and all others are filled. 

		if currentColNumber%2 == 0:
			# For even numbered columns
			if currentRowNumber%2 == 0:
				# This is an even row in an even column, for instance B2, which should be blank
				blankTest = True
			else:
				# This is an odd row in an even column, for instance B1, which should have something
				# in it
				blankTest = False
		else:
			# For odd numbered columns
			if currentRowNumber%2 == 0:
				# This is an even row in an odd column, for instance A2, which should have something
				# in it
				blankTest = False
			else:
				# This is an odd row in an odd column, for instance A1, which should be blank
				blankTest = True
		
	elif pattern == 'all':
		# If every well is filled, then no well should be blank. 
		blankTest = False
	else:
		print('Currently only doing re-array into all wells or checkerboard.')
		return
	
	return blankTest
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
def IncrementWellInColumnsOrder(currentPlateNumber, currentRowNumber, currentColNumber, rows, cols):
# What I mean here is that the current well is incremented along columns in the plate. For example:
# A1, B1, C1, D1, .... which is I think the most natural way to fill a plate. 		
	
	if currentRowNumber == len(rows):
		nextRowNumber = 1
		nextColNumber = currentColNumber+1
	else:
		nextRowNumber = currentRowNumber + 1
		nextColNumber = currentColNumber
	
	if nextColNumber == len(cols) + 1:
		nextPlateNumber = currentPlateNumber + 1
		nextColNumber = 1
	else:
		nextPlateNumber = currentPlateNumber
			
	return [nextPlateNumber, nextRowNumber, nextColNumber]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def IncrementWellInRowsOrder(currentPlateNumber, currentRowNumber, currentColNumber, rows, cols):
# What I mean here is that the current well is incremented along rows in the plate. For example:
# A1, A2, A3, A4, .... which is I think is a slightly unnatural way to fill a plate. 		
	
	if currentColNumber == len(cols):
		nextColNumber = 1
		nextRowNumber = currentRowNumber + 1
	else:
		nextRowNumber = currentRowNumber
		nextColNumber = currentColNumber + 1
	
	if nextRowNumber == len(rows) + 1:
		nextPlateNumber = currentPlateNumber + 1
		nextRowNumber = 1
	else:
		nextPlateNumber = currentPlateNumber
		
	return [nextPlateNumber, nextRowNumber, nextColNumber]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def IncrementWell(currentPlateNumber, currentRowNumber, currentColNumber, rows, cols, fillOrder):
	if fillOrder == 'rows':
		[nextPlateNumber, nextRowNumber, nextColNumber] = \
		IncrementWellInRowsOrder(currentPlateNumber, currentRowNumber, currentColNumber, rows96, \
		columns96)
	elif fillOrder == 'columns':
		[nextPlateNumber, nextRowNumber, nextColNumber] = \
		IncrementWellInColumnsOrder(currentPlateNumber, currentRowNumber, currentColNumber, \
		rows96, columns96)

	return [nextPlateNumber, nextRowNumber, nextColNumber]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Some useful global variables for use in condensation
blankSudokuCoord = SudokuGenomicCoord(-1, 'NA', 'NA', 'NA')
blankSudokuCoord.locusTag.append('BLANK')
blankSudokuCoord.featureName.append('BLANK')
blankSudokuCoord.distanceFromFeatureTranslationStart.append('NA')
blankSudokuCoord.fracDistanceFromFeatureTranslationStart.append('NA')

blankSudokuGridEntrySummary = {\
'addressDict': {'progenitor': {'col': 'NA','plateCol': 'NA','plateName': 'NA', 'plateRow': 'NA',\
'row': 'NA'}},\
'col': 'NA', 'distFromTranslationStart': 'NA', 'fracDistFromTranslationStart': 'NA', \
'locatability': 'NA', 'multipleOccupancy': 'NA','plateCol': 'NA', \
'plateName': 'NA','plateRow': 'NA','readAlignmentCoord': -1,'row': 'NA', 'wellOccupants': 0, \
'condensationType':'NA'}

columns96 = [1,2,3,4,5,6,7,8,9,10,11,12]
rows96 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
wellsPer96WellPlate = len(columns96)*len(rows96)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindSingleRepresentativeForColonyPurifiedMutantsInCondensedCollection(cdSudokuWellArray, \
rearrayPickOrder='columns'):
# Used in the kosudoku-qc program
	
# How is this supposed to work
	return

# ------------------------------------------------------------------------------------------------ #
