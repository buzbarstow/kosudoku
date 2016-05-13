# kosudoku-condense2.py

from numpy import ceil
from kosudoku.grid import SudokuWell, SudokuGenomicCoord
from kosudoku.condense import rows96, columns96, SudokuCondensedWell, \
TestIfCondensedCollectionPlateWellIsBlank, blankSudokuCoord, blankSudokuGridEntrySummary, \
IncrementWell
from operator import itemgetter
from copy import deepcopy
import pdb





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
def WriteSudokuGridCondensationInstructions(condensationFileName, condensedWellArray):
	
	from kosudoku.xml import ExportListForXML
	import pdb
	
	headers = [\
	'Condensed Collection Plate', 'Condensed Collection Row', \
	'Condensed Collection Column', 'Condensed Collection Well', \
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
		condCollectionWell, progCollectionPlate, progCollectionRow, progCollectionCol, \
		progCollectionWell, condensationType, desiredTransposonCoord, desiredFeature, \
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









featuresAndBestWellArrayCopy = deepcopy(featuresAndBestWellArray)

# ----------------------------------------------------------------------------------------------- #
# Separate out the representatives into ones that are singly occupied and ones that are multiply
# occupied.
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
		singlyOccupiedRepresenativeWells[-1][1]['sudokuGridEntry']['condensationType'] = 'Rearray'
		totalWellsForRearrayOfSinglyOccupiedWells += coloniesToPick

totalWellsForRearrayOfSinglyOccupiedWells = int(totalWellsForRearrayOfSinglyOccupiedWells)
totalWellsForRearrayOfMultiplyOccupiedWells = int(totalWellsForRearrayOfMultiplyOccupiedWells)
# ----------------------------------------------------------------------------------------------- #

# Let's figure out how many plates we need. In future versions we can actually calculate these 
# numbers algorithmically, but for now, we will calculate them with human power. 


# Sort the singly occupied wells array
singlyOccupiedRepresenativeWells = \
sorted(singlyOccupiedRepresenativeWells, key = lambda x: (x[1]['sudokuGridEntry']['plateName'], \
x[1]['sudokuGridEntry']['col'], x[1]['sudokuGridEntry']['row']))

# Sort the multiply occupied wells array
multiplyOccupiedRepresentativeWells = \
sorted(multiplyOccupiedRepresentativeWells, key = lambda x: (x[1]['sudokuGridEntry']['plateName'], \
x[1]['sudokuGridEntry']['col'], x[1]['sudokuGridEntry']['row']))

# Construct a flattened well list

fillOrder = 'columns'
fillPattern = 'checkerboard'



# Assign the colony purified mutants from multiply occupied wells
[condensedWellArrayForSingles, currentPlateNumber, currentRowNumber, currentColNumber] = \
AssignProgenitorWellsToCondensedCollection(singlyOccupiedRepresenativeWells, 1, 1, 1, rows96, \
columns96, fillOrder=fillOrder, fillPattern=fillPattern, fillLastPlateWithBlanks=True)

# Assign the colony purified mutants from multiply occupied wells
[condensedWellArrayForMultiples, currentPlateNumber, currentRowNumber, currentColNumber] = \
AssignProgenitorWellsToCondensedCollection(multiplyOccupiedRepresentativeWells, currentPlateNumber,\
currentRowNumber, currentColNumber, rows96, \
columns96, fillOrder=fillOrder, fillPattern=fillPattern, fillLastPlateWithBlanks=True)


condensedWellArray = condensedWellArrayForSingles + condensedWellArrayForMultiples


WriteSudokuGridCondensationInstructions(rearrayInstructionsFileName, condensedWellArray)