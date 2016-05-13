# from kosudoku.xml import ExportListForXML



# ------------------------------------------------------------------------------------------------ #
class SudokuPlate:

	def __init__(self, plateName, plateRow, plateCol, rowPools, colPools):
		
		self.plateName = plateName
		self.rowPools = rowPools
		self.colPools = colPools
		self.plateRow = plateRow
		self.plateCol = plateCol
		
		self.wellGrid = self.generateWellGrid(plateName, plateRow, plateCol, rowPools, colPools)
	
	def generateWellGrid(self, plateName, plateRow, plateCol, rowPools, colPools):
		
		wellGrid = {}
		
		i = 0
		while i < len(rowPools):
			wellGrid[rowPools[i]] = {}
			j = 0
			while j < len(colPools):
				wellGrid[rowPools[i]][colPools[j]] = SudokuWell(plateName, plateRow, plateCol, \
				rowPools[i], colPools[j])
				j += 1
			i += 1
		
		return wellGrid
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class SudokuWell:
	def __init__(self, plateName, plateRow, plateCol, row, col, OD=1.0, addressSystem='progenitor'):
		self.plateName = plateName
		self.row = row
		self.col = col
		self.libraryAddress = str(plateRow) + '_' + str(plateCol) + '_' +  str(row) + '_' + str(col)
		self.OD = OD
		self.readAlignmentCoords = []
		self.plateRow = plateRow
		self.plateCol = plateCol
		
		self.addressDict = {}
		self.addressDict[addressSystem] = {'plateName':plateName, 'plateRow':plateRow, \
		'plateCol':plateCol, 'row':row, 'col':col}

 
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
class SudokuGenomicCoord:
	def __init__(self, coord, locatability, locatabilityScore, readCount):
		self.coord = coord
		self.locatability = locatability
		self.locatabilityScore = locatabilityScore
		self.readCount = readCount
		self.featureName = []
		self.distanceFromFeatureTranslationStart = []
		self.fracDistanceFromFeatureTranslationStart = []
		self.locusTag = []

# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ImportSudokuGridLayout(sudokuGridLayoutFile, rowPools, colPools):
	
	fileHandle = open(sudokuGridLayoutFile, 'r')
	
	data = fileHandle.readlines()
	
	pcPools = data[0].strip().split(',')
	pcPools.pop(0)

	
	prPools = []
	plateNameGrid = []
	
	i = 1
	while i < len(data):
		
		dataLine = data[i].strip().split(',')
		
		prPools.append(dataLine[0])
		
		plateNameGrid.append([])
		
		j = 1
		while j < len(dataLine):
			 plateNameGrid[i-1].append(dataLine[j])
			 j += 1
			 
		i += 1
		
	# The sudokuGridLookupDict is very similar to the plateRowLookupDict that was hard-coded into
	# earlier versions of the sudoku code. It is initially indexed by the plate row, followed by 
	# plate columns
	
	sudokuGridLookupDict = {}
	
	i = 0
	while i < len(prPools):
		sudokuGridLookupDict[prPools[i]] = {}
		j = 0
		while j < len(pcPools):
		
			plateName = plateNameGrid[i][j]
			
			sudokuGridLookupDict[prPools[i]][pcPools[j]] = SudokuPlate(plateName, \
			prPools[i], pcPools[j], rowPools, colPools)

			j += 1
		i += 1

	
	return [sudokuGridLookupDict, prPools, pcPools]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def InitializeEmptySudokuGridLookupDict(prPools, pcPools, rowPools, colPools):
	
	sudokuGridLookupDict = {}
	
	i = 0
	while i < len(prPools):
		sudokuGridLookupDict[prPools[i]] = {}
		j = 0
		while j < len(pcPools):
		
			sudokuGridLookupDict[prPools[i]][pcPools[j]] = SudokuPlate('Blank', \
			prPools[i], pcPools[j], rowPools, colPools)

			j += 1
		i += 1
	
	return sudokuGridLookupDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def AddAddressToSudokuGrid(coord, locatability, possibleAddressesAndScoresEntry, \
sudokuGridLookupDict):
	
	addressCoords = possibleAddressesAndScoresEntry[0].split('_')
	locatabilityScore = possibleAddressesAndScoresEntry[2]
	readCount = possibleAddressesAndScoresEntry[3]
					 					
	row = addressCoords[0]
	col = addressCoords[1]
	pr = addressCoords[2]
	pc = addressCoords[3]
		
	sudokuCoord = SudokuGenomicCoord(coord, locatability, locatabilityScore, readCount)
	
	sudokuGridLookupDict[pr][pc].wellGrid[row][col].readAlignmentCoords.append(sudokuCoord)
	
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, threshold, \
maxEntriesInSingleAddressPoolAxis, maxTotalCoords):
	
	import pdb
	
	[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
	CalculatePoolCoordsForLine(\
	poolPresenceTableLine, rowPools, colPools, prPools, pcPools, controlPools, \
	threshold=threshold)
	
	[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
	CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
	addresses_pr, addresses_pc)
	
	nTotalCoords = nRowCoords + nColCoords + nPRCoords + nPCCoords
	
	[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
	CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, \
	addresses_pc, addresses_control)

	nAddressPoolAxesWithMoreThanOneEntry = \
	CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, \
	addresses_pr, addresses_pc, addresses_control)
	
	maxSinglePoolCoordNumber = \
	CalculateMaxNumberOfEntriesInSinglePoolAxis(addresses_r, addresses_c, \
	addresses_pr, addresses_pc)

	# Decide on the locatability of the genomic coordinate
	if nAddressPoolAxesWithEntries < 3:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None
		
	
	elif maxSinglePoolCoordNumber > maxEntriesInSingleAddressPoolAxis or \
	nTotalCoords > maxTotalCoords:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None
	
	elif nAddressPoolAxesWithEntries == 3:
		if nAddressPoolAxesWithMoreThanOneEntry > 1:
			locatability = 'unlocatable'
			possibleAddressesAndScores = None
		elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
			locatability = 'guessable'
			possibleAddressesAndScores = None
	
	elif nAddressPoolAxesWithEntries == 4:
		
		possibleAddressesAndScores = \
		CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, fitDict, \
		areaDict, threshold)
		
		if len(possibleAddressesAndScores) == 0:
			pdb.set_trace()
	
		if nAddressPoolAxesWithMoreThanOneEntry > 1:
			locatability = 'ambiguous'
		elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
			locatability = 'unambiguous'
	
	else:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None

	return [possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
locatability):
	j = 0
	while j < len(possibleAddressesAndScores):
		AddAddressToSudokuGrid(coord, locatability, possibleAddressesAndScores[j], \
		sudokuGridLookupDict)
		j += 1

	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def AssignAmbiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
locatability, maxEntriesInSingleAddressPoolAxis):

	import operator
	import pdb
	
	sortedPossibleAddressesAndScores = sorted(possibleAddressesAndScores, \
	key=operator.itemgetter(2), reverse=True)
	
	# if  coord==5038711:
# 		pdb.set_trace()
	
	j = 0
	continueAddingAddressesToSudokuGrid = True
	
	while continueAddingAddressesToSudokuGrid == True:
		
		try:
			addressAndScore = sortedPossibleAddressesAndScores[j]
		except IndexError:
			pdb.set_trace()
		
		
		AddAddressToSudokuGrid(coord, locatability, addressAndScore, sudokuGridLookupDict)
		
		j += 1
		
		if (j < maxEntriesInSingleAddressPoolAxis) \
		and j < len(sortedPossibleAddressesAndScores):
			continueAddingAddressesToSudokuGrid = True
		else:
			continueAddingAddressesToSudokuGrid = False

	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignGuessableAddresses(sudokuGridLookupDict, coord, poolPresenceTableLine, \
rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, \
maxSinglePoolCoordNumber, maxTotalCoords):
	
	temporaryThreshold = 1
	validAddressesFound = False
	
	while validAddressesFound == False and temporaryThreshold <= readCountThreshold:
		
		[possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber] = \
		CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, temporaryThreshold, \
		maxSinglePoolCoordNumber, maxTotalCoords)
		
		if locatability == 'unambiguous':
			validAddressesFound = True
		else:
			temporaryThreshold += 1
	
	if validAddressesFound == True:
		AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
		'unambiguous')
		
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PopulateSudokuGrid3(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, scoreThreshold, \
maxEntriesInSingleAddressPoolAxis, maxTotalCoords):

	import pdb
	
	i = 0
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
		
# 		pdb.set_trace()
		
		coord = poolPresenceTableLine['readAlignmentCoord']
		
		[possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber] = \
		CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, \
		maxEntriesInSingleAddressPoolAxis, maxTotalCoords)
		
		# if  (coord==5038710) or (coord==5038711) or (coord==5038712):
# 			pdb.set_trace()
	
		if locatability == 'unambiguous':
			AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
			locatability)
			
		elif locatability == 'ambiguous':
			AssignAmbiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
			locatability, maxSinglePoolCoordNumber)
		
		elif locatability == 'guessable':
			AssignGuessableAddresses(sudokuGridLookupDict, coord, poolPresenceTableLine, \
			rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, \
			readCountThreshold, maxEntriesInSingleAddressPoolAxis, maxTotalCoords)
							
		i += 1
	
	return 
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
class GridTypeError(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ConvertPlateNumberToPlateRowAndPlateColCoords(plateNumber, sudokuGridLookupDict, prPools, \
pcPools):
	import re	

	plateNameRe = re.compile(r'\d+')
	
	matchingPlateFound = False
	matchingPlateCol = None
	matchingPlateRow = None
	
	i = 0
	while (i < len(prPools)) and (matchingPlateFound == False):
		plateRow = prPools[i]
		j = 0
		while (j < len(pcPools)) and (matchingPlateFound == False):
			plateCol = pcPools[j]
			plate =  sudokuGridLookupDict[plateRow][plateCol].plateName
			plateNameMatch = plateNameRe.match(plate)

			if plateNameMatch != None:
				if int(plate) == int(plateNumber):
					matchingPlateFound = True
					matchingPlateCol = plateCol
					matchingPlateRow = plateRow
			j += 1
		i += 1
	
	return [matchingPlateRow, matchingPlateCol]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindSudokuWellByPlateName(sudokuGridLookupDict, prPools, pcPools, rowPools, colPools, \
sourcePlate, sourceRow, sourceCol, copy=True):
	
	# Return a deep copy of sudoku well object after finding it using a plate name reference
	
	import pdb
	
	from copy import deepcopy
	
	plateFound = False
	
	i = 0
	while i < len(prPools) and plateFound == False:
		prPool = prPools[i]
		
		j = 0
		while j < len(pcPools) and plateFound == False:
			pcPool = pcPools[j]
			plateName = sudokuGridLookupDict[prPool][pcPool].plateName
			
			if str(plateName).zfill(0) == str(sourcePlate).zfill(0):
# 				print('Search Plate: ' + str(sourcePlate) + ', Found Plate: ' + str(plateName) \
# 				+ ', Well: ' + str(sourceRow) + str(sourceCol))
				wellGrid = sudokuGridLookupDict[prPool][pcPool].wellGrid
				
				if copy == True:
					sudokuWell = deepcopy(wellGrid[sourceRow][str(sourceCol).zfill(0)])
				else:
					sudokuWell = wellGrid[sourceRow][str(sourceCol).zfill(0)]
			j += 1
		i += 1
	
	return sudokuWell
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def AnnotateSudokuGridWithFeautureNames(sudokuGridLookupDict, featureArray, prPools, pcPools, \
rowPools, colPools):

	import pdb

	for prPool in prPools:
		for pcPool in pcPools:
			plate = sudokuGridLookupDict[prPool][pcPool]
			wellGrid = plate.wellGrid
			
			for row in rowPools:
				for col in colPools:
					well = wellGrid[row][col]
					readAlignmentCoords = well.readAlignmentCoords
					
					for readAlignmentCoord in readAlignmentCoords:
						coord = readAlignmentCoord.coord
						
						for feature in featureArray:
							startCoord = feature.startCoord
							endCoord = feature.endCoord
														
							if startCoord < coord < endCoord:
								
								try:
									locusTag = feature.tagDict['locus_tag'][0]
								except:
									locusTag = ''
								
								featureName = feature.featureName
								
								startTranslation = feature.startTranslation
								endTranslation = feature.endTranslation
								translationLength = abs(endTranslation-startTranslation)
							
								distanceFromFeatureTranslationStart = abs(coord - startTranslation)
								
								fracDistanceFromFeatureTranslationStart = \
								distanceFromFeatureTranslationStart/translationLength
								
								readAlignmentCoord.featureName.append(featureName)
								
								try:
									readAlignmentCoord.distanceFromFeatureTranslationStart.append(\
									"{:.0f}".format(int(distanceFromFeatureTranslationStart)))
								except:
									pdb.set_trace()
								
								readAlignmentCoord.fracDistanceFromFeatureTranslationStart.append(\
								"{:.2f}".format(float(fracDistanceFromFeatureTranslationStart)))
								
								readAlignmentCoord.locusTag.append(locusTag)
							
						lenLocusTag = len(readAlignmentCoord.locusTag)
						lenFeatureName = len(readAlignmentCoord.featureName)
						lenFracDist = len(\
						readAlignmentCoord.fracDistanceFromFeatureTranslationStart)
						lenDist = len(readAlignmentCoord.distanceFromFeatureTranslationStart)
						
						if lenLocusTag == 0:
							readAlignmentCoord.locusTag.append('None')
						if lenFeatureName == 0:
							readAlignmentCoord.featureName.append('None')
						if lenFracDist == 0:
							readAlignmentCoord.fracDistanceFromFeatureTranslationStart.append('NA')
						if lenDist == 0:
							readAlignmentCoord.distanceFromFeatureTranslationStart.append('NA')
	
								
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, rowPools, colPools):
	
	flattenedWellList = []
	
	for prPool in prPools:
		for pcPool in pcPools:
			for rowPool in rowPools:
				for colPool in colPools:
					well = sudokuGridLookupDict[prPool][pcPool].wellGrid[rowPool][colPool]
					flattenedWellList.append(well)

	return flattenedWellList
# ------------------------------------------------------------------------------------------------ #





# ----------------------------------------------------------------------------------------------- #
def WriteSudokuGridSummaryTable(sudokuGridSummaryFileName, sudokuGridLookupDict, \
prPools, pcPools, rowPools, colPools):
# Modification of WriteSudokuGridSummary that allows output of feature names as well
	import pdb
	from .xml import ExportListForXML
	
	fileHandle = open(sudokuGridSummaryFileName, 'w')
	outputStr = ''
	outputStr += 'Plate,Plate Row,Plate Column,Well,Transposon Coordinates,' \
	+ 'Transposon Locatabilities,' \
	+ 'Features Disrupted,Distances of Transposon from Translation Starts,' \
	+ 'Fractional Distances of Transposon from Translation Starts,Locus Tags\n'
	
	for prPool in prPools:
		for pcPool in pcPools:
			plate = sudokuGridLookupDict[prPool][pcPool]
			
			plateName = plate.plateName
			wellGrid = plate.wellGrid
			
			for rowPool in rowPools:
				for colPool in colPools:
					well = wellGrid[rowPool][colPool]
					readAlignmentCoords = well.readAlignmentCoords
					
					
					coords = []
					locatabilities = []
					featureNameStrings = []
					distsFromFeatureTranslationStartStrings = []
					fracDistsFromFeatureTranslationStartStrings = []
					locusTagsStrings = []
					
					i = 0
					while i < len(readAlignmentCoords):
						coord = readAlignmentCoords[i]
						coords.append(coord.coord)
						locatabilities.append(coord.locatability)
						
						featureNames = coord.featureName
						featureNamesStr = ExportListForXML(featureNames, delimeter='\\')
						featureNameStrings.append(featureNamesStr)
						
						
						distsFromFeatureTranslationStart = \
						coord.distanceFromFeatureTranslationStart
						
						distsFromFeatureTranslationStartStr = \
						ExportListForXML(distsFromFeatureTranslationStart, delimeter='\\')
						
						distsFromFeatureTranslationStartStrings.append(\
						distsFromFeatureTranslationStartStr)
						
						
						fracDistsFromFeatureTranslationStart = \
						coord.fracDistanceFromFeatureTranslationStart
						
						fracDistsFromFeatureTranslationStartStr = \
						ExportListForXML(fracDistsFromFeatureTranslationStart, delimeter='\\')
						
						fracDistsFromFeatureTranslationStartStrings.append(\
						fracDistsFromFeatureTranslationStartStr)
						
						
						locusTags = coord.locusTag
						locusTagsStr = ExportListForXML(locusTags, delimeter='\\')
						locusTagsStrings.append(locusTagsStr)
						
						
						i += 1
							
					featureNameStrings = '"' + ExportListForXML(featureNameStrings, delimeter=', ')\
					+ '"'
					
					coordsStr = '"' + ExportListForXML(coords, delimeter=', ') + '"'
					locatabilitiesStr = '"' + ExportListForXML(locatabilities, delimeter=', ') + '"'
					
					distsStr = '"' + ExportListForXML(distsFromFeatureTranslationStartStrings, \
					delimeter=',')  + '"'
					
					fracDistsStr = '"' + ExportListForXML(\
					fracDistsFromFeatureTranslationStartStrings, delimeter=', ')  + '"'
					
					locusTagStr = '"' + ExportListForXML(locusTagsStrings, delimeter=', ')  + '"'
					
					
					outputStr += plateName + ',' + prPool + ',' + pcPool +  ',' + rowPool \
					+ colPool + ',' + coordsStr + ',' \
					+ locatabilitiesStr + ',' + featureNameStrings + ',' + distsStr + ',' \
					+ fracDistsStr + ',' + locusTagStr + '\n'
	
	fileHandle.write(outputStr)
	fileHandle.close()
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateSudokuGridOccupancyTaxonomy(sudokuGridLookupDict, rowPools, colPools, \
prPools, pcPools):

	from scipy import unique
	import pdb

	taxonomyDict = {'emptyWells':0, 'singlyOccupiedWells':0, 'multiplyOccupiedWells':0, \
	'totalMutants':0, 'mutantsInSinglyOccupiedWells':0, 'mutantsInMultiplyOccupiedWells':0}

	collectedGenomicCoords = []
	mutantsInSinglyOccupiedWells = []
	mutantsInMultiplyOccupiedWells = []

	i = 0
	while i < len(prPools):
		j = 0
		while j < len(pcPools):
		
			k = 0
			while k < len(rowPools):
				l = 0
				while l < len(colPools):
				
					pr = prPools[i]
					pc = pcPools[j]
					row = rowPools[k]
					col = colPools[l]
					
					well = sudokuGridLookupDict[pr][pc].wellGrid[row][col]
					
					readAlignmentCoords = well.readAlignmentCoords
					
					nGenomicCoords = len(readAlignmentCoords)
					
					for coord in readAlignmentCoords:
# 						pdb.set_trace()
						collectedGenomicCoords.append(int(coord.coord))
					
					if nGenomicCoords == 0:
						taxonomyDict['emptyWells'] += 1
					elif nGenomicCoords == 1:
						taxonomyDict['singlyOccupiedWells'] += 1
						taxonomyDict['mutantsInSinglyOccupiedWells'] += 1
						
						for coord in readAlignmentCoords:
							mutantsInSinglyOccupiedWells.append(int(coord.coord))
						
					elif nGenomicCoords > 1:
						taxonomyDict['multiplyOccupiedWells'] += 1
						
						for coord in readAlignmentCoords:
							mutantsInMultiplyOccupiedWells.append(int(coord.coord))
					
					l += 1
				k += 1
			j += 1
		i += 1
		
		uniqueTotalCoords = unique(collectedGenomicCoords)
		uniqueCoordsInMultiplyOccupiedWells = unique(mutantsInMultiplyOccupiedWells)
		uniqueCoordsInSinglyOccupiedWells = unique(mutantsInSinglyOccupiedWells)
		
		taxonomyDict['totalMutants'] = len(uniqueTotalCoords)
		taxonomyDict['mutantsInSinglyOccupiedWells'] = len(uniqueCoordsInSinglyOccupiedWells)
		taxonomyDict['mutantsInMultiplyOccupiedWells'] = len(uniqueCoordsInMultiplyOccupiedWells)


	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PrintSudokuGridOccupancyTaxonomy(taxonomyDict):
	
	keysInPrintOrder = [\
	'totalMutants',\
	'mutantsInSinglyOccupiedWells',\
	'mutantsInMultiplyOccupiedWells',\
	'emptyWells',\
	'singlyOccupiedWells',\
	'multiplyOccupiedWells']
	
	outputStr = ''
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'
	
	print(outputStr)
	return


# ------------------------------------------------------------------------------------------------ #

