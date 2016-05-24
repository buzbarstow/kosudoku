from kosudoku.grid import SudokuWell, SudokuGenomicCoord


# ------------------------------------------------------------------------------------------------ #
def ConvertSudokuColonyPurifiedWellBackToSudokuWell(well):
	
	from copy import deepcopy
	
	newWell = deepcopy(well)
	
	newWell.__class__ = SudokuWell
	
	if newWell.hasPredictionForContents and newWell.predictionCorrect == True :
	
		newWell.readAlignmentCoords = []
		newWell.simplifiedReadAlignmentCoords = []
		for coord in newWell.likelyReadAlignmentCoords:
			newWell.readAlignmentCoords.append(coord)
	
	elif newWell.hasPredictionForContents and newWell.predictionCorrect == False:
		newWell.readAlignmentCoords = []
		newCoord = SudokuGenomicCoord(-2, 'Incorrect', 0, 0)
		newCoord.featureName = 'Incorrect'
		newCoord.locusTag = ['Incorrect']
		newWell.readAlignmentCoords.append(newCoord)
	
	elif newWell.hasPredictionForContents == False:
		newWell.readAlignmentCoords = []
		newCoord = SudokuGenomicCoord(-1, 'NA', 0, 0)
		newCoord.featureName = 'Blank'
		newCoord.locusTag = ['Blank']
		newWell.readAlignmentCoords.append(newCoord)
		
	return newWell
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertSudokuGridLookupDictBackToSudokuWells(sudokuGridLookupDict, rowPools, colPools, \
prPools, pcPools):
	
	from copy import deepcopy

	newSudokuGrid = deepcopy(sudokuGridLookupDict)
	
	for pr in prPools:
		for pc in pcPools:
			plate = newSudokuGrid[pr][pc]
			for row in rowPools:
				for col in colPools:
					well = plate.wellGrid[row][col]
					plate.wellGrid[row][col] = ConvertSudokuColonyPurifiedWellBackToSudokuWell(well)
	
	return newSudokuGrid
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ExpandOutUniqueCoords(updatedColonyPurificationDict):
	from numpy import arange
	from scipy import unique
	
	cpKeys = updatedColonyPurificationDict.keys()
	coordsPutInExpanded = []

	for key in cpKeys:
		progenitorContents = updatedColonyPurificationDict[key]['progenitorContents']
		for content in progenitorContents:
			coord = content[0]
			coordExpanded = arange(coord - 1, coord + 2, 1, int)
			coordsPutInExpanded += list(coordExpanded)
	
	uniqueCoordsPutInExpanded = unique(coordsPutInExpanded)
	return uniqueCoordsPutInExpanded
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ReducePoolPresenceTable(poolPresenceTable, uniqueCoordsPutInExpanded):
	import numpy	
	reducedPoolPresenceTable = []
	linesToDelete = []
	i = 0
	while i < len(poolPresenceTable):
		coord = poolPresenceTable[i]['readAlignmentCoord']
		if coord not in uniqueCoordsPutInExpanded:
			linesToDelete.append(i)
		i += 1

	reducedPoolPresenceTable = numpy.delete(poolPresenceTable, linesToDelete)
	
	return reducedPoolPresenceTable
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class SudokuColonyPurifiedWell(SudokuWell):

	def __init__(self, sourcePlateName, sourcePlateRow, sourcePlateCol, sourceRow, \
	sourceCol, OD=1.0):
				
		SudokuWell.__init__(self, sourcePlateName, sourcePlateRow, sourcePlateCol, sourceRow, \
		sourceCol, OD=1.0)

		self.hasPredictionForContents = False
		# This variable is supposed to hold the predicted progenitor contents that have been 
		# expanded out to allow for slop in transposon position identification
		self.predictionsForContents = []
		self.predictionCorrect = False
		self.hopedForPresent = False
		self.hopedForCoord = -1
		self.simplifiedReadAlignmentCoords = []
		self.simplifiedLikelyReadAlignmentCoords = []
		self.likelyReadAlignmentCoords = []
		# This variable is supposed to hold the predicted progenitor contents
		self.progenitorContents = []
		self.progenitorLocatabilities = []
		self.condensationType = 'NA'
	
	def updateSimplifiedReadAlignmentCoords(self):
		i = 0
		self.simplifiedReadAlignmentCoords = []
		while i < len(self.readAlignmentCoords):
			simplifiedCoord = self.readAlignmentCoords[i].coord
			self.simplifiedReadAlignmentCoords.append(simplifiedCoord)
			i += 1
			
	def updateLikelyReadAlignmentCoords(self):
		self.likelyReadAlignmentCoords = []
		i = 0
		while i < len(self.simplifiedLikelyReadAlignmentCoords):
			coord = self.simplifiedLikelyReadAlignmentCoords[i]
			newCoord = SudokuGenomicCoord(coord, 'direct', 0, 0)
			self.likelyReadAlignmentCoords.append(newCoord)
			i += 1

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def UpdateSudokuGridLookupDictToSudokuColonyPurifiedWell(sudokuGridLookupDict, rowPools, colPools,\
prPools, pcPools):
	from copy import deepcopy

	for prPool in prPools:
		for pcPool in pcPools:
			for rowPool in rowPools:
				for colPool in colPools:
					well = sudokuGridLookupDict[prPool][pcPool].wellGrid[rowPool][colPool]
					well.__class__ = SudokuColonyPurifiedWell	
					
					well.hasPredictionForContents = False
					well.predictionsForContents = []
					well.predictionCorrect = False
					well.hopedForPresent = False
					well.simplifiedReadAlignmentCoords = []
					well.simplifiedLikelyReadAlignmentCoords = []
					well.likelyReadAlignmentCoords = []
					well.hopedForCoord = -1
					well.progenitorContents = []
					well.progenitorLocatabilities = []
					well.condensationType = 'NA'
					
					cpAddressData = deepcopy(well.addressDict['progenitor'])
					well.addressDict['colonyPurified'] = cpAddressData
					well.addressDict['progenitor'] = {}
					
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GroupGenomicCoords(coordArray, maxGap=1):
	
	import pdb
	import numpy
	import operator
	
	i = 0
	expect = None
	run = []
	result = [run]
	currentLineMatchesPrevious = True
	
	while i < len(coordArray):
		if currentLineMatchesPrevious:
			run.append(coordArray[i])
		else:
			run = [coordArray[i]]
			result.append(run)
			
		currentLineMatchesPrevious = False
		
		if i < len(coordArray) - 1:
			currentCoord = coordArray[i]
			nextCoord = coordArray[i+1]
			expect = currentCoord + maxGap
			if nextCoord <= expect:
				currentLineMatchesPrevious = True
		i += 1
		
	
	i = 0
	medianResults = []
	while i < len(result):
		medianResults.append(int(numpy.median(result[i])))
		i += 1
		
	return medianResults
	
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressesForPoolPresenceTableLine(poolPresenceTableLine, \
rowPools, colPools, prPools, pcPools, threshold):
	
	import pdb
	from kosudoku.pool import FindAddressCoords2
	
	line = poolPresenceTableLine
	
	[addresses_r, rowScoreDict] = FindAddressCoords2(line, rowPools, threshold=threshold)
	[addresses_c, colScoreDict] = FindAddressCoords2(line, colPools, threshold=threshold)
	[addresses_pr, prScoreDict] = FindAddressCoords2(line, prPools, threshold=threshold)
	[addresses_pc, pcScoreDict] = FindAddressCoords2(line, pcPools, threshold=threshold)		
	
	nAddresses_r = len(addresses_r)
	nAddresses_c = len(addresses_c)
	nAddresses_pr = len(addresses_pr)
	nAddresses_pc = len(addresses_pc)	
	
	nTotalCoords = nAddresses_r + nAddresses_c + nAddresses_pr + nAddresses_pc
	
# 	pdb.set_trace()
	
	possibleAddresses = []
	
	for address_pr in addresses_pr:
		for address_pc in addresses_pc:
			for address_c in addresses_c:
				for address_r in addresses_r:
					row = address_r
					col = address_c
					pr = address_pr
					pc = address_pc
					
					rowReads = rowScoreDict[row]
					colReads = colScoreDict[col]
					prReads = prScoreDict[pr]
					pcReads = pcScoreDict[pc]
					
					reads = rowReads + colReads + prReads + pcReads
					
					possibleAddress = row + '_' + col + '_' + pr + '_' + pc
					possibleAddresses.append([possibleAddress, reads])

	return possibleAddresses
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def PopulateColonyPurifiedSudokuGrid(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
	prPools, pcPools, readCountThreshold):

	import pdb
	from kosudoku.grid import SudokuGenomicCoord
	
	i = 0
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
				
		coord = poolPresenceTableLine['readAlignmentCoord']
		
		possibleAddresses = \
		CalculateLibraryAddressesForPoolPresenceTableLine(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, readCountThreshold)
		
		for address in possibleAddresses:
			[row, col, pr, pc] = address[0].split('_')
			readCount = address[1]
			newCoord = SudokuGenomicCoord(coord, 'direct', 0.0, readCount)
			sudokuGridLookupDict[pr][pc].wellGrid[row][col].readAlignmentCoords.append(newCoord)
		
		i += 1
	
	return 
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SolvePoolPresenceTableForColonyPurifiedCollection(sudokuGridLookupDict, poolPresenceTable, \
rowPools, colPools, prPools, pcPools, readCountThreshold):

	import pdb

	# pdb.set_trace()
	
	PopulateColonyPurifiedSudokuGrid(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
	prPools, pcPools, readCountThreshold)	
	
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignLikelyContentsToSudokuGridLookupDictAndCheckPredictions(sudokuGridLookupDict, \
updatedColonyPurificationDict):
	
	from numpy import arange
	from scipy import intersect1d, unique
	import pdb
	from copy import deepcopy
	
	cpKeys = list(updatedColonyPurificationDict.keys())

	totalHopedForCorrect = 0
	totalPredictedWells = 0
	totalPredictionsCorrect = 0
	hopedForTransposonsInCollection = []

	# Go through the colony purification dict
	for key in cpKeys:
		line = updatedColonyPurificationDict[key]
		rowPool = line['rowPool']
		colPool = str(int(line['colPool']))
		pcPool = line['pcPool']
		prPool = line['prPool']
		
		# Takes the progenitor contents and expands them a little to allow for slop in 
		# transposon position
		progenitorContentsExpanded = []
		progenitorContentsUnexpanded = []
		progenitorContentsUnexpandedLocatabilities = []
		i = 0

		while i < len(line['progenitorContents']):
			progCoord = line['progenitorContents'][i][0]
			progLocatability = line['progenitorContents'][i][1]
			progCoordExpanded = list(arange(progCoord - 1, progCoord + 2, 1, int))
			progenitorContentsExpanded += progCoordExpanded
			progenitorContentsUnexpanded.append(progCoord)
			progenitorContentsUnexpandedLocatabilities.append(progLocatability)
			i += 1
	
		hopedForCoord = int(line['hopedForTransposonCoord'])
		hopedForTransposonCoords = list(arange(hopedForCoord - 1, hopedForCoord + 2, 1, int))
	
		sudokuWell = sudokuGridLookupDict[prPool][pcPool].wellGrid[rowPool][colPool]
		sudokuWell.updateSimplifiedReadAlignmentCoords()
		sudokuWell.hasPredictionForContents = True
		sudokuWell.predictionsForContents = progenitorContentsExpanded
		sudokuWell.hopedForCoord = hopedForCoord
		sudokuWell.progenitorContents = progenitorContentsUnexpanded
		sudokuWell.progenitorLocatabilities = progenitorContentsUnexpandedLocatabilities
		sudokuWell.condensationType = line['condensationType']
		
		progenitorCol = line['progenitorCol']
		progenitorPlate = line['progenitorPlate']
		progenitorRow = line['progenitorRow']
		
		sudokuWell.addressDict['progenitor'] = {'plateName':progenitorPlate, \
		'row':progenitorRow, 'col':progenitorCol}
	
		# Intersect the predicted and hoped for transposon coords with the coords that intersect
		# with the wel location
	
		intersectPredictionsContents = intersect1d(sudokuWell.predictionsForContents, \
		sudokuWell.simplifiedReadAlignmentCoords)
	
		intersectHopedForContents = intersect1d(hopedForTransposonCoords, \
		sudokuWell.simplifiedReadAlignmentCoords)
		
		
		# Make some new read alignment coords
			
		# Check to see if the predicted and hoped for coords are there
	
		if len(intersectPredictionsContents) > 0:
			sudokuWell.predictionCorrect = True
			
			groupedIntersectPredictionsContents = \
			GroupGenomicCoords(intersectPredictionsContents, maxGap=1)
			
			sudokuWell.simplifiedLikelyReadAlignmentCoords = groupedIntersectPredictionsContents
						
			totalPredictionsCorrect += 1
		
		if len(intersectHopedForContents) > 0:
			sudokuWell.hopedForPresent = True
			sudokuWell.predictionCorrect = True
			
			groupedIntersectHopedForContents = \
			GroupGenomicCoords(intersectHopedForContents, maxGap=1)
			
			sudokuWell.simplifiedLikelyReadAlignmentCoords = groupedIntersectHopedForContents
			
			hopedForTransposonsInCollection.append(hopedForCoord)
			totalHopedForCorrect += 1
			
		totalPredictedWells += 1
		
		newReadAlignmentCoords = []
		
		if sudokuWell.predictionCorrect == True:
			for likelyCoord in sudokuWell.simplifiedLikelyReadAlignmentCoords:
				for readAlignmentCoord in sudokuWell.readAlignmentCoords:
					coord = readAlignmentCoord.coord
					if likelyCoord - 3 <= coord <= likelyCoord + 3:
						newReadAlignmentCoords.append(deepcopy(readAlignmentCoord))
		
		if sudokuWell.predictionCorrect and len(newReadAlignmentCoords) == 0:
			print('No coordinates found for well with correct prediction')
			pdb.set_trace()
		
		sudokuWell.readAlignmentCoords = newReadAlignmentCoords

	hopedForTransposonsInCollection = unique(hopedForTransposonsInCollection)
	print('Total Predicted Wells: ' + str(totalPredictedWells))
	print('Total Hoped for Correct: ' + str(totalHopedForCorrect))
	print('Total Predictions Correct: ' + str(totalPredictionsCorrect))
	print('Total Unique Hoped For Transposons: ' + str(len(hopedForTransposonsInCollection)))
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def UpdateLikelyReadAlignmentCoords(flattenedSudokuWellList):
	
	import pdb
	
	for well in flattenedSudokuWellList:
		well.updateLikelyReadAlignmentCoords()
	
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def AnnotateCondensedCollectionSudokuGridWithFeautureNames(flattenedSudokuWellList, featureArray):

	import pdb
	
	for well in flattenedSudokuWellList:

		readAlignmentCoords = well.likelyReadAlignmentCoords
	
		for readAlignmentCoord in readAlignmentCoords:
			coord = readAlignmentCoord.coord
		
			for feature in featureArray:
				startCoord = feature.startCoord
				endCoord = feature.endCoord
										
				if startCoord <= coord <= endCoord:
				
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
	
	return							
# ------------------------------------------------------------------------------------------------ #




# ----------------------------------------------------------------------------------------------- #
def WriteCondensedCollectionSummaryTable(rearrayedSudokuWellList, \
cpSudokuGridSummaryFileName, \
cpCollectionAddressDictKey='colonyPurified', \
progenitorCollectionAddressDictKey='progenitor', printBlankEntries=True):
# Modification of WriteSudokuGridSummary that allows output of feature names as well
	import pdb
	from kosudoku.xml import ExportListForXML
	
	fileHandle = open(cpSudokuGridSummaryFileName, 'w')
	
	outputStr = ''
	
	for well in rearrayedSudokuWellList:
		
# 		pdb.set_trace()
		
		coords = []
		featureNameStrings = []
		distsFromFeatureTranslationStartStrings = []
		fracDistsFromFeatureTranslationStartStrings = []
		locusTagsStrings = []
		
		cpAddressDict = well.addressDict[cpCollectionAddressDictKey]
		cpPlate = cpAddressDict['plateName']
		cpWell = cpAddressDict['row'] + cpAddressDict['col'].zfill(2)
		cpPlateRow = well.plateRow
		cpPlateCol = well.plateCol
		
		try:
			hopedForCoord = well.hopedForCoord
		except:
			pdb.set_trace()
		
		progAddressDict = well.addressDict[progenitorCollectionAddressDictKey]
		
		try:
			progPlate = progAddressDict['plateName']
		except:
			progPlate = ''
	
		try:
			progWell = progAddressDict['row'] + progAddressDict['col'].zfill(2)
		except:
			progWell = ''
		
		if well.predictionCorrect:
			predictionCorrectStr = 'Yes'
		else:
			predictionCorrectStr = 'No'
		
		if well.hasPredictionForContents and well.predictionCorrect:
		
			likelyCoords = well.likelyReadAlignmentCoords
		
			i = 0
			while i < len(likelyCoords):
				coord = likelyCoords[i]
				coords.append(coord.coord)
			
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
			
			progenitorContents = well.progenitorContents
			progenitorContentsStr = '"' + ExportListForXML(progenitorContents, delimeter=', ') + '"'
			
			progenitorLocatabilities = well.progenitorLocatabilities
			progenitorLocStr = '"' + ExportListForXML(progenitorLocatabilities, delimeter=', ') + '"'
			
			
			featureNameStrings = '"' + ExportListForXML(featureNameStrings, delimeter=', ') + '"'
		
			coordsStr = '"' + ExportListForXML(coords, delimeter=', ') + '"'

			distsStr = '"' + ExportListForXML(distsFromFeatureTranslationStartStrings, \
			delimeter=',')  + '"'

			fracDistsStr = '"' + ExportListForXML(\
			fracDistsFromFeatureTranslationStartStrings, delimeter=', ')  + '"'
		
			locusTagStr = '"' + ExportListForXML(locusTagsStrings, delimeter=', ')  + '"'
			
			hopedForCoordStr = str(hopedForCoord)
			
			if well.hopedForPresent == True:
				hopedForStr = 'Yes'
			else:
				hopedForStr = 'No'
			
			
			
		elif well.hasPredictionForContents and well.predictionCorrect == False:
			coordsStr = '-2'
			distsStr = 'NA'
			fracDistsStr = 'NA'
			locusTagStr = 'Prediction Incorrect'
			featureNameStrings = 'NA'
			hopedForStr = 'No'
			hopedForCoordStr = str(hopedForCoord)
			progenitorContents = well.progenitorContents
			progenitorContentsStr = '"' + ExportListForXML(progenitorContents, delimeter=', ') + '"'
			progenitorLocatabilities = well.progenitorLocatabilities
			progenitorLocStr = '"' + ExportListForXML(progenitorLocatabilities, delimeter=', ') + '"'
			
		
		elif well.hasPredictionForContents == False:
			coordsStr = '-1'
			distsStr = 'NA'
			fracDistsStr = 'NA'
			locusTagStr = 'Blank'
			featureNameStrings = 'NA'
			hopedForStr = 'NA'
			hopedForCoordStr = 'NA'
			progenitorContentsStr = 'NA'
			progenitorLocStr = 'NA'
			
		if printBlankEntries:
			try:
				outputStr += cpPlate + ',' + cpPlateRow + ',' + cpPlateCol + ',' + cpWell + ',' \
				+ progPlate + ',' + progWell + ',' \
				+ progenitorContentsStr + ',' + progenitorLocStr + ',' \
				+ hopedForCoordStr + ',' \
				+ coordsStr + ',' + featureNameStrings + ',' + distsStr + ',' + fracDistsStr + ',' \
				+ locusTagStr + ',' \
				+ predictionCorrectStr + ',' + hopedForStr + '\n'
			except:
				pdb.set_trace()
	
	headerStr = 'QC Plate,QC Plate Row,QC Plate Column,QC Well,' \
	+ 'Progenitor Plate,Progenitor Well,' \
	+ 'Predicted Progenitor Transposon Coordinates,Progenitor Locatabilities,' \
	+ 'Hoped For Transposon,' \
	+ 'Confirmed Predicted Transposon Coordinates,'\
	+ 'Features Disrupted,' \
	+ 'Distances of Transposon from Translation Starts,' \
	+ 'Fractional Distances of Transposon from Translation Starts,Locus Tags,' \
	+ 'Prediction Correct,Got What I Wanted\n'
	
	outputStr = headerStr + outputStr
	
	fileHandle.write(outputStr)
	fileHandle.close()
	return
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def UpdateCondensedCollectionCatalogWithProgenitorCollection(colonyPurificationDict, \
sudokuGridLookupDict, prPools, pcPools, rowPools, colPools):

	cpKeys = sorted(colonyPurificationDict.keys())
	for key in cpKeys:
		progenitorPlate = colonyPurificationDict[key]['progenitorPlate']
		progenitorRow = colonyPurificationDict[key]['progenitorRow']
		progenitorCol = colonyPurificationDict[key]['progenitorCol']
	
		for prPool in prPools:
			for pcPool in pcPools:
				plateName = sudokuGridLookupDict[prPool][pcPool].plateName
				wellGrid = sudokuGridLookupDict[prPool][pcPool].wellGrid
				if plateName == progenitorPlate:
					well = wellGrid[progenitorRow][progenitorCol]
					readAlignmentCoords = well.readAlignmentCoords
				
					for coord in readAlignmentCoords:
						coordInt = coord.coord
						coordLocatability = coord.locatability
						colonyPurificationDict[key]['progenitorContents'].append(\
						[coordInt, coordLocatability])

	return colonyPurificationDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def UpdateCondensedCollectionCatalogWithInternalAddressFormats(colonyPurificationDict, \
sudokuGridLookupDict, prPools, pcPools, prPoolKey='prPool', pcPoolKey='pcPool', \
rowPoolKey='rowPool', colPoolKey='colPool'):

	from kosudoku.utils import ConvertWellIDToRowAndColumn
	import pdb
	
	cpKeys = sorted(colonyPurificationDict.keys())
# 	pdb.set_trace()
	
	for key in cpKeys:
		matchFound = False
		[plate, well] = key.split('_')
		[row, col] = ConvertWellIDToRowAndColumn(well)
	
		j = 0
		while j < len(prPools) and matchFound == False:
			prPool = prPools[j]
			for pcPool in pcPools:
				plateName = sudokuGridLookupDict[prPool][pcPool].plateName
								
				if plateName == plate:
					colonyPurificationDict[key][prPoolKey] = prPool
					colonyPurificationDict[key][pcPoolKey] = pcPool
					colonyPurificationDict[key][rowPoolKey] = row
					colonyPurificationDict[key][colPoolKey] = col
					matchFound = True
					break
			j += 1
		if matchFound == False:
# 			pdb.set_trace()
			print('Match not for found for: ' + key)


	return colonyPurificationDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FillCoordinatePools(poolPresenceTable, rowPools, colPools, prPools, pcPools, threshold=1, \
thresholdForAddressCount=10):
	
	import pdb
	from scipy import intersect1d
	from kosudoku.pool import FindAddressCoords2
	
	rowPoolArrays = {}
	colPoolArrays = {}
	prPoolArrays = {}
	pcPoolArrays = {}
	
	for rowPool in rowPools:
		rowPoolArrays[rowPool] = []
		
	for colPool in colPools:
		colPoolArrays[colPool] = []
	
	for prPool in prPools:
		prPoolArrays[prPool] = []
		
	for pcPool in pcPools:
		pcPoolArrays[pcPool] = []
	
	i = 0 
	while i < len(poolPresenceTable):
		
		line = poolPresenceTable[i]
		coord = line['readAlignmentCoord']
		
		
		[addresses_r, rowScoreDict] = FindAddressCoords2(line, rowPools, threshold=thresholdForAddressCount)
		[addresses_c, colScoreDict] = FindAddressCoords2(line, colPools, threshold=thresholdForAddressCount)
		[addresses_pr, prScoreDict] = FindAddressCoords2(line, prPools, threshold=thresholdForAddressCount)
		[addresses_pc, pcScoreDict] = FindAddressCoords2(line, pcPools, threshold=thresholdForAddressCount)		
		
		nAddresses_r = len(addresses_r)
		nAddresses_c = len(addresses_c)
		nAddresses_pr = len(addresses_pr)
		nAddresses_pc = len(addresses_pc)
		
		coordNumberThreshold = 4
		
		addCoord = True
		
		if nAddresses_r > 8 or nAddresses_c > 5 or nAddresses_pr > 2 or nAddresses_pr > 2:
# 			print(str(coord) + " has lots of adresses coords")
			addCoord = False

# 		pdb.set_trace()
		if addCoord == True:
			for rowPool in rowPools:
				if line[rowPool] >= threshold:
					rowPoolArrays[rowPool].append(coord)
	# 				pdb.set_trace()
		
			for colPool in colPools:
				if line[colPool] >= threshold:
					colPoolArrays[colPool].append(coord)
	# 				pdb.set_trace()
				
			for prPool in prPools:
				if line[prPool] >= threshold:
					prPoolArrays[prPool].append(coord)
	# 				pdb.set_trace()
				
			for pcPool in pcPools:
				if line[pcPool] >= threshold:
					pcPoolArrays[pcPool].append(coord)

		
		
# 		pdb.set_trace()
		i += 1
	
	return [rowPoolArrays, colPoolArrays, prPoolArrays, pcPoolArrays]
# ------------------------------------------------------------------------------------------------ #


# ----------------------------------------------------------------------------------------------- #
def CalculateIntersectionReadCount(intersection, poolPresenceTable, pr, pc, row, col):
	
	i = 0
	intersectionFound = False
	while i < len(poolPresenceTable) and intersectionFound == False:
		coord = poolPresenceTable[i]['readAlignmentCoord']
		if intersection == coord:
			intersectionFound = True
			intersectionSize = poolPresenceTable[i][row]
			intersectionSize += poolPresenceTable[i][col]
			intersectionSize += poolPresenceTable[i][pr]
			intersectionSize += poolPresenceTable[i][pc]
		
		i += 1

	if intersectionFound == False:
		print('Intersection Not Found!')
	
	return intersectionSize
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def CalculateIntersectionReadCounts(intersections, poolPresenceTable, pr, pc, row, col):
	i = 0
	intersectionSizes = []
	while i < len(intersections):
		intersection = intersections[i]
		intersectionSize = \
		CalculateIntersectionReadCount(intersection, poolPresenceTable, pr, pc, row, col)
		intersectionSizes.append(intersectionSize)
		i += 1
	return intersectionSizes
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GetAddressCoordsForGenomicCoord(queryCoord, poolPresenceDict):
# Get the address coords for and around a requested genomic coordinate
	
	from numpy import arange
	import pdb
	
	i = 0
	
# 	pdb.set_trace()
	
	addressCoordsForQueryGenomicCoord = []
	
	queryRange = arange(queryCoord - 4,  queryCoord + 5, 1, int)
# 	pdb.set_trace()
	i = 0 
	while i < len(queryRange):
		if queryRange[i] in poolPresenceDict.keys():
			addressCoordsForQueryGenomicCoord += poolPresenceDict[queryRange[i]]
	
		i += 1
	
	return addressCoordsForQueryGenomicCoord
# ----------------------------------------------------------------------------------------------- #
		

# ----------------------------------------------------------------------------------------------- #
def UpdateCondensedCollectionCatalogWithPoolPresenceTableData(updatedColonyPurificationDict, \
poolPresenceDict):
# Gets the pool coordinate for each entry in the colony purification dict. Tells you the coordinates
# of the hoped for and possible contents. 
	

	
	cpKeys = updatedColonyPurificationDict.keys()
	
	for cpKey in cpKeys:
		hopedForCoord = updatedColonyPurificationDict[cpKey]['desiredTransposonCoord']
		
		progenitorContents = updatedColonyPurificationDict[cpKey]['progenitorContents']
		
		row = updatedColonyPurificationDict[cpKey]['rowPool']
		col = updatedColonyPurificationDict[cpKey]['colPool']
		pc = updatedColonyPurificationDict[cpKey]['pcPool']
		pr = updatedColonyPurificationDict[cpKey]['prPool']
		
		
		addressCoordsForHopedForCoord = GetAddressCoordsForGenomicCoord(hopedForCoord, \
		poolPresenceDict)
		
		i = 0
		addressCoordsForProgenitorContents = []
		while i < len(progenitorContents):
			progenitorCoord = progenitorContents[i][0]
			addressCoords = GetAddressCoordsForGenomicCoord(progenitorCoord, poolPresenceDict)
			addressCoordsForProgenitorContents.append(addressCoords)
			i += 1
		
		updatedColonyPurificationDict[cpKey]['addressCoordsForHopedForCoord'] = \
		addressCoordsForHopedForCoord
		
		updatedColonyPurificationDict[cpKey]['addressCoordsForProgenitorContents'] = \
		addressCoordsForProgenitorContents
		
		
	return updatedColonyPurificationDict
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GeneratePoolPresenceDict(poolPresenceTable, addressPools):
	from kosudoku.pool import FindAddressCoords
	i = 0
	poolPresenceDict = {}
	while i < len(poolPresenceTable):
		coord = poolPresenceTable[i]['readAlignmentCoord']
	
		poolPresenceDict[coord] = FindAddressCoords(poolPresenceTable[i], addressPools)
	
		i += 1
	return poolPresenceDict
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
def FindSingleRepresentativeForColonyPurifiedMutantsInCondensedCollection(cdSudokuWellArray):
	
	from copy import deepcopy
	from scipy import unique
	
	# Figure out which of the colony purified mutants has the correct transposon
	colonyPurifiedWellsWithCorrectTransposons = []
	isolatedHopedForTransposonCoords = []
	i = 0
	while i < len(cdSudokuWellArray):
		well = cdSudokuWellArray[i]
		
		if well.condensationType == 'Colony Purification' and well.hopedForPresent == True:
			colonyPurifiedWellsWithCorrectTransposons.append(deepcopy(well))
			isolatedHopedForTransposonCoords.append(well.hopedForCoord)
		
		i += 1
		
		
	# Now, we need to figure out a representative for all of these
	uniqTransposons = unique(isolatedHopedForTransposonCoords)
	
	cpWellDict = {}
	
	for transposon in uniqTransposons:
		cpWellDict[transposon] = []
		
	for well in colonyPurifiedWellsWithCorrectTransposons:
		coord = well.hopedForCoord
		cpWellDict[coord].append(well)
	
	
	uniqueRepresentativeWells = []
	
	coordKeys = cpWellDict.keys()
	
	for key in coordKeys:
		uniqueRepresentativeWells.append(cpWellDict[key][0])
	
	
	return uniqueRepresentativeWells
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignRepresentativeColonyPurifiedWellsToLaterInCondensedCollection(wellsToBeAssignedInput, \
startPlateNumber, \
startRowNumber, startColNumber, rows, columns, fillOrder='columns', fillPattern='checkerboard', \
plateRowForCondensedWells=1, fillLastPlateWithBlanks=True):

# Note on input parameters: rows and columns are ordered listings of all of the rows and columns
# on a plate. 

	from copy import deepcopy
	from .condense import TestIfCondensedCollectionPlateWellIsBlank, IncrementWell, \
	blankSudokuCoord, blankSudokuGridEntrySummary
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
	
		newWell = SudokuColonyPurifiedWell(currentPlateNumber, plateRowForCondensedWells, \
		currentPlateNumber, currentRow, currentCol)
		
		newWell.addressDict['new_condensed'] = \
		{'col': str(currentColNumber), 'row': str(currentRowNumber),\
		'plateName':str(currentPlateNumber), 'plateCol':'PC'+str(currentPlateNumber).zfill(2), \
		'plateRow': 'PR'+str(plateRowForCondensedWells).zfill(2)}
	
		# Test if this is supposed to be a blank well
		blank = TestIfCondensedCollectionPlateWellIsBlank(currentPlateNumber, currentRowNumber, \
		currentColNumber, fillPattern)
	
		if blank:
		# If the well is supposed to be blank, set readalignment coords to blank
			newWell.readAlignmentCoords.append(deepcopy(blankSudokuCoord))
			newWell.addressDict['progenitor'] = \
			blankSudokuGridEntrySummary['addressDict']['progenitor']		
		else:
		# If the well isn't supposed to be blank, copy over information from wells to be 
		# assigned
			wellToBeAssigned = wellsToBeAssigned.pop()
			addressDictData = deepcopy(wellToBeAssigned.addressDict['colonyPurified'])
			readAlignmentCoords = deepcopy(wellToBeAssigned.readAlignmentCoords)
			
			newWell.addressDict['progenitor'] = addressDictData
			newWell.readAlignmentCoords = readAlignmentCoords
			
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
	
			newWell = SudokuColonyPurifiedWell(currentPlateNumber, plateRowForCondensedWells, \
			currentPlateNumber, currentRow, currentCol)
			
			newWell.addressDict['new_condensed'] = \
			{'col': str(currentColNumber), 'row': str(currentRowNumber),\
			'plateName':str(currentPlateNumber), 'plateCol':'PC'+str(currentPlateNumber).zfill(2), \
			'plateRow': 'PR'+str(plateRowForCondensedWells).zfill(2)}
			
			newWell.addressDict['progenitor'] = \
			blankSudokuGridEntrySummary['addressDict']['progenitor']
			
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
def WriteColonyPurifiedMutantCondensationInstructions(condensationFileName, condensedWellArray):
	
	from kosudoku.xml import ExportListForXML
	import pdb
	
	headers = [\
	'#New Condensed Collection Plate', 'New Condensed Collection Row', \
	'New Condensed Collection Column', 'New Condensed Collection Well', \
	'New Condensed Collection Plate Row', 'New Condensed Collection Plate Column', \
	'Old Condensed Collection Plate', 'Old Condensed Collection Row', \
	'Old Condensed Collection Column', 'Old Condensed Collection Well', \
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
			condAddressDict = addressDict['new_condensed']
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
			condensationType = 'Colony Purification Rearray'
		except:
			pdb.set_trace()
		
		readAlignmentCoords = condensedWell.readAlignmentCoords
		desiredTransposonCoord = str(readAlignmentCoords[0].coord)
		
		desiredFeature = str(readAlignmentCoords[0].featureName[0])
		desiredLocusTag = str(readAlignmentCoords[0].locusTag[0])

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

