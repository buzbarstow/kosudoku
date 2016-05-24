
# ------------------------------------------------------------------------------------------------ #
class Feature3:
# An updated and simplified version of the feature2 class.
	def __init__(self, featureType, coordinates):
		self.coordinates = coordinates
		self.featureType = featureType
		self.tagDict = {}
		
		self.sudokuGridEntries = []
		
		self.featureName = ''
	
		[self.startCoord, self.endCoord, self.startTranslation, self.endTranslation] \
		= self.calculateStartAndEndCoords(coordinates)
		

	
	def updateTagDict(self, tag, tagData):
		
		tagDictKeys = list(self.tagDict.keys())
		
		if tag not in tagDictKeys:
			self.tagDict[tag] = []
			self.tagDict[tag].append(tagData)
		else:
			self.tagDict[tag].append(tagData)
		
		return
		
	def updateName(self):
		featureTagDictKeys = list(self.tagDict.keys())
	
		if 'gene' in featureTagDictKeys:
			geneName = self.tagDict['gene'][0]
		elif 'note' in featureTagDictKeys:
			geneName = self.tagDict['note'][0]
		elif 'locus_tag' in featureTagDictKeys:
			geneName = self.tagDict['locus_tag'][0]
		else:
			geneName = ''
	
		self.featureName = geneName
		
		return
	
	def calculateStartAndEndCoords(self, coordinates):
		import re
		coordsRe = re.compile(r'(\d+)\.\.(\d+)')
		complementRe = re.compile(r'complement\(')
		joinRe = re.compile(r'join\(')
		joinedCoordsRe = re.compile(r'(\d+)\.\.(\d+)\,(\d+)\.\.(\d+)')
		
		coordsMatch = coordsRe.search(coordinates)
		complementMatch = complementRe.search(coordinates)
		joinMatch = joinRe.search(coordinates)
		joinedCoordsMatch = joinedCoordsRe.search(coordinates)
		
		if (coordsMatch == None) and (joinedCoordsMatch == None):
			print('Error in matching coordinates!')
		else:
			if joinMatch != None and complementMatch == None:
				startCoord = int(joinedCoordsMatch.group(1))
				endCoord = int(joinedCoordsMatch.group(4))
				startTranslation = int(joinedCoordsMatch.group(1))
				endTranslation = int(joinedCoordsMatch.group(4))
			elif joinMatch != None and complementMatch != None:
				startCoord = int(joinedCoordsMatch.group(1))
				endCoord = int(joinedCoordsMatch.group(4))
				startTranslation = int(joinedCoordsMatch.group(4))
				endTranslation = int(joinedCoordsMatch.group(1))
			elif complementMatch == None and joinMatch == None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(1))
				endTranslation = int(coordsMatch.group(2))
			elif complementMatch != None and joinMatch == None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(2))
				endTranslation = int(coordsMatch.group(1))
		return [startCoord, endCoord, startTranslation, endTranslation]
	
# ------------------------------------------------------------------------------------------------ #





		
# ------------------------------------------------------------------------------------------------ #
class SudokuGridEntrySummaryLine:

	def __init__(self, plateName, row, col, od, rearrayedPlateName, rearrayedRow, \
	rearrayedCol, cellDensityScore, noCoordsNotInFeature, coordsInFeatureArray, \
	totalScore, totalProbabilityOfPicking, gridChoice, referringFeature):
		self.row = row
		self.col = col
		self.plateName = plateName
		self.rearrayedRow = rearrayedRow
		self.rearrayedCol = rearrayedCol
		self.rearrayedPlateName = rearrayedPlateName
		self.gridChoice = gridChoice
		self.cellDensityScore = cellDensityScore
		self.noCoordsNotInFeature = noCoordsNotInFeature
		self.noCoordsInFeature = len(coordsInFeatureArray)
		self.od = od
		self.coordsInFeatureArray = coordsInFeatureArray 
		self.totalScore = totalScore
		self.referringFeature = referringFeature
		self.totalProbabilityOfPicking = totalProbabilityOfPicking
		self.colonyPurificationWells = None
		
		return
	
	def calculateMinimumColoniesToPick(self, mode='inverseProbability',pickProbability=0.8):
		
		from numpy import log
		import pdb
		import math
		from math import ceil, floor
		
		
		
		if mode == 'inverseProbability':
			pA = self.totalProbabilityOfPicking
		elif mode == 'inverseFeatureCount':
			pA = (self.noCoordsInFeature/(self.noCoordsNotInFeature + self.noCoordsInFeature))
		
		if pA == 1.0:
			minimumColoniesToPick = 1.0
		else:
			minimumColoniesToPick = floor(log(1-pickProbability)/log(1-pA))
			
		
		return minimumColoniesToPick
	
	def getPlateName(self, fillNumber=3):
		
		if self.gridChoice == 'sudokuGrid':
			plateNameToReturn = str(self.plateName).zfill(fillNumber)
		elif self.gridChoice == 'rearrayedGrid':
			plateNameToReturn = str(self.rearrayedPlateName).zfill(fillNumber)
		
		return plateNameToReturn
	
	def getRow(self, fillNumber=3):
		
		if self.gridChoice == 'sudokuGrid':
			rowToReturn = str(self.row).zfill(fillNumber)
		elif self.gridChoice == 'rearrayedGrid':
			rowToReturn = str(self.rearrayedRow).zfill(fillNumber)
		
		return rowToReturn
	
	def getCol(self, fillNumber=3):
		
		if self.gridChoice == 'sudokuGrid':
			colToReturn = str(self.col).zfill(fillNumber)
		elif self.gridChoice == 'rearrayedGrid':
			colToReturn = str(self.rearrayedCol).zfill(fillNumber)
		
		return colToReturn
	
	def getGrid(self):
		
		if self.gridChoice == 'sudokuGrid':
			gridToReturn = 'sudoku'
		elif self.gridChoice == 'rearrayedGrid':
			gridToReturn = 'rearrayed'
		
		return gridToReturn

	
	def printSummaryLine(self, gridReportingChoice='sudokuGrid', gridReportingOn=False, \
	printHeaderLine=False, directOutput=False, minimumPickMode='inverseFeatureCount'):
		# Grid choice here is the grid that we want to report on
		
		if gridReportingChoice == 'sudokuGrid':
			row = self.row
			col = self.col
			plateName = self.plateName
		
		elif gridReportingChoice == 'rearrayedGrid':
			row = self.rearrayedRow
			col = self.rearrayedCol
			plateName = self.rearrayedPlateName
			
			
		minimumColoniesToPick = self.calculateMinimumColoniesToPick(mode=minimumPickMode)
		
		dataStr = str(plateName) + '_' + str(row) + str(col) + ': '
		
		dataStr += WriteSummaryOfGridEntryGenomicCoordsArray2(self.coordsInFeatureArray)
		
		dataStr += ', ' + str(self.noCoordsInFeature) + ', ' + str(self.noCoordsNotInFeature) \
		+ ', ' + str(self.od) + ', ' + str(self.totalScore) + ', ' \
		+ str(self.totalProbabilityOfPicking)  + ', ' + str(minimumColoniesToPick)
		
		if gridReportingOn == True:
			dataStr += ', ' + gridReportingChoice
		
		return dataStr
	
	
	def printSummaryLineForDirectCSVOutput(self, gridReportingChoice='sudokuGrid', \
	gridReportingOn=False, printHeaderLine=False, minimumPickMode='inverseFeatureCount', \
	includeCarriageReturn=True):
	
		import pdb
		
		# Grid choice here is the grid that we want to report on
		
		if gridReportingChoice == 'sudokuGrid':
			row = self.row
			col = self.col
			plateName = self.plateName
		
		elif gridReportingChoice == 'rearrayedGrid':
			row = self.rearrayedRow
			col = self.rearrayedCol
			plateName = self.rearrayedPlateName
		
		
		if self.gridChoice == 'sudokuGrid':
			plateToPickFrom = self.plateName
			rowToPickFrom = self.row
			colToPickFrom = self.col
		elif self.gridChoice == 'rearrayedGrid':
			plateToPickFrom = self.rearrayedPlateName
			rowToPickFrom = self.rearrayedRow
			colToPickFrom = self.rearrayedCol
		
		minimumColoniesToPick = self.calculateMinimumColoniesToPick(mode=minimumPickMode)
		
		headerStr = 'Plate,Row,Column,Rearrayed Plate,Rearrayed Row,Rearrayed Col,' \
		+ '"Summary of Genomic Coords (Coords, Distance To Translation Start, Locatability, ' \
		+ 'Pick Probability)",No Coords in Feature,No Coords not in Feature,OD,Total Score,' \
		+ 'Total Probability of Picking,Minimum Colonies To Pick,Colony Purification Well,' \
		+ 'Grid Choice,Plate To Pick From,Row to Pick From,Col to Pick From'
		
		dataStr = str(plateName) + ',' + str(row) + ',' + str(col) + ','
		dataStr += str(self.rearrayedPlateName) + ',' + str(self.rearrayedRow) + ',' \
		+ str(self.rearrayedCol) + ','
		
		dataStr += '"' + WriteSummaryOfGridEntryGenomicCoordsArray3(self.coordsInFeatureArray) + '"'
		
		dataStr += ',' + str(self.noCoordsInFeature) + ',' + str(self.noCoordsNotInFeature) \
		+ ',' + str(self.od) + ',' + str(self.totalScore) + ',' \
		+ str(self.totalProbabilityOfPicking) + ',' + str(minimumColoniesToPick) + ',' \
		+ '"' + str(self.colonyPurificationWells) + '"'
		
		dataStr += ','  + str(self.gridChoice) + ',' + str(plateToPickFrom) \
		+ ',' + str(rowToPickFrom) + ',' + str(colToPickFrom)
		
# 		pdb.set_trace()
		
		if includeCarriageReturn == True:
			EOLChar = '\n'
		else:
			EOLChar = ''
		
		if gridReportingOn == True:
			dataStr += ',' + gridReportingChoice + EOLChar
			headerStr += ',Grid Reporting Choice' + '\n'
		else:
			headerStr += EOLChar
			dataStr += EOLChar
		
		if printHeaderLine == True:
			outputStr = headerStr + dataStr
		else:
			outputStr = dataStr
		
		return outputStr
	
	
	
	def printPickingInstructions(self, minimumPickMode='inverseFeatureCount'):
	
		import pdb
				
		if self.gridChoice == 'sudokuGrid':
			plateToPickFrom = self.plateName
			rowToPickFrom = self.row
			colToPickFrom = self.col
			gridChoice = 'S'
		elif self.gridChoice == 'rearrayedGrid':
			plateToPickFrom = self.rearrayedPlateName
			rowToPickFrom = self.rearrayedRow
			colToPickFrom = self.rearrayedCol
			gridChoice = 'R'	
		
		minimumColoniesToPick = self.calculateMinimumColoniesToPick(mode=minimumPickMode)
		
		colonyPurificationWellsStr = ''
		
		if len(self.colonyPurificationWells) == 1:
			well = self.colonyPurificationWells[0]
			colonyPurificationWellsStr = str(well[1]) + '_' + str(well[0])
		
		elif len(self.colonyPurificationWells) == 2:
			well0 = self.colonyPurificationWells[0]
			well1 = self.colonyPurificationWells[1]
			colonyPurificationWellsStr = str(well0[1]) + '_' + str(well0[0])
			colonyPurificationWellsStr += ', ' + str(well1[1]) + '_' + str(well1[0])
		
		elif len(self.colonyPurificationWells) > 2:
			firstWell = self.colonyPurificationWells[0]
			lastWell = self.colonyPurificationWells[-1]			
			colonyPurificationWellsStr = str(firstWell[1]) + '_' + str(firstWell[0])
			colonyPurificationWellsStr += ' - ' + str(lastWell[1]) + '_' + str(lastWell[0])
		
		
		dataStr = '"' + str(gridChoice) + '_' + str(plateToPickFrom) + '_' \
		+ str(rowToPickFrom) + str(colToPickFrom) + '\n' + colonyPurificationWellsStr + '"'
				
		return dataStr

# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def FindWellsWithDisruptionBetweenStartAndEndCoords(sudokuWellList, startCoord, endCoord, \
startTranslation, endTranslation):
	
	i = 0
	
	wellsToReturn = []
	
	translationLength = abs(endTranslation - startTranslation)
	
	for well in sudokuWellList:
		
		row = well.row
		col = well.col
		pr = well.plateRow
		pc = well.plateCol
		addressDict = well.addressDict
		
		readAlignmentCoordsArray = well.readAlignmentCoords
		tempWellsToReturn = []
		wellOccupants = len(readAlignmentCoordsArray)
		
		for sudokuCoord in readAlignmentCoordsArray:
			
			coord = sudokuCoord.coord
			locatability = sudokuCoord.locatability
			fracDistFromTranslationStart = abs(coord - startTranslation)/translationLength
			distFromTranslationStart = abs(coord - startTranslation)
			
			if startCoord < coord < endCoord:
				tempWellsToReturn.append({'plateName':well.plateName, 'multipleOccupancy':False, \
				'readAlignmentCoord':coord, 'row':row, 'col':col, 'plateRow':pr, 'plateCol':pc, \
				'locatability':locatability, \
				'fracDistFromTranslationStart':fracDistFromTranslationStart, \
				'addressDict':addressDict, 'distFromTranslationStart':distFromTranslationStart, \
				'wellOccupants':wellOccupants})
		
		if len(readAlignmentCoordsArray) > 1:
			for tempWell in tempWellsToReturn:
				tempWell['multipleOccupancy'] = True
		
		for tempWell in tempWellsToReturn:
			wellsToReturn.append(tempWell)

	
	return wellsToReturn
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def UpdateFeatureArrayWithSudokuWells(featureArray, sudokuWellList):

	import pdb

	i = 0 
	while i < len(featureArray):
		feature = featureArray[i]		
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		startTranslation = feature.startTranslation
		endTranslation = feature.endTranslation
		
		wellsWithDisruptionsMatchingFeature = \
		FindWellsWithDisruptionBetweenStartAndEndCoords(sudokuWellList, startCoord, endCoord, \
		startTranslation, endTranslation)
		
		for well in wellsWithDisruptionsMatchingFeature:
			feature.sudokuGridEntries.append(well)
		
		i += 1
			
		
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateFeatureArrayTaxonomy(featureArray):
	
	import pdb
	
	keysInPrintOrder = [\
	'totalFeatures',\
	'featuresWithAtLeastOneDisruption',\
	'featuresWithNoDisruptions',\
	'featuresWithNoDisruptionsInSinglyOccupiedWells', \
	'featuresWithDisruptionsInSinglyOccupiedWells', \
	'featuresWithAtLeastOneUnambiguouslyMappingDisruption', \
	'featuresWithAtLeastOneDisruptionInFirstQuarter', \
	'featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstQuarter', \
	'featuresWithAtLeastOneSingleOccupanyDisruptionInFirstQuarter', \
	'featuresWithAtLeastOneDisruptionInFirstTenth', \
	'featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstTenth', \
	'featuresWithAtLeastOneSingleOccupanyDisruptionInFirstTenth']
	
	taxonomyDict = {}
	
	for key in keysInPrintOrder:
		taxonomyDict[key] = 0
	
	for feature in featureArray:
		taxonomyDict['totalFeatures'] += 1
		unambiguousEntryFound = False
		entryInFirstQuarterOfCDSFound = False
		unambiguousEntryInFirstQuarterOfCDSFound = False
		singleOccupancyEntryInFirstQuarterFound = False
		entryInFirstTenthOfCDSFound = False
		unambiguousEntryInFirstTenthOfCDSFound = False
		singleOccupancyEntryInFirstTenthFound = False
		
		if len(feature.sudokuGridEntries) == 0:
			taxonomyDict['featuresWithNoDisruptions'] += 1
		
		if len(feature.sudokuGridEntries) > 0:
			taxonomyDict['featuresWithAtLeastOneDisruption'] += 1
			
			entriesInSinglyOccupiedWells = 0
			
			for entry in feature.sudokuGridEntries:
				if entry['multipleOccupancy'] == False:
					entriesInSinglyOccupiedWells += 1
					
				if entry['fracDistFromTranslationStart'] <= 0.25:
					entryInFirstQuarterOfCDSFound = True
				
				if entry['fracDistFromTranslationStart'] <= 0.1:
					entryInFirstTenthOfCDSFound = True
				
				if entry['locatability'] == 'unambiguous':
					unambiguousEntryFound = True
				
				if entry['locatability'] == 'unambiguous' and \
				entry['fracDistFromTranslationStart'] <= 0.25:
					unambiguousEntryInFirstQuarterOfCDSFound = True
				
				if entry['locatability'] == 'unambiguous' and \
				entry['fracDistFromTranslationStart'] <= 0.1:
					unambiguousEntryInFirstTenthOfCDSFound = True
					
				if entry['multipleOccupancy'] == False and \
				entry['fracDistFromTranslationStart'] <= 0.25:
					singleOccupancyEntryInFirstQuarterFound = True
				
				if entry['multipleOccupancy'] == False and \
				entry['fracDistFromTranslationStart'] <= 0.1:
					singleOccupancyEntryInFirstTenthFound = True
				
				
			if entriesInSinglyOccupiedWells == 0:
				taxonomyDict['featuresWithNoDisruptionsInSinglyOccupiedWells'] += 1
			elif entriesInSinglyOccupiedWells > 0:
				taxonomyDict['featuresWithDisruptionsInSinglyOccupiedWells'] += 1
			
			if unambiguousEntryFound:
				taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruption'] += 1
			
			if unambiguousEntryInFirstQuarterOfCDSFound:
				taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstQuarter']\
				+= 1
			
			if entryInFirstQuarterOfCDSFound:
				taxonomyDict['featuresWithAtLeastOneDisruptionInFirstQuarter'] += 1
			
			if singleOccupancyEntryInFirstQuarterFound:
				taxonomyDict['featuresWithAtLeastOneSingleOccupanyDisruptionInFirstQuarter'] += 1
			
			if unambiguousEntryInFirstTenthOfCDSFound:
				taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstTenth']\
				+= 1
			
			if entryInFirstTenthOfCDSFound:
				taxonomyDict['featuresWithAtLeastOneDisruptionInFirstTenth'] += 1
			
			if singleOccupancyEntryInFirstTenthFound:
				taxonomyDict['featuresWithAtLeastOneSingleOccupanyDisruptionInFirstTenth'] += 1
			
			
	
	
	PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder)
	
	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder):
	
	outputStr = ''
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'
	
	print(outputStr)
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ImportFeatureArrayFromGenBank(genBankFile):
# A heavily updated version of FindFeaturesInGenome. Can deal with all sorts of tags, rather than
# just genes and can deal with joined CDSs. 

	import re
	import pdb
	
	featureNameRe = re.compile(r'\s+/gene=[\'|"](\w+)[\'|"]')
	
	featureStartRe = re.compile(r'\s+(\w+)\s+([complement\(|join\(]*\d+\.\.\d+[\,\d+\.\.\d+]*[\)]*)')
	
	featureTagRe = \
	re.compile(r'\s+/([a-zA-Z0-9_\-]+)=([\'|"]*)([\[\]a-zA-Z0-9_\-\.\,\:?!><\- \t\&\*\(\)/;\+\/`\']+)([\'|"]*)')
# 	re.compile(r'\s+/([a-zA-Z0-9_\-]+)=([\'|"]*)([\[a\]-zA-Z0-9_\-\.\,\:\?\! \t\&\*\(\)/;\+\/`\']+)([\'|"]*)')

	
	featureTagRe2 = \
	re.compile(r'\s+/([a-zA-Z0-9_\-]+)')
	
	sequenceSectionStartRe = re.compile(r'\s*ORIGIN')
	featureSectionStartRe = re.compile(r'\s*FEATURES')

	
	fileHandle = open(genBankFile, 'r')
	
	data = fileHandle.readlines()
	fileHandle.close()
	
	
	featureArray = []
	
	i = 0
	currentCoordinates = None
	currentFeature = None
	
	inFeatureEntry = False
	
	inFeatureSection = False
	inSequenceSection = False
	
	continuationNeeded = False
	
	
	while i < len(data):
			
		line = data[i]
		
		# Start looking for the start of a gene record
		
		featureSectionStartMatch = featureSectionStartRe.match(line)
		sequenceSectionStartMatch = sequenceSectionStartRe.match(line)
		
		if featureSectionStartMatch != None:
			inFeatureSection = True
			inSequenceSection = False
			inFeatureEntry = False
			
		if sequenceSectionStartMatch != None:
			inFeatureSection = False
			inSequenceSection = True
			inFeatureEntry = False
			
			if currentFeature != None:
				featureArray.append(currentFeature)
			
		
		if inFeatureSection and not inSequenceSection:
			
			featureStartMatch = featureStartRe.search(line)
					
			if (featureStartMatch != None):
				if currentFeature != None:
					featureArray.append(currentFeature)
				
				# Now, start getting the new feature
				inInterestingRecord = True
				currentFeatureType = featureStartMatch.group(1)
				currentCoordinates = featureStartMatch.group(2)
			
				# Initialize a new current feature
				currentFeature = Feature3(currentFeatureType, currentCoordinates)
				inFeatureEntry = True
				
			else:
				if inFeatureEntry == True:
					
					featureTagMatch = featureTagRe.search(line)
					featureTagMatch2 = featureTagRe2.match(line)
					
					if (featureTagMatch != None):
						
						if continuationNeeded == True:
							print("continuation error step 1!")
							pdb.set_trace()
							return
						
						
						currentTag = featureTagMatch.group(1)
						currentTagData = featureTagMatch.group(3).strip()
						
						currentStartQuote = featureTagMatch.group(2).strip()
						currentEndQuote = featureTagMatch.group(4).strip()
						
						# If the tag data starts with a quote but doesn't end with a quote, 
						# then continue onto the next line
						if ((currentStartQuote == '"') or (currentStartQuote == '\'')) \
						and (currentEndQuote == ''):
							continuationNeeded = True
					
					if (featureTagMatch2 != None) and (featureTagMatch == None) \
					and (continuationNeeded == False):
						# print(featureTagMatch2.group(1))
# 						print(currentCoordinates)
						currentTag = featureTagMatch2.group(1)
						currentTagData = True
					
							
					if continuationNeeded == True and featureTagMatch == None:
						currentTagData += line.strip()
						currentTagEnding = currentTagData[-1]
						
						if ((currentTagEnding == '"') or (currentTagEnding == '\'')):
							currentTagData = currentTagData[0:-1]
							continuationNeeded = False
						
					if continuationNeeded == False:
# 						print(currentTag + str(currentTagData))
						currentFeature.updateTagDict(currentTag, currentTagData)
				
		i += 1
	
	for feature in featureArray:
		feature.updateName()
		
	
	return featureArray
# ------------------------------------------------------------------------------------------------ #

