# ------------------------------------------------------------------------------------------------ #
def CalculateFeatureTaxonomyOfRearrayedCollection(featureArray):
	
	import pdb
	from sudokuutils13.feature import PrintFeatureArrayTaxonomyDict
	
	keysInPrintOrder = [\
	'totalFeatures',\
	'featuresWithOneDisruption',\
	'featuresWithTwoDisruptions',\
	'featuresWithThreeDisruptions', \
	'featuresWithFourOrMoreDisruptions', \
	'featuresWithAtLeastOneDisruption', \
	'featuresWithAtLeastOneDisruptionInSinglyOccupiedWell', \
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
		
		hasDisruptionInSinglyOccupiedWell = False
		unambiguousEntryFound = False
		entryInFirstQuarterOfCDSFound = False
		unambiguousEntryInFirstQuarterOfCDSFound = False
		singleOccupancyEntryInFirstQuarterFound = False
		entryInFirstTenthOfCDSFound = False
		unambiguousEntryInFirstTenthOfCDSFound = False
		singleOccupancyEntryInFirstTenthFound = False
		
		taxonomyDict['totalFeatures'] += 1
		
		if len(feature.sudokuGridEntries) > 0:
			taxonomyDict['featuresWithAtLeastOneDisruption'] += 1
			
		if len(feature.sudokuGridEntries) == 1:
			taxonomyDict['featuresWithOneDisruption'] += 1	
		elif len(feature.sudokuGridEntries) == 2:
			taxonomyDict['featuresWithTwoDisruptions'] += 1
		elif len(feature.sudokuGridEntries) == 3:
			taxonomyDict['featuresWithThreeDisruptions'] += 1
		elif len(feature.sudokuGridEntries) >= 4:
			taxonomyDict['featuresWithFourOrMoreDisruptions'] += 1			
		
		for entry in feature.sudokuGridEntries:
			if entry['multipleOccupancy'] == False:
				hasDisruptionInSinglyOccupiedWell = True
								
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
			entry['fracDistFromTranslationStart'] <= 0.1:
				singleOccupancyEntryInFirstTenthFound = True
		
			if entry['multipleOccupancy'] == False and \
			entry['fracDistFromTranslationStart'] <= 0.25:
				singleOccupancyEntryInFirstQuarterFound = True
			
		
		if hasDisruptionInSinglyOccupiedWell:
			taxonomyDict['featuresWithAtLeastOneDisruptionInSinglyOccupiedWell'] += 1
		
		if unambiguousEntryFound:
			taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruption'] += 1
	
		if entryInFirstQuarterOfCDSFound:
			taxonomyDict['featuresWithAtLeastOneDisruptionInFirstQuarter'] += 1
		
		if unambiguousEntryInFirstQuarterOfCDSFound:
			taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstQuarter']\
			+= 1
		
		if singleOccupancyEntryInFirstQuarterFound:
			taxonomyDict['featuresWithAtLeastOneSingleOccupanyDisruptionInFirstQuarter'] += 1
		
		if entryInFirstTenthOfCDSFound:
			taxonomyDict['featuresWithAtLeastOneDisruptionInFirstTenth'] += 1
			
		if unambiguousEntryInFirstTenthOfCDSFound:
			taxonomyDict['featuresWithAtLeastOneUnambiguouslyMappingDisruptionInFirstTenth']\
			+= 1
		
		if singleOccupancyEntryInFirstTenthFound:
			taxonomyDict['featuresWithAtLeastOneSingleOccupanyDisruptionInFirstTenth'] += 1
		
	
	PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder)
	
	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ImportSimplifiedCatalogOfIsaosLibrary(filename, sourcePlateCol=0, sourceWellCol=1, \
rearrayedPlateCol=2, rearrayedWellCol=3):
	
	import re
	import pdb
	from numpy import zeros
	
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	noGrowthReportedRe = re.compile(r'No Growth', re.IGNORECASE)
	disruptionCoordRe = re.compile('\d+')
	
	fileHandle = open(filename, 'r')
	data = fileHandle.readlines()
	
	
	dtypeArray = [('sourcePlate', '|S10'), ('sourceRow', '|S10'), ('sourceCol', 'i8'), \
	('rearrayedPlate', '|S10'), ('rearrayedRow', '|S10'), ('rearrayedCol', 'i8')]
	
	
	rearrayedWellList = []
	nonCommentLines = 0
	i = 0
	while i < len(data):
		if data[i][0] != '#':
			dataLine = data[i].strip().split(',')
			
			sourcePlate = str(dataLine[sourcePlateCol])
			sourceWell = dataLine[sourceWellCol]
			rearrayedPlate = dataLine[rearrayedPlateCol]		
			rearrayedWell = dataLine[rearrayedWellCol]
		
			sourceWellMatch = wellIDRe.search(sourceWell)
		
			if sourceWellMatch != None:
				sourceRow = sourceWellMatch.group(1)
				sourceCol = sourceWellMatch.group(2)
			else:
				pdb.set_trace()

			rearrayedWellMatch = wellIDRe.search(rearrayedWell)
			rearrayedRow = rearrayedWellMatch.group(1)
			rearrayedCol = rearrayedWellMatch.group(2)
		
			
			rearrayedWellList.append([sourcePlate, sourceRow, sourceCol, rearrayedPlate, \
			rearrayedRow, rearrayedCol])
			
			nonCommentLines += 1
			
		i += 1
	
	rearrayedWellArray = zeros(nonCommentLines, dtype=dtypeArray)
	i = 0
	while i < len(rearrayedWellArray):
		rearrayedWellArray[i]['sourcePlate'] = rearrayedWellList[i][0]
		rearrayedWellArray[i]['sourceRow'] = rearrayedWellList[i][1]
		rearrayedWellArray[i]['sourceCol'] = rearrayedWellList[i][2]
		rearrayedWellArray[i]['rearrayedPlate'] = rearrayedWellList[i][3]
		rearrayedWellArray[i]['rearrayedRow'] = rearrayedWellList[i][4]
		rearrayedWellArray[i]['rearrayedCol'] = rearrayedWellList[i][5]
		
		i += 1
		
	return rearrayedWellArray
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def UpdateSudokuGridWithRearrayedData(rearrayedWellArray, sudokuGridLookupDict, prPools, pcPools, \
rowPools, colPools):
	
	import pdb
	import re
	from copy import deepcopy
	from sudokuutils13.grid import FindSudokuWellByPlateName
	
		
	rearrayedSudokuWellList = []
	
	i = 0
	while i < len(rearrayedWellArray):
		
		sourceRow = rearrayedWellArray[i]['sourceRow'].decode('utf-8')
		sourceCol = rearrayedWellArray[i]['sourceCol']
		sourcePlate = rearrayedWellArray[i]['sourcePlate'].decode('utf-8')
		rearrayedRow = rearrayedWellArray[i]['rearrayedRow'].decode('utf-8')
		rearrayedCol = str(rearrayedWellArray[i]['rearrayedCol'])
		rearrayedPlate = rearrayedWellArray[i]['rearrayedPlate'].decode('utf-8')
		
		# Find the corresponding well in the sudokuGridLookupDict
			
		rearrayedSudokuWell = FindSudokuWellByPlateName(sudokuGridLookupDict, prPools, pcPools, \
		rowPools, colPools, sourcePlate, sourceRow, sourceCol)
		
		rearrayedSudokuWell.addressDict['isaosCollection'] = \
		{'plateName':rearrayedPlate,'row':rearrayedRow, 'col':rearrayedCol}
		 
		rearrayedSudokuWellList.append(rearrayedSudokuWell)
		
		i += 1
	

	return rearrayedSudokuWellList
# ------------------------------------------------------------------------------------------------ #



# ----------------------------------------------------------------------------------------------- #
def WriteRearrayedCollectionSummaryTable(sudokuGridSummaryFileName, rearrayedSudokuWellList, \
rearrayedCollectionAddressDictKey='isaosCollection', \
progenitorCollectionAddressDictKey='progenitor'):
# Modification of WriteSudokuGridSummary that allows output of feature names as well
	import pdb
	from sudokuutils13.xml import ExportListForXML
	
	fileHandle = open(sudokuGridSummaryFileName, 'w')
	outputStr = ''
	outputStr += 'Plate,Plate Row,Plate Column,Well,' \
	+ 'Rearrayed Plate,Rearrayed Well,' \
	+ 'Transposon Coordinates,' \
	+ 'Transposon Locatabilities,' \
	+ 'Features Disrupted,Distances of Transposon from Translation Starts,' \
	+ 'Fractional Distances of Transposon from Translation Starts,Locus Tags\n'
	
	for well in rearrayedSudokuWellList:
	
		readAlignmentCoords = well.readAlignmentCoords
		
		coords = []
		locatabilities = []
		featureNameStrings = []
		distsFromFeatureTranslationStartStrings = []
		fracDistsFromFeatureTranslationStartStrings = []
		locusTagsStrings = []
		
		rearrayAddressDict = well.addressDict[rearrayedCollectionAddressDictKey]
		
		try:
			rearrayPlate = rearrayAddressDict['plateName']
		except:
			pdb.set_trace()
		
		rearrayWell = rearrayAddressDict['row'] + rearrayAddressDict['col']
		
		plateRow = well.plateRow
		plateCol = well.plateCol
		plateName = well.plateName
		row = well.row
		col = well.col
		wellID = row+col
		
		
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
		
		
		outputStr += plateName + ',' + plateRow + ',' + plateCol +  ',' + wellID + ',' \
		+ rearrayPlate + ',' + rearrayWell + ',' + coordsStr + ',' \
		+ locatabilitiesStr + ',' + featureNameStrings + ',' + distsStr + ',' \
		+ fracDistsStr + ',' + locusTagStr + '\n'
	
	fileHandle.write(outputStr)
	fileHandle.close()
	return
# ------------------------------------------------------------------------------------------------ #


