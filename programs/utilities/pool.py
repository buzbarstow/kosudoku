
####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Functions for handling pool presence table analysis
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

from kosudoku.utils import gSeparatorString, UpdateLogFileData
from kosudoku.output import generateOutputMatrixWithHeaders, writeOutputMatrix
from kosudoku.grid import SudokuGenomicCoord, CalculateSudokuGridOccupancyTaxonomy, \
PrintSudokuGridOccupancyTaxonomy



# ------------------------------------------------------------------------------------------------ #
def GenerateDTypeArrayForPoolPresenceTableImport(poolColumns):
	
	import pdb
	from numpy import int32
	
	dtypeArray = []
	
	i = 0
	while i < len(poolColumns):
		dtypeArray.append((poolColumns[i], int32))
		i += 1
	
	return dtypeArray
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateDTypeArrayForPoolPresenceDict(indexLookupTable):
	
	import pdb
	from numpy import int32
	
	poolNames = indexLookupTable.values()
	poolNames = sorted(poolNames)
		
	dtypeArray = []
	dtypeArray.append(('readAlignmentCoord', int32))
	
	# Add in columns for pools
	i = 0
	while i < len(poolNames):
		dtypeArray.append((poolNames[i], int32))
		i += 1
	
	
	return dtypeArray
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GeneratePoolPresenceTable(uniqueCoords, sortedValidGenomeArray, indexLookupTable, dtypeDict):
	import numpy
	import pdb
	
	try:
		poolPresenceTable = numpy.zeros(len(uniqueCoords), dtype=dtypeDict)
	except:
		pdb.set_trace()
	
	i = 0
	while i < len(uniqueCoords):
		poolPresenceTable[i] = uniqueCoords[i]
		i += 1
	
	
	i = 0
	j = 0
	
	while j < len(uniqueCoords):
		while(i < len(sortedValidGenomeArray) and \
		sortedValidGenomeArray[i]['readAlignmentCoord'] == poolPresenceTable[j]['readAlignmentCoord']):
		
			index = str(sortedValidGenomeArray[i]['index'])
			column = indexLookupTable[index]
			poolPresenceTable[j][column] += 1
		
			i += 1
		
		j += 1
		
	return poolPresenceTable
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportPoolPresenceTable(poolPresenceTableFileName, poolColumns):
	
	import numpy
	import gc
	from pdb import set_trace
	import csv
	
	
	fHandle = open(poolPresenceTableFileName, 'r')
	
	poolColumnToHeaderIndexDict = {}
	i = 0
	datareader = csv.reader(fHandle)
	for row in datareader:
		if i == 0:
			headers = row
			poolColumnToHeaderIndexDict = {}
			for col in poolColumns:
				poolColumnToHeaderIndexDict[col] = headers.index(col)
		i += 1
	fHandle.close()
	dtypeArray = GenerateDTypeArrayForPoolPresenceTableImport(poolColumns)
	poolPresenceTable = numpy.zeros(i-1, dtype=dtypeArray)
	
	
	i = 0
	colKeys = list(poolColumnToHeaderIndexDict.keys())
	colIndices = []
	while i < len(colKeys):
		colIndices.append(poolColumnToHeaderIndexDict[colKeys[i]])
		i += 1
	
	
	
	fHandle = open(poolPresenceTableFileName, 'r')
	i = 0
	datareader = csv.reader(fHandle)
	for row in datareader:
		if i > 0:
			j = 0
			while j < len(colKeys):
				poolPresenceTable[i-1][colKeys[j]] = row[colIndices[j]]
				j += 1
		i += 1
	fHandle.close()
	
	return poolPresenceTable

# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def WritePoolPresenceTable3(filename, poolPresenceTable, poolColumns):

	fhandle = open(filename, 'w')
	
			
	# Write out the header line
	
	headerLine = ''
	
	i = 0
	while i < len(poolColumns):
		headerLine += poolColumns[i]
		if i < len(poolColumns) - 1:
			headerLine += ','
		i += 1
	
	headerLine += '\n'
	
	totalStr = ''
	
	i = 0
	while i < poolPresenceTable.shape[0]:
		outputStr = ''
		j = 0
		while j < len(poolColumns):
			outputStr += str(poolPresenceTable[i][poolColumns[j]]) 
			if j < len(poolColumns) - 1:
				outputStr += ','
			j += 1
		
		outputStr += '\n'
		totalStr += outputStr
		i += 1
		
	writeStr = headerLine + totalStr
	
	fhandle.write(writeStr)
	fhandle.close()

	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def GenerateBlankPhysicalAddressDict(poolPresenceTable):

	uniqueCoords = poolPresenceTable['readAlignmentCoord']

	physicalAddressDict = {}

	i = 0

	while i < len(uniqueCoords):
		physicalAddressDict[int(uniqueCoords[i])] = []
		i += 1

	return physicalAddressDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindAddressCoords(poolPresenceTableLine, axisPools, threshold=5):
	
	axisAddresses = []
	
	i = 0
	while i < len(axisPools):
		if poolPresenceTableLine[axisPools[i]] >= threshold:
			axisAddresses.append([axisPools[i], poolPresenceTableLine[axisPools[i]]])
		i += 1
	

	return axisAddresses
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def FindAddressCoords2(poolPresenceTableLine, axisPools, threshold=5):
	# Very much the same as FindAddressCoords, but returns a list of just axis addresses and a 
	# dict of the read counts for each of these axis addresses, rather than combining this dict
	# into the axisAddresses array
	
	axisAddresses = []
	axisAddressScoreDict = {}
	
	i = 0
	while i < len(axisPools):
		if poolPresenceTableLine[axisPools[i]] >= threshold:
			axisAddresses.append(axisPools[i])
			axisAddressScoreDict[axisPools[i]] = poolPresenceTableLine[axisPools[i]]
		i += 1
	

	return [axisAddresses, axisAddressScoreDict]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePoolCoordsForLine(poolPresenceTableLine, rowPools, colPools, prPools, \
pcPools, controlPools, threshold=5):
# Very much like CalculatePhysicalAddressCoordsForLine but also reports contents of control pools
# as well.  
# Used in generation of pool presence table taxonomy

	addresses_r = FindAddressCoords2(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords2(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords2(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords2(poolPresenceTableLine, pcPools, threshold=threshold)
	
	if controlPools != None:
		addresses_control = FindAddressCoords2(poolPresenceTableLine, controlPools, \
		threshold=threshold)
	else:
		addresses_control = None
	
	# Remember, each line in the addresses array is a 2 element list, the first containing the pool
	# name, and the second containing the number of reads associated with it. 
	
	return [addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, \
addresses_pr, addresses_pc, addresses_control):
	
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	nAddressPoolAxesWithMoreThanOneEntry = 0
	if nRowPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nColPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nPRPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nPCPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	

	return nAddressPoolAxesWithMoreThanOneEntry
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
addresses_control):
# Used in generation of pool presence table taxonomy
# Used to calculate how many lines can be used to calculate library addresses
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	nAddressPoolAxesWithEntries = 0
	if nRowPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nColPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nPRPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nPCPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
		
	
	if addresses_control != None:
		nControlPoolsWithEntries = len(addresses_control[0])
	else:
		nControlPoolsWithEntries = 0

	return [nAddressPoolAxesWithEntries, nControlPoolsWithEntries]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, addresses_pc):
# Calculate the unambiguous library addresses that that a read alignment coord maps to
# The score is the total number of reads that are associated with the library address assignment
	possibleAddresses = []
	
	for address_r in addresses_r:
		for address_c in addresses_c:
			for address_pr in addresses_pr:
				for address_pc in addresses_pc:
					
					row = address_r[0]
					row_score = address_r[1]
					col = address_c[0]
					col_score = address_c[1]
					pr = address_pr[0]
					pr_score = address_pr[1]
					pc = address_pc[0]
					pc_score = address_pc[1]
					
					possibleAddress = row + '_' + col + '_' + pr + '_' + pc
					possibleAddressScore = row_score + col_score + pr_score + pc_score
					
					possibleAddresses.append([possibleAddress, possibleAddressScore])

	return possibleAddresses
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressesForPoolPresenceTableLine(poolPresenceTableLine, rowPools, colPools, \
prPools, pcPools, threshold=5):
# Used in generation of pool presence table taxonomy
# Used to calculate possible addresses for pool presence table line
	
	coord = int(poolPresenceTableLine['readAlignmentCoord'])
	
	addresses_r = FindAddressCoords(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords(poolPresenceTableLine, pcPools, threshold=threshold)
	
	possibleAddresses = CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, \
	addresses_pc)
	
	return possibleAddresses
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, rowPools, colPools, \
prPools, pcPools, logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict, \
threshold):

# Used in generation of pool presence table taxonomy
# Used to calculate possible addresses for pool presence table line
	
	coord = int(poolPresenceTableLine['readAlignmentCoord'])
	
	addresses_r = FindAddressCoords(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords(poolPresenceTableLine, pcPools, threshold=threshold)
	
	possibleAddressesAndScores = \
	CalculatePossibleAddresses2(addresses_r, addresses_c, addresses_pr, addresses_pc, \
	logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict)
	
	return possibleAddressesAndScores
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
addresses_pr, addresses_pc):
	
	import pdb
	
	nRowCoords = len(addresses_r[0])
	nColCoords = len(addresses_c[0])
	nPRCoords = len(addresses_pr[0])
	nPCCoords = len(addresses_pc[0])
	
# 	pdb.set_trace()
	
	return [nRowCoords, nColCoords, nPRCoords, nPCCoords]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateMaxNumberOfEntriesInSinglePoolAxis(addresses_r, addresses_c, \
addresses_pr, addresses_pc):
# This function finds the pool axis with the most coordinate entries and reports this number
# This number is important for guessing how many library addresses an ambiguous line might map to
# It is the minimum number of addresses that the line could map to. 
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	poolAxisEntries = [nRowPoolAxisEntries, nColPoolAxisEntries, nPRPoolAxisEntries, \
	nPCPoolAxisEntries]
	
	maxEntriesInSingleAddressPoolAxis = max(poolAxisEntries)
	

	return maxEntriesInSingleAddressPoolAxis
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Functions for calculating a Voigt function to a read ratio histogram
def voigtFunction(x, p):
	
	from scipy.special import erfc
	from numpy import exp
	from numpy import sqrt
	from numpy import pi, float64

	a = p[0]
	c = p[1]
	delta = p[2]
	sigma = p[3]
	
	firstArg = ((-1.j)*(-c + x) + delta)/(sqrt(2)*sigma)
	secondArg = ((1.j)*(-c + x) + delta)/(sqrt(2)*sigma)

	voigtEquation = a\
	*(exp(firstArg**2)*erfc(firstArg) \
	+ exp(secondArg**2)*erfc(secondArg) ) \
	/ (2*sqrt(2*pi)*sigma)
	
	voigtEquation = float64(voigtEquation)
	
	return voigtEquation

	
def voigtResiduals(p, y, x): 
	err = y - voigtFunction(x,p) 
	return err

def voigtFit(x,y, p0):
	import scipy
	from scipy.optimize import leastsq
	
	plsq = leastsq(voigtResiduals, p0, args=(y, x), maxfev=2000)
	
	return [plsq[0], voigtFunction(x, plsq[0])]

# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def CalculateVoigtScore(logNRatio, plsq, normalizationFactor):
	
	voigtScore = voigtFunction(logNRatio, plsq)/normalizationFactor

	return voigtScore
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddresses2(addresses_r, addresses_c, addresses_pr, addresses_pc, plsqDict, \
normalizationFactorDict):
# Calculate the unambiguous library addresses that that a read alignment coord maps to
# Same as CalculatePossibleAddresses
# However, the reported score is slightly different. It is a collection of the ratios of pool axes
# entries. 

# Bookmark: this is where I'll add the Voigt function
	
	import numpy
	import pdb

	possibleAddressesAndScores = []
	
	for address_r in addresses_r:
		for address_c in addresses_c:
			for address_pr in addresses_pr:
				for address_pc in addresses_pc:
					
					row = address_r[0]
					col = address_c[0]
					pr = address_pr[0]
					pc = address_pc[0]
					
					possibleAddress = row + '_' + col + '_' + pr + '_' + pc
					
					nRowReads = address_r[1]
					nColReads = address_c[1]
					nPRReads = address_pr[1]
					nPCReads = address_pc[1]
					
					totalReads = nRowReads + nColReads + nPRReads + nPCReads
					
					nr2nc = nRowReads/nColReads
					nr2npr = nRowReads/nPRReads
					nr2npc = nRowReads/nPCReads
					nc2npr = nColReads/nPRReads
					nc2npc = nColReads/nPCReads
					npr2npc = nPRReads/nPCReads
					
					logNr2nc = numpy.log(nr2nc)
					logNr2npr = numpy.log(nr2npr)
					logNr2npc = numpy.log(nr2npc)
					logNc2npr = numpy.log(nc2npr)
					logNc2npc = numpy.log(nc2npc)
					logNpr2npc = numpy.log(npr2npc)
										
					voigtScoreNr2nc = CalculateVoigtScore(logNr2nc, plsqDict['nr2nc'], \
					normalizationFactorDict['nr2nc'])
					
					voigtScoreNr2npr = CalculateVoigtScore(logNr2npr, plsqDict['nr2npr'], \
					normalizationFactorDict['nr2npr'])
					
					voigtScoreNr2npc = CalculateVoigtScore(logNr2npc, plsqDict['nr2npc'], \
					normalizationFactorDict['nr2npc'])
					
					voigtScoreNc2npr = CalculateVoigtScore(logNc2npr, plsqDict['nc2npr'], \
					normalizationFactorDict['nc2npr'])
					
					voigtScoreNc2npc = CalculateVoigtScore(logNc2npc, plsqDict['nc2npc'], \
					normalizationFactorDict['nc2npc'])
					
					voigtScoreNpr2npc = CalculateVoigtScore(logNpr2npc, plsqDict['npr2npc'], \
					normalizationFactorDict['npr2npc'])

					scoreDict = {'nr2nc':voigtScoreNr2nc, 'nr2npr':voigtScoreNr2npr, \
					'nr2npc':voigtScoreNr2npc, 'nc2npr':voigtScoreNc2npr, \
					'nc2npc':voigtScoreNc2npc, 'npr2npc':voigtScoreNpr2npc}
					
					logReadCountRatioDict = {'logNr2nc':logNr2nc, 'logNr2npr':logNr2npr, \
					'logNr2npc':logNr2npc, 'logNc2npr':logNc2npr, 'logNc2npc':logNc2npc, \
					'logNpr2npc':logNpr2npc}
					
					score = voigtScoreNr2nc * voigtScoreNr2npr * voigtScoreNr2npc \
					* voigtScoreNc2npr * voigtScoreNc2npc * voigtScoreNpr2npc

					possibleAddressesAndScores.append([possibleAddress, scoreDict, score, \
					totalReads, logReadCountRatioDict])

	return possibleAddressesAndScores
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def CalculatePhysicalAddresses(poolPresenceTable, rowPools, colPools, prPools, pcPools, \
threshold=5):
	
	physicalAddressDict = GenerateBlankPhysicalAddressDict(poolPresenceTable)
	
	i = 0
	while i < len(poolPresenceTable):
		entry = poolPresenceTable[i]
		coord = int(poolPresenceTable[i]['readAlignmentCoord'])
		
		addresses_r = FindAddressCoords(poolPresenceTable[i], rowPools, threshold=threshold)
		addresses_c = FindAddressCoords(poolPresenceTable[i], colPools, threshold=threshold)
		addresses_pr = FindAddressCoords(poolPresenceTable[i], prPools, threshold=threshold)
		addresses_pc = FindAddressCoords(poolPresenceTable[i], pcPools, threshold=threshold)
		
		possibleAddresses = CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, \
		addresses_pc)
		
		physicalAddressDict[coord] = possibleAddresses
		
		i += 1
	
	return physicalAddressDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePhysicalAddressCoordsForLine(poolPresenceTableLine, rowPools, colPools, prPools, \
pcPools, threshold=5):

	addresses_r = FindAddressCoords2(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords2(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords2(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords2(poolPresenceTableLine, pcPools, threshold=threshold)
	
	# Remember, each line in the addresses array is a 2 element list, the first containing the pool
	# name, and the second containing the number of reads associated with it. 
	
	return [addresses_r, addresses_c, addresses_pr, addresses_pc]
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def CalculateAverageReadAlignmentCoordinate(group, rowPools, colPools, prPools, pcPools, \
averagingType='median'):
	
	import numpy
	from pdb import set_trace
	import scipy.stats
	
	# Calculate the median and mean read alignment coordinate
	[readAlignmentCoords, totalReads] = \
	GenerateHistogramOfReadsVersusReadAlignmentCoord(group, rowPools, colPools, prPools, pcPools)
	
	readAlignmentCoordList = []
	i = 0
	while i < len(totalReads):
		j = 0
		while j < totalReads[i]:
			readAlignmentCoordList.append(readAlignmentCoords[i])
			j += 1
		i += 1
	
# 	set_trace()
	
# 	print(str(readAlignmentCoordList))
	
	if len(readAlignmentCoordList) == 0:
		averageReadAlignmentCoord = 0
		includeInSummedTable = False
	elif averagingType == 'median':
		averageReadAlignmentCoord = int(numpy.median(readAlignmentCoordList))
		includeInSummedTable = True
	elif averagingType == 'mode':
		averageReadAlignmentCoord = int(scipy.stats.mode(readAlignmentCoordList))
		includeInSummedTable = True
	elif averagingType == 'mean':
		averageReadAlignmentCoord = int(numpy.mean(readAlignmentCoordList))
		includeInSummedTable = True
	else:
		averageReadAlignmentCoord = int(numpy.median(readAlignmentCoordList))
		includeInSummedTable = True
	
	return [averageReadAlignmentCoord, includeInSummedTable]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def GenerateHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
rowPools, colPools, prPools, pcPools):

	i = 0
	
	readAlignmentCoords = []
	totalReads = []
	
	while i < len(groupedPoolPresenceTableGroup):
		
		readAlignmentCoord = groupedPoolPresenceTableGroup[i]['readAlignmentCoord']
		
		readAlignmentCoords.append(readAlignmentCoord)
		
		readCount = \
		CountReadsAssociatedWithCoordinateThatAreInLocationPools(groupedPoolPresenceTableGroup[i],\
		rowPools, colPools, prPools, pcPools)
		
		totalReads.append(readCount)
	
		i += 1
	return [readAlignmentCoords, totalReads]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CountReadsAssociatedWithCoordinateThatAreInLocationPools(groupedPoolPresenceTableGroupLine,\
rowPools, colPools, prPools, pcPools):
	totalReads = 0
	
	i = 0
	while i < len(rowPools):
		totalReads += groupedPoolPresenceTableGroupLine[rowPools[i]]
		i += 1
	
	i = 0
	while i < len(colPools):
		totalReads += groupedPoolPresenceTableGroupLine[colPools[i]]
		i += 1
	
	i = 0
	while i < len(prPools):
		totalReads += groupedPoolPresenceTableGroupLine[prPools[i]]
		i += 1
		
	i = 0
	while i < len(pcPools):
		totalReads += groupedPoolPresenceTableGroupLine[pcPools[i]]
		i += 1

	return totalReads
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CountReadsAssociatedWithLineBooleanOperation(axisCoordIntersections, nextLineCoords, \
currentLineCoords):
	
	readsCount = 0
	
	for coord in axisCoordIntersections:
		if coord in nextLineCoords[0]:
			readsCount += nextLineCoords[1][coord]
		if coord in currentLineCoords[0]:
			readsCount += currentLineCoords[1][coord]
	
	return readsCount
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def CalculateTotalReadsInLine(poolPresenceTableLine, rowPools, colPools, prPools, pcPools, \
controlPools):

	i = 0
	totalReads = 0
	while i < len(rowPools):
		totalReads += poolPresenceTableLine[rowPools[i]]
		i +=1
	
	i = 0
	while i < len(colPools):
		totalReads += poolPresenceTableLine[colPools[i]]
		i +=1
	
	i = 0
	while i < len(prPools):
		totalReads += poolPresenceTableLine[prPools[i]]
		i +=1

	
	i = 0
	while i < len(pcPools):
		totalReads += poolPresenceTableLine[pcPools[i]]
		i +=1	

	if controlPools != None:
		i = 0
		while i < len(controlPools):
			totalReads += poolPresenceTableLine[controlPools[i]]
			i +=1
	
	return totalReads
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def CalculatePoolPresenceTableTaxonomyDict(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):

	keysInPrintOrder = [\
	'totalLines',\
	'linesThatMapToLibraryAddresses',\
	'linesThatMapToSingleLibraryAddresses',\
	'linesThatMapToMultipleLibraryAddresses',\
	'linesThatMapToUnambiguousLibraryAddresses',\
	'linesThatMapToAmbiguousLibraryAddresses', \
	'linesThatDoNotMapToLibraryAddresses',\
	'linesThatHaveNoReadsAboveThresholdInAnyPool',\
	'linesThatMapToControlIndexesOnly',\
	'linesThatHaveCoordinatesInNoPoolAxis',\
	'linesThatHaveCoordinatesInOnlyOnePoolAxis',\
	'linesThatHaveCoordinatesInOnlyTwoPoolAxes',\
	'linesThatHaveCoordinatesInOnlyThreePoolAxes',\
	'totalReads']

	# Define the pool presence line taxonomy dict 
	poolPresenceTableTaxonomyDict = {}

	for key in keysInPrintOrder:
		poolPresenceTableTaxonomyDict[key] = 0

	poolPresenceTablePoolsCoordsList = {}
	possibleAddressesDict = {}
	numberPoolAxesWithEntriesForLine = {}
	numberPoolAxesWithMoreThanOneEntry = {}

	i = 0
	while i < len(poolPresenceTable):
		
		readsInLine = CalculateTotalReadsInLine(poolPresenceTable[i], \
		rowPools, colPools, prPools, pcPools, controlPools)
		
		poolPresenceTableTaxonomyDict['totalReads'] += readsInLine
		
		coord = int(poolPresenceTable[i]['readAlignmentCoord'])
		
		poolPresenceTableTaxonomyDict['totalLines'] += 1
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		poolCoordsForLine = [addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] 
		
		[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
		CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
		addresses_control)
	
		
		# Calculate if there will be an ambiguous library address calculation by calculating the number
		# of pool axes with more than one entry. One axis with multiple (even if possible coords are
		# filled) is fine, but more axis with more than one entry leads to cross terms that are 
		# ambiguous. 
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
	

		if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries == 0):
			poolPresenceTableTaxonomyDict['linesThatHaveNoReadsAboveThresholdInAnyPool'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries > 0):
			poolPresenceTableTaxonomyDict['linesThatMapToControlIndexesOnly'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 1):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyOnePoolAxis'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 2):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyTwoPoolAxes'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 3):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyThreePoolAxes'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
	
		if (nAddressPoolAxesWithEntries == 0):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInNoPoolAxis'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord))
	
	
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			poolPresenceTableTaxonomyDict['linesThatMapToAmbiguousLibraryAddresses'] += 1


		if nAddressPoolAxesWithMoreThanOneEntry <= 1 and lenPossibleAddresses >= 1:
			poolPresenceTableTaxonomyDict['linesThatMapToUnambiguousLibraryAddresses'] += 1
	
		
	
		if lenPossibleAddresses == 0:
			poolPresenceTableTaxonomyDict['linesThatDoNotMapToLibraryAddresses'] += 1
		elif lenPossibleAddresses == 1:
			poolPresenceTableTaxonomyDict['linesThatMapToSingleLibraryAddresses'] += 1
			poolPresenceTableTaxonomyDict['linesThatMapToLibraryAddresses'] += 1
		elif lenPossibleAddresses > 1:
			poolPresenceTableTaxonomyDict['linesThatMapToMultipleLibraryAddresses'] += 1
			poolPresenceTableTaxonomyDict['linesThatMapToLibraryAddresses'] += 1
		
		i += 1
	
	PrintPoolPresenceTableTaxonomyDict(threshold, poolPresenceTableTaxonomyDict, keysInPrintOrder)
	
	return poolPresenceTableTaxonomyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PrintPoolPresenceTableTaxonomyDict(threshold, poolPresenceTableTaxonomyDict, keysInPrintOrder):
	
	outputStr = ''
	
	outputStr += 'threshold: ' + str(threshold) + '\n'
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(poolPresenceTableTaxonomyDict[key]) + '\n'
	
	print(outputStr)
	return
# ------------------------------------------------------------------------------------------------ #






####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 6: Functions for initially populating the pool presence table
# ------------------------------------------------------------------------------------------------ #
####################################################################################################


# ------------------------------------------------------------------------------------------------ #
def GenerateUniqueCoordsList(genomeArray):
	from scipy import unique
	
	coords = genomeArray['readAlignmentCoord']
	
	uniqueCoords = unique(coords)
	
	return uniqueCoords
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def GenerateValidCoordsList(genomeArray):
	
	import numpy
	from numpy import int32
	
	validCoords = 0
	
	i = 0
	maxReadIDLength = 15
	
	while i < len(genomeArray):
		
		coord = genomeArray[i]
		
		if coord['alignmentQuality'] > 1 and coord['alignmentFound'] == 1 \
		and coord['multipleAlignmentsFound'] == 0 and coord['himarRecognized'] == 1 \
		and coord['index'] > 0:
			validCoords += 1
			readIDLength = len(coord[0])
			if readIDLength > maxReadIDLength:
				maxReadIDLength = readIDLength
			
		i += 1
	
	readIDFieldCode = 'a' + str(maxReadIDLength+2)
	
	validGenomeArray = numpy.zeros(validCoords, \
	dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'index'], \
	'formats':[readIDFieldCode, int32, int32, int32, int32]})
	
	i = 0
	j = 0
	
	while i < len(genomeArray) and j < validCoords:
		
		coord = genomeArray[i]
		
		if coord['alignmentQuality'] > 1 and coord['alignmentFound'] == 1 \
		and coord['multipleAlignmentsFound'] == 0 and coord['himarRecognized'] == 1 \
		and coord['index'] > 0 and coord['strangeFlagsSum'] == 0:
			validGenomeArray[j]['readID'] = coord['readID']
			validGenomeArray[j]['readAlignmentCoord'] = coord['readAlignmentCoord']
			validGenomeArray[j]['alignmentQuality'] = coord['alignmentQuality']
			validGenomeArray[j]['index'] = coord['index']
			j += 1
		
		i += 1
	
	return validGenomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GeneratePoolNameToPoolCodeLookupTable(barcodeFile):
# Note that this
	barcodeFileHandle = open(barcodeFile, 'r')
	barcodeFileData = barcodeFileHandle.readlines()
	barcodeFileHandle.close()

	indexLookupTable = {}

	for line in barcodeFileData:
		if line[0] != '#':
			lineData = line.strip().split(',')
			poolName = lineData[0]
			forwardSeq = lineData[1]
			revCompl = lineData[2]
			barcodeNumber = lineData[3]
			indexLookupTable[str(barcodeNumber)] = poolName
		
	
	return indexLookupTable
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def BuildInitialPoolPresenceTable(genomeArray, outputLog, barcodeFile, poolPresenceTableFileName):
	
	import numpy
	import pdb
	
	outputStr = gSeparatorString
	outputStr += 'Building Initial Pool Presence Table\n'
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	
	# Compile the valid reads
	outputStr = "Making Valid Genome Array\n"
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	validGenomeArray = GenerateValidCoordsList(genomeArray)
	validGenomeArray = numpy.sort(validGenomeArray, order='readAlignmentCoord')
	
	# Generate the unique coordinates list
	outputStr = "Generating Unique Coordinates List\n"
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	uniqueCoords = GenerateUniqueCoordsList(validGenomeArray)

	# Make the first round of the pool presence table
	indexLookupTable = GeneratePoolNameToPoolCodeLookupTable(barcodeFile)

	dtypeArray = GenerateDTypeArrayForPoolPresenceDict(indexLookupTable)
	
	poolKeys = sorted(indexLookupTable.keys())
	poolColumns = ['readAlignmentCoord']
	i = 0
	while i < len(poolKeys):
		poolColumns.append(indexLookupTable[poolKeys[i]])
		i += 1
	
	print("Generating Pool Presence Table")
	poolPresenceTable = GeneratePoolPresenceTable(uniqueCoords, validGenomeArray, \
	indexLookupTable, dtypeArray)
		
	WritePoolPresenceTable3(poolPresenceTableFileName, poolPresenceTable, poolColumns)

	return poolPresenceTable
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 7: Functions for analyzing the pool presence table
# ------------------------------------------------------------------------------------------------ #
####################################################################################################



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
def SolvePoolPresenceTable(sudokuGridLookupDict, poolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools, logReadNumberRatioHistogramFitDict, \
logReadNumberRatioHistogramIntegralDict, \
readCountThreshold, voigtScoreThreshold, maxSinglePoolCoordNumber, maxTotalCoords, \
maxGapForCoordGrouping):

	import pdb

	# pdb.set_trace()
	
	PopulateSudokuGrid3(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
	prPools, pcPools, controlPools, logReadNumberRatioHistogramFitDict, \
	logReadNumberRatioHistogramIntegralDict, readCountThreshold, voigtScoreThreshold, \
	maxSinglePoolCoordNumber, maxTotalCoords)
	
	sudokuGridTaxonomyDictPreGrouping = CalculateSudokuGridOccupancyTaxonomy(sudokuGridLookupDict, \
	rowPools, colPools, prPools, pcPools)
	
	GroupReadAlignmentCoordsInSudokuGrid(sudokuGridLookupDict, prPools, pcPools, rowPools, \
	colPools, maxGap=maxGapForCoordGrouping)
	
	sudokuGridTaxonomyDict = CalculateSudokuGridOccupancyTaxonomy(sudokuGridLookupDict, rowPools, \
	colPools, prPools, pcPools)
	
	print('Sudoku Taxonomy Pre-Grouping')
	PrintSudokuGridOccupancyTaxonomy(sudokuGridTaxonomyDictPreGrouping)
	print('Sudoku Taxonomy Post-Grouping')
	PrintSudokuGridOccupancyTaxonomy(sudokuGridTaxonomyDict)
	
	
	
	return sudokuGridTaxonomyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GroupSudokuGenomicCoords(coordArray, maxGap=1):
	
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
			currentCoord = coordArray[i].coord
			nextCoord = coordArray[i+1].coord
			expect = currentCoord + maxGap
			if nextCoord <= expect:
				currentLineMatchesPrevious = True
		i += 1
		
		
	groupedCoords = []
	
	i = 0
	while i < len(result):
		
		if len(result[i]) == 1:
# 			pdb.set_trace()
			groupedCoords.append(result[i][0])
		
		elif len(result[i]) > 1:
			coords = sorted(result[i], key=operator.attrgetter('readCount'), reverse=True)
			
			j = 0
			readCountList = []
			while j < len(coords):
				k = 0
				while k < coords[j].readCount:
					readCountList.append(coords[j].coord)
					k += 1
				j += 1
			
			representativeCoord = numpy.median(readCountList)
			
			j = 0 
			totalReadCount = 0
			locatabilityArray = []
			while j < len(coords):
				totalReadCount += coords[j].readCount
				locatabilityArray.append(coords[j].locatability)
				j +=1
						
			if 'unambiguous' in locatabilityArray:
				locatability = 'unambiguous'
			elif 'ambiguous' in locatabilityArray:
				locatability = 'ambiguous'
			elif 'guessable' in locatabilityArray:
				locatability = 'guessable'
			else:
				locatability = 'merged'
			
			
			locatabilityScore = 'merged'
			
			groupedCoords.append(SudokuGenomicCoord(representativeCoord, locatability, \
			locatabilityScore, totalReadCount))
		
		i += 1
	
	return groupedCoords
	
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GroupReadAlignmentCoordsInSudokuGrid(sudokuGridLookupDict, prPools, pcPools, rowPools, \
colPools, maxGap=4):
	
	import operator
	
	for prPool in prPools:
		
		for pcPool in pcPools:
		
			sudokuPlate = None
			try:
				sudokuPlate = sudokuGridLookupDict[prPool][pcPool]
			except IndexError:
				print('No plate at: ' + rowPool + '_' + colPool)
				pass
			
		
			if sudokuPlate != None:
				plateName = sudokuPlate.plateName
		
				for colPool in colPools:
					for rowPool in rowPools:
						sudokuWell = sudokuPlate.wellGrid[rowPool][colPool]
						readAlignmentCoords = sudokuWell.readAlignmentCoords
						
						readAlignmentCoords = sorted(readAlignmentCoords, \
						key=operator.attrgetter('coord'))
						
						groupedReadAlignmentCoords = GroupSudokuGenomicCoords(readAlignmentCoords, \
						maxGap=maxGap)
						sudokuWell.readAlignmentCoords = groupedReadAlignmentCoords
						
						
# ------------------------------------------------------------------------------------------------ #

