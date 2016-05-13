# ----------------------------------------------------------------------------------------------- #
def SummarizePoolPresenceTableLine(line, rowPools, colPools, prPools, pcPools, ctrlPools, \
calculateAddresses=False, fitDict=None, areaDict=None, readCountThreshold=5, \
maxSinglePoolCoordNumber=8, maxTotalCoords=8):
	
	from sudokuutils13.pool import CalculateLibraryAddressLocatabilityForPoolPresenceTableLine
	import operator
	import pdb
	
	readAlignmentCoord = line['readAlignmentCoord']
	
	rowPoolSummary = []
	for rowPool in rowPools:
		if line[rowPool] > 0:
			rowPoolSummary.append(rowPool + ':' + str(line[rowPool]))
	
	colPoolSummary = []
	for colPool in colPools:
		if line[colPool] > 0:
			colPoolSummary.append(colPool + ':' + str(line[colPool]))
	
	prPoolSummary = []
	for prPool in prPools:
		if line[prPool] > 0:
			prPoolSummary.append(prPool + ':' + str(line[prPool]))
	
	pcPoolSummary = []
	for pcPool in pcPools:
		if line[pcPool] > 0:
			pcPoolSummary.append(pcPool + ':' + str(line[pcPool]))
	
	
	ctrlPoolSummary = []
	for ctrlPool in ctrlPools:
		if line[ctrlPool] > 0:
			ctrlPoolSummary.append(ctrlPool + ':' + str(line[ctrlPool]))
	
	
	[possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber] = \
	CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(line, rowPools, colPools, \
	prPools, pcPools, ctrlPools, fitDict, areaDict, readCountThreshold, \
	maxSinglePoolCoordNumber, maxTotalCoords)
		
	possibleAddressSummary = []
	i = 0
	if possibleAddressesAndScores != None:
		while i < len(possibleAddressesAndScores):
			summaryLine = [possibleAddressesAndScores[i][0], \
			possibleAddressesAndScores[i][1], possibleAddressesAndScores[i][2], \
			possibleAddressesAndScores[i][4]]		
			possibleAddressSummary.append(summaryLine)
			i += 1
	
		possibleAddressSummary = sorted(possibleAddressSummary, key=operator.itemgetter(2), \
		reverse=True)
		
		return [readAlignmentCoord, rowPoolSummary, colPoolSummary, prPoolSummary, pcPoolSummary, \
		ctrlPoolSummary, locatability, maxSinglePoolCoordNumber, possibleAddressSummary]
	else:
		return [readAlignmentCoord, rowPoolSummary, colPoolSummary, prPoolSummary, pcPoolSummary, \
		ctrlPoolSummary]

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GeneratePoolPresenceTableLineSummary(lineSummaryArray, delimeter=','):
	
	import pdb
	
	from sudokuutils13.xml import ExportListForXML
	
	headerStr = 'Feature Name,Read Alignment Coord,Row Pool Summary,Column Pool Summary,' \
	+ 'Plate Row Pool Summary,'\
	+ 'Plate Column Summary,Locatability,' \
	+ 'Maximum Single Pool Coordinate Number,Possible Address Summary\n'
	
	summaryStr = ''
	
	i = 0
	while i < len(lineSummaryArray):
		outputStr = ''
		
		summaryDict = lineSummaryArray[i]
		
		
		outputStr += '"' + ExportListForXML(summaryDict['featureNames']) + '"' + delimeter
		outputStr += str(summaryDict['coord']) + delimeter 
		outputStr += '"' + ExportListForXML(summaryDict['rowCoords']) + '"' + delimeter
		outputStr += '"' + ExportListForXML(summaryDict['colCoords']) + '"' + delimeter
		outputStr += '"' + ExportListForXML(summaryDict['prCoords']) + '"' + delimeter
		outputStr += '"' + ExportListForXML(summaryDict['pcCoords']) + '"' + delimeter
		
		if 'locatability' in list(summaryDict.keys()):
			outputStr += str(summaryDict['locatability']) + delimeter
		if 'maxSinglePoolCoordNumber' in list(summaryDict.keys()):
			outputStr += str(summaryDict['maxSinglePoolCoordNumber']) + delimeter
		
		if 'possibleAddressSummary' in list(summaryDict.keys()):
			scoresAndAddresses = summaryDict['possibleAddressSummary']
			j = 0
			lineStr = '"'
			while j < len(scoresAndAddresses):
				lineStr += scoresAndAddresses[j][0] + ': ' + str(scoresAndAddresses[j][2])
				if j < len(scoresAndAddresses) - 1:
					lineStr += '\n'
				j += 1
			lineStr += '"'
			outputStr += lineStr
			
		
		outputStr += '\n'
		summaryStr += outputStr
		i += 1
	
	returnStr = headerStr + summaryStr
	
	return returnStr
# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
def FindFeatureAssociatedWithInsertionCoord(insertionCoord, featureArray):
	
	from sudokuutils13.xml import ExportListForXML
	
	featureNames = []
	
	for feature in featureArray:
		startCoord = feature.startCoord
		endCoord = feature.endCoord
									
		if startCoord <= insertionCoord <= endCoord:
				
			featureName = feature.featureName
			
			featureNames.append(featureName)
		
	
	return featureNames
# ----------------------------------------------------------------------------------------------- #



