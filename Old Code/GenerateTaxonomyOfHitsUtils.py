# ------------------------------------------------------------------------------------------------ #
def ImportSimplifiedCatalogOfHits(filename, sourcePlateCol=0, sourceWellCol=1):
	
	import re
	import pdb
	from numpy import zeros
	
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	fileHandle = open(filename, 'r')
	data = fileHandle.readlines()
	
	
	dtypeArray = [('sourcePlate', '|S10'), ('sourceRow', '|S10'), ('sourceCol', 'i8')]
	
	rearrayedWellList = []
	nonCommentLines = 0
	i = 0
	while i < len(data):
		if data[i][0] != '#':
			dataLine = data[i].strip().split(',')
			
			sourcePlate = str(dataLine[sourcePlateCol])
			sourceWell = dataLine[sourceWellCol]
			
			sourceWellMatch = wellIDRe.search(sourceWell)
		
			if sourceWellMatch != None:
				sourceRow = sourceWellMatch.group(1)
				sourceCol = sourceWellMatch.group(2)
			else:
				pdb.set_trace()


			rearrayedWellList.append([sourcePlate, sourceRow, sourceCol])
			
			nonCommentLines += 1
			
		i += 1
	
	rearrayedWellArray = zeros(nonCommentLines, dtype=dtypeArray)
	i = 0
	while i < len(rearrayedWellArray):
		rearrayedWellArray[i]['sourcePlate'] = rearrayedWellList[i][0]
		rearrayedWellArray[i]['sourceRow'] = rearrayedWellList[i][1]
		rearrayedWellArray[i]['sourceCol'] = rearrayedWellList[i][2]
		
		i += 1
		
	return rearrayedWellArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SortFeatureArrayAndAssignNumber(featureArray):

	import operator
	
	featureArray = sorted(featureArray, key=operator.attrgetter('startCoord'))
	
	locusNumber = 1
	for feature in featureArray:
		feature.tagDict['locus_number'] = locusNumber
		locusNumber += 1

	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SelectFeaturesWithHits(geneFeatureArray):
	
	featureArrayWithHits = []
	for feature in geneFeatureArray:
		if len(feature.sudokuGridEntries) > 0:
			featureArrayWithHits.append(feature)
	
	return featureArrayWithHits
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SelectWellsWithHits(sudokuWells, hits):
	
	wellsWithHits = []
	
	for well in sudokuWells:
		plateName = well.addressDict['colonyPurified']['plateName']
		row = well.addressDict['colonyPurified']['row']
		col = well.addressDict['colonyPurified']['col']

		i = 0
		while i < len(hits):
			hit = hits[i]
			hitPlateName = hit['sourcePlate'].decode('utf-8')
			hitRow =  hit['sourceRow'].decode('utf-8')
			hitCol =  str(hit['sourceCol'])
			
			if hitPlateName == plateName and hitRow == row and hitCol == col:
				wellsWithHits.append(well)
			i += 1

	return wellsWithHits


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GroupHitsByGenomicLocation(hitsFeatureArray, maxGap=5):
	
	i = 0
	expect = None
	run = []
	result = [run]
	currentLineMatchesPrevious = True
	
	while i < len(hitsFeatureArray):
		if currentLineMatchesPrevious:
			run.append(hitsFeatureArray[i])
		else:
			run = [hitsFeatureArray[i]]
			result.append(run)
			
		currentLineMatchesPrevious = False
		
		if i < len(hitsFeatureArray) - 1:
			currentCoord = hitsFeatureArray[i].tagDict['locus_number']
			nextCoord = hitsFeatureArray[i+1].tagDict['locus_number']
			expect = currentCoord + maxGap
			if nextCoord <= expect:
				currentLineMatchesPrevious = True
		i += 1
	
	

	return result
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindFeatureNameCorrespondingToStartAndEndCoords(queryStartCoord, queryEndCoord, featureArray):

	import pdb
	from numpy import mean

	queryMidCoord = int(mean([queryStartCoord, queryEndCoord]))

	featureNamesToReturn = []
	for feature in featureArray:
		featureStartCoord = feature.startCoord
		featureEndCoord = feature.endCoord
								
		if featureStartCoord <= queryMidCoord <= featureEndCoord:
			featureNamesToReturn.append(feature.featureName)
	
	return featureNamesToReturn						
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def MakeGeneFeatureArray(featureArray):
	from copy import deepcopy
	
	geneFeatureArray = []
	for feature in featureArray:
		if feature.featureType == 'gene':
			geneFeatureArray.append(deepcopy(feature))


	return geneFeatureArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def BlastFeatureAgainstSubjectGenome(queryFeature, referenceGenomeFile, referenceFeatureArray, \
removeTempFiles=True, tempFileDir='./', evalue=0.1, task='blastn', perc_identity=50):
	
	from Bio.Blast.Applications import NcbiblastnCommandline
	from Bio.Blast import NCBIXML
	import os
	
	featureName = queryFeature.featureName
	querySequence = queryFeature.tagDict['coding_sequence']
	
	querySequenceFileName = tempFileDir + '/' + featureName + '.seq'
	blastResultFileName = tempFileDir + '/' + featureName + '.xml'
	
	fileHandle = open(querySequenceFileName, 'w')
	fileHandle.write('>' + featureName + '\n')
	fileHandle.write(querySequence)
	fileHandle.close()
	

	blastn_cline = NcbiblastnCommandline(query=querySequenceFileName, subject=referenceGenomeFile, \
	evalue=evalue, outfmt=5, out=blastResultFileName, task=task, perc_identity=perc_identity)

	blastn_cline()

	result_handle = open(blastResultFileName)

	xmlresults = NCBIXML.parse(result_handle)
	xmlresults = list(xmlresults)

	parsedResults = []

	hsps = []
	for result in xmlresults:
		alignments = []
		for alignment in result.alignments:
			hsps = []
			for hsp in alignment.hsps:
				start = int(hsp.sbjct_start)
				end = int(hsp.sbjct_end)
				evalue = hsp.expect
				
				featureNames = \
				FindFeatureNameCorrespondingToStartAndEndCoords(start, end, referenceFeatureArray)
				
				hsps.append({'start':start, 'end':end, 'evalue':evalue, \
				'featureNames':featureNames})
				
	if removeTempFiles:
		os.remove(querySequenceFileName)
		os.remove(blastResultFileName)

	return hsps
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def MakeTaxonomyOfGroupedHits(groupedHits):
	singleMemberGroups = 0
	multiMemberGroups = 0
	for group in groupedHits:
		if len(group) > 1:
			multiMemberGroups += 1
		else:
			singleMemberGroups += 1
	print('Single member groups: ' + str(singleMemberGroups))
	print('Multiple member groups: ' + str(multiMemberGroups))

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def MakeSummaryTableOfFeaturesWithHits(fileName, featuresWithHits, delimeter=',', \
printOnlyOneCoord=True):
	
	outputDictKeys = ['featureName','locus_tag','locus_number','startCoord','endCoord',\
	'readAlignmentCoord','fracDistFromTranslationStart','coding_sequence']
	
	outputHandle = open(fileName, 'w')
	
	headerStr = 'Feature Name,Locus Tag,Locus Number,Start Coord,End Coord,Transposon Location,' \
	+ 'Fractional Distance From Translation Start,Coding Sequence\n'
	
	outputStr = ''
	
	for feature in featuresWithHits:
		
		if printOnlyOneCoord == True:
			maxLength = 1
		else:
			maxLength = len(feature.sudokuGridEntries)
		
		i = 0
		while i < maxLength:
			coord = feature.sudokuGridEntries[i]		
			outputDict = {}
			outputDict['endCoord'] = feature.endCoord
			outputDict['startCoord'] = feature.startCoord
			outputDict['fracDistFromTranslationStart'] = coord['fracDistFromTranslationStart']
			outputDict['readAlignmentCoord'] = coord['readAlignmentCoord']
			outputDict['featureName'] = feature.featureName
			outputDict['locus_number'] = feature.tagDict['locus_number']
			outputDict['locus_tag'] = feature.tagDict['locus_tag'][0]
			outputDict['coding_sequence'] = feature.tagDict['coding_sequence']
		
			lineStr = ''
			j = 0
			while j < len(outputDictKeys):
				lineStr += str(outputDict[outputDictKeys[j]])
				if j < len(outputDictKeys) - 1:
					lineStr += delimeter
				else:
					lineStr += '\n'
				j += 1
				
			outputStr += lineStr
			i += 1
	
	totalStr = headerStr + outputStr
	outputHandle.write(totalStr)
	outputHandle.close()
	return totalStr
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# def MakeSummaryTableOfSudokuWellsWithHits(fileName, flattenedCPSudokuGrid):
# 	
# 	outputDictKeys = ['featureName','locus_tag','locus_number','startCoord','endCoord',\
# 	'readAlignmentCoord','fracDistFromTranslationStart','coding_sequence']
# 	
# 	outputHandle = open(fileName, 'w')
# 	
# 	headerStr = 'Feature Name,Locus Tag,Locus Number,Plate,Well,Start Coord,End Coord,' \
# 	+ 'Transposon Location,' \
# 	+ 'Fractional Distance From Translation Start,Coding Sequence\n'
# 	
# 	outputStr = ''
# 	
# 	
# 	
# 	totalStr = headerStr + outputStr
# 	outputHandle.write(totalStr)
# 	outputHandle.close()
# 
# 	return totalStr
# 
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def AnnotateSudokuWellListWithFeautureNames(flattenedSudokuWellList, featureArray):

	import pdb
	
	for well in flattenedSudokuWellList:

		readAlignmentCoords = well.readAlignmentCoords
	
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
