# Sanger





####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Sanger verification code

# ------------------------------------------------------------------------------------------------ #
class BlastAlignment:

	def __init__(self, score, expectation):
		
		import re
		
		self.score = score
		self.expectation = expectation
		self.strandOrientation = ''
		self.alignmentLines = []
		self.readAlignmentCoord = 0
		self.strandOrientationRe = re.compile(r"(Plus|Minus)/(Minus|Plus)")
		self.alignmentSubjectLineRe = re.compile(r'Sbjct\s+(\d+)\s+\w+\s+(\d+)')
	
	def CalculateReadAlignmentCoord(self):
		
		import pdb
		
# 		print(self.strandOrientation)
		strandOrientationMatch = self.strandOrientationRe.search(self.strandOrientation)
		
		strandOrientationQuery = strandOrientationMatch.group(1)
		strandOrientationSubject = strandOrientationMatch.group(2)
		
		if strandOrientationQuery == 'Plus':
			if strandOrientationSubject == 'Minus':
# 				pdb.set_trace()
				readAlignmentCoord = self.alignmentSubjectLineRe.search(self.alignmentLines[0]).group(1)
			elif strandOrientationSubject == 'Plus':
				readAlignmentCoord = self.alignmentSubjectLineRe.search(self.alignmentLines[0]).group(1)
		elif strandOrientationQuery == 'Minus':
			print("Error!")
			readAlignmentCoord = 0
			
		self.readAlignmentCoord = readAlignmentCoord
		
		return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ParseBlastOutput3(alignmentLines):
	
	import re
	import pdb
	
	strandOrientationRe = re.compile(r"Strand=((Plus|Minus)/(Minus|Plus))")
	alignmentSubjectLineRe = re.compile(r'Sbjct\s+(\d+)\s+\w+\s+(\d+)')
	noHitsRe = re.compile(r'\*\*\*\*\* No hits found \*\*\*\*\*')
	newAlignmentRe = re.compile(r'Score\s+=\s+(\d+)\s+bits\s+\(\d+\),\s+Expect\s+=\s+(\d+e*\.*[+|-]*\d+)')
	
	alignments = []
	noHitsFound = False
	currentAlignment = None
	
	for line in alignmentLines:
		if noHitsRe.search(line) != None:
			noHitsFound = True
			readAlignmentCoords = None
	
	
	if noHitsFound == False:
		for line in alignmentLines:
		
			newAlignmentMatch = newAlignmentRe.search(line)
			strandOrientationMatch = strandOrientationRe.search(line)
			alignmentSubjectLineMatch = alignmentSubjectLineRe.search(line)
		
			if newAlignmentMatch != None:
			
				if currentAlignment != None:
# 					print(currentAlignment.score)
					currentAlignment.CalculateReadAlignmentCoord()
					alignments.append(currentAlignment)
			
				score = newAlignmentMatch.group(1)
				expectation = newAlignmentMatch.group(2)
				currentAlignment = BlastAlignment(score, expectation)
		
			elif (newAlignmentMatch == None) and (currentAlignment != None):
				if strandOrientationMatch != None:
					currentAlignment.strandOrientation = strandOrientationMatch.group(1)
				if alignmentSubjectLineMatch != None:
					currentAlignment.alignmentLines.append(line)
			
	if noHitsFound == False:
# 		pdb.set_trace()
# 		print(currentAlignment.score)	
		currentAlignment.CalculateReadAlignmentCoord()
		alignments.append(currentAlignment)	
	

	return alignments
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AlignSequenceToReferenceByBlast(referenceFile, queryFilePath):

	import subprocess
	
	alignmentData = \
	subprocess.check_output(['blastn', '-subject', referenceFile, '-query', queryFilePath])

	alignmentStr = alignmentData.decode("utf-8")
	alignmentLines = alignmentStr.split('\n')
	
	# Find read alignment coordinates
	readAlignmentCoords = ParseBlastOutput3(alignmentLines)

	return readAlignmentCoords

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportVerificationDataFile(sangerVerificationFile):
	
	fileHandle = open(sangerVerificationFile, 'r')
	
	data = fileHandle.readlines()
	
	verificationDict = {}
	
	i = 1
	while i < len(data):
		
		dataLine = data[i].strip().split(',')
		
		sequenceName = dataLine[0]
		fileDir = dataLine[1]
		fileName = dataLine[2]
		well = dataLine[3]
		plate = dataLine[4]
		predictionRound = dataLine[5]
		
		
		verificationDict[sequenceName] = \
		{'fileDir':fileDir, 'fileName':fileName, 'well':well, 'plate':plate, \
		'predictionRound':predictionRound}
		
		i += 1
	
	return verificationDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateSudokuReadAlignmentOutputStrings(sudokuReadAlignmentCoords):
	
	import pdb
	
	sudokuReadAlignmentCoordsArray = []
	
	j = 0
	while j < len(sudokuReadAlignmentCoords):
	
# 		pdb.set_trace()
		
		coord = sudokuReadAlignmentCoords[j].coord
		sudokuReadAlignmentCoordsArray.append(coord)
		j += 1
	
	sudokuReadAlignmentCoordsStr = '"'
	
	j = 0
	
	while j < len(sudokuReadAlignmentCoordsArray):
		
		sudokuReadAlignmentCoordsStr += str(sudokuReadAlignmentCoordsArray[j])
		
		if j < len(sudokuReadAlignmentCoordsArray) - 1:
			sudokuReadAlignmentCoordsStr += ","
	
		j += 1
	
	
	sudokuReadAlignmentCoordsStr += '"'

	return sudokuReadAlignmentCoordsStr

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateSangerReadAlignmentOutputStrings(sangerReadAlignmentCoords):
	
	sangerReadAlignmentCoordsArray = []
	sangerScoresArray = []
	sangerExpectationsArray = []
	sangerOrientationsArray = []
	
	j = 0
	while j < len(sangerReadAlignmentCoords):
		
		coord = sangerReadAlignmentCoords[j].readAlignmentCoord
		expectation = sangerReadAlignmentCoords[j].expectation
		score = sangerReadAlignmentCoords[j].score
		orientation = sangerReadAlignmentCoords[j].strandOrientation
		
		sangerReadAlignmentCoordsArray.append(coord)
		sangerScoresArray.append(expectation)
		sangerExpectationsArray.append(score)
		sangerOrientationsArray.append(orientation)
		
		j += 1
	
	sangerReadAlignmentCoordsStr = '"'
	sangerScoresStr = '"'
	sangerExpectationsStr = '"'
	sangerOrientationsStr = '"'
	
	j = 0
	
	while j < len(sangerReadAlignmentCoordsArray):
		
		sangerReadAlignmentCoordsStr += str(sangerReadAlignmentCoordsArray[j])
		sangerScoresStr += str(sangerScoresArray[j])
		sangerExpectationsStr += str(sangerExpectationsArray[j])
		sangerOrientationsStr += str(sangerOrientationsArray[j])
		
		if j < len(sangerReadAlignmentCoordsArray) - 1:
			sangerReadAlignmentCoordsStr += ","
			sangerScoresStr += ","
			sangerExpectationsStr += ","
			sangerOrientationsStr += ","
	
		j += 1
	
	
	sangerReadAlignmentCoordsStr += '"'
	sangerScoresStr += '"'
	sangerExpectationsStr += '"'
	sangerOrientationsStr += '"'

	return [sangerReadAlignmentCoordsStr, sangerScoresStr, sangerExpectationsStr, \
	sangerOrientationsStr]

# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def OutputVerificationDict(verificationDict, outputLog):
	
	import pdb
	
	fileHandle = open(outputLog, 'w')

	verificationKeys = sorted(list(verificationDict.keys()))
	
	headerStr = "Nominal Genetic Locus,File Directory,File Name,Well,Plate,Prediction Round," \
	+ "Sanger Read Alignment Coordinates,Sanger Read Alignment Score," \
	+ "Sanger Read Alignment Expectation,Sanger Read Alignment Orientation," \
	+ "Sudoku Read Alignment Coordinates\n"
	
	fileHandle.write(headerStr)
	
	outputStr = ""
	
	i = 0
	while i < len(verificationKeys):
		key = verificationKeys[i]
		
		verificationDataDict = verificationDict[key]
		
		fileDir = verificationDataDict['fileDir']
		fileName = verificationDataDict['fileName']
		well = verificationDataDict['well']
		plate = verificationDataDict['plate']
		predictionRound = verificationDataDict['predictionRound']
		sangerReadAlignmentCoords = verificationDataDict['sangerReadAlignmentCoords']
		sudokuReadAlignmentCoords = verificationDataDict['sudokuReadAlignmentCoords']
		
		[sangerReadAlignmentCoordsStr, sangerScoresStr, sangerExpectationsStr, \
		sangerOrientationsStr] = GenerateSangerReadAlignmentOutputStrings(sangerReadAlignmentCoords)
		
		sudokuReadAlignmentOutputStr = \
		GenerateSudokuReadAlignmentOutputStrings(sudokuReadAlignmentCoords)
		
		
		outputStr += key + ',' + fileDir + ',' + fileName + ',' + well + ',' + plate \
		+ ',' + predictionRound + ',' + sangerReadAlignmentCoordsStr + ',' + sangerScoresStr \
		+ ',' + sangerExpectationsStr + ',' + sangerOrientationsStr + ',' \
		+ sudokuReadAlignmentOutputStr + '\n'
		
		
		i += 1
	
	fileHandle.write(outputStr)

	fileHandle.close()

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportCatalogOfSangerReads(filename):
	
	import re
	import pdb
	from numpy import zeros
	
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	fileHandle = open(filename, 'r')
	
	data = fileHandle.readlines()
	
	
	catalogArray = []
	i = 0
	while i < len(data):
		if data[i][0] != '#':
			dataLine = data[i].strip().split(',')
			
			sourcePlate = str(dataLine[0])
			sourceWell = dataLine[1]
			baseDir = dataLine[2]		
			fileName = dataLine[3]
		
			sourceWellMatch = wellIDRe.search(sourceWell)
		
			if sourceWellMatch != None:
				sourceRow = sourceWellMatch.group(1)
				sourceCol = sourceWellMatch.group(2)
			else:
				pdb.set_trace()

			catalogEntry = {'progPlate':sourcePlate, 'progRow':sourceRow, 'progCol':sourceCol, \
			'baseDir':baseDir, 'fileName':fileName}
		
			catalogArray.append(catalogEntry)
			
		i += 1
	
	return catalogArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def BlastSangerFileAgainstSubjectGenome(queryFileName, referenceGenomeFile, \
removeTempFiles=True, tempFileDir='./', evalue=0.1, task='blastn', perc_identity=50):
	
	from Bio.Blast.Applications import NcbiblastnCommandline
	from Bio.Blast import NCBIXML
	import os
	from kosudoku.utils import ProcessFileNameExtension
	
	[queryFileNameBase, queryFileExt] = ProcessFileNameExtension(queryFileName)
	
	blastResultFileName = tempFileDir + '/' + queryFileNameBase + '.xml'

	blastn_cline = NcbiblastnCommandline(query=queryFileName, subject=referenceGenomeFile, \
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
								
				hsps.append({'start':start, 'end':end, 'evalue':evalue})
				
	if removeTempFiles:
# 		os.remove(querySequenceFileName)
		os.remove(blastResultFileName)

	return hsps
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindFeatureNameCorrespondingToCoord(queryCoord, featureArray):

	import pdb
	from numpy import mean

	featureNamesToReturn = []
	for feature in featureArray:
		featureStartCoord = feature.startCoord
		featureEndCoord = feature.endCoord
								
		if featureStartCoord <= queryCoord <= featureEndCoord:
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

