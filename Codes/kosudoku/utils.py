

gSeparatorString = '# ----------------------------------------------------------------------------------------------- #\n'


# ------------------------------------------------------------------------------------------------ #
def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateFileNamesFromFileListAndBaseDir(baseDir, files):
	
	fullFileNames = []
	
	for file in files:
		fullFileNames.append(baseDir + '/' + file)
		
	return fullFileNames
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ProcessFileNameExtension(fileName):
	import re
	import os
	import pdb
	
# 	pdb.set_trace()
	
	baseName = os.path.basename(fileName)
	[fileName, fileExtension] = os.path.splitext(baseName)
		
	return [fileName, fileExtension]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def WriteFastaSequence(sequence, name, fileName, width=100):
	
	fileHandle = open(fileName, 'w')
	
	fileHandle.write('>'+name+'\n')
	
	i = 0
	while i < len(sequence):
		
		if i > 0:
			if i%width == 0:
				fileHandle.write('\n')
		
		fileHandle.write(sequence[i])
		
		i += 1
	
	fileHandle.close()
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportFastaSequence(fastaFile):
	import re
	from pdb import set_trace
	
	
	fileHandle = open(fastaFile, 'r')
	fileData = fileHandle.readlines()
	
	commentRe = re.compile('>(\w+)')
	sequence = ''
	
	i = 0
	while i < len(fileData):
		line = fileData[i].strip()
		commentMatch = commentRe.match(line)
		
		if commentMatch != None:
			sequenceName = commentMatch.group(1)
		else:
			sequence += line
		
		i += 1
	
	return sequence
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportGenBankSequence(genBankFile):
	import re
	from pdb import set_trace
	
	originRe = re.compile(r'ORIGIN')
	endRe = re.compile(r'//')
	sequenceRe = re.compile(r'(\d+)\s+((\w+\s*)+)')
	
	fileHandle = open(genBankFile, 'r')
	fileData = fileHandle.readlines()
	
	i = 0
	
	inSequence = False
	sequence = ''

	sequenceLineCount = 0
	
	while i < len(fileData):
		
		line = fileData[i].strip()
		
		originMatch = originRe.match(line)
		
		if originMatch != None:
			inSequence = True
		
		endMatch = endRe.search(line)
		
		if endMatch != None:
			inSequence = False
		
		if (originMatch == None) and (inSequence == True):
			sequenceMatch = sequenceRe.search(line)
			
			if sequenceMatch != None:
				sequenceLineCount += 1				
				sequenceText = sequenceMatch.group(2)				
				sequence += sequenceText
			else:
				print("Line not recognized")
	
		i += 1
		

	sequenceNoSpaces = sequence.replace(" ", "")
	
	return sequenceNoSpaces
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindATandTAPositions(genomeFile, format='genbank'):
# Does the same thing as FindATandTAPositions but can work with a GenBank or a Fasta file, \
# so you only  need one file format
	
	import re
	from pdb import set_trace
	
	if format == 'genbank':
		sequence = ImportGenBankSequence(genomeFile)
	elif format == 'fasta':
		sequence = ImportFastaSequence(genomeFile)
	
	ATandTAPositions = []
	
	atRegex = re.compile('(at|ta)', re.IGNORECASE)
	
# 	set_trace()
	
	i = 0
	while i < len(sequence) - 1:
		atMatch = atRegex.match(sequence[i:i+2])
		
		if atMatch != None:
			ATandTAPositions.append(i+1)
		
		i += 1
	
	
	
	return [ATandTAPositions, sequence]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertWellIDToRowAndColumn(well):

	import re

	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	wellMatch = wellIDRe.match(well)
	
	row = wellMatch.group(1)
	col = wellMatch.group(2)
	
	return [row, col]
# ------------------------------------------------------------------------------------------------ #



####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 0: Initialize folders
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ----------------------------------------------------------------------------------------------- #
def InitializeOutputFoldersAndFiles(outputLog, himarRecognitionBaseDir, indexSummaryBaseDir, \
trimmedSequencesBaseDir, genomeAlignmentBaseDir, genomeArrayFileName, \
poolPresenceTableFileName):
	
	import os.path
	
	import pdb
	
	ensure_dir(himarRecognitionBaseDir)
	ensure_dir(indexSummaryBaseDir)
	ensure_dir(genomeAlignmentBaseDir)
	ensure_dir(trimmedSequencesBaseDir)
	
	outputLogDir = os.path.dirname(outputLog)
	genomeArrayDir = os.path.dirname(genomeArrayFileName)
	poolPresenceTableDir = os.path.dirname(poolPresenceTableFileName)
	
# 	pdb.set_trace()
	
	ensure_dir(outputLogDir + '/')
	ensure_dir(genomeArrayDir + '/')
	ensure_dir(poolPresenceTableDir + '/')
	
	fileHandle = open(outputLog, 'w')
	fileHandle.close()

# ----------------------------------------------------------------------------------------------- #




# ------------------------------------------------------------------------------------------------ #
def ParseCSVLine(line):
	
	import pdb
	
	line = line.strip()
	
	i = 0
	columns = []
	currentColumn = ''
	inQuotedBlock = False
	
	while i < len(line):
		letter = line[i]
		
		if letter == '"':
			if inQuotedBlock == True:
				inQuotedBlock = False
			elif inQuotedBlock == False:
				inQuotedBlock = True
		
		
		elif letter == ',' and inQuotedBlock == False:
			columns.append(currentColumn)
			currentColumn = ''
		
		else:
			currentColumn += letter
		
		i += 1
	
	columns.append(currentColumn)
	
# 	pdb.set_trace()
		
	
	return columns
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def ReverseComplement(sequence):
	
	import pdb
	
	revComplement = ''
	
	i = len(sequence) - 1
	while i >= 0:
		
		sequenceLetter = sequence[i]
		
		if sequenceLetter == 'a' or sequenceLetter == 'A':
			revComplementLetter = 'T'
		elif sequenceLetter == 't' or sequenceLetter == 'T':
			revComplementLetter = 'A'
		elif sequenceLetter == 'c' or sequenceLetter == 'C':
			revComplementLetter = 'G'
		elif sequenceLetter == 'g' or sequenceLetter == 'G':
			revComplementLetter = 'C'
		
		revComplement += revComplementLetter
		
		
		i -= 1

	return revComplement
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def UpdateLogFileData(outputLog, outputStr):
	
	logFileHandle = open(outputLog, 'a')
	logFileHandle.write(outputStr)
	logFileHandle.close()
# ------------------------------------------------------------------------------------------------ #

