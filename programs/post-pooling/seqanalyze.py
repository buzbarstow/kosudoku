#!/opt/local/bin/python3.7
# Change this to your preferred python distribution.
# Make sure that kosudoku package is in your python path. 

# ----------------------------------------------------------------------------------------------- #
# kosudoku-seqanalyze
# Program for reducing Illumina sequencing dataset to a set of transposon genomic insertion 
# coordinates and pool coordinates. 
# All code necessary to build pool presence table from Sudoku reads. 


# Created by Buz Barstow 2015-6-14
# Last updated by Buz Barstow 2019-04-20
# ----------------------------------------------------------------------------------------------- #



import pdb
import gc
import re
import subprocess
import sys
import datetime

from os import remove


from utilities.input import get_input


from utilities.utils import GenerateFileNamesFromFileListAndBaseDir, \
InitializeOutputFoldersAndFiles3, ensure_dir

from utilities.align import CheckHimarContentAndTrimSequences, SummarizeIndexReads, \
AlignReadsToGenome, CompileAlignmentsWithIndexAndHimarSummaries3

from utilities.pool import InitializeEmptyPoolPresenceTable, WritePoolPresenceTable3, \
UpdatePoolPresenceTable




# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = ['outputLog', 'bowtieAlignmentMode', \
'referenceIndexPrefix', 'himarSequence', 'genomeAlignmentBaseDir', 'genomeAlignmentFilePrefix', \
'barcodeFile', 'barcodePostfix', \
'poolPresenceTableFileName', \
'trimmedSequencesBaseDir', 'trimmedSequencesFilePrefix', \
'indexSummaryBaseDir', 'indexSummaryFilePrefix', \
'himarRecognitionBaseDir', \
'himarRecognitionFilePrefix', \
'fastqFileBase', 'indexFastqFileBase', 'fastqFileNumDigits', 'fastqFileFirst', 'fastqFileLast', \
'deleteTempFiles', \
'allowedNsInIndexSeq', \
'indexMismatchesAllowed', \
'indexSequenceSampleLength', \
'maxTrialsForIndexMatch']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'seqanalyze.inp'
inputParameterValues = get_input(file_name, inputParameters)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data

# Input data for himar recognition
fastqFiles = []
indexFastqFiles = []

fastqFileBase = inputParameterValues['fastqFileBase']
indexFastqFileBase = inputParameterValues['indexFastqFileBase']
fastqFileNumDigits = inputParameterValues['fastqFileNumDigits']
firstFastqFile = int(inputParameterValues['fastqFileFirst'])
lastFastqFile = int(inputParameterValues['fastqFileLast'])

for i in range(firstFastqFile,lastFastqFile + 1):
    fastqFiles.append(fastqFileBase + ("{:0" + fastqFileNumDigits + "d}").format(i) + '.fastq')
    indexFastqFiles.append(indexFastqFileBase + ("{:0" + fastqFileNumDigits + "d}").format(i) + '.fastq')

outputLog = file_intro + inputParameterValues['outputLog']
himarSequence = inputParameterValues['himarSequence']
barcodeFile = file_intro + inputParameterValues['barcodeFile']
barcodePostfix = inputParameterValues['barcodePostfix']
himarRecognitionFilePrefix = inputParameterValues['himarRecognitionFilePrefix']
himarRecognitionBaseDir = file_intro + inputParameterValues['himarRecognitionBaseDir']
trimmedSequencesFilePrefix = inputParameterValues['trimmedSequencesFilePrefix']
trimmedSequencesBaseDir = file_intro + inputParameterValues['trimmedSequencesBaseDir']
indexSummaryFilePrefix = inputParameterValues['indexSummaryFilePrefix']
indexSummaryBaseDir = file_intro + inputParameterValues['indexSummaryBaseDir']
genomeAlignmentFilePrefix = inputParameterValues['genomeAlignmentFilePrefix']
genomeAlignmentBaseDir = file_intro + inputParameterValues['genomeAlignmentBaseDir']
bowtieAlignmentMode = inputParameterValues['bowtieAlignmentMode']
referenceIndexPrefix = file_intro + inputParameterValues['referenceIndexPrefix']
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
deleteTempFiles = inputParameterValues['deleteTempFiles']
allowedNsInIndexSeq = int(inputParameterValues['allowedNsInIndexSeq'])
indexMismatchesAllowed = int(inputParameterValues['indexMismatchesAllowed'])
indexSequenceSampleLength = int(inputParameterValues['indexSequenceSampleLength'])
maxTrialsForIndexMatch = int(inputParameterValues['maxTrialsForIndexMatch'])
# ----------------------------------------------------------------------------------------------- #






# ----------------------------------------------------------------------------------------------- #
# Non-Looped Part 
# Initialize the log file and output folders
InitializeOutputFoldersAndFiles3(outputLog, himarRecognitionBaseDir, indexSummaryBaseDir, \
trimmedSequencesBaseDir, genomeAlignmentBaseDir, poolPresenceTableFileName)

# Import reads, check Himar content and write out trimmed reads and himar recognition summary
[himarRecognitionFiles, trimmedSequenceFiles] = CheckHimarContentAndTrimSequences(fastqFiles, \
himarSequence, himarRecognitionFilePrefix, himarRecognitionBaseDir, \
trimmedSequencesFilePrefix, trimmedSequencesBaseDir, outputLog)

# Import index reads, summarize with a code where possible, if not indicate an unreadable index
indexAlignmentFiles = SummarizeIndexReads(indexFastqFiles, barcodeFile, indexSummaryFilePrefix, \
indexSummaryBaseDir, outputLog, allowedNsInIndexSeq=allowedNsInIndexSeq, \
indexMismatchesAllowed=indexMismatchesAllowed, indexSequenceSampleLength=indexSequenceSampleLength,\
maxTrialsForIndexMatch=maxTrialsForIndexMatch, barcodePrefix='', barcodePostfix=barcodePostfix)


# Import the mfa files, align against the genome and then write summaries to log file
bowtieMode = '--'+bowtieAlignmentMode

genomeAlignmentFiles = AlignReadsToGenome(trimmedSequenceFiles, \
outputLog, genomeAlignmentBaseDir, genomeAlignmentFilePrefix, referenceIndexPrefix, \
bowtieMode=bowtieMode)

# Initialize the pool presence table here. 
poolPresenceTable, indexLookupTable, dtypeArray, poolColumns = \
InitializeEmptyPoolPresenceTable(barcodeFile, outputLog)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Looped part of algorithm
# Update the pool presence table one genome alignment file at a time to save memory

lengthGenomeAlignmentFiles = len(genomeAlignmentFiles)
lengthIndexAlignmentFiles = len(indexAlignmentFiles)
lengthHimarRecognitionFiles = len(himarRecognitionFiles)

if lengthGenomeAlignmentFiles != lengthIndexAlignmentFiles \
and lengthGenomeAlignmentFiles != lengthHimarRecognitionFiles:
    sys.exit("Error. Number of genome alignment, index alignment, and himar recognition" \
    + " files are not the same.")

i = 0
while i < lengthGenomeAlignmentFiles:

    genomeArray = CompileAlignmentsWithIndexAndHimarSummaries3(genomeAlignmentFiles[i], \
    indexAlignmentFiles[i], himarRecognitionFiles[i], outputLog)


    poolPresenceTable = \
    UpdatePoolPresenceTable(poolPresenceTable, genomeArray, indexLookupTable, dtypeArray, outputLog)

    i += 1

# Write out the pool presence table
WritePoolPresenceTable3(poolPresenceTableFileName, poolPresenceTable, poolColumns)

gc.collect()
# ----------------------------------------------------------------------------------------------- #

