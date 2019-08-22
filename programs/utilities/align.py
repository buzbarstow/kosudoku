# ------------------------------------------------------------------------------------------------ #
# Sudoku code for alignment of massively parallel sequencing reads

# Last modified by Buz Barstow 2018-11-27
# Added new index matching functions
# ------------------------------------------------------------------------------------------------ #

from kosudoku.utils import gSeparatorString, UpdateLogFileData


# ------------------------------------------------------------------------------------------------ #
class Alignment:

    def __init__(self, readID, readAlignmentCoord, alignmentQuality, alignmentFound, \
    multipleAlignmentsFound, strangeFlagsSum):
        self.readID = readID
        self.readAlignmentCoord = readAlignmentCoord
        self.alignmentQuality = alignmentQuality
        self.alignmentFound = alignmentFound
        self.multipleAlignmentsFound = multipleAlignmentsFound
        self.strangeFlagsSum = strangeFlagsSum
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def ParseBowtie2Output(alignmentLines, reverseCorrection=-1):
# The key update to version 2 of the AlignSequencesToGenome code is in here. 
# We take a look at the second field of the alignment line

# Note, I've added the option to add 1 bp to the readAlignmentCoord if the read aligns to the 
# reverse strand. I can easily take this out in the future if I redesign the primers to avoid 
# binding to the ITR at the opposite end of the transposon. 

    import re
    import pdb

    xsRe = re.compile('XS:i:')

    lowScoreCount = 0
    noAlignCount = 0
    perfectAlignmentCount = 0
    multipleAlignmentCount = 0
    xsCount = 0
    xsIndices = []
    lowScoreIndices = []
    noAlignIndices = []
    readsWithStrangeFlagsSums = 0

    genomeAlignmentDict = {}

    i = 0
    while i < len(alignmentLines):
        line = alignmentLines[i]
        if len(line) > 0:
            if line[0] != '@':
                if xsRe.search(line) != None:
                    xsCount += 1
                    xsIndices.append(i)
                    multipleAlignmentsFound = True
                else:
                    multipleAlignmentsFound = False


                lineData = line.split('\t')
                readID = lineData[0]
                flagsSum = int(lineData[1])
                readAlignmentCoord = int(lineData[3])
                alignmentQuality = int(lineData[4])
                alignmentSequenceLength = len(lineData[9])

# 				pdb.set_trace()


                if flagsSum == 0 or flagsSum == 16:
                    # Flags sum should be zero as we are analyzing non paired reads, if
                    # read aligns to forward strand.
                    if flagsSum == 0:
                        readAlignmentCoord = readAlignmentCoord
                        strangeFlagsSum = False
                    # Flags sum should be 16 as we are analyzing non paired reads, if
                    # read aligns to reverse strand.
                    elif flagsSum == 16:
                        readAlignmentCoord = readAlignmentCoord + alignmentSequenceLength \
                        + reverseCorrection
                        # Note, I've added the option to add 1 bp here to normalize the result with
                        # the Sanger read.
                        strangeFlagsSum = False
                else:
                    # This really shouldn't be needed, but in case it is, we have a way to
                    # report a strange flags sum.
                    readAlignmentCoord = readAlignmentCoord
                    readsWithStrangeFlagsSums += 1
                    strangeFlagsSum = True

                if alignmentQuality <= 1:
                    lowScoreCount += 1
                    lowScoreIndices.append(i)

                if readAlignmentCoord == 0:
                    noAlignIndices.append(i)
                    noAlignCount += 1
                    alignmentFound = False
                else:
                    alignmentFound = True
                    if multipleAlignmentsFound == False:
                        perfectAlignmentCount += 1
                    elif multipleAlignmentsFound == True:
                        multipleAlignmentCount += 1


                genomeAlignmentDict[readID] = Alignment(readID, readAlignmentCoord, \
                alignmentQuality, alignmentFound, multipleAlignmentsFound, strangeFlagsSum)

        i += 1
    return [genomeAlignmentDict, perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, \
    noAlignCount, readsWithStrangeFlagsSums]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def OutputGenomeAlignmentSummary(outputFileName, genomeAlignmentDict):
# Updated to include strange flags sum
    keys = genomeAlignmentDict.keys()

    outputHandle = open(outputFileName, 'w')

    for key in keys:

        readAlignmentCoord = str(genomeAlignmentDict[key].readAlignmentCoord)
        alignmentQuality = str(genomeAlignmentDict[key].alignmentQuality)
        alignmentFound = str(genomeAlignmentDict[key].alignmentFound)
        multipleAlignmentsFound = str(genomeAlignmentDict[key].multipleAlignmentsFound)
        strangeFlagsSum = str(genomeAlignmentDict[key].strangeFlagsSum)

        outputString = str(key) + ',' + readAlignmentCoord + ',' + alignmentQuality + ',' \
        + alignmentFound + ',' + multipleAlignmentsFound + ',' + strangeFlagsSum + '\n'

        outputHandle.write(outputString)

    outputHandle.close()
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateLogStringForIndividualMfaFile(startTime, endTime, currentFileNumber, totalFiles, \
perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, noAlignCount, strangeFlagsSumCount):
# Updated to include strange flags sum
    duration = (endTime - startTime).total_seconds()

    summary = "Reads with Perfect Genome Alignments: " + str(perfectAlignmentCount) + '\n'
    summary += "Reads with Multiple Genome Alignnments: " + str(multipleAlignmentCount) + '\n'
    summary += "Reads with Low Alignment Scores: " + str(lowScoreCount) + '\n'
    summary += "Reads with No Alignment Scores: " + str(noAlignCount) + '\n'
    summary += "Reads with Strange Flags Sum: " + str(strangeFlagsSumCount) + '\n'

    lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'

    timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
    timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  + '\n'

    logStr = summary + lastOperationDurationStr + timeRemainingStr

    return logStr
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateLogFileSummary(totalPerfectAlignmentCount, totalMultipleAlignmentCount, \
totalLowScoreCount, totalNoAlignCount, totalStrangeFlagsSumCount, initialStartTime, completeEndTime):
# Updated to include strange flags sum
    jobDuration = (completeEndTime - initialStartTime).total_seconds()

    summary = '# ----------------------------------------------------------------------------------------------- #\n'
    summary = "Total Reads with Perfect Genome Alignments: " + str(totalPerfectAlignmentCount) + '\n'
    summary += "Total Reads with Multiple Genome Alignnments: " + str(totalMultipleAlignmentCount) + '\n'
    summary += "Total Reads with Low Alignment Scores: " + str(totalLowScoreCount) + '\n'
    summary += "Total Reads with No Alignment Scores: " + str(totalNoAlignCount) + '\n'
    summary += "Total Reads with Strange Flags Sum: " + str(totalStrangeFlagsSumCount) + '\n'

    summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
    summary += '# ----------------------------------------------------------------------------------------------- #\n'

    return summary
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CountEntriesInSummaryFiles(alignmentFiles, calculateMaxReadIDLength=False):
    import gc

    alignmentFileEntries = 0

    maxReadIDLength = 10

    for file in alignmentFiles:
# 		print("Processing file: " + file)

        with open(file) as infile:
            for line in infile:
                alignmentFileEntries += 1
                if calculateMaxReadIDLength == True:
                    lineArray = line.strip().split(',')
                    readIDLength = len(lineArray[0])
                    if readIDLength > maxReadIDLength:
                        maxReadIDLength = readIDLength

        gc.collect()

    if calculateMaxReadIDLength:
        return [alignmentFileEntries, maxReadIDLength]
    else:
        return alignmentFileEntries
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CountEntriesInSummaryFile3(alignmentFile, calculateMaxReadIDLength=False):
    import gc

    alignmentFileEntries = 0

    maxReadIDLength = 10

    with open(alignmentFile) as infile:
        for line in infile:
            alignmentFileEntries += 1
            if calculateMaxReadIDLength == True:
                lineArray = line.strip().split(',')
                readIDLength = len(lineArray[0])
                if readIDLength > maxReadIDLength:
                    maxReadIDLength = readIDLength

    gc.collect()

    if calculateMaxReadIDLength:
        return [alignmentFileEntries, maxReadIDLength]
    else:
        return alignmentFileEntries
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CountEntriesInAllSummaryFiles(genomeAlignmentFiles, himarRecognitionFiles, indexAlignmentFiles):

    print("Counting entries in genome alignment files and calculating maximum length of read ids")
    [genomeAlignmentEntries, maxReadIDLength] = CountEntriesInSummaryFiles(genomeAlignmentFiles, \
    calculateMaxReadIDLength=True)
    print("Counting entries in himar alignment files")
    himarRecognitionEntries = CountEntriesInSummaryFiles(himarRecognitionFiles)
    print("Counting entries in index alignment files")
    indexAlignmentEntries = CountEntriesInSummaryFiles(indexAlignmentFiles)

    return [genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CountEntriesInAllSummaryFiles3(genomeAlignmentFile, himarRecognitionFile, indexAlignmentFile):

    print("Counting entries in genome alignment file and calculating maximum length of read ids")
    [genomeAlignmentEntries, maxReadIDLength] = CountEntriesInSummaryFile3(genomeAlignmentFile, \
    calculateMaxReadIDLength=True)
    print("Counting entries in himar alignment file")
    himarRecognitionEntries = CountEntriesInSummaryFile3(himarRecognitionFile)
    print("Counting entries in index alignment file")
    indexAlignmentEntries = CountEntriesInSummaryFile3(indexAlignmentFile)

    return [genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength]
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def ParseGenomeAlignmentFiles(readIDFieldCode, genomeAlignmentFiles, genomeAlignmentEntries):
# Updated to import flag that indicates if bowtie alignment flags sum was wrong
    import numpy
    import ast
    from numpy import int32

    genomeArray = numpy.zeros(genomeAlignmentEntries, \
    dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'alignmentFound', \
    'multipleAlignmentsFound', 'himarRecognized', 'index', 'strangeFlagsSum'], \
    'formats':[readIDFieldCode, int32, int32, int32, int32, int32, int32, int32]})

    i = 0
    for file in genomeAlignmentFiles:
        print("Processing file: " + file)

        with open(file) as infile:
            for line in infile:
                lineArray = line.strip().split(',')

                genomeArray[i]['readID'] = lineArray[0]
                genomeArray[i]['readAlignmentCoord'] = int(lineArray[1])
                genomeArray[i]['alignmentQuality'] = int(lineArray[2])

                if ast.literal_eval(lineArray[3]):
                    genomeArray[i]['alignmentFound'] = 1
                else:
                    genomeArray[i]['alignmentFound'] = 0

                if ast.literal_eval(lineArray[4]):
                    genomeArray[i]['multipleAlignmentsFound'] = 1
                else:
                    genomeArray[i]['multipleAlignmentsFound'] = 0

                if ast.literal_eval(lineArray[5]):
                    genomeArray[i]['strangeFlagsSum'] = 1
                else:
                    genomeArray[i]['strangeFlagsSum'] = 0

                i += 1

    return genomeArray
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ParseGenomeAlignmentFile3(readIDFieldCode, genomeAlignmentFile, genomeAlignmentEntries):
# Updated to import flag that indicates if bowtie alignment flags sum was wrong
# Also updated to deal with with only one alignment file at once
    import numpy
    import ast
    from numpy import int32

    genomeArray = numpy.zeros(genomeAlignmentEntries, \
    dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'alignmentFound', \
    'multipleAlignmentsFound', 'himarRecognized', 'index', 'strangeFlagsSum'], \
    'formats':[readIDFieldCode, int32, int32, int32, int32, int32, int32, int32]})

    i = 0
    print("Processing file: " + genomeAlignmentFile)

    with open(genomeAlignmentFile) as infile:
        for line in infile:
            lineArray = line.strip().split(',')

            genomeArray[i]['readID'] = lineArray[0]
            genomeArray[i]['readAlignmentCoord'] = int(lineArray[1])
            genomeArray[i]['alignmentQuality'] = int(lineArray[2])

            if ast.literal_eval(lineArray[3]):
                genomeArray[i]['alignmentFound'] = 1
            else:
                genomeArray[i]['alignmentFound'] = 0

            if ast.literal_eval(lineArray[4]):
                genomeArray[i]['multipleAlignmentsFound'] = 1
            else:
                genomeArray[i]['multipleAlignmentsFound'] = 0

            if ast.literal_eval(lineArray[5]):
                genomeArray[i]['strangeFlagsSum'] = 1
            else:
                genomeArray[i]['strangeFlagsSum'] = 0

            i += 1

    return genomeArray
# ------------------------------------------------------------------------------------------------ #











# ------------------------------------------------------------------------------------------------ #
def ParseIndexAlignmentFiles(readIDFieldCode, indexAlignmentFiles, indexAlignmentEntries):

    import numpy
    import ast
    from numpy import int32

    indexArray = numpy.zeros(indexAlignmentEntries, dtype={'names':['readID', 'index'], \
    'formats':[readIDFieldCode, int32]})


    i = 0
    for file in indexAlignmentFiles:
        print("Processing file: " + file)

        with open(file) as infile:
            for line in infile:
                lineArray = line.strip().split(',')
                indexArray[i]['readID'] = lineArray[0]

                # If the index is not recognized, set the index to zero
                if ast.literal_eval(lineArray[2]) == False:
                    index = 0
                else:
                    indexArray[i]['index'] = lineArray[1]

                i += 1

    return indexArray
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
# Updated to deal with only one file at once
def ParseIndexAlignmentFile3(readIDFieldCode, indexAlignmentFile, indexAlignmentEntries):

    import numpy
    import ast
    from numpy import int32
    import pdb

    indexArray = numpy.zeros(indexAlignmentEntries, dtype={'names':['readID', 'index'], \
    'formats':[readIDFieldCode, int32]})


    i = 0
    print("Processing file: " + indexAlignmentFile)

    infile = open(indexAlignmentFile, 'r')
    data = infile.readlines()

    for line in data:

        lineArray = line.strip().split(',')
        try:
            indexArray[i]['readID'] = lineArray[0]
        except:
            pdb.set_trace()

        # If the index is not recognized, set the index to zero
        if ast.literal_eval(lineArray[2]) == False:
            index = 0
        else:
            indexArray[i]['index'] = lineArray[1]

        i += 1

    return indexArray
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def ParseHimarAlignmentFiles(readIDFieldCode, himarRecognitionFiles, himarRecognitionEntries):
    import numpy
    import ast
    from numpy import int32

    himarArray = numpy.zeros(himarRecognitionEntries, dtype={'names':['readID', 'himarRecognized'], \
    'formats':[readIDFieldCode, int32]})


    i = 0
    for file in himarRecognitionFiles:
        print("Processing file: " + file)

        with open(file) as infile:
            for line in infile:
                lineArray = line.strip().split(',')
                readID = lineArray[0]
                himarRecognized = ast.literal_eval(lineArray[1])

                himarArray[i]['readID'] = readID

                if himarRecognized:
                    himarArray[i]['himarRecognized'] = 1
                else:
                    himarArray[i]['himarRecognized'] = 0

                i += 1
    return himarArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Updated to deal with only one file at once
def ParseHimarAlignmentFile3(readIDFieldCode, himarRecognitionFile, himarRecognitionEntries):
    import numpy
    import ast
    from numpy import int32

    himarArray = numpy.zeros(himarRecognitionEntries, dtype={'names':['readID', 'himarRecognized'], \
    'formats':[readIDFieldCode, int32]})

    print("Processing file: " + himarRecognitionFile)

    i = 0
    with open(himarRecognitionFile) as infile:
        for line in infile:
            lineArray = line.strip().split(',')
            readID = lineArray[0]
            himarRecognized = ast.literal_eval(lineArray[1])

            himarArray[i]['readID'] = readID

            if himarRecognized:
                himarArray[i]['himarRecognized'] = 1
            else:
                himarArray[i]['himarRecognized'] = 0

            i += 1

    return himarArray
# ------------------------------------------------------------------------------------------------ #














# ------------------------------------------------------------------------------------------------ #
def UpdateGenomeArrayWithHimarAndIndex(genomeArray, indexArray, himarArray):

    from pdb import set_trace
    import numpy

    sortedGenomeArray = numpy.sort(genomeArray, order='readID')
    sortedHimarArray = numpy.sort(himarArray, order='readID')
    sortedIndexArray = numpy.sort(indexArray, order='readID')

# 	set_trace()

    i = 0
    while i < len(sortedGenomeArray):
        if sortedGenomeArray[i]['readID'] == sortedHimarArray[i]['readID'] \
        and sortedGenomeArray[i]['readID'] == sortedIndexArray[i]['readID']:

            sortedGenomeArray[i]['himarRecognized'] = sortedHimarArray[i]['himarRecognized']
            sortedGenomeArray[i]['index'] = sortedIndexArray[i]['index']
        else:
            print("Read ids don't match")
        i += 1

    return sortedGenomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WriteCompiledGenomeArray(outputFile, sortedGenomeArray):
# Updated to print out strangeFlagsSum flag

    fileHandle = open(outputFile, 'w')

    i = 0
    while i < len(sortedGenomeArray):

        outputStr = str(sortedGenomeArray[i]['readID'].tostring().decode('utf-8')) + ','
        outputStr += str(sortedGenomeArray[i]['readAlignmentCoord']) + ','
        outputStr += str(sortedGenomeArray[i]['alignmentQuality']) + ','
        outputStr += str(sortedGenomeArray[i]['alignmentFound']) + ','
        outputStr += str(sortedGenomeArray[i]['multipleAlignmentsFound']) + ','
        outputStr += str(sortedGenomeArray[i]['himarRecognized']) + ','
        outputStr += str(sortedGenomeArray[i]['index']) + ','
        outputStr += str(sortedGenomeArray[i]['strangeFlagsSum']) + '\n'

        fileHandle.write(outputStr)

        i += 1
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateSudokuReadTaxonomy(genomeArray, outputLog, operationTimingInterval=100000):
# Updated to include strange flags flag in read taxonomy	

    import datetime

    logStr = gSeparatorString
    logStr += 'Generating Sudoku Read Taxonomy'

    print(logStr)

    UpdateLogFileData(outputLog, logStr)

    taxDict = {'totalReads':0, 'totalReadsWithHimar':0,	'totalReadsWithIndex':0,\
    'totalReadsWithNoStrangeFlags':0, 'totalReadsWithNoHimar':0, 'totalReadsWithNoIndex':0,\
    'totalReadsWithStrangeFlags':0, 'totalReadsWithPerfectGenome':0, \
    'totalReadsWithMultipleAlignment':0, 'totalReadsWithNoAlignment':0,	\
    'totalReadsWithLowScoreAlignment':0, \
    'totalReadsWithIndexAndHimar':0, 'totalReadsWithIndexAndHimarAndPerfectGenome':0, \
    'totalReadsWithIndexAndHimarAndMultipleAlignment':0, \
    'totalReadsWithIndexAndHimarButNoAlignment':0,	\
    'totalReadsWithIndexAndHimarAndLowScoreAlignment':0,\
    'totalReadsWithIndexAndHimarAndPerfectGenomeAndNoStrangeFlags':0,
    'totalReadsWithIndexAndNoHimar':0, 'totalReadsWithIndexNoHimarAndPerfectGenome':0,\
    'totalReadsWithIndexNoHimarAndMultipleAlignment':0,	\
    'totalReadsWithIndexNoHimarButNoAlignment':0, \
    'totalReadsWithIndexNoHimarAndLowScoreAlignment':0, 'totalReadsWithNoIndexAndHimar':0, \
    'totalReadsWithNoIndexAndHimarAndPerfectGenome':0,\
    'totalReadsWithNoIndexAndHimarAndMultipleAlignment':0,\
    'totalReadsWithNoIndexAndHimarButNoAlignment':0,\
    'totalReadsWithNoIndexAndHimarAndLowScoreAlignment':0, 'totalReadsWithNoIndexAndNoHimar':0,\
    'totalReadsWithNoIndexNoHimarAndPerfectGenome':0,\
    'totalReadsWithNoIndexNoHimarAndMultipleAlignment':0,\
    'totalReadsWithNoIndexNoHimarButNoAlignment':0,	\
    'totalReadsWithNoIndexNoHimarAndLowScoreAlignment':0}


    i = 0

    startTime = datetime.datetime.now()

    while i < len(genomeArray):

        sudokuRead = genomeArray[i]

        taxDict['totalReads'] += 1

        if sudokuRead['himarRecognized'] == 1:
            taxDict['totalReadsWithHimar'] += 1
        else:
            taxDict['totalReadsWithNoHimar'] += 1

        if sudokuRead['index'] != 0:
            taxDict['totalReadsWithIndex'] += 1
        else:
            taxDict['totalReadsWithNoIndex'] += 1

        if sudokuRead['strangeFlagsSum'] == 0:
            taxDict['totalReadsWithNoStrangeFlags'] += 1
        else:
            taxDict['totalReadsWithStrangeFlags'] += 1

        if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
            taxDict['totalReadsWithPerfectGenome'] += 1
        elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] and sudokuRead['alignmentQuality'] > 1:
            taxDict['totalReadsWithMultipleAlignment'] += 1
        elif sudokuRead['alignmentFound'] == 0:
            taxDict['totalReadsWithNoAlignment'] += 1
        elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
            taxDict['totalReadsWithLowScoreAlignment'] += 1


        if sudokuRead['himarRecognized'] == 1 and sudokuRead['index'] != 0:
            taxDict['totalReadsWithIndexAndHimar'] += 1
            if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithIndexAndHimarAndPerfectGenome'] += 1
                if sudokuRead['strangeFlagsSum'] == 0:
                    taxDict['totalReadsWithIndexAndHimarAndPerfectGenomeAndNoStrangeFlags'] += 1
            elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithIndexAndHimarAndMultipleAlignment'] += 1
            elif sudokuRead['alignmentFound'] == 0:
                taxDict['totalReadsWithIndexAndHimarButNoAlignment'] += 1
            elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
                taxDict['totalReadsWithIndexAndHimarAndLowScoreAlignment'] += 1


        elif sudokuRead['himarRecognized'] == 0 and sudokuRead['index'] != 0:
            taxDict['totalReadsWithIndexAndNoHimar'] += 1
            if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithIndexNoHimarAndPerfectGenome'] += 1
            elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithIndexNoHimarAndMultipleAlignment'] += 1
            elif sudokuRead['alignmentFound'] == 0:
                taxDict['totalReadsWithIndexNoHimarButNoAlignment'] += 1
            elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
                taxDict['totalReadsWithIndexNoHimarAndLowScoreAlignment'] += 1


        elif sudokuRead['himarRecognized'] and sudokuRead['index'] == 0:
            taxDict['totalReadsWithNoIndexAndHimar'] += 1
            if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithNoIndexAndHimarAndPerfectGenome'] += 1
            elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithNoIndexAndHimarAndMultipleAlignment'] += 1
            elif sudokuRead['alignmentFound'] == 0:
                taxDict['totalReadsWithNoIndexAndHimarButNoAlignment'] += 1
            elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
                taxDict['totalReadsWithNoIndexAndHimarAndLowScoreAlignment'] += 1


        elif sudokuRead['himarRecognized'] == 0 and sudokuRead['index'] == 0:
            taxDict['totalReadsWithNoIndexAndNoHimar'] += 1
            if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithNoIndexNoHimarAndPerfectGenome'] += 1
            elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
                taxDict['totalReadsWithNoIndexNoHimarAndMultipleAlignment'] += 1
            elif sudokuRead['alignmentFound'] == 0:
                taxDict['totalReadsWithNoIndexNoHimarButNoAlignment'] += 1
            elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
                taxDict['totalReadsWithNoIndexNoHimarAndLowScoreAlignment'] += 1

        if i % operationTimingInterval == 0:
            endTime = datetime.datetime.now()
            duration = (endTime - startTime).total_seconds()
            timeRemainingEstimate = (duration/operationTimingInterval)*(len(genomeArray)-i)/60

            logStr = "Time for last " + str(operationTimingInterval) + " operations: " \
            + str(duration) + " seconds\n"
            logStr += "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes\n"

            print(logStr)
            UpdateLogFileData(outputLog, logStr)

            startTime = datetime.datetime.now()

        i += 1

    return taxDict
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def OutputReadTaxonomy(taxonomyDict, outputLog):

    outputLogHandle  = open(outputLog, 'w')

    taxonomyDictKeys = taxonomyDict.keys()

    outputStr = ''

    for key in taxonomyDictKeys:
        outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'

    print(outputStr)

    UpdateLogFileData(outputLog, outputStr)

    return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportGenomeCompilationFile(genomeCompilationFile, includesb=False):
    import gc
    import numpy
    import pdb

    readIDRe = re.compile('\d+:\d+:\d+:\d+')

    # Determine length of genomeCompilationFile

    fHandle = open(genomeCompilationFile)
    compilationEntries = 0
    maxReadIDLength = 15

    with open(genomeCompilationFile) as infile:
        for line in infile:
            compilationEntries += 1
            lineArray = line.strip().split(',')
            readIDLength = len(lineArray[0])
            if readIDLength > maxReadIDLength:
                maxReadIDLength = readIDLength

    gc.collect()
# 	pdb.set_trace()

    readIDFieldCode = 'a' + str(maxReadIDLength+2)

    genomeArray = numpy.zeros(compilationEntries, \
    dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'alignmentFound', \
    'multipleAlignmentsFound', 'himarRecognized', 'index', 'strangeFlagsSum'], \
    'formats':[readIDFieldCode, 'i32', 'i16', 'i8', 'i8', 'i8', 'i8', 'i8']})

    i = 0
    with open(genomeCompilationFile) as infile:
        for line in infile:
            lineArray = line.strip().split(',')
            if includesb:
                genomeArray[i]['readID'] = readIDRe.search(lineArray[0]).group()
            else:
                genomeArray[i]['readID'] = lineArray[0]
            genomeArray[i]['readAlignmentCoord'] = int(lineArray[1])
            genomeArray[i]['alignmentQuality'] = int(lineArray[2])
            genomeArray[i]['alignmentFound'] = int(lineArray[3])
            genomeArray[i]['multipleAlignmentsFound'] = int(lineArray[4])
            genomeArray[i]['himarRecognized'] = int(lineArray[5])
            genomeArray[i]['index'] = int(lineArray[6])
            genomeArray[i]['strangeFlagsSum'] = int(lineArray[7])


            i += 1

    return genomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GenerateBarcodeLookupTable(barcodeFile, barcodePrefix='', barcodePostfix=''):
# Modified with a kludge to include extra Illumina index sequence that gets read in with the
# Sudoku barcode as part of the index. Adds the prefix and postfix to the reverse complement

    barcodeFileHandle = open(barcodeFile, 'r')
    barcodeFileData = barcodeFileHandle.readlines()
    barcodeFileHandle.close()

    barcodeLookupTable = {}

    for line in barcodeFileData:
        if line[0] != '#':
            lineData = line.strip().split(',')
            revCompl = barcodePrefix + lineData[2] + barcodePostfix
            barcodeNumber = lineData[3]
            barcodeLookupTable[revCompl] = int(barcodeNumber)

    barcodes = barcodeLookupTable.keys()
    return [barcodeLookupTable, barcodes]
# ------------------------------------------------------------------------------------------------ #

# Creates dictionary of every reasonable mismatch and what barcode they refer to
# ------------------------------------------------------------------------------------------------ #
def GenerateBarcodeCorrectionMapping(barcodeList, maxMismatches):
    import itertools

    correctionMapping = {}
    noGo = {}
    idxs = list(itertools.combinations(range(len(barcodeList[0])), maxMismatches))
    baseCombos = list(itertools.product('ACTGN', repeat=maxMismatches))
    for barcode in barcodeList:
        for idx in idxs:
            for baseCombo in baseCombos:
                temp = list(barcode)
                for i,base in zip(idx, baseCombo):
                    temp[i] = base
                newBarcode = ''.join(temp)
                dist = hamming_distance(newBarcode, barcode)
                if newBarcode not in correctionMapping or correctionMapping[newBarcode][0] > dist:
                    correctionMapping[newBarcode] = (dist, barcode)
                    if newBarcode in noGo:
                        del noGo[newBarcode]
                elif newBarcode in correctionMapping and correctionMapping[newBarcode][0] == dist and correctionMapping[newBarcode][1] != barcode:
                    noGo[newBarcode] = dist


    # delete all mismatched barcodes that are equally close to at least two real barcodes
    for barcode in noGo:
        del correctionMapping[barcode]
    return correctionMapping

# ------------------------------------------------------------------------------------------------ #

# depreciated, will be removed later
# ------------------------------------------------------------------------------------------------ #
def GenerateNestedDecisionDict(barcodeList):

    nestedDecisionDict = {}
    j = 0

    while j < len(barcodeList):
        firstLetter = barcodeList[j][0]
        secondLetter = barcodeList[j][1]

        nestedDecisionDictKeys = nestedDecisionDict.keys()

        if firstLetter not in nestedDecisionDictKeys:
            nestedDecisionDict[firstLetter] = {secondLetter:[barcodeList[j]]}
        else:
            secondLevelKeys = nestedDecisionDict[firstLetter]

            if secondLetter not in secondLevelKeys:
                nestedDecisionDict[firstLetter][secondLetter] = [barcodeList[j]]
            else:
                nestedDecisionDict[firstLetter][secondLetter].append(barcodeList[j])

        j += 1

    return nestedDecisionDict
# ------------------------------------------------------------------------------------------------ #

# depreciated will be removed later
# ------------------------------------------------------------------------------------------------ #
def GenerateMultipleEntryDict(barcodeList, indexSequenceSampleLength):
    # Generate the search dictionary for sequence indexes
    barcodeLen = len(barcodeList[0])

    multipleEntryPointDecisionDict = {}
    i = 0
    while i < barcodeLen - indexSequenceSampleLength:
        j = 0
        while j < len(barcodeList):
            letters = barcodeList[j][i:i+indexSequenceSampleLength+1]
            code = str(i) + '_' + letters

            multipleEntryPointDecisionDictKeys = multipleEntryPointDecisionDict.keys()

            if code not in multipleEntryPointDecisionDictKeys:
                multipleEntryPointDecisionDict[code] = [barcodeList[j]]
            else:
                multipleEntryPointDecisionDict[code].append(barcodeList[j])

            j += 1
        i += 1

    multipleEntryPointDecisionDictLengths = {}
    multipleEntryPointDecisionDictKeys = multipleEntryPointDecisionDict.keys()

    for key in multipleEntryPointDecisionDictKeys:
        multipleEntryPointDecisionDictLengths[key] = len(multipleEntryPointDecisionDict[key])

    return multipleEntryPointDecisionDict
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def ParseGenomeArrayToPoolFiles(genomeArray, barcodeLookupTable, poolFileBaseDir, poolFilePrefix):	
    import re
    import gc
    import pdb

    # Generate pool file handles
    poolFileHandleDict = {}
    poolFiles = []

    indices = barcodeLookupTable.values()

    for index in indices:
        poolFileName = poolFileBaseDir + poolFilePrefix + str(index) + '.csv'
        poolFileHandleDict[index] = open(poolFileName, 'w')
        poolFiles.append(poolFileName)

    # Run through genome array and parse good reads into the pool files
    i = 0
    while i < len(genomeArray):
        read = genomeArray[i]
        readID = genomeArray[i]['readID'].tostring().decode('utf-8')
        readAlignmentCoord = genomeArray[i]['readAlignmentCoord']
        alignmentQuality = genomeArray[i]['alignmentQuality']
        alignmentFound = genomeArray[i]['alignmentFound']
        multipleAlignmentsFound= genomeArray[i]['multipleAlignmentsFound']
        himarRecognized = genomeArray[i]['himarRecognized']
        index = genomeArray[i]['index']
        strangeFlagsSum = genomeArray[i]['strangeFlagsSum']

        if alignmentQuality > 1 and multipleAlignmentsFound == 0 and himarRecognized == 1 \
        and index != 0 and strangeFlagsSum == 0:
            outputStr = str(readID) + ',' + str(readAlignmentCoord) + '\n'
            try:
                poolFileHandleDict[index].write(outputStr)
            except:
                pdb.set_trace()
        i += 1

    # Close up all of the pool file handles
    for index in indices:
        poolFileHandleDict[index].close()

    return poolFiles
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def GetIndexCodeFromPoolFileName(poolFileName, poolFileNamePrefix='pool_'):

    import re

    import pdb

    poolNameRe = re.compile(poolFileNamePrefix + '(\d+)')

# 	pdb.set_trace()
    [fileName, fileExtension] = ProcessFileNameExtension(poolFileName)

    poolNameMatch = poolNameRe.match(fileName)

    if poolNameRe.match(fileName) != None:
        indexCode = poolNameMatch.group(1)
    else:
        print("Error. Pool file name does not match regular expression: " + poolFileName)

    return indexCode
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 1: Functions for Checking Himar Content and Trimming Sequence Reads
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ----------------------------------------------------------------------------------------------- #
def Mismatch(Search,n):
    import itertools
    SearchL = list(Search)
    List     = [] # hold output
    # print list of indices to replace with '$'
    idxs = itertools.combinations(range(len(SearchL)),n)
    # for each combination `idx` in idxs, replace str[idx] with '$':
    for idx in idxs:
        str = SearchL[:] # make a copy
        for i in idx:
            str[i]='[ACTGN]'
        List.append( ''.join(str) ) # convert back to string
    return List
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateHimarMismatchRegexes(himarSequence, mismatchNumber=2):
    import re

    himarMismatchRegexList = []
    fullSequenceHimarMismatchRegexList = []

    searchList = Mismatch(himarSequence,mismatchNumber)
    searchListRegexes = []
    for searchTerm in searchList:
        himarMismatchRegexList.append(re.compile(r'' + searchTerm))
        fullSequenceHimarMismatchRegexList.append(re.compile(r'([NACTG]{0,8})' + r'(' + searchTerm + r')' + r'([ACTGN]+)'))

    return [himarMismatchRegexList, fullSequenceHimarMismatchRegexList]

# ----------------------------------------------------------------------------------------------- #

# Allows for at least 2 mismatches (other code confirms that there are no more than 2)
#------------------------------------------------------------------------------------------------ #
def GenerateHimarMismatchRegexes2(himarSequence):
    import re

    himarMismatchRegexList = []
    fullSequenceHimarMismatchRegexList = []
    half_len = len(himarSequence) // 2

    firstTermVars = Mismatch(himarSequence[1:half_len + 1], 1)

    searchTerm = himarSequence[0] + '[NACTG]{' +  str(half_len) + '}' + himarSequence[half_len + 1:]
    himarMismatchRegexList.append(re.compile(r'' + searchTerm))
    fullSequenceHimarMismatchRegexList.append(
        re.compile(r'([NACTG]{0,8})' + r'(' + searchTerm + r')' + r'([ACTGN]+)'))

    for term in firstTermVars:
        searchTerm = '[NACTG]' + term + '[NACTG]{' + str(half_len - 1) + '}'
        himarMismatchRegexList.append(re.compile(r'' + searchTerm))
        fullSequenceHimarMismatchRegexList.append(
            re.compile(r'([NACTG]{0,8})' + r'(' + searchTerm + r')' + r'([ACTGN]+)'))

    return [himarMismatchRegexList, fullSequenceHimarMismatchRegexList]

#------------------------------------------------------------------------------------------------ #

# ----------------------------------------------------------------------------------------------- #
def CheckAndTrimHimarReadSequence(readSeq, idString, himarRe, himarMismatchRegexList, \
fullSequenceHimarMismatchRegexList,himarSequence):

    from pdb import set_trace

    himarReMatch = himarRe.match(readSeq)


    if himarReMatch != None:
        himarAssignment = True
        trimmedSequence = himarReMatch.group(3)
        assignmentType = 'perfect'
    else:
        # if we don't find a perfect match for the himar sequence, check for mismatch
        himarMismatchMatched = False

        trimmedSequence = ''

        i = 0
        while i < len(himarMismatchRegexList):
            regex = himarMismatchRegexList[i]
            regexRes = regex.findall(readSeq)
            if regexRes:
                missed = 0
                for match in regexRes:
                    missed = 0
                    for letter in range(len(himarSequence)):
                        if himarSequence[letter] != match[letter]:
                            missed += 1
                    if missed <= 2:
                        break
                if missed > 2:
                    i += 1
                    continue
                himarMismatchMatched = True
                fullSequenceRe = fullSequenceHimarMismatchRegexList[i].search(readSeq)
                if fullSequenceRe == None:
# 					set_trace()
                    trimmedSequence = ''
                else:
                    trimmedSequence = fullSequenceRe.group(3)
                break
            i += 1

        if trimmedSequence == '':
            himarAssignment = False
            assignmentType = 'unmatchable'
            trimmedSequence = readSeq
        elif himarMismatchMatched:
            himarAssignment = True
            assignmentType = 'matchable'
        else:
            himarAssignment = False
            assignmentType = 'unmatchable'
            trimmedSequence = readSeq


    return [himarAssignment, assignmentType, trimmedSequence]
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def ParseHimarReadFile(file, himarMismatchRegexList, fullSequenceHimarMismatchRegexList, himarRe, himarSequence):
    # Import the fastq files and summarize

    from pdb import set_trace
    import re

    # We'll assume the following format for the read id

    # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    # EAS139 is the unique instrument name
    # 136 is the run id
    # FC706VJ is the flowcell id
    # 2	is the flowcell lane
    # 2104 is the tile number within the flowcell lane
    # 15343	is the 'x'-coordinate of the cluster within the tile
    # 197393 is the 'y'-coordinate of the cluster within the tile
    # 1	is the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
    # Y	is Y if the read is filtered, N otherwise
    # 18i is 0 when none of the control bits are on, otherwise it is an even number
    # ATCACG index sequence

    # '@JLK5VL1:704:H3M72BCXX:2:1101:1209:2029 2:N:0:'
    identifierLineRe = re.compile('@(\w+):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)')

    # In these dicts, the key is the read number, and the associated value is true or false for
    # himar recognition

    himarRecognitionDict = {}
    perfectHimarReads = 0
    imperfectButMatchableHimarReads = 0
    unMatchableHimarReads = 0
    himarTrimmedSequenceDict = {}


    fileHandle = open(file, 'r')
    fileData = fileHandle.readlines()

    i = 0
    while i < len(fileData):
        line = fileData[i]

        identifierMatch = identifierLineRe.match(line)

        if identifierMatch != None:

            lane = identifierMatch.group(4)
            tileNumber = identifierMatch.group(5)
            xCoord = identifierMatch.group(6)
            yCoord = identifierMatch.group(7)
            idString = lane + ':' + tileNumber + ':' + xCoord + ':' + yCoord
            readSeq = fileData[i+1].strip()


            [himarAssignment, assignmentType, trimmedSequence] = \
            CheckAndTrimHimarReadSequence(readSeq, idString, himarRe, himarMismatchRegexList, \
            fullSequenceHimarMismatchRegexList, himarSequence)

            himarTrimmedSequenceDict[idString] = trimmedSequence
            himarRecognitionDict[idString] = himarAssignment

            if assignmentType == 'perfect':
                perfectHimarReads += 1
            elif assignmentType == 'matchable':
                imperfectButMatchableHimarReads += 1
            elif assignmentType == 'unmatchable':
                unMatchableHimarReads += 1

        i += 1

    return [himarRecognitionDict, himarTrimmedSequenceDict, perfectHimarReads, \
    imperfectButMatchableHimarReads, unMatchableHimarReads]
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def WriteHimarSequenceSummaries(himarRecognitionOutputFileName, himarRecognitionDict):

    summaryFileHandle = open(himarRecognitionOutputFileName, 'w')

    himarRecognitionDictKeys = himarRecognitionDict.keys()

    for key in himarRecognitionDictKeys:
        summaryFileHandle.write(str(key) + ',' + str(himarRecognitionDict[key]) + '\n')

    summaryFileHandle.close()
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def WriteTrimmedSequences(trimmedFileName, himarTrimmedSequenceDict):

    mfastaHandle = open(trimmedFileName, 'w')

    himarSequenceDictKeys = himarTrimmedSequenceDict.keys()

    for key in himarSequenceDictKeys:
        mfastaHandle.write('>' + str(key) + '\n')
        mfastaHandle.write(str(himarTrimmedSequenceDict[key]) + '\n')

    mfastaHandle.close()
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateLogStringForIndividualFastQFileForHimarRecognition(startTime, endTime, \
currentFileNumber, totalFiles, perfectIndexReads, mismatchedButMatchableIndexReads, \
unMatchableIndexReads ):

    duration = (endTime - startTime).total_seconds()

    indexReadSummary = "Perfect Himar Reads: " + str(perfectIndexReads) + '\n'
    indexReadSummary += "Mismatched But Readable Himar Reads: " + \
    str(mismatchedButMatchableIndexReads) + '\n'
    indexReadSummary += "Unmatchable Himar Reads: " + str(unMatchableIndexReads) + '\n'

    lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'

    timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
    timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  \
    + '\n'

    logStr = indexReadSummary + lastOperationDurationStr + timeRemainingStr

    return logStr
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateLogFileSummaryForHimarRecognition(totalPerfectHimarReads, \
totalImperfectButMatchableHimarReads, totalUnMatchableHimarReads, initialStartTime, \
completeEndTime):

    jobDuration = (completeEndTime - initialStartTime).total_seconds()

    summary = gSeparatorString
    summary += "Total Perfect Himar Reads: " + str(totalPerfectHimarReads) + '\n'
    summary += "Total Mismatched But Readable Himar Reads: " + \
    str(totalImperfectButMatchableHimarReads) + '\n'
    summary += "Total Unmatchable Himar Reads: " + str(totalUnMatchableHimarReads) + '\n'

    summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
    summary += gSeparatorString

    return summary
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
def CheckHimarContentAndTrimSequences(fastqFiles, \
himarSequence, himarRecognitionFilePrefix, himarRecognitionBaseDir, \
trimmedSequencesFilePrefix, trimmedSequencesBaseDir, outputLog):


# Import the fastq files, check the himar and summarize the himar content and trim the sequence
# to only the genomic material

# In these dicts, the key is the read number, and the associated value is true or false for himar
# recognition

# Note, as of the current revision this is a bit of kludgy hold over from earlier versions of the 
# code. This function writes out the summary files and returns the file names, only for the 
# next function in the analysis chain to read them back in. 

    import datetime
    import re

    # Generate himar regular expressions and mismatch regular expressions
    himarRe = re.compile(r'([NACTG]{0,8})' + r'(' + himarSequence + r')' + r'([ACTGN]+)')

    [himarMismatchRegexList, fullSequenceHimarMismatchRegexList] = \
    GenerateHimarMismatchRegexes2(himarSequence)


    logFileHandle = open(outputLog, 'a')
    initialStartTime = datetime.datetime.now()

    j = 1
    totalPerfectHimarReads = 0
    totalImperfectButMatchableHimarReads = 0
    totalUnMatchableHimarReads = 0

    himarRecognitionFiles = []
    trimmedMfaFiles = []

    for file in fastqFiles:

        startTime = datetime.datetime.now()

        # Parse the individual himar read fastq file
        processStr = 'Processing ' + str(j) + ' of ' + str(len(fastqFiles)) + ': ' + file + '\n'
        print(processStr)

        [himarRecognitionDict, himarTrimmedSequenceDict, perfectHimarReads, \
        imperfectButMatchableHimarReads, unMatchableHimarReads] \
        = ParseHimarReadFile(file, himarMismatchRegexList, fullSequenceHimarMismatchRegexList, \
        himarRe, himarSequence)


        totalPerfectHimarReads += perfectHimarReads
        totalImperfectButMatchableHimarReads += imperfectButMatchableHimarReads
        totalUnMatchableHimarReads += unMatchableHimarReads

        # Write out the himar sequence summaries
        himarRecognitionOutputFileName = himarRecognitionBaseDir + '/' \
        + himarRecognitionFilePrefix + '_' + str(j) + '.csv'

        WriteHimarSequenceSummaries(himarRecognitionOutputFileName, himarRecognitionDict)

        himarRecognitionFiles.append(himarRecognitionOutputFileName)

        # Write out the trimmed sequences
        trimmedFileName = trimmedSequencesBaseDir + '/' \
        + trimmedSequencesFilePrefix + '_' + str(j) + '.mfa'

        WriteTrimmedSequences(trimmedFileName, himarTrimmedSequenceDict)

        trimmedMfaFiles.append(trimmedFileName)

        endTime = datetime.datetime.now()

        summaryLogStr = GenerateLogStringForIndividualFastQFileForHimarRecognition(startTime, \
        endTime, j, len(fastqFiles), perfectHimarReads, imperfectButMatchableHimarReads, \
        unMatchableHimarReads)

        print(summaryLogStr)

        # Write data to log file
        UpdateLogFileData(outputLog, processStr+summaryLogStr)


        j += 1

    completeEndTime = datetime.datetime.now()

    summaryStr = GenerateLogFileSummaryForHimarRecognition(totalPerfectHimarReads, \
    totalImperfectButMatchableHimarReads, \
    totalUnMatchableHimarReads, initialStartTime, completeEndTime)

    print(summaryStr)

    UpdateLogFileData(outputLog, summaryStr)

    return [himarRecognitionFiles, trimmedMfaFiles]
# ------------------------------------------------------------------------------------------------ #





####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 2: Functions for Index Summarization
# ------------------------------------------------------------------------------------------------ #
####################################################################################################


# ----------------------------------------------------------------------------------------------- #
def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    
    import numpy
    import pdb
    
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
        
        
    zipseq = zip(s1, s2)
        
    sum = 0
    
    for el1, el2 in zipseq:
        if el1 != el2:
            sum += 1
        else:
            sum += 0
    
    return sum
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def CountNs(sequence):

    nCount = 0
    i = 0
    while i < len(sequence):
        if sequence[i] == 'N':
            nCount += 1
        i += 1

    return nCount
# ----------------------------------------------------------------------------------------------- #


# depreciated will be removed later
# ----------------------------------------------------------------------------------------------- #
def AssignBarcodeReadToPool(indexSeq, barcodeLookupTable, barcodeList, nestedDecisionDict, \
multipleEntryPointDecisionDict, indexIndexArray, allowedNs, mismatchesAllowed, \
sampleLength, maxTrials, minBestCandidatesForNonPerfectMatch):

    # New version of function to parse a single index read using a nested decision dictionary

    from scipy import unique
    from pdb import set_trace

    nCount = CountNs(indexSeq)

    if nCount >= allowedNs:
        # If we have more than the allowed number of Ns in the sequence, reject it immediately.
        barcodeAssignment = [indexSeq, False]
        assignmentType = 'unmatchable'
    else:
        # Try to find the index sequence, or something like it using the nested decision dict.
        # Note: we might not be able to, that's why we have the exception
        firstLetter = indexSeq[0]
        secondLetter = indexSeq[1]

        try:
            candidateSequences = nestedDecisionDict[firstLetter][secondLetter]

            if indexSeq in candidateSequences:
                # If we can find can find the first two letters of the index sequence with
                # the nested decision dict, it's likely to be perfectly matched.
                assignmentType = 'perfect'
                barcodeAssignment = [barcodeLookupTable[indexSeq], True]
            else:
                # But, it might not be.
                [barcodeAssignment, assignmentType] = \
                MatchNonPerfectIndexSequence(indexSeq, indexIndexArray, barcodeList, \
                barcodeLookupTable, multipleEntryPointDecisionDict, maxTrials, sampleLength,
                mismatchesAllowed, minBestCandidatesForNonPerfectMatch)

        except:
            # In the event that we can't find the index sequence with the nested decision dict
            # go straight to the non-perfect matching algorithm.
            [barcodeAssignment, assignmentType] = \
            MatchNonPerfectIndexSequence(indexSeq, indexIndexArray, barcodeList, \
            barcodeLookupTable, multipleEntryPointDecisionDict, maxTrials, sampleLength,
            mismatchesAllowed, minBestCandidatesForNonPerfectMatch)

    return [barcodeAssignment, assignmentType]

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def AssignBarcodeRealToPool2(indexSeq, barcodeLookupTable, barcodeCorrectionMapping, allowedNs):

    import itertools
    nCount = CountNs(indexSeq)

    if nCount >= allowedNs:
        # If we have more than the allowed number of Ns in the sequence, reject it immediately.
        barcodeAssignment = [indexSeq, False]
        assignmentType = 'unmatchable'
    else:
        # Try to find the index sequence, or something like it using the nested decision dict.
        # Note: we might not be able to, that's why we have the exception
        if indexSeq in barcodeCorrectionMapping:
            if barcodeCorrectionMapping[indexSeq][0] == 0:
                assignmentType = 'perfect'
                barcodeAssignment = [barcodeLookupTable[indexSeq], True]
            else:
                assignmentType = 'matchable'
                barcodeAssignment = [barcodeLookupTable[barcodeCorrectionMapping[indexSeq][1]], True]
        else:
            barcodeAssignment = [indexSeq, False]
            assignmentType = 'unmatchable'

    return [barcodeAssignment, assignmentType]

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def MatchNonPerfectIndexSequence(indexSeq, indexIndexArray, barcodeList, barcodeLookupTable, \
multipleEntryPointDecisionDict, \
maxTrials, sampleLength, mismatchesAllowed, minBestCandidates):

    from numpy.random import random, choice
    import pdb
    from collections import Counter

    trial = 1
    bestCandidateMatches = []


    while trial <= maxTrials:
        # Select maxTrials (usually 3) random positions in the index sequence and pull out a
        # short sample. Use this to interrogate the multipleEntryPointDecisionDict to get some
        # candidate barcode match sequences.
        # Sometimes the index will be mismatched, so we pull out less than maxTrials sets of
        # barcodes.
        indexChoice = choice(indexIndexArray)
        code = str(indexChoice) + '_' + indexSeq[indexChoice:indexChoice+sampleLength]

        try:
            candidateSequences = multipleEntryPointDecisionDict[code]
            hammingDistances = []

            for candidateSequence in candidateSequences:
                # For each of the candidate sequences that we pull out, calculate its Hamming
                # distance from the index sequence.
                hammingDistances.append([hamming_distance(candidateSequence, indexSeq), \
                candidateSequence])

            # Pick the candidate sequence with the smallest Hamming distance from the index
            sortedHammingDistances = sorted(hammingDistances, key = lambda x: x[0])

            pdb.set_trace()

            bestCandidateMatches.append(sortedHammingDistances[0][1])
            trial += 1
        except:
            trial += 1
            continue

    # Now, we vote on the best candidate matches that we get.
    tempHammingDistances = []
    if len(bestCandidateMatches) <= minBestCandidates:
        # If we don't get sufficient candidates through the multipleEntryPointDecisionDict,
        # brute force it
        # and try to match the sequence against every possible barcode


        for barcode in barcodeList:
            tempHammingDistances.append([hamming_distance(indexSeq, barcode), barcode])

    else:
        for candidate in bestCandidateMatches:
            tempHammingDistances.append([hamming_distance(indexSeq, candidate), candidate])

    sortedTempHammingDistances = sorted(tempHammingDistances, key = lambda x: x[0])

    lowestHammingDistance = sortedTempHammingDistances[0][0]
    bestMatchingBarcode = sortedTempHammingDistances[0][1]

    if lowestHammingDistance >= mismatchesAllowed:
        matchType = 'unmatchable'
        barcodeAssignment = [indexSeq, False]
    else:
        matchType = 'matchable'
        barcodeAssignment = [barcodeLookupTable[bestMatchingBarcode], True]


    return [barcodeAssignment, matchType]
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateLogFileSummaryForIndexSummary(totalPerfectIndexReads, \
totalMismatchedButMatchableIndexReads, \
totalUnMatchableIndexReads, initialStartTime, completeEndTime):

    jobDuration = (completeEndTime - initialStartTime).total_seconds()

    summary = gSeparatorString
    summary += "Total Perfect Index Reads: " + str(totalPerfectIndexReads) + '\n'
    summary += "Total Mismatched But Readable Index Reads: " + \
    str(totalMismatchedButMatchableIndexReads) + '\n'
    summary += "Total Unmatchable Index Reads: " + str(totalUnMatchableIndexReads) + '\n'

    summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
    summary += gSeparatorString

    return summary
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateLogStringForIndividualFastQFileForIndexSummary(startTime, endTime, currentFileNumber, \
totalFiles, perfectIndexReads, mismatchedButMatchableIndexReads, unMatchableIndexReads ):

    duration = (endTime - startTime).total_seconds()

    indexReadSummary = "Perfect Index Reads: " + str(perfectIndexReads) + '\n'
    indexReadSummary += "Mismatched But Readable Index Reads: " + \
    str(mismatchedButMatchableIndexReads) + '\n'
    indexReadSummary += "Unmatchable Index Reads: " + str(unMatchableIndexReads) + '\n'

    lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'

    timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
    timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  + '\n'

    logStr = indexReadSummary + lastOperationDurationStr + timeRemainingStr

    return logStr
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def ParseIndexReadFile(file, barcodeLookupTable, barcodeCorrectionMapping, allowedNsInIndexSeq):

    # New version of function to parse a single index read using a nested decision dictionary

    import re
    from numpy import arange

    # We'll assume the following format for the read id

    # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    # EAS139 is the unique instrument name
    # 136 is the run id
    # FC706VJ is the flowcell id
    # 2	is the flowcell lane
    # 2104 is the tile number within the flowcell lane
    # 15343	is the 'x'-coordinate of the cluster within the tile
    # 197393 is the 'y'-coordinate of the cluster within the tile
    # 1	is the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
    # Y	is Y if the read is filtered, N otherwise
    # 18i is 0 when none of the control bits are on, otherwise it is an even number
    # ATCACG index sequence

    # '@JLK5VL1:704:H3M72BCXX:2:1101:1209:2029 2:N:0:'
    identifierLineRe = re.compile('@(\w+):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)')

    barcodeRecognitionDict = {}

    fileHandle = open(file, 'r')
    fileData = fileHandle.readlines()

    perfectIndexReads = 0
    mismatchedButMatchableIndexReads = 0
    unMatchableIndexReads = 0

    # Import index sequences

    i = 0

    while i < len(fileData):
        line = fileData[i]

        identifierMatch = identifierLineRe.match(line)

        if identifierMatch != None:

            lane = identifierMatch.group(4)
            tileNumber = identifierMatch.group(5)
            xCoord = identifierMatch.group(6)
            yCoord = identifierMatch.group(7)
            idString = lane + ':' + tileNumber + ':' + xCoord + ':' + yCoord

            indexSeq = fileData[i+1].strip()

            if i == 0:
                indexIndexArray = arange(0, len(indexSeq) - 2, 1, dtype=int)

            [barcodeAssignment, assignmentType] = AssignBarcodeRealToPool2(indexSeq, barcodeLookupTable,
                                                                           barcodeCorrectionMapping, allowedNsInIndexSeq)


            barcodeRecognitionDict[idString] = barcodeAssignment

            if assignmentType == 'perfect':
                perfectIndexReads += 1
            elif assignmentType == 'matchable':
                mismatchedButMatchableIndexReads += 1
            elif assignmentType == 'unmatchable':
                unMatchableIndexReads += 1

        i += 1

    return [barcodeRecognitionDict, perfectIndexReads, mismatchedButMatchableIndexReads, \
    unMatchableIndexReads]
# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
def OutputBarcodeRecognitionTable(barcodeRecognitionDict, barcodeRecognitionOutputFileName):

    barcodeRecognitionKeys = list(barcodeRecognitionDict.keys())

    fileHandle = open(barcodeRecognitionOutputFileName, 'w')

    for key in barcodeRecognitionKeys:

        outputStr = str(key) + ',' + str(barcodeRecognitionDict[key][0]) \
        + ',' + str(barcodeRecognitionDict[key][1]) + '\n'

        fileHandle.write(outputStr)

    fileHandle.close()
# ----------------------------------------------------------------------------------------------- #




# ------------------------------------------------------------------------------------------------ #
def SummarizeIndexReads(indexFastqFiles, barcodeFile, indexSummaryFilePrefix, indexSummaryBaseDir, \
outputLog, allowedNsInIndexSeq=3, indexMismatchesAllowed=3, barcodePrefix='', barcodePostfix='A'):
# Import the fastq files and summarize

# In these dicts, the key is the read number, and the associated value is the code for the pool
# barcode

    import datetime


    # Generate the barcode lookup table and the mismatch regexes
    [barcodeLookupTable, barcodes] = GenerateBarcodeLookupTable(barcodeFile, \
    barcodePrefix=barcodePrefix, barcodePostfix=barcodePostfix)
    barcodeList = list(barcodes)
    barcodeCorrectionMapping = GenerateBarcodeCorrectionMapping(barcodeList, indexMismatchesAllowed)

    j = 1
    totalPerfectIndexReads = 0
    totalMismatchedButMatchableIndexReads = 0
    totalUnMatchableIndexReads = 0

    indexAlignmentFiles = []

    logFileHandle = open(outputLog, 'a')

    initialStartTime = datetime.datetime.now()

    for file in indexFastqFiles:

        startTime = datetime.datetime.now()

        # Parse the individual index read fastq file
        processStr = 'Processing ' + str(j) + ' of ' + str(len(indexFastqFiles)) + ': ' \
        + file + '\n'

        print(processStr)

        [barcodeRecognitionDict, perfectIndexReads, mismatchedButMatchableIndexReads, \
        unMatchableIndexReads] = \
        ParseIndexReadFile(file, barcodeLookupTable, barcodeCorrectionMapping, allowedNsInIndexSeq)

        totalPerfectIndexReads += perfectIndexReads
        totalMismatchedButMatchableIndexReads += mismatchedButMatchableIndexReads
        totalUnMatchableIndexReads += unMatchableIndexReads

        # Output the barcode recognition table
        barcodeRecognitionOutputFileName = indexSummaryBaseDir + '/' \
        + indexSummaryFilePrefix + '_' + str(j) + '.csv'


        OutputBarcodeRecognitionTable(barcodeRecognitionDict, barcodeRecognitionOutputFileName)

        indexAlignmentFiles.append(barcodeRecognitionOutputFileName)

        endTime = datetime.datetime.now()

        summaryLogStr = GenerateLogStringForIndividualFastQFileForIndexSummary(startTime, \
        endTime, j, len(indexFastqFiles), perfectIndexReads, mismatchedButMatchableIndexReads, \
        unMatchableIndexReads)

        print(summaryLogStr)

        # Write data to log file
        UpdateLogFileData(outputLog, processStr+summaryLogStr)

        j += 1


    completeEndTime = datetime.datetime.now()

    summaryStr = GenerateLogFileSummaryForIndexSummary(totalPerfectIndexReads, \
    totalMismatchedButMatchableIndexReads, totalUnMatchableIndexReads, initialStartTime, \
    completeEndTime)

    print(summaryStr)

    UpdateLogFileData(outputLog, summaryStr)


    return indexAlignmentFiles
# ------------------------------------------------------------------------------------------------ #




####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 3: Functions for Genome Alignment
# ------------------------------------------------------------------------------------------------ #
####################################################################################################


# ------------------------------------------------------------------------------------------------ #
def AlignReadsToGenome(mfaFiles, outputLog, \
genomeAlignmentBaseDir, genomeAlignmentFilePrefix, referenceIndexPrefix, \
bowtieMode='--end-to-end'):
# Code to align illumina reads to reference genome and return alignment point

    import gc
    import subprocess
    import datetime
    from pdb import set_trace
    import re

    totalLowScoreCount = 0
    totalNoAlignCount = 0
    totalMultipleAlignmentCount = 0
    totalPerfectAlignmentCount = 0
    totalStrangeFlagsSumCount = 0

    genomeAlignmentFiles = []

    xsRe = re.compile('XS:i:')

    logFileHandle = open(outputLog, 'a')
    initialStartTime = datetime.datetime.now()

    j = 1

    for file in mfaFiles:

        startTime = datetime.datetime.now()

        processStr = 'Processing ' + str(j) + ' of ' + str(len(mfaFiles)) + ': ' + file + '\n'
        print(processStr)

# 		fullFilePath = mfaBaseDir + '/' + file

# 		set_trace()

        alignmentData = \
        subprocess.check_output(['bowtie2', '-x', referenceIndexPrefix, '-f', file, \
        bowtieMode, '--threads', '4'])

        alignmentStr = alignmentData.decode("utf-8")
        alignmentLines = alignmentStr.split('\n')

        [genomeAlignmentDict, perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, \
        noAlignCount, strangeFlagsSumCount] = ParseBowtie2Output(alignmentLines)

        # Update the read summary counts
        totalLowScoreCount += lowScoreCount
        totalNoAlignCount += noAlignCount
        totalMultipleAlignmentCount += multipleAlignmentCount
        totalPerfectAlignmentCount += perfectAlignmentCount
        totalStrangeFlagsSumCount += strangeFlagsSumCount

        # Write out the genome alignment files
        outputFileName = \
        genomeAlignmentBaseDir + '/' + genomeAlignmentFilePrefix + '_' + str(j) + '.csv'
        OutputGenomeAlignmentSummary(outputFileName, genomeAlignmentDict)

        genomeAlignmentFiles.append(outputFileName)

        # End timing
        endTime = datetime.datetime.now()

        summaryLogStr = GenerateLogStringForIndividualMfaFile(startTime, endTime, j, \
        len(mfaFiles), perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, noAlignCount, \
        strangeFlagsSumCount)

        print(summaryLogStr)

        # Write data to log file
        UpdateLogFileData(outputLog, processStr+summaryLogStr)

        j += 1

    completeEndTime = datetime.datetime.now()

    summaryStr = GenerateLogFileSummary(totalPerfectAlignmentCount, totalMultipleAlignmentCount, \
    totalLowScoreCount, totalNoAlignCount, totalStrangeFlagsSumCount, initialStartTime, \
    completeEndTime)

    print(summaryStr)

    UpdateLogFileData(outputLog, summaryStr)

    gc.collect()

    return genomeAlignmentFiles
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def StatisticsForCoordinateGroupsWithMoreThanOneMember(\
groupedPoolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools, \
dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0,\
hiAddEntThresh_r=None, hiAddEntThresh_c=None, hiAddEntThresh_pr=None, hiAddEntThresh_pc=None):

# For the groups with more than one member, count how many can be assigned to library addresses
# Note to be include, address entries must be greater than loAddEntThresh, not equal to or greater
# than.

    groupsWithMoreThanOneEntry = 0
    groupsWithMoreThanOneEntryNotToIncludeInSummedTable = 0
    groupsWithMoreThanOneEntryToIncludeInSummedTable = 0
    groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress = 0

    if hiAddEntThresh_r == None:
        hiAddEntThresh_r = len(rowPools)
    if hiAddEntThresh_c == None:
        hiAddEntThresh_c = len(colPools)
    if hiAddEntThresh_pr == None:
        hiAddEntThresh_pr = len(prPools)
    if hiAddEntThresh_pc == None:
        hiAddEntThresh_pc = len(pcPools)

    i = 0

    while i < len(groupedPoolPresenceTable):

        if len(groupedPoolPresenceTable[i]) > 1:

            groupsWithMoreThanOneEntry += 1

            [summedPoolPresenceTableGroup, includeInSummedTable] = \
            SumPoolPresenceTableGroup(groupedPoolPresenceTable[i], \
            dTypeDictForPoolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools)

            if includeInSummedTable == True:

                groupsWithMoreThanOneEntryToIncludeInSummedTable += 1

                addresses_r = FindAddressCoords(summedPoolPresenceTableGroup, rowPools, threshold=threshold)
                addresses_c = FindAddressCoords(summedPoolPresenceTableGroup, colPools, threshold=threshold)
                addresses_pr = FindAddressCoords(summedPoolPresenceTableGroup, prPools, threshold=threshold)
                addresses_pc = FindAddressCoords(summedPoolPresenceTableGroup, pcPools, threshold=threshold)

                nEntries_r = len(addresses_r)
                nEntries_c = len(addresses_c)
                nEntries_pr = len(addresses_pr)
                nEntries_pc = len(addresses_pc)

                if (loAddEntThresh < nEntries_r < hiAddEntThresh_r) \
                and (loAddEntThresh < nEntries_c < hiAddEntThresh_c) \
                and (loAddEntThresh < nEntries_pr < hiAddEntThresh_pr) \
                and (loAddEntThresh < nEntries_pc < hiAddEntThresh_pc):
                    groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress += 1
            else:
                groupsWithMoreThanOneEntryNotToIncludeInSummedTable += 1



        i += 1

    return [groupsWithMoreThanOneEntry, groupsWithMoreThanOneEntryToIncludeInSummedTable, \
    groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress, \
    groupsWithMoreThanOneEntryNotToIncludeInSummedTable]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def StatisticsForCoordinateGroupsWithOnlyOneMember(\
groupedPoolPresenceTable, rowPools, colPools, prPools, pcPools, \
dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0,\
hiAddEntThresh_r=None, hiAddEntThresh_c=None, hiAddEntThresh_pr=None, hiAddEntThresh_pc=None):


    groupsWithOnlyOneEntry = 0
    groupsWithOnlyOneEntryWithValidLibraryAddress = 0

    i = 0

    if hiAddEntThresh_r == None:
        hiAddEntThresh_r = len(rowPools)
    if hiAddEntThresh_c == None:
        hiAddEntThresh_c = len(colPools)
    if hiAddEntThresh_pr == None:
        hiAddEntThresh_pr = len(prPools)
    if hiAddEntThresh_pc == None:
        hiAddEntThresh_pc = len(pcPools)


    while i < len(groupedPoolPresenceTable):

        if len(groupedPoolPresenceTable[i]) <= 1:
            groupsWithOnlyOneEntry += 1

            addresses_r = FindAddressCoords(groupedPoolPresenceTable[i][0], rowPools, threshold=threshold)
            addresses_c = FindAddressCoords(groupedPoolPresenceTable[i][0], colPools, threshold=threshold)
            addresses_pr = FindAddressCoords(groupedPoolPresenceTable[i][0], prPools, threshold=threshold)
            addresses_pc = FindAddressCoords(groupedPoolPresenceTable[i][0], pcPools, threshold=threshold)

            nEntries_r = len(addresses_r)
            nEntries_c = len(addresses_c)
            nEntries_pr = len(addresses_pr)
            nEntries_pc = len(addresses_pc)

            if (loAddEntThresh < nEntries_r < hiAddEntThresh_r) \
            and (loAddEntThresh < nEntries_c < hiAddEntThresh_c) \
            and (loAddEntThresh < nEntries_pr < hiAddEntThresh_pr) \
            and (loAddEntThresh < nEntries_pc < hiAddEntThresh_pc):
                groupsWithOnlyOneEntryWithValidLibraryAddress += 1



        i += 1

    return [groupsWithOnlyOneEntry, groupsWithOnlyOneEntryWithValidLibraryAddress]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def AnalyzeGroupedPoolPresenceTable(groupedPoolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0):

    [groupsWithOnlyOneEntry, groupsWithOnlyOneEntryWithValidLibraryAddress] = \
    StatisticsForCoordinateGroupsWithOnlyOneMember(groupedPoolPresenceTable, \
    rowPools, colPools, prPools, pcPools, dTypeDictForPoolPresenceTable, \
    threshold=threshold, loAddEntThresh=loAddEntThresh)


    [groupsWithMoreThanOneEntry, groupsWithMoreThanOneEntryToIncludeInSummedTable, \
    groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress, \
    groupsWithMoreThanOneEntryNotToIncludeInSummedTable] = \
    StatisticsForCoordinateGroupsWithMoreThanOneMember(groupedPoolPresenceTable, \
    rowPools, colPools, prPools, pcPools, controlPools, dTypeDictForPoolPresenceTable, \
    threshold=threshold, loAddEntThresh=loAddEntThresh)


    outputStr = "groupsWithOnlyOneEntry: " + str(groupsWithOnlyOneEntry) + '\n'

    outputStr += "groupsWithOnlyOneEntryWithValidLibraryAddress: " \
    + str(groupsWithOnlyOneEntryWithValidLibraryAddress) + '\n'

    outputStr += "groupsWithMoreThanOneEntry" + str(groupsWithMoreThanOneEntry) + '\n'

    outputStr += "groupsWithMoreThanOneEntryNotToIncludeInSummedTable: " \
    + str(groupsWithMoreThanOneEntryNotToIncludeInSummedTable) + '\n'

    outputStr += "groupsWithMoreThanOneEntryToIncludeInSummedTable: " \
    + str(groupsWithMoreThanOneEntryToIncludeInSummedTable) + '\n'

    outputStr += "groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress: " \
    + str(groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress)

    groupedPoolPresenceAnalysisDict = {}
    groupedPoolPresenceAnalysisDict['groupsWithOnlyOneEntry'] = groupsWithOnlyOneEntry
    groupedPoolPresenceAnalysisDict['groupsWithOnlyOneEntryWithValidLibraryAddress'] = groupsWithOnlyOneEntryWithValidLibraryAddress
    groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntry: '] = groupsWithMoreThanOneEntry
    groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryNotToIncludeInSummedTable'] = groupsWithMoreThanOneEntryNotToIncludeInSummedTable
    groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryToIncludeInSummedTable'] = groupsWithMoreThanOneEntryToIncludeInSummedTable
    groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress'] = groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress



    return [outputStr, groupedPoolPresenceAnalysisDict]

# ------------------------------------------------------------------------------------------------ #






####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 4: Functions for Alignment Compilation
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def CompileAlignmentsWithIndexAndHimarSummaries(genomeAlignmentFiles, indexAlignmentFiles, \
himarRecognitionFiles, outputLog, genomeArrayFileName):

    import gc
    import pdb

    outputStr = gSeparatorString
    outputStr += 'Compiling Genome Alignments \n'
    UpdateLogFileData(outputLog, outputStr)

    # Count lines in the summary files
    [genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength] \
    = CountEntriesInAllSummaryFiles(genomeAlignmentFiles, indexAlignmentFiles, \
    himarRecognitionFiles)

    genomeHimarIndexLinesSummaryEntry = "Genome alignment entries: " + str(genomeAlignmentEntries) \
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Index alignment entries: " + str(indexAlignmentEntries)\
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Himar alignment entries: " + str(himarRecognitionEntries)\
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Maximum read id length: " + str(maxReadIDLength) + '\n'
    print(genomeHimarIndexLinesSummaryEntry)


    UpdateLogFileData(outputLog, genomeHimarIndexLinesSummaryEntry)


    # Allocate arrays for import of data
    readIDFieldCode = 'a' + str(maxReadIDLength+2)

    genomeArray = ParseGenomeAlignmentFiles(readIDFieldCode, genomeAlignmentFiles, \
    genomeAlignmentEntries)
    indexArray = ParseIndexAlignmentFiles(readIDFieldCode, indexAlignmentFiles, \
    indexAlignmentEntries)
    himarArray =  ParseHimarAlignmentFiles(readIDFieldCode, himarRecognitionFiles, \
    himarRecognitionEntries)

    genomeArray = UpdateGenomeArrayWithHimarAndIndex(genomeArray, indexArray, himarArray)

    # Write compiled genome array
    WriteCompiledGenomeArray(genomeArrayFileName, genomeArray)


    # Generate the read taxonomy
    taxonomyDict = GenerateSudokuReadTaxonomy(genomeArray, outputLog)
    OutputReadTaxonomy(taxonomyDict, outputLog)

    UpdateLogFileData(outputLog, gSeparatorString)
    print(gSeparatorString)

    gc.collect()

    return genomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# This variant of CompileAlignmentsWithIndexAndHimarSummaries works on only one set of alignment 
# files at once, rather than trying to deal with all of them at the same time. 
def CompileAlignmentsWithIndexAndHimarSummaries3(genomeAlignmentFile, indexAlignmentFile, \
himarRecognitionFile, outputLog, genomeArrayFileName=None):

    import gc
    import pdb

    outputStr = gSeparatorString
    outputStr += 'Compiling Genome Alignments \n'
    UpdateLogFileData(outputLog, outputStr)

    # Count lines in the summary files
    [genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength] \
    = CountEntriesInAllSummaryFiles3(genomeAlignmentFile, indexAlignmentFile, \
    himarRecognitionFile)

    genomeHimarIndexLinesSummaryEntry = "Genome alignment entries: " + str(genomeAlignmentEntries) \
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Index alignment entries: " + str(indexAlignmentEntries)\
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Himar alignment entries: " + str(himarRecognitionEntries)\
    + '\n'
    genomeHimarIndexLinesSummaryEntry += "Maximum read id length: " + str(maxReadIDLength) + '\n'
    print(genomeHimarIndexLinesSummaryEntry)


    UpdateLogFileData(outputLog, genomeHimarIndexLinesSummaryEntry)


    # Allocate arrays for import of data
    readIDFieldCode = 'a' + str(maxReadIDLength+2)

    genomeArray = ParseGenomeAlignmentFile3(readIDFieldCode, genomeAlignmentFile, \
    genomeAlignmentEntries)
    indexArray = ParseIndexAlignmentFile3(readIDFieldCode, indexAlignmentFile, \
    indexAlignmentEntries)
    himarArray =  ParseHimarAlignmentFile3(readIDFieldCode, himarRecognitionFile, \
    himarRecognitionEntries)

    genomeArray = UpdateGenomeArrayWithHimarAndIndex(genomeArray, indexArray, himarArray)

    # Write compiled genome array
    if genomeArrayFileName != None:
        WriteCompiledGenomeArray(genomeArrayFileName, genomeArray)

    # Generate the read taxonomy
    taxonomyDict = GenerateSudokuReadTaxonomy(genomeArray, outputLog)
    OutputReadTaxonomy(taxonomyDict, outputLog)

    UpdateLogFileData(outputLog, gSeparatorString)
    print(gSeparatorString)

    gc.collect()

    return genomeArray
# ------------------------------------------------------------------------------------------------ #

