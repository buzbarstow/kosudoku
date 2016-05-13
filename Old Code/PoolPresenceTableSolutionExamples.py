#!/sw/bin/python3.4

# ----------------------------------------------------------------------------------------------- #
# PoolPresenceTableSolve.py
# Created by Buz Barstow 2015-6-14
# Last updated by Buz Barstow 2015-8-26

# All code necessary to build pool presence table from Sudoku reads, except for code that performs
# summary of himar section of reads, and summarizes index reads. 
# Now uses sudokuutils version 7.
# ----------------------------------------------------------------------------------------------- #


from PoolPresenceTableSolutionExamplesUtils import SummarizePoolPresenceTableLine, \
GeneratePoolPresenceTableLineSummary, FindFeatureAssociatedWithInsertionCoord

from sudokuutils13.input import get_input

from sudokuutils13.grid import ImportSudokuGridLayout, WriteSudokuGridSummaryTable, \
AnnotateSudokuGridWithFeautureNames, FlattenSudokuGridLookupDict

from sudokuutils13.pool import ImportPoolPresenceTable, SolvePoolPresenceTable, \
CalculateLibraryAddressLocatabilityForPoolPresenceTableLine

from sudokuutils13.xml import ExportSudokuGridToXML, ImportSudokuGridFromXML

from sudokuutils13.utils import ensure_dir

from sudokuutils13.feature import ImportFeatureArrayFromGenBank

import sys
import pdb
import os

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'controlPools', 'poolPresenceTableFileName', \
'sudokuGridLayout', 'outputLog', 'maxGapForCoordGrouping', \
'readCountThreshold', 'voigtScoreThreshold', 'maxTotalCoords', 'maxSinglePoolCoordNumber', \
'sudokuGridXMLFile', 'sudokuGridSummaryFileName', 'shewyGenomeGenBankFileName']

inputParameterValues = get_input(['',\
'/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/' \
+ 'PoolPresenceTableSolutionExamples.inp'], \
inputParameters)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
controlPools = inputParameterValues['controlPools'].split(',')
poolPresenceTableFileName = inputParameterValues['poolPresenceTableFileName']

outputLog = inputParameterValues['outputLog']
sudokuGridLayout = inputParameterValues['sudokuGridLayout']
maxGapForCoordGrouping = int(inputParameterValues['maxGapForCoordGrouping'])
readCountThreshold = int(inputParameterValues['readCountThreshold'])
voigtScoreThreshold = float(inputParameterValues['voigtScoreThreshold'])
maxTotalCoords = int(inputParameterValues['maxTotalCoords'])
maxSinglePoolCoordNumber = int(inputParameterValues['maxSinglePoolCoordNumber'])
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
outputLogDir = os.path.dirname(outputLog)
	
ensure_dir(outputLogDir + '/')
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, colPools)

sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools + controlPools

poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
logReadNumberRatioHistogramFitDict = {\
'nr2nc': array([  4.40167575e+03,  -6.14416490e-01,   2.31441504e-01, 6.10917305e-01]), \
'npr2npc': array([  4.37639039e+03,  -1.93463028e-01,   2.49560416e-01, 5.69134210e-01]), \
'nc2npc': array([  4.26554263e+03,  -1.28866329e+00,   1.02679161e-01, 8.49570121e-01]), \
'nr2npr': array([  4.29226485e+03,  -1.77483793e+00,   1.67029253e-01, 8.50852254e-01]), \
'nc2npr': array([  4.33491451e+03,  -9.78105660e-01,   2.26399445e-01, 7.66314686e-01]), \
'nr2npc': array([  4.24814935e+03,  -2.11653643e+00,   7.17920853e-02, 9.20773485e-01]) }
          

logReadNumberRatioHistogramIntegralDict = {\
'nr2nc': 4335.5755955204122, 'npr2npc': 4305.8760588073219, 'nc2npc': 4236.5746634922598, \
'nr2npr': 4244.0194718705397, 'nc2npr': 4270.6449466057529, 'nr2npc': 4227.2659786258046 }
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)
geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(feature)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Select the pool presence table lines that I like
# These ones are for the mtr operons

useAmbiguousLinesOnly = False

poolPresenceTableShort = []

mtrStartPos = 1856172
mtrEndPos = 1869219

ccmStartPos = 260541
ccmEndPos = 268692

hypStartPos = 2185668
hypEndPos = 2198009

menStartPos = 4768718
menEndPos = 4772146

otherCoords = [2377953, 5092640, 3495167, 2691769, 2811903, 30488, 4184454, 1637582, 2931139, \
3554361, 2292305, 3828101, 4606388, 4606388, 2756427, 2756427, 966582, 3140568, 3140568, 641373, \
641373, 2903204, 2903204, 1610209, 291316, 1697670, 1697670, 1813478, 1813478, 932464, 932464, \
2184691, 1898865, 1898865, 1193232, 1193232, 2128603, 2128603, 3910492, 3910492, 3005695, 3005695, \
749857, 2480469, 2480469, 3835250, 2102435, 2102435, 1703729, 1703729, 829647, 2085296, 332908, \
332908, 52044, 52044, 718445, 718445, 940467, 2930899, 2930899, 2416708, 2416708, 4917161, 4917161,\
3964600, 3964600, 2741866, 4069459, 4069459, 1642454, 1642454, 3202030, 3202030, 1637272, 2079928,\
1656686, 1168256, 847292, 3011186, 1096597, 549236, 293625, 1583953, 3850567, 2653801, 4454863, \
1799131, 4067947, 38388, 3414719, 3503675, 2673533, 3920296, 2723477, 4767186, 1004803, 3063049, \
1762609, 1696532, 4052160, 2393110, 103338, 1896287, 4369798, 2190757, 2766492, 1735720, 3014637, \
2597384, 3096631, 141334, 1442097, 3300969, 3447165, 1487683, 3344725, 2041378, 3423455, 912465, \
3652497, 106914, 1696307, 500688, 1912454, 623863, 3474541, 383681, 656270, 1902621, 4120453, \
178622, 3783872, 2208013, 509664, 1648184, 2369719, 543186, 1543625, 1656734]


i = 0
while i < len(poolPresenceTable):
	readAlignmentCoord = poolPresenceTable[i]['readAlignmentCoord']
	
	if (mtrStartPos <= readAlignmentCoord <= mtrEndPos) \
	or (ccmStartPos <= readAlignmentCoord <= ccmEndPos) \
	or (hypStartPos <= readAlignmentCoord <= hypEndPos) \
	or (menStartPos <= readAlignmentCoord <= menEndPos) \
	or readAlignmentCoord in otherCoords:
		poolPresenceTableShort.append(poolPresenceTable[i])
	i += 1
	

i = 0
summaryArray = []
while i < len(poolPresenceTableShort):
	
	poolPresenceTableLine = poolPresenceTableShort[i]
	
	summary = SummarizePoolPresenceTableLine(poolPresenceTableLine, rowPools, colPools, prPools, \
	pcPools, controlPools, calculateAddresses=True, fitDict=logReadNumberRatioHistogramFitDict, \
	areaDict=logReadNumberRatioHistogramIntegralDict, readCountThreshold=5, \
	maxSinglePoolCoordNumber=maxSinglePoolCoordNumber, maxTotalCoords=maxTotalCoords)
	
	readAlignmentCoord = poolPresenceTableLine['readAlignmentCoord']
	
	featureNames = FindFeatureAssociatedWithInsertionCoord(readAlignmentCoord, geneFeatureArray)
	
	summary.append(featureNames)
	
	if len(summary) > 7:
		summarize=False
		summaryDict = {}
		
		if useAmbiguousLinesOnly:
			if summary[6] == 'ambiguous':
				summarize=True
		else:
			summarize=True
			
		if summarize == True:
			summaryDict['coord'] = summary[0]
			summaryDict['rowCoords'] = summary[1]
			summaryDict['colCoords'] = summary[2]
			summaryDict['prCoords'] = summary[3]
			summaryDict['pcCoords'] = summary[4]
			summaryDict['controlPools'] = summary[5]
			summaryDict['locatability'] = summary[6]
			summaryDict['maxSinglePoolCoordNumber'] = summary[7]
			summaryDict['possibleAddressSummary'] = summary[8]
			summaryDict['featureNames'] = summary[9]
			
			summaryArray.append(summaryDict)

# 	summaryArray.append(summary)

	
	i += 1

summaryStr = GeneratePoolPresenceTableLineSummary(summaryArray)
outputHandle = open(outputLog, 'w')
outputHandle.write(summaryStr)
outputHandle.close()

# ------------------------------------------------------------------------------------------------ #

