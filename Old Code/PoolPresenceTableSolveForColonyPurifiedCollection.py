#!/sw/bin/python3.4


# ----------------------------------------------------------------------------------------------- #
# PoolPresenceTableSolveForColonyPurifiedCollection.py
# Created by Buz Barstow 2015-6-14
# ----------------------------------------------------------------------------------------------- #



from IsItInThereUtils import GenerateColonyPurificationDict, \
ImportDetailedColonyPurificationCatalog, UpdateColonyPurificationDictWithProgenitorCollection, \
UpdateColonyPurificationDictWithInternalAddressFormats, FillCoordinatePools, \
UpdateColonyPurificationDictWithPoolPresenceTableData, GeneratePoolPresenceDict, \
CalculateIntersectionReadCounts, ExpandOutUniqueCoords, ReducePoolPresenceTable, \
SolvePoolPresenceTableForColonyPurifiedCollection, \
UpdateSudokuGridLookupDictToSudokuColonyPurifiedWell, \
AssignLikelyContentsToSudokuGridLookupDictAndCheckPredictions, UpdateLikelyReadAlignmentCoords, \
AnnotateColonyPurifiedSudokuGridWithFeautureNames, WriteColonyPurifiedCollectionSummaryTable, \
ConvertSudokuGridLookupDictBackToSudokuWells


from sudokuutils13.input import get_input

from sudokuutils13.utils import ensure_dir

from sudokuutils13.grid import ImportSudokuGridLayout, WriteSudokuGridSummaryTable, \
AnnotateSudokuGridWithFeautureNames, FlattenSudokuGridLookupDict

from sudokuutils13.pool import ImportPoolPresenceTable, SolvePoolPresenceTable

from sudokuutils13.xml import ExportSudokuGridToXML, ImportSudokuGridFromXML

from sudokuutils13.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
CalculateFeatureArrayTaxonomy


import sys
import pdb
import os

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'poolPresenceTableFileName', \
'sudokuGridLayout', 'outputLog', 'maxGapForCoordGrouping', \
'readCountThreshold', 'maxTotalCoords', 'maxSinglePoolCoordNumber', \
'pgsudokuGridXMLFile', 'shewyGenomeGenBankFileName', 'colonyPurificationDataFileName', \
'cpSudokuGridSummaryFileName', 'cpSudokuGridFileName']


inputParameterFile = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/SudokuAnalyze13_PoolPresenceSolveForColonyPurifiedCollection.inp'

inputParameterValues = get_input(['', inputParameterFile], inputParameters)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
poolPresenceTableFileName = inputParameterValues['poolPresenceTableFileName']

# outputLog = inputParameterValues['outputLog']
sudokuGridLayout = inputParameterValues['sudokuGridLayout']
maxGapForCoordGrouping = int(inputParameterValues['maxGapForCoordGrouping'])
readCountThreshold = int(inputParameterValues['readCountThreshold'])
maxTotalCoords = int(inputParameterValues['maxTotalCoords'])
maxSinglePoolCoordNumber = int(inputParameterValues['maxSinglePoolCoordNumber'])
pgsudokuGridXMLFile = inputParameterValues['pgsudokuGridXMLFile']
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
colonyPurificationDataFileName = inputParameterValues['colonyPurificationDataFileName']
cpSudokuGridSummaryFileName = inputParameterValues['cpSudokuGridSummaryFileName']
cpSudokuGridFileName = inputParameterValues['cpSudokuGridFileName']
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
cpSudokuGridSummaryFileNameDir = os.path.dirname(cpSudokuGridSummaryFileName)
ensure_dir(cpSudokuGridSummaryFileNameDir + '/')
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
colPools)

sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools

poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Read the detailed colony purification data file in. This contains the progenitor collection 
# locations for the colony purified bugs

detailedColonyPurificationData = \
ImportDetailedColonyPurificationCatalog(colonyPurificationDataFileName)

colonyPurificationDict = GenerateColonyPurificationDict(detailedColonyPurificationData)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the progenitor collection Sudoku Grid
pSudokuGridLookupDict, pgprPools, pgpcPools, pgrowPools, pgcolPools = \
ImportSudokuGridFromXML(pgsudokuGridXMLFile)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Make a list of transposons that we added into the colony purified collection
flatPSudokuGrid = FlattenSudokuGridLookupDict(pSudokuGridLookupDict, pgprPools, pgpcPools, \
pgrowPools, pgcolPools)

coordsPutIn = []

for entry in flatPSudokuGrid:
	for coord in entry.readAlignmentCoords:
		coordsPutIn.append(coord.coord)

uniqueCoordsPutIn = unique(coordsPutIn)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Go through the progenitor sudoku grid and add data from it to the colonyPurificationDict
updatedColonyPurificationDict = \
UpdateColonyPurificationDictWithProgenitorCollection(colonyPurificationDict, pSudokuGridLookupDict,\
pgprPools, pgpcPools, pgrowPools, pgcolPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Expand out the unique coords put in so that there is some wiggle room for alignment error. 
# Expand it out by about 4 bp either side.
uniqueCoordsPutInExpanded = ExpandOutUniqueCoords(updatedColonyPurificationDict)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Reduce down the pool presence table
reducedPoolPresenceTable = ReducePoolPresenceTable(poolPresenceTable, uniqueCoordsPutInExpanded)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Update plate name addresses in updatedColonyPurificationDict with internal names
updatedColonyPurificationDict = \
UpdateColonyPurificationDictWithInternalAddressFormats(updatedColonyPurificationDict, \
sudokuGridLookupDict, prPools, pcPools)			
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Solve the pool presence table
# This is a much simpler solution than that used for the progenitor collection, so I'm using
# a special function for this case. 
SolvePoolPresenceTableForColonyPurifiedCollection(sudokuGridLookupDict, reducedPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, 5)

# Update the sudoku well objects in the sudokugridlookupdict to a sub-class (colony purified well)
# that has a little extra data in it. 
UpdateSudokuGridLookupDictToSudokuColonyPurifiedWell(sudokuGridLookupDict, rowPools, colPools,\
prPools, pcPools)


# Check the sudoku grid lookup dict against the predicted catalog and assign likely contents
AssignLikelyContentsToSudokuGridLookupDictAndCheckPredictions(sudokuGridLookupDict, \
updatedColonyPurificationDict)

flattenedSudokuWellList = FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, \
rowPools, colPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Update the flattened sudoku well list to make a list of sudoku genomic coord objects from the 
# simplified likely read alignment coords

UpdateLikelyReadAlignmentCoords(flattenedSudokuWellList)


# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)

geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(feature)

# UpdateFeatureArrayWithSudokuWells(geneFeatureArray, flattenedSudokuGrid)
# featureArrayTaxonomyDict = CalculateFeatureArrayTaxonomy(geneFeatureArray)

AnnotateColonyPurifiedSudokuGridWithFeautureNames(flattenedSudokuWellList, geneFeatureArray)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Conver the sudoku grid back to regular wells so that it can be written out with regular functions
newSudokuGrid = ConvertSudokuGridLookupDictBackToSudokuWells(sudokuGridLookupDict, rowPools, \
colPools, prPools, pcPools)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Write out the colony purified collection summary
WriteColonyPurifiedCollectionSummaryTable(flattenedSudokuWellList, cpSudokuGridSummaryFileName, \
cpCollectionAddressDictKey='colonyPurified', progenitorCollectionAddressDictKey='progenitor')
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Write out the new sudoku grid as xml
ExportSudokuGridToXML(newSudokuGrid, prPools, pcPools, rowPools, colPools, cpSudokuGridFileName)


# ------------------------------------------------------------------------------------------------ #
