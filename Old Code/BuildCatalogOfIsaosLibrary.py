from sudokuutils13.utils import ensure_dir
from sudokuutils13.input import get_input
from sudokuutils13.xml import ImportSudokuGridFromXML, ExportListForXML
from sudokuutils13.grid import FlattenSudokuGridLookupDict, AnnotateSudokuGridWithFeautureNames
from sudokuutils13.feature import ImportFeatureArrayFromGenBank, PrintFeatureArrayTaxonomyDict, \
UpdateFeatureArrayWithSudokuWells, CalculateFeatureArrayTaxonomy

import os

from copy import deepcopy

from BuildCatalogOfIsaosLibraryUtils import ImportSimplifiedCatalogOfIsaosLibrary, \
UpdateSudokuGridWithRearrayedData, CalculateFeatureTaxonomyOfRearrayedCollection, \
WriteRearrayedCollectionSummaryTable


# ------------------------------------------------------------------------------------------------ #
inputParameters = [\
'shewyGenomeGenBankFileName', 'sudokuGridXMLFile', 'cherryPickedFileName', \
'rearrayedCollectionSummaryTableFileName']

inputParameterValues = get_input(['','/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/SudokuAnalyze13_BuildCatalogOfIsaosLibrary.inp'], \
inputParameters)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
sudokuGridXMLFile = inputParameterValues['sudokuGridXMLFile']
cherryPickedFileName = inputParameterValues['cherryPickedFileName']
rearrayedCollectionSummaryTableFileName = inputParameterValues['rearrayedCollectionSummaryTableFileName']
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
rearrayedCollectionSummaryTableFileNameDir = os.path.dirname(rearrayedCollectionSummaryTableFileName)
ensure_dir(rearrayedCollectionSummaryTableFileNameDir + '/')
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import Sudoku grid and build gene feature array
sudokuGridLookupDict, prPools, pcPools, rowPools, colPools = \
ImportSudokuGridFromXML(sudokuGridXMLFile)

flattenedSudokuGrid = FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, rowPools,\
colPools)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Build the gene feature array
featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)

geneFeatureArray = []
progenitorGeneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(deepcopy(feature))
		progenitorGeneFeatureArray.append(deepcopy(feature))

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Import the cherry picked list
rearrayedWellArray = ImportSimplifiedCatalogOfIsaosLibrary(cherryPickedFileName)

rearrayedSudokuWellList = UpdateSudokuGridWithRearrayedData(rearrayedWellArray, \
sudokuGridLookupDict, prPools, pcPools, rowPools, colPools)

UpdateFeatureArrayWithSudokuWells(geneFeatureArray, rearrayedSudokuWellList)

AnnotateSudokuGridWithFeautureNames(sudokuGridLookupDict, geneFeatureArray, prPools, pcPools, \
rowPools, colPools)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Calculate the feature taxonomy of Isao's collection
taxonomyDict = CalculateFeatureTaxonomyOfRearrayedCollection(geneFeatureArray)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
WriteRearrayedCollectionSummaryTable(rearrayedCollectionSummaryTableFileName, \
rearrayedSudokuWellList, rearrayedCollectionAddressDictKey='isaosCollection', \
progenitorCollectionAddressDictKey='progenitor')
# ------------------------------------------------------------------------------------------------ #

