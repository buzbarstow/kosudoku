#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# kosudoku-poolsolve
# Created by Buz Barstow 2015-06-14
# Last updated by Buz Barstow 2016-05-03

# Code needed for Bayesian inference solution of pool presence table
# ----------------------------------------------------------------------------------------------- #

from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.grid import ImportSudokuGridLayout, FlattenSudokuGridLookupDict, \
    AnnotateSudokuGridWithFeautureNames, WriteSudokuGridSummaryTable
from utilities.pool import ImportPoolPresenceTable, SolvePoolPresenceTable
from utilities.xml import ExportSudokuGridToXML
from utilities.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
    CalculateFeatureArrayTaxonomy
from utilities.poolanalysis import ReadBayesianInferenceParameters

import sys
import pdb
import os
import warnings

# ----------------------------------------------------------------------------------------------- #
warnings.simplefilter("ignore")
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'poolsolve.inp'



inputParameters = [ \
    'rowPools', 'colPools', 'controlPools', 'sudokuGridLayout', 'genBankFileName', \
    'poolPresenceTableFileName', \
    'readCountThreshold', 'maxGapForCoordGrouping', 'voigtScoreThreshold', 'maxTotalCoords', \
    'maxSinglePoolCoordNumber', \
    'sudokuGridXMLFile', 'sudokuGridSummaryFileName', \
    'bayesianParameterFileName']

inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
controlPools = inputParameterValues['controlPools'].split(',')
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
genBankFileName = file_intro + inputParameterValues['genBankFileName']
sudokuGridLayout = file_intro + inputParameterValues['sudokuGridLayout']
maxGapForCoordGrouping = int(inputParameterValues['maxGapForCoordGrouping'])
readCountThreshold = int(inputParameterValues['readCountThreshold'])
voigtScoreThreshold = float(inputParameterValues['voigtScoreThreshold'])
maxTotalCoords = int(inputParameterValues['maxTotalCoords'])
maxSinglePoolCoordNumber = int(inputParameterValues['maxSinglePoolCoordNumber'])
sudokuGridXMLFile = file_intro + inputParameterValues['sudokuGridXMLFile']
sudokuGridSummaryFileName = file_intro + inputParameterValues['sudokuGridSummaryFileName']
bayesianParameterFileName = file_intro + inputParameterValues['bayesianParameterFileName']
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
                                                                  colPools)

if controlPools[0].lower() == 'none':
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools
    controlPools = None
else:
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools \
                        + controlPools




poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Read in the Bayesian inference parameters

[logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict] = \
    ReadBayesianInferenceParameters(bayesianParameterFileName)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
ensure_dir(sudokuGridXMLFile)
ensure_dir(sudokuGridSummaryFileName)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
sudokuGridTaxonomyDict = SolvePoolPresenceTable(sudokuGridLookupDict, poolPresenceTable, \
                                                rowPools, colPools, prPools, pcPools, controlPools, logReadNumberRatioHistogramFitDict, \
                                                logReadNumberRatioHistogramIntegralDict, \
                                                readCountThreshold, voigtScoreThreshold, maxSinglePoolCoordNumber, maxTotalCoords, \
                                                maxGapForCoordGrouping)


featureArray = ImportFeatureArrayFromGenBank(genBankFileName)

geneFeatureArray = []
for feature in featureArray:
    if feature.featureType == 'gene':
        geneFeatureArray.append(feature)

flattenedSudokuGrid = FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, rowPools, \
                                                  colPools)

UpdateFeatureArrayWithSudokuWells(geneFeatureArray, flattenedSudokuGrid)

AnnotateSudokuGridWithFeautureNames(sudokuGridLookupDict, geneFeatureArray, prPools, pcPools, \
                                    rowPools, colPools)

featureArrayTaxonomyDict = CalculateFeatureArrayTaxonomy(geneFeatureArray)

ExportSudokuGridToXML(sudokuGridLookupDict, prPools, pcPools, rowPools, colPools, sudokuGridXMLFile)

WriteSudokuGridSummaryTable(sudokuGridSummaryFileName, sudokuGridLookupDict, prPools, pcPools, \
                            rowPools, colPools)

# ------------------------------------------------------------------------------------------------ #
