#!/sw/bin/python3.5


# ----------------------------------------------------------------------------------------------- #
# poolfit
# Program for fitting pool presence table read count ratio histograms to find Bayesian inference
# parameters for solution of pool presence table

# Created by Buz Barstow 2015-6-14
# Last updated by Sean Medin 2019-08-23
# ----------------------------------------------------------------------------------------------- #

import sys
import warnings
from matplotlib.pyplot import ion, ioff, show
import os

from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.grid import ImportSudokuGridLayout
from utilities.pool import ImportPoolPresenceTable
from utilities.poolanalysis import CalculateRatioOfReadNumbersForSingleAddressLines, \
    WriteBayesianInferenceParameters

# ----------------------------------------------------------------------------------------------- #
warnings.simplefilter("ignore")
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Turn the matplotlib interactive features on to allow multiple plots to be shown from the command
# line
ion()

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'poolfit.inp'


# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [ \
    'rowPools', 'colPools', 'controlPools', \
    'poolPresenceTableFileName', 'sudokuGridLayout', \
    'readCountThreshold', \
    'logReadNumberRatioHistogramsFileName', 'logReadNumberRatioFitFileName', \
    'bayesianParameterFileName']

inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data

# Some input parameters that we don't need right now, but might be useful

controlPools = inputParameterValues['controlPools'].split(',')
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
sudokuGridLayout = file_intro + inputParameterValues['sudokuGridLayout']
readCountThreshold = int(inputParameterValues['readCountThreshold'])
logReadNumberRatioHistogramsFileName = file_intro + inputParameterValues['logReadNumberRatioHistogramsFileName']
logReadNumberRatioFitFileName = file_intro + inputParameterValues['logReadNumberRatioFitFileName']
bayesianParameterFileName = file_intro + inputParameterValues['bayesianParameterFileName']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Initialize the output files
ensure_dir(logReadNumberRatioHistogramsFileName)
ensure_dir(logReadNumberRatioFitFileName)
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Get the sudoku grid plate layout and the plate row and plate column pools from the sudoku
# grid definition file
[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
                                                                  colPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the sudoku pool columns

if controlPools[0].lower() == 'none':
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools
    controlPools = None
else:
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools \
                        + controlPools
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table
poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Calculate the ratios of read numbers

[logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict] = \
    CalculateRatioOfReadNumbersForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, pcPools,
                                                     controlPools, readCountThreshold, logReadNumberRatioHistogramsFileName,
                                                     logReadNumberRatioFitFileName)

# Output the fit parameters for Bayesian inference
WriteBayesianInferenceParameters(logReadNumberRatioHistogramFitDict, 
                                 logReadNumberRatioHistogramIntegralDict, bayesianParameterFileName)
print('Written Bayesian inference parameters: ' + os.path.abspath(bayesianParameterFileName))
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Show the plots generated in the script
# Turn the interactive mode off to allow plots to stay on screen. 
ioff()
show()
# ------------------------------------------------------------------------------------------------ #
