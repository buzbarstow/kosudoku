# ------------------------------------------------------------------------------------------------ #

# Finish of Shewanella Sudoku Collection

# Import input data from the input file

from sudokuutils13.grid import FlattenSudokuGridLookupDict, AnnotateSudokuGridWithFeautureNames
from sudokuutils13.input import get_input
from sudokuutils13.utils import ensure_dir
from sudokuutils13.xml import ExportSudokuGridToXML, ImportSudokuGridFromXML

from sudokuutils13.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
CalculateFeatureArrayTaxonomy

from FinishQCSudokuCollectionUtils import FindFeaturesWithDisruptionsInProgenitorButNotInQC, \
PickBestSudokuWellToDisruptFeature, OutputPickingInstructions

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'pgsudokuGridXMLFile', 'cpSudokuGridFileName', \
'shewyGenomeGenBankFileName', 'pickingInstructionsFileName']


inputParameterFile = \
'/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/' \
+ 'FinishQCSudokuCollection.inp'

inputParameterValues = get_input(['', inputParameterFile], inputParameters)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
pgsudokuGridXMLFile = inputParameterValues['pgsudokuGridXMLFile']
cpSudokuGridFileName = inputParameterValues['cpSudokuGridFileName']
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
pickingInstructionsFileName = inputParameterValues['pickingInstructionsFileName']

ensure_dir(pickingInstructionsFileName)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import the quality controlled collection catalog
qcSudokuGridLookupDict, qcprPools, qcpcPools, qcrowPools, qccolPools = \
ImportSudokuGridFromXML(cpSudokuGridFileName)

# Flatten the qcSudokuGridLookupDict
qcSudokuWellArray = FlattenSudokuGridLookupDict(qcSudokuGridLookupDict, qcprPools, qcpcPools, \
qcrowPools, qccolPools)


# Import the progenitor collection catalog
pgSudokuGridLookupDict, pgprPools, pgpcPools, pgrowPools, pgcolPools = \
ImportSudokuGridFromXML(pgsudokuGridXMLFile)

# Flatten the progenitor collection catalog
pgSudokuWellArray = FlattenSudokuGridLookupDict(pgSudokuGridLookupDict, pgprPools, pgpcPools, \
pgrowPools, pgcolPools)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import a feature array
featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)

# Make a gene feature array
geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(feature)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Update the gene feature array with data from the progenitor collection
UpdateFeatureArrayWithSudokuWells(geneFeatureArray, qcSudokuWellArray)
UpdateFeatureArrayWithSudokuWells(geneFeatureArray, pgSudokuWellArray)

CalculateFeatureArrayTaxonomy(geneFeatureArray)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Calculate the features in the genome that have representatives in the progenitor only, and in the
# progenitor and QC collection
[featuresWithDisruptionsInQCAndProgenitor, featuresWithDisruptionsInProgenitorOnly] = \
FindFeaturesWithDisruptionsInProgenitorButNotInQC(geneFeatureArray)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Pick the best sudoku well objects from the features that only appear in the progenitor collection

featuresAndBestWellArray = PickBestSudokuWellToDisruptFeature(\
featuresWithDisruptionsInProgenitorOnly, pgSudokuGridLookupDict)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
OutputPickingInstructions(pickingInstructionsFileName, featuresAndBestWellArray)
# ----------------------------------------------------------------------------------------------- #
