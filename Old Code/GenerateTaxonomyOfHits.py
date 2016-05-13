# ------------------------------------------------------------------------------------------------ #
def AttachGeneSequenceDataToFeatures(sequence, features):
	
	for feature in features:
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		
		feature.tagDict['coding_sequence'] = sequence[startCoord-1:endCoord]
	
	

# ------------------------------------------------------------------------------------------------ #


from sudokuutils13.utils import ensure_dir, ImportGenBankSequence
from sudokuutils13.input import get_input
from sudokuutils13.xml import ImportSudokuGridFromXML, ExportListForXML
from sudokuutils13.grid import FlattenSudokuGridLookupDict, AnnotateSudokuGridWithFeautureNames
from sudokuutils13.feature import ImportFeatureArrayFromGenBank, PrintFeatureArrayTaxonomyDict, \
UpdateFeatureArrayWithSudokuWells, CalculateFeatureArrayTaxonomy

import os

from copy import deepcopy

from GenerateTaxonomyOfHitsUtils import ImportSimplifiedCatalogOfHits, \
SortFeatureArrayAndAssignNumber, SelectFeaturesWithHits, SelectWellsWithHits, \
GroupHitsByGenomicLocation, BlastFeatureAgainstSubjectGenome, MakeGeneFeatureArray, \
MakeSummaryTableOfFeaturesWithHits, AnnotateSudokuWellListWithFeautureNames


# ------------------------------------------------------------------------------------------------ #
inputParameters = [\
'shewyGenomeGenBankFileName', 'cpSudokuGridFileName', 'hitsSummaryFileName', \
'hitsSummaryTableFileName']

inputFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/' \
+ 'Version 13/GenerateTaxonomyOfHits.inp'

inputParameterValues = get_input(['', inputFileName], \
inputParameters)

shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
cpSudokuGridFileName = inputParameterValues['cpSudokuGridFileName']
hitsSummaryFileName = inputParameterValues['hitsSummaryFileName']
hitsSummaryTableFileName = inputParameterValues['hitsSummaryTableFileName']


# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
# hitsSummaryTableFileNameDir = os.path.dirname(hitsSummaryTableFileName)
# ensure_dir(hitsSummaryTableFileNameDir + '/')
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Import Sudoku grid and build gene feature array
cpSudokuGridLookupDict, prPools, pcPools, rowPools, colPools = \
ImportSudokuGridFromXML(cpSudokuGridFileName)

flattenedCPSudokuGrid = FlattenSudokuGridLookupDict(cpSudokuGridLookupDict, prPools, pcPools, \
rowPools, colPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import a gene feature array
featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)

geneFeatureArray = []
progenitorGeneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(deepcopy(feature))
		progenitorGeneFeatureArray.append(deepcopy(feature))
		
shewySequence = ImportGenBankSequence(shewyGenomeGenBankFileName)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Import the cherry picked list
hitsWellArray = ImportSimplifiedCatalogOfHits(hitsSummaryFileName)

# Select sudoku wells with hits
wellsWithHits = SelectWellsWithHits(flattenedCPSudokuGrid, hitsWellArray)

# Give all of the features a number like the locus tag, but is an integer representing
# the position in the genome

SortFeatureArrayAndAssignNumber(geneFeatureArray)

UpdateFeatureArrayWithSudokuWells(geneFeatureArray, wellsWithHits)

featuresWithHits = SelectFeaturesWithHits(geneFeatureArray)

# Group the hits together by locus number and calculate the number of groups
groupedHits = GroupHitsByGenomicLocation(featuresWithHits, maxGap=5)

# Figure out how many groups there are
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Get sequence data for hits
AttachGeneSequenceDataToFeatures(shewySequence, featuresWithHits)
AnnotateSudokuWellListWithFeautureNames(flattenedCPSudokuGrid, geneFeatureArray)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
MakeSummaryTableOfFeaturesWithHits(hitsSummaryTableFileName, featuresWithHits, delimeter=',', \
printOnlyOneCoord=True)

# MakeSummaryTableOfSudokuWellsWithHits(flattenedCPSudokuGrid)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
tempBlastDir = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Analysis/Version 13/temp'

eColiFastaFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/U00096.3.fa'
eColiGenBankFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/U00096.3.gb'

eColiFeatureArray = ImportFeatureArrayFromGenBank(eColiGenBankFileName)
eColiGeneFeatureArray = MakeGeneFeatureArray(eColiFeatureArray)


featureAndHSPSArray = []

for feature in featuresWithHits:
	
	hsps = BlastFeatureAgainstSubjectGenome(feature, eColiFastaFileName, eColiGeneFeatureArray, \
	removeTempFiles=True, tempFileDir=tempBlastDir, evalue=0.1, task='blastn', perc_identity=30)
	
	featureAndHSPSArray.append([feature, hsps])
	
# ------------------------------------------------------------------------------------------------ #
	
	
	






