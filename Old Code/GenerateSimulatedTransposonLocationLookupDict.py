#!/sw/bin/python3.4


# ------------------------------------------------------------------------------------------------ #
def ImportEssentialityData(fileName):
	
	fileHandle = open(fileName, 'r')
	
	data = fileHandle.readlines()
	
	headers = data[0].strip().split(',')
	
	dataArray = []
		
	# Ignore the header line
	i = 1
	while i < len(data):
		j = 0
		dataLine = data[i].strip().split(',')
		
		dataDict = {}
		
		while j < len(headers):
			dataDict[headers[j]] = dataLine[j]
			j += 1
		
		dataArray.append(dataDict)
		i += 1

	return [headers, dataArray]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ConvertColumnInEssentialityDataDictsToArray(essentialityDataArray, header):
	
	column = []
	
	i = 0
	while i < len(essentialityDataArray):
		column.append(essentialityDataArray[i][header])
		i += 1
	
	return column
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Build essentiality data dict that is keyed by locus tag
def BuildEssentialityDictThatIsKeyedByLocusTag(headers, dataArray):
	
	essentialityDict = {}
	
	locusTags = []
	
	headersWithoutSysName = []
	
	i = 0
	while i < len(headers):
		if headers[i] != 'sysName':
			headersWithoutSysName.append(headers[i])
		i += 1
	
	dataDict = {}
	
	for line in dataArray:
		
		
		dataDict[line['sysName']] = {}
		
		for header in headersWithoutSysName:
			dataDict[line['sysName']][header] = line[header]
	
	return dataDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def BuildLocusToGeneDictThatIsKeyedByLocusTag(features):
	i = 0
	locusToGeneDict = {}
	while i < len(features):
		
		locusTags = features[i].tagDict['locus_tag']
		
		j = 0
		while j < len(locusTags):
			locusToGeneDict[locusTags[j]] = features[i]
			j += 1
		i += 1
	
	return locusToGeneDict
# ------------------------------------------------------------------------------------------------ #




# ----------------------------------------------------------------------------------------------- #
# GenerateTransposonLocationLookupDict.py
# Created by Buz Barstow 2015-6-29
# Last updated by Buz Barstow 2016-2-17

# Code needed to generate a lookup dictionary connecting randomly generated transposon locations to 
# features in the genome and their essentiality. Speeds up Monte Carlo code enormously
# ----------------------------------------------------------------------------------------------- #


from sudokuutils13.input import get_input
from sudokuutils13.utils import FindATandTAPositions2
from sudokuutils13.feature import ImportFeatureArrayFromGenBank
from sudokuutils13.xml import ExportListForXML


import sys
import pdb
import re
from scipy import unique, setdiff1d, union1d, in1d, intersect1d
import numpy

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'shewyGenomeGenBankFileName', 'shewyGenomeFastaFileName', 'essentialityDataFileName', \
'outputFileName']


inputFileName = \
'/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' + \
'2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/' + \
'ShewyMonteCarlo.inp'


inputParameterValues = get_input(['', inputFileName], inputParameters)

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
shewyGenomeFastaFileName = inputParameterValues['shewyGenomeFastaFileName']
essentialityDataFileName = inputParameterValues['essentialityDataFileName']
outputFileName = inputParameterValues['outputFileName']

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
[ATandTAPositions, sequence] = FindATandTAPositions2(shewyGenomeGenBankFileName, format='genbank')

[essentialityHeaders, essentialityDataArray] = ImportEssentialityData(essentialityDataFileName)

featureArray = ImportFeatureArrayFromGenBank(shewyGenomeGenBankFileName)
geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'gene':
		geneFeatureArray.append(feature)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
features = []
locusTags = []


for feature in geneFeatureArray:
	
	tagDictKeys = list(feature.tagDict.keys())
	
	if 'locus_tag' in tagDictKeys:
		featureLocusTags = feature.tagDict['locus_tag']
			
		for tag in featureLocusTags:
			locusTags.append(tag)
			


# Figure out the feature types in the genome
uniqueLocusTags = unique(locusTags)


# Get all of the locus tags in the Deutschbauer data
locusTagsEssentiality = ConvertColumnInEssentialityDataDictsToArray(essentialityDataArray, 'sysName')

locusTagsEssentialityUnique = unique(locusTagsEssentiality)

essentialLocusTagsNotInGenBank = setdiff1d(locusTagsEssentialityUnique, uniqueLocusTags)
genbankLocusTagsNotInEssential = setdiff1d(uniqueLocusTags, locusTagsEssentialityUnique)


essentialityDict = \
BuildEssentialityDictThatIsKeyedByLocusTag(essentialityHeaders, essentialityDataArray)

cdsDict = \
BuildLocusToGeneDictThatIsKeyedByLocusTag(geneFeatureArray)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Go through the gene feature array and find essentiality data for that locus

essentialityTags = ['ADP_Essential', 'PEC Essential', 'Gerdes Essential', 'Essential Class']

essentialityDictKeys = list(essentialityDict.keys())

geneFeaturesWithEssentialityData = 0
locusTagsWithNoEssentialityData = []
essentialityClassifications = []

i = 0
while i < len(geneFeatureArray):
	locusTags = geneFeatureArray[i].tagDict['locus_tag']
	
	for locusTag in locusTags:
	
		if locusTag in essentialityDictKeys:
			# If there is essentiality data for the locus, copy it over
			essentialityData = essentialityDict[locusTag]
	
			for tag in essentialityTags:
				geneFeatureArray[i].updateTagDict(tag, essentialityData[tag])
		
			essentialityClassifications.append(essentialityData['Essential Class'])
			
			geneFeaturesWithEssentialityData += 1
		else:
			# If there is no essentiality data for the locus, guess that it is dispensable
			locusTagsWithNoEssentialityData.append(locusTag)
			geneFeatureArray[i].updateTagDict('Essential Class', 'GuessDispensable')
	
	i += 1
	
uniqueEssentialClasses = unique(essentialityClassifications)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Build a simplified gene locus array that lists that start and end of a locus and its 
# essentiality
i = 0
geneLocusDict = {}

while i < len(geneFeatureArray):
	locusTags = geneFeatureArray[i].tagDict['locus_tag']
	
	j = 0
	while j < len(locusTags):
	
		locusTag = locusTags[j]
	
		startCoord = geneFeatureArray[i].startCoord
		endCoord = geneFeatureArray[i].endCoord
	
		if endCoord < startCoord:
			print("Start coord is earlier than end coord in " + str(locusTag))
	
		essentialityTag = geneFeatureArray[i].tagDict['Essential Class'][j]
	
		geneLocusDict[locusTag] = \
		{'startCoord':startCoord, 'endCoord':endCoord, 'essentialityTag':essentialityTag}
		
		j += 1
		
	i += 1

# Note, here the locusKeys are generated before we add the non-coding part!
locusKeys = list(geneLocusDict.keys())
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Next, tie the AT and TA positions to gene loci
# Use this stage of the algorithm to calculate what is hittable and what isn't. 
# This algorithm allows for a coordinate to be inside multiple loci

transposonCoordToFeatureDict = {}
ATandTAPositionsThatAreHittable = []

i = 0
while i < len(ATandTAPositions):
	
	coord = ATandTAPositions[i]
	
	transposonCoordToFeatureDict[coord] = {'locusTags':[], 'essentialityTags':[]}
	
	j = 0
	coordFoundInsideLoci = False
	while j < len(locusKeys):
		key = locusKeys[j]
		startCoord = geneLocusDict[key]['startCoord']
		endCoord = geneLocusDict[key]['endCoord']
		essentiality = geneLocusDict[key]['essentialityTag']
		
		if (startCoord < coord < endCoord):
			transposonCoordToFeatureDict[coord]['locusTags'].append(key)
			transposonCoordToFeatureDict[coord]['essentialityTags'].append(essentiality)
			coordFoundInsideLoci = True
		j += 1
	
	if coordFoundInsideLoci == False:
		transposonCoordToFeatureDict[coord]['locusTags'].append('None')
		transposonCoordToFeatureDict[coord]['essentialityTags'].append('Dispensable')
		
	i += 1

transposonCoordToFeatureDictFileHandle = open(outputFileName, 'w')

for key in transposonCoordToFeatureDict.keys():
	data = transposonCoordToFeatureDict[key]
	
	outputStr = str(key) + ',' \
	+ '"' + ExportListForXML(data['locusTags'], delimeter=',') + '"' + ',' \
	+ '"' + ExportListForXML(data['essentialityTags'], delimeter=',') + '"' + '\n'
	
	transposonCoordToFeatureDictFileHandle.write(outputStr)

transposonCoordToFeatureDictFileHandle.close()
# ----------------------------------------------------------------------------------------------- #


