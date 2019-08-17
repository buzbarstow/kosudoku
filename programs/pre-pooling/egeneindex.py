#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# kosudoku-egeneindex
# Created by Buz Barstow 2016-05-07
# Last modified by Sean Medin 2019-08-16
#
# Code that uses gene essentiality data to build a lookup table connecting the location
# of each possible transposon insertion site in a genome to a locus and its essentiality.
# At this point, we expect that the essentiality data and the genome sequence file have been
# conditioned so that there is an essentiality tag for each locus.
# ----------------------------------------------------------------------------------------------- #

from utilities.utils import ensure_dir, FindATandTAPositions
from utilities.input import get_input
from utilities.feature import ImportFeatureArrayFromGenBank
from utilities.montecarlo import ImportEssentialityData
from utilities.xml import indent
import xml.etree.ElementTree as ET
import pdb
import sys

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file
inputParameters = [ \
    'essentialityDataFileName', 'referenceGenomeGenBankFileName', 'outputFileName', \
    'featureTypeToIndex']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/pre-pool-files/'
file_name = file_intro + 'egeneindex.inp'
inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
referenceGenomeGenBankFileName = file_intro + inputParameterValues['referenceGenomeGenBankFileName']
essentialityDataFileName = file_intro + inputParameterValues['essentialityDataFileName']
outputFileName = file_intro + inputParameterValues['outputFileName']
featureTypeToIndex = inputParameterValues['featureTypeToIndex']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
ensure_dir(outputFileName)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
[ATandTAPositions, sequence] = FindATandTAPositions(referenceGenomeGenBankFileName, \
                                                    format='genbank')

essentialityDict = ImportEssentialityData(essentialityDataFileName)
featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Truncate the AT and TA position list to make the computation faster
# pdb.set_trace()
# ATandTAPositions = ATandTAPositions[0:101]
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Go through all of the features in the imported genome and extract the ones that refer to genes
geneFeatures = []

if featureTypeToIndex == 'CDS':
    for feature in featureArray:
        if feature.featureType == 'CDS':
            geneFeatures.append(feature)
elif featureTypeToIndex == 'gene':
    for feature in featureArray:
        if feature.featureType == 'gene':
            geneFeatures.append(feature)
else:
    for feature in featureArray:
        if feature.featureType == 'CDS':
            geneFeatures.append(feature)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Go through the gene feature array (which is where we have gene coordinate information) and
# see if we can find a matching locus in the essentiality data. If so, mark the gene with that data,
# if not, mark the gene as unknown which implies non-essential (giving the smallest number of
# essential genes and the largest number of non-essentials).

# Makes the assumption that there is only one locus tag per gene

essentialityLocusTags = essentialityDict.keys()

geneLocusDict = {}

for feature in geneFeatures:
    locusTag = feature.tagDict['locus_tag'][0]
    startCoord = feature.startCoord
    endCoord = feature.endCoord

    if locusTag in essentialityLocusTags:
        feature.tagDict['Essentiality'] = essentialityDict[locusTag][1]
    else:
        feature.tagDict['Essentiality'] = 'Dispensable'

    essentialityTag = feature.tagDict['Essentiality']

    geneLocusDict[startCoord] = \
        {'locusTag': locusTag, 'endCoord': endCoord, 'essentialityTag': essentialityTag}

startCoords = sorted(list(geneLocusDict.keys()))
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Next, tie the AT and TA positions to gene loci
# Use this stage of the algorithm to calculate what is hittable and what isn't.

transposonCoordToFeatureDict = {}

i = 0
j = 0
transposonCoordsWithMoreThanOneFeature = 0

while i < len(ATandTAPositions):

    coord = ATandTAPositions[i]
    transposonCoordToFeatureDict[coord] = []
    while j < len(startCoords) and geneLocusDict[startCoords[j]]['endCoord'] < coord:
        j += 1
    k = 0
    while j + k < len(startCoords) and startCoords[j + k] <= coord:
        endCoord = geneLocusDict[startCoords[j + k]]['endCoord']
        if endCoord < coord:
            del startCoords[j+k]
            continue
        locuskey = geneLocusDict[startCoords[j+k]]['locusTag']
        essentiality = geneLocusDict[startCoords[j+k]]['essentialityTag']
        if coord <= endCoord:
            transposonCoordToFeatureDict[coord].append([locuskey, essentiality])
        k += 1

    if len(transposonCoordToFeatureDict[coord]) > 1:
        transposonCoordsWithMoreThanOneFeature += 1
    elif len(transposonCoordToFeatureDict[coord]) == 0:
        transposonCoordToFeatureDict[coord] = [['non-coding', 'Dispensable']]

    i += 1
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Write out the coordinate to feature dict as an XML file
coordListRoot = ET.Element('coordList')

coordKeys = transposonCoordToFeatureDict.keys()

for coordKey in coordKeys:
    coordSubElement = ET.SubElement(coordListRoot, 'coord')
    coordSubElement.set('coord', str(coordKey))

    locusArray = transposonCoordToFeatureDict[coordKey]

    for locus in locusArray:
        locusSubElement = ET.SubElement(coordSubElement, 'locus')
        locusSubElement.set('locus', locus[0])
        locusSubElement.set('essentiality', locus[1])

indent(coordListRoot)
coordListTree = ET.ElementTree(coordListRoot)
coordListTree.write(outputFileName)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
print("Transposons with more than one corresponding Feature: " \
      + str(transposonCoordsWithMoreThanOneFeature))
print("Total transposon insertion sites: " + str(len(ATandTAPositions)))
# ----------------------------------------------------------------------------------------------- #