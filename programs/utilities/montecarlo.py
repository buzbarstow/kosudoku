# ------------------------------------------------------------------------------------------------ #
def ImportEssentialityData(fileName):
# Not yet ready for prime time
# Import a defined format essentiality data file
# Assumes that data is in the format: locus tag, gene name, essentiality	

	from .utils import ParseCSVLine

	fileHandle = open(fileName, 'r')
	data = fileHandle.readlines()
	dataDict = {}
	
	i = 0
	while i < len(data):	
		# Ignore comment lines
		if data[i][0] != '#':
			dataLine = ParseCSVLine(data[i])
			dataDict[dataLine[0]] = [dataLine[1], dataLine[2]]
		i += 1

	return dataDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def BuildEssentialityDictThatIsKeyedByLocusTag(dataArray):
# Not yet ready for prime time
# Build essentiality data dict that is keyed by locus tag
	
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
def BuildCDSDictThatIsKeyedByLocusTag(cdsFeatures):
# Not yet ready for prime time

	i = 0
	cdsDict = {}
	while i < len(cdsFeatures):
		locusTag = cdsFeatures[i].tagDict['locus_tag'][0]
		cdsDict[locusTag] = cdsFeatures[i]
		i += 1
	
	return cdsDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SimulatePicking(hittableFeatures, notHittableFeatures, hittableTransposonCoords, \
transposonCoordToFeatureDict, maxMutants):
	
	from numpy.random import choice
	import pdb
	
	nonEssentialGeneCount = len(hittableFeatures)
	
	
	featureHitCountDict = {}
	for feature in hittableFeatures:
		featureHitCountDict[feature] = 0
	
	featuresHitAtLeastOnce = 0
	
	featuresHitAtLeastOnceVersusMutant = []
	
	i = 1
	while i <= maxMutants:
		randomCoord = int(choice(hittableTransposonCoords))
				
		featuresToBeHit = transposonCoordToFeatureDict[randomCoord]
		
		isAnyFeatureIncludingThisCoordNotHittable = False
		
		for featureToBeHit in featuresToBeHit:
			if featureToBeHit in notHittableFeatures:
				isAnyFeatureIncludingThisCoordNotHittable = True
		
		if isAnyFeatureIncludingThisCoordNotHittable == False:
			for featureToBeHit in featuresToBeHit:
				try:
					featureHitCountDict[featureToBeHit] += 1
				except:
					pdb.set_trace()
		
				if featureHitCountDict[featureToBeHit] == 1:
					featuresHitAtLeastOnce += 1
		
		featuresHitAtLeastOnceVersusMutant.append(featuresHitAtLeastOnce)
		
		i += 1
	
	return featuresHitAtLeastOnceVersusMutant
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SimulateMultiplePickings(transposonCoordToFeatureDictFile, numberOfTrials, maxMutants):
	
	from scipy import unique, intersect1d
	from numpy import mean, std, arange
	import xml.etree.ElementTree as ET
	import pdb
	
	
	transposonCoordToFeatureDictFileHandle = open(transposonCoordToFeatureDictFile, 'r')


	transposonCoordToFeatureDict = {}
	hittableFeatures = []
	hittableTransposonCoords = []
	notHittableTransposonCoords = []
	notHittableFeatures = []
	otherFeatures = []


	tree = ET.parse(transposonCoordToFeatureDictFile)
	root = tree.getroot()
	importedCoordsList = root.findall('coord')

	for coord in importedCoordsList:
		
		coordinate = int(coord.attrib['coord'])
		loci = coord.findall('locus')
	
		importedCoordsKeys = transposonCoordToFeatureDict.keys()
	
		if coordinate not in importedCoordsKeys:
			transposonCoordToFeatureDict[coordinate] = []
	
		for locus in loci:
			locusName = locus.attrib['locus']
			essentiality = locus.attrib['essentiality']
			transposonCoordToFeatureDict[coordinate].append(locusName)

			if essentiality == 'Dispensable':
				hittableTransposonCoords.append(coordinate)
				hittableFeatures.append(locusName)
			elif essentiality == 'Essential':
				notHittableFeatures.append(locusName)
				notHittableTransposonCoords.append(coordinate)
			else:
				otherFeatures.append(locusName)
				print(locusName)


	hittableFeatures = unique(hittableFeatures)
	hittableTransposonCoords = unique(hittableTransposonCoords)
	notHittableFeatures = unique(notHittableFeatures)
	otherFeatures = unique(otherFeatures)
	
	intersection = intersect1d(hittableFeatures, notHittableFeatures)
	
		
	# Simulate a number of picking runs
	featuresHitAtLeastOnceTrialsArray = []
	i = 0
	while i < numberOfTrials:
		featuresHitAtLeastOnceVersusMutant = \
		SimulatePicking(hittableFeatures, notHittableFeatures, hittableTransposonCoords, \
		transposonCoordToFeatureDict, maxMutants)
	
		featuresHitAtLeastOnceTrialsArray.append(featuresHitAtLeastOnceVersusMutant)	
		i += 1
		
	# Collect together then data from the picking runs for calculation of mean and standard 
	# deviation of number of hits picked

	i = 0
	collectedFeatureHitCountArray = []
	while i < len(featuresHitAtLeastOnceTrialsArray[0]):
		collectedFeatureHitCountArray.append([])
		i += 1

	i = 0
	while i < len(collectedFeatureHitCountArray):
		j = 0
		while j < len(featuresHitAtLeastOnceTrialsArray):
			collectedFeatureHitCountArray[i].append(featuresHitAtLeastOnceTrialsArray[j][i])
			j += 1
		i += 1

	averageFeatureHitCount = []
	sdFeatureHitCount = []
	featureHitCountUpperBound = []
	featureHitCountLowerBound = []

	# Calculate the mean and standard deviation of the number of unique features hit at each pick
	# from the trials

	i = 0
	while i < len(collectedFeatureHitCountArray):
		averageFeatureHitCount.append(mean(collectedFeatureHitCountArray[i]))
		sdFeatureHitCount.append(std(collectedFeatureHitCountArray[i]))

		featureHitCountUpperBound.append(averageFeatureHitCount[i] + sdFeatureHitCount[i])
		featureHitCountLowerBound.append(averageFeatureHitCount[i] - sdFeatureHitCount[i])

		i += 1
	
	
	# Prepare an x axis (the number of mutants picked) for the output
	iAxis = arange(1, maxMutants+1, 1)
	
	noUniqHittableFeatures = len(hittableFeatures)
	
	return [iAxis, averageFeatureHitCount, sdFeatureHitCount, featureHitCountUpperBound, \
	featureHitCountLowerBound, noUniqHittableFeatures ]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PoissonEstimateOfGenesHit(iAxis, noUniqHittableFeatures):
	
	from numpy import exp, array, float
	
	uniqueGenesHit = []
	i = 0
	while i < len(iAxis):
		ans = noUniqHittableFeatures*(1-exp(-iAxis[i]/noUniqHittableFeatures))
		uniqueGenesHit.append(ans)
		i += 1
		
	uniqueGenesHit = array(uniqueGenesHit, float)
	
	return uniqueGenesHit
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def FindATandTAPositions2(genomeFile, format='genbank'):
# Does the same thing as FindATandTAPositions but can work with a GenBank or a Fasta file, \
# so you only  need one file format
	
	import re
	from pdb import set_trace
	
	if format == 'genbank':
		sequence = ImportGenBankSequence(genomeFile)
	elif format == 'fasta':
		sequence = ImportFastaSequence(genomeFile)
	
	ATandTAPositions = []
	
	atRegex = re.compile('(at|ta)', re.IGNORECASE)
	
# 	set_trace()
	
	i = 0
	while i < len(sequence) - 1:
		atMatch = atRegex.match(sequence[i:i+2])
		
		if atMatch != None:
			ATandTAPositions.append(i+1)
		
		i += 1
	
	
	
	return [ATandTAPositions, sequence]
# ------------------------------------------------------------------------------------------------ #


