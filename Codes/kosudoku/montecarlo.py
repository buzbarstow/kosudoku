# ------------------------------------------------------------------------------------------------ #
def ImportEssentialityData(fileName):
# Not yet ready for prime time
# Import a defined format essentiality data file
# Assumes that data is in the format: locus tag, gene name, essentiality	
	fileHandle = open(fileName, 'r')
	data = fileHandle.readlines()
	dataDict = {}
	
	i = 0
	while i < len(data):	
		# Ignore comment lines
		if data[i][0] != '#':
			dataLine = data[i].strip().split(',')
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
def SimulatePicking(hittableFeatures, hittableTransposonCoords, transposonCoordToFeatureDict, \
maxMutants):
	
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
		
# 		pdb.set_trace()
		
		featureHitCountDict[transposonCoordToFeatureDict[randomCoord]] += 1
		
		if featureHitCountDict[transposonCoordToFeatureDict[randomCoord]] == 1:
			featuresHitAtLeastOnce += 1
		
		featuresHitAtLeastOnceVersusMutant.append(featuresHitAtLeastOnce)
		
		i += 1
	
	return featuresHitAtLeastOnceVersusMutant
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SimulateMultiplePickings(transposonCoordToFeatureDictFile, numberOfTrials, maxMutants):
	
	from scipy import unique
	from numpy import mean, std, arange
	import pdb
	
	transposonCoordToFeatureDictFileHandle = open(transposonCoordToFeatureDictFile, 'r')

	transposonCoordData = transposonCoordToFeatureDictFileHandle.readlines()

	transposonCoordToFeatureDict = {}
	hittableFeatures = []
	hittableTransposonCoords = []

	# Decide what is an essential gene and what isn't

	i = 0
	while i < len(transposonCoordData):
		line = transposonCoordData[i].strip().split(',')
		coord = int(line[0])
		locus = line[1]
		essentiality = line[2]
	
		if essentiality == 'Dispensable':
			transposonCoordToFeatureDict[coord] = locus
			hittableTransposonCoords.append(coord)
			hittableFeatures.append(locus)
	
		i += 1
	
	# Come up with a list of unique hittable features

	hittableFeatures = unique(hittableFeatures)
	hittableTransposonCoords = unique(hittableTransposonCoords)
	
	
	# Simulate a number of picking runs
	featuresHitAtLeastOnceTrialsArray = []
	i = 0
	while i < numberOfTrials:
		featuresHitAtLeastOnceVersusMutant = \
		SimulatePicking(hittableFeatures, hittableTransposonCoords, transposonCoordToFeatureDict, \
		maxMutants)
	
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


