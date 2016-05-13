# ------------------------------------------------------------------------------------------------ #
def FindFeaturesWithDisruptionsInProgenitorButNotInQC(featureArray): 
# Take an updated feature array and check it for features that have representative mutants in the 
# progenitor but not in the QC collection

	i = 0

	featuresWithDisruptionsInQCAndProgenitor = []
	featuresWithDisruptionsInProgenitorOnly = []
	

	while i < len(featureArray):
		if len(featureArray[i].sudokuGridEntries) > 0:
			
			featureHasAtLeastOneRepInQCCollection = False
			featureHasAtLeastOneRepInProgenitorCollection = False
			
			sudokuGridEntries = featureArray[i].sudokuGridEntries
			j = 0
			while j < len(sudokuGridEntries):
				
				addressDictKeys = list(sudokuGridEntries[j]['addressDict'].keys())
				
				if 'colonyPurified' in addressDictKeys:
					featureHasAtLeastOneRepInQCCollection = True
				
				if 'progenitor' in addressDictKeys:
					featureHasAtLeastOneRepInProgenitorCollection = True
				
				j += 1
				
			if featureHasAtLeastOneRepInProgenitorCollection == True \
			and featureHasAtLeastOneRepInQCCollection == True:
				featuresWithDisruptionsInQCAndProgenitor.append(featureArray[i])
			elif featureHasAtLeastOneRepInProgenitorCollection == True \
			and featureHasAtLeastOneRepInQCCollection == False:
				featuresWithDisruptionsInProgenitorOnly.append(featureArray[i])

		
		i += 1

	return [featuresWithDisruptionsInQCAndProgenitor, featuresWithDisruptionsInProgenitorOnly]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Pick Best Mutant To Disrupt Feature
def PickBestSudokuWellToDisruptFeature(featureArray, sudokuGridLookupDict, pickProbability=0.95):
	
	import pdb
	from operator import itemgetter
	from math import floor, ceil
	from numpy import log
	
	featuresAndBestWellArray = []
	
	i = 0
	while i < len(featureArray):
		
		feature = featureArray[i]
		featureLocusTag = feature.tagDict['locus_tag']

		sudokuGridEntries = feature.sudokuGridEntries
		sudokuWellObjectsToScore = []
		
		j = 0
		while j < len(sudokuGridEntries) and len(sudokuGridEntries) > 0:
				
			addressDictKeys = list(sudokuGridEntries[j]['addressDict'].keys())
			
			if 'progenitor' in addressDictKeys:
				plateCol = sudokuGridEntries[j]['addressDict']['progenitor']['plateCol']
				plateRow = sudokuGridEntries[j]['addressDict']['progenitor']['plateRow']
				row = sudokuGridEntries[j]['addressDict']['progenitor']['row']
				col = sudokuGridEntries[j]['addressDict']['progenitor']['col']
				wellOccupants = sudokuGridEntries[j]['wellOccupants']
				fracDistFromTranslationStart = sudokuGridEntries[j]['fracDistFromTranslationStart']
				locatability = sudokuGridEntries[j]['locatability']
				
				sudokuWell = sudokuGridLookupDict[plateRow][plateCol].wellGrid[row][col]
				
				if locatability == 'unambiguous':
					locatabilityDiscountFactor = 1.0
				else:
					locatabilityDiscountFactor = 0.5
				
				purifiability = \
				(1.0 - fracDistFromTranslationStart)*locatabilityDiscountFactor*(1.0/wellOccupants)
				
				pA = (1.0/(float(wellOccupants)))
		
				if pA == 1.0:
					minimumColoniesToPick = 1.0
				else:
					minimumColoniesToPick = floor(log(1-pickProbability)/log(1-pA))
				
				sudokuWellObjectsToScore.append({'well':sudokuWell, 'purifiability':purifiability,\
				'sudokuGridEntry':sudokuGridEntries[j], \
				'minimumColoniesToPick':minimumColoniesToPick})
				
			j += 1
				
		
		sudokuWellObjectsToScore = sorted(sudokuWellObjectsToScore, reverse=True, \
		key=itemgetter('purifiability'))
		bestWell = sudokuWellObjectsToScore[0]

		featuresAndBestWellArray.append([feature, bestWell])
		
		i += 1

	return featuresAndBestWellArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def OutputPickingInstructions(outputFile, featuresAndBestWellArray):
	headerStr = 'Feature Name,Locus Tag,Transposon Coord,Plate Col,PlateRow,Plate,Well,' \
	+ 'Multiple Occupancy,Well Occupants,Locatability,Distance From Translation Start,' \
	+ 'Fractional Distance From Translation Start,Purifiability,Minimum Colonies To Pick\n'

	outputHandle = open(outputFile, 'w')
	outputHandle.write(headerStr)

	i = 0
	while i < len(featuresAndBestWellArray):
		
		
		
		featureName = featuresAndBestWellArray[i][0].featureName
		locusTag = featuresAndBestWellArray[i][0].tagDict['locus_tag'][0]
		transposonCoord = featuresAndBestWellArray[i][1]['sudokuGridEntry']['readAlignmentCoord']
	
		plateCol = \
		featuresAndBestWellArray[i][1]['sudokuGridEntry']['addressDict']['progenitor']['plateCol']
	
		plateRow = \
		featuresAndBestWellArray[i][1]['sudokuGridEntry']['addressDict']['progenitor']['plateRow']
		col = featuresAndBestWellArray[i][1]['sudokuGridEntry']['addressDict']['progenitor']['col']
		row = featuresAndBestWellArray[i][1]['sudokuGridEntry']['addressDict']['progenitor']['row']
	
		plateName = \
		featuresAndBestWellArray[i][1]['sudokuGridEntry']['addressDict']['progenitor']['plateName']
	
		multipleOccupancy = featuresAndBestWellArray[i][1]['sudokuGridEntry']['multipleOccupancy']
		distFromTranslationStart = \
		featuresAndBestWellArray[i][1]['sudokuGridEntry']['distFromTranslationStart']
	
		fracDistFromTranslationStart = \
		featuresAndBestWellArray[i][1]['sudokuGridEntry']['fracDistFromTranslationStart']
	
		locatability = featuresAndBestWellArray[i][1]['sudokuGridEntry']['locatability']
		readAlignmentCoord = featuresAndBestWellArray[i][1]['sudokuGridEntry']['readAlignmentCoord']
		wellOccupants = featuresAndBestWellArray[i][1]['sudokuGridEntry']['wellOccupants']
		
		minimumColoniesToPick = featuresAndBestWellArray[i][1]['minimumColoniesToPick']
		
		purifiability = featuresAndBestWellArray[i][1]['purifiability']
		
	
		outputStr = str(featureName) + ','
		outputStr += str(locusTag) + ','
		outputStr += str(transposonCoord) + ','
		outputStr += str(plateCol) + ','
		outputStr += str(plateRow) + ','
		outputStr += str(plateName) + ','
		outputStr += str(row + col.zfill(2)) + ','
		outputStr += str(multipleOccupancy) + ','
		outputStr += str(wellOccupants) + ','
		outputStr += str(locatability) + ','
		outputStr += str(distFromTranslationStart) + ','
		outputStr += str(fracDistFromTranslationStart) + ','
		outputStr += str(purifiability) + ','
		outputStr += str(minimumColoniesToPick)
		outputStr += '\n'

		outputHandle.write(outputStr)
	
		i += 1

	outputHandle.close()

# ------------------------------------------------------------------------------------------------ #
