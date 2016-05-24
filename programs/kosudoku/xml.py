# xml.py


# ------------------------------------------------------------------------------------------------ #
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
class SudokuGridXMLFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ParsePoolListFromXML(inputPoolList, delimeter=','):
	outputList = inputPoolList.strip().split(delimeter)
	
	return outputList
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CheckPoolDefinitionsInXMLFile(prPoolElements, pcPoolElements, rowPoolElements, colPoolElements):
	
	if len(prPoolElements) != 1:
		ex = SudokuGridXMLFailure('Number of prPool entries wrong')
		raise ex
	else:
		prPools = ParsePoolListFromXML(prPoolElements[0].text)

	if len(pcPoolElements) != 1:
		ex = SudokuGridXMLFailure('Number of pcPool entries wrong')
		raise ex
	else:
		pcPools = ParsePoolListFromXML(pcPoolElements[0].text)


	if len(rowPoolElements) != 1:
		ex = SudokuGridXMLFailure('Number of rowPool entries wrong')
		raise ex
	else:
		rowPools = ParsePoolListFromXML(rowPoolElements[0].text)


	if len(colPoolElements) != 1:
		ex = SudokuGridXMLFailure('Number of colPool entries wrong')
		raise ex
	else:
		colPools = ParsePoolListFromXML(colPoolElements[0].text)

	return prPools, pcPools, rowPools, colPools
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportSudokuCoordFromXML(xmlCoordObject):
	
	from .grid import SudokuGenomicCoord
	import pdb
	
	coord = float(xmlCoordObject.attrib['coord'])
	locatability = str(xmlCoordObject.attrib['locatability'])
	readCount = float(xmlCoordObject.attrib['readCount'])

	try:
		locatabilityScore = float(xmlCoordObject.attrib['locatabilityScore'])
	except:
		locatabilityScore = str(xmlCoordObject.attrib['locatabilityScore'])

	try:
		featureName = str(xmlCoordObject.attrib['featureName'])
		featureName = featureName.split('//')
	except:
		featureName = []
	
	try:
		dists = str(xmlCoordObject.attrib['distanceFromFeatureTranslationStart'])
		dists = dists.split(',')
		for dist in dists:
			dist = int(dist)
	except:
		dists = []
	
	try:
		fracDists = str(xmlCoordObject.attrib['fracDistanceFromFeatureTranslationStart'])
		fracDists = fracDists.split(',')
		for dist in fracDists:
			dist = float(dist)
	except:
		fracDists = []
	
	
	try:
		locusTags = str(xmlCoordObject.attrib['locusTag'])
		locusTags = locusTags.split(',')
	except:
		locusTags = []
	
	
	gCoord = SudokuGenomicCoord(int(coord), locatability, locatabilityScore, \
	int(readCount))
	
	gCoord.featureName = featureName
	gCoord.distanceFromFeatureTranslationStart = dists
	gCoord.fracDistanceFromFeatureTranslationStart = fracDists
	gCoord.locusTag = locusTags
	
	return gCoord
# ------------------------------------------------------------------------------------------------ #
	

# ------------------------------------------------------------------------------------------------ #
def ImportSudokuWellFromXML(xmlWellObject):
	
	from .grid import SudokuWell
	
	plateName = xmlWellObject.attrib['plateName']
	row = xmlWellObject.attrib['row']
	col = xmlWellObject.attrib['col']
	libraryAddress = xmlWellObject.attrib['libraryAddress']
	od = float(xmlWellObject.attrib['od'])
	plateRow = xmlWellObject.attrib['plateRow']
	plateCol = xmlWellObject.attrib['plateCol']
	
	sudokuWell = SudokuWell(plateName, plateRow, plateCol, row, col, OD=od)
	sudokuWell.libraryAddress = libraryAddress
	
	gCoords = xmlWellObject.findall('genomicCoord')
	
	for gCoord in gCoords:
		sudokuCoord = ImportSudokuCoordFromXML(gCoord)
		sudokuWell.readAlignmentCoords.append(sudokuCoord)
	
	addressSystems = xmlWellObject.findall('addressSystem')
	
	for addressSystem in addressSystems:
		name = addressSystem.attrib['name']
		coordsElements = addressSystem.findall('coords')
		if len(coordsElements) != 1:
			ex = SudokuGridXMLFailure('Too many coords!')
			raise ex
		else:
			coordsDict = {}
			for key in coordsElements[0].attrib.keys():
				coordsDict[key] = coordsElements[0].attrib[key]
			sudokuWell.addressDict[name] = coordsDict
	
	return sudokuWell
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportSudokuColonyPurifiedWellFromXML(xmlWellObject):
	
	from .isitinthere import SudokuColonyPurifiedWell
	from ast import literal_eval
	
	
	plateName = xmlWellObject.attrib['plateName']
	row = xmlWellObject.attrib['row']
	col = xmlWellObject.attrib['col']
	libraryAddress = xmlWellObject.attrib['libraryAddress']
	od = float(xmlWellObject.attrib['od'])
	plateRow = xmlWellObject.attrib['plateRow']
	plateCol = xmlWellObject.attrib['plateCol']
	
	# Attributes specific to SudokuColonyPurifiedWell
	condensationType = str(xmlWellObject.attrib['condensationType'])
	hopedForCoord = int(xmlWellObject.attrib['hopedForCoord'])
	hopedForPresent = literal_eval(xmlWellObject.attrib['hopedForPresent'])
	
	sudokuWell = SudokuColonyPurifiedWell(plateName, plateRow, plateCol, row, col, OD=od)
	sudokuWell.libraryAddress = libraryAddress
	sudokuWell.condensationType = condensationType
	sudokuWell.hopedForCoord = hopedForCoord
	sudokuWell.hopedForPresent = hopedForPresent
	
	gCoords = xmlWellObject.findall('genomicCoord')
	
	for gCoord in gCoords:
		sudokuCoord = ImportSudokuCoordFromXML(gCoord)
		sudokuWell.readAlignmentCoords.append(sudokuCoord)
	
	addressSystems = xmlWellObject.findall('addressSystem')
	
	for addressSystem in addressSystems:
		name = addressSystem.attrib['name']
		coordsElements = addressSystem.findall('coords')
		if len(coordsElements) != 1:
			ex = SudokuGridXMLFailure('Too many coords!')
			raise ex
		else:
			coordsDict = {}
			for key in coordsElements[0].attrib.keys():
				coordsDict[key] = coordsElements[0].attrib[key]
			sudokuWell.addressDict[name] = coordsDict
	
	return sudokuWell
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ImportSudokuPlateFromXML(plateElement, expectedRows=[], expectedCols=[], \
expectedPlateName=None, expectedPlateRow=None, expectedPlateCol=None, \
useSudokuColonyPurifiedWells=False):

	import pdb
	from .grid import SudokuPlate
	
	if expectedRows == []:
		rows = plateElement.attrib['rows'].split(',')
	else:
		rows = expectedRows
		
	if expectedCols == []:
		cols = plateElement.attrib['cols'].split(',')
	else:
		cols = expectedCols
	
	if expectedPlateName == None:
		plateName = plateElement.attrib['plateName']
	else:
		plateName = expectedPlateName
		
	if expectedPlateRow == None:
		plateRow = plateElement.attrib['plateRow']
	else:
		plateRow = expectedPlateRow
		
	if expectedPlateCol == None:
		plateCol = plateElement.attrib['plateCol']
	else:
		plateCol = expectedPlateCol
	
	
	sudokuPlate = SudokuPlate(plateName, plateRow, plateCol, rows, cols)
	
	rowElements = plateElement.findall('row')
	
	i = 0
	while i < len(rowElements):
		colElements = rowElements[i].findall('col')
		currentRow = rowElements[i].attrib['row']
			
		for colElement in colElements:
		
			currentCol = colElement.attrib['col']
# 			print(plateName + '_' + currentRow + currentCol)
			wellElements = colElement.findall('well')
			
			if len(wellElements) != 1:
				ex = SudokuGridXMLFailure('Number of well entries is wrong!')
				raise ex
			else:
				if useSudokuColonyPurifiedWells == True:
					well = ImportSudokuColonyPurifiedWellFromXML(wellElements[0])
				else:
					well = ImportSudokuWellFromXML(wellElements[0])
				
				sudokuPlate.wellGrid[currentRow][currentCol] = well
		i += 1
		
	return sudokuPlate
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ImportSudokuGridFromXML(fileName, useSudokuColonyPurifiedWells=False):
	
	import xml.etree.ElementTree as ET
	from .grid import InitializeEmptySudokuGridLookupDict
	
	tree = ET.parse(fileName)
	root = tree.getroot()

	tempSudokuGrid = {}

	# Import the grid layout, including plate rows, plate columns, rows and columns

	prPoolElements = root.findall('prPools')
	pcPoolElements = root.findall('pcPools')
	rowPoolElements = root.findall('rowPools')
	colPoolElements = root.findall('colPools')

	prPools, pcPools, rowPools, colPools = \
	CheckPoolDefinitionsInXMLFile(prPoolElements, pcPoolElements, rowPoolElements, colPoolElements)

	sudokuGridLookupDict = InitializeEmptySudokuGridLookupDict(prPools, pcPools, rowPools, colPools)
	
	# Start at the level of plate rows

	plateRowElements = root.findall('PR')

	for pr in plateRowElements:
		plateRow = str(pr.attrib['plateRow'])
		plateColElements = pr.findall('PC')
	
		for pc in plateColElements:
			plateCol = str(pc.attrib['plateCol'])
		
			plateElements = pc.findall('plate')
			
			if len(plateElements) != 1:
				ex = SudokuGridXMLFailure('Number of plate entries is greater than 1')
				raise ex
			else:
				plate = ImportSudokuPlateFromXML(plateElements[0], \
				useSudokuColonyPurifiedWells=useSudokuColonyPurifiedWells)
				
				sudokuGridLookupDict[plateRow][plateCol] = plate
	
	
	return sudokuGridLookupDict, prPools, pcPools, rowPools, colPools
# ------------------------------------------------------------------------------------------------ #
				
				
# ------------------------------------------------------------------------------------------------ #
def ExportListForXML(listToExport, delimeter=','):
	
	outputStr = ''
	
	i = 0
	while i < len(listToExport):
		outputStr += str(listToExport[i])
		if i < len(listToExport) - 1:
			outputStr += delimeter
		i += 1
	
	return outputStr
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def ExportSudokuGridToXML(sudokuGridLookupDict, plateRowPools, plateColPools, rowPools, colPools, \
fileName, useSudokuColonyPurifiedWells=False):

	import xml.etree.ElementTree as ET

	sudokuGridTreeRoot = ET.Element('grid')
	
	prPoolSubElement = ET.SubElement(sudokuGridTreeRoot, 'prPools')
	prPoolSubElement.text = ExportListForXML(plateRowPools)
	
	pcPoolSubElement = ET.SubElement(sudokuGridTreeRoot, 'pcPools')
	pcPoolSubElement.text = ExportListForXML(plateColPools)
	
	rowPoolSubElement = ET.SubElement(sudokuGridTreeRoot, 'rowPools')
	rowPoolSubElement.text = ExportListForXML(rowPools)
	
	colPoolSubElement = ET.SubElement(sudokuGridTreeRoot, 'colPools')
	colPoolSubElement.text = ExportListForXML(colPools)
	

	for prKey in plateRowPools:
	
		PRSubElement = ET.SubElement(sudokuGridTreeRoot, 'PR')
		PRSubElement.set('plateRow', prKey)
	
		pcPoolKeys = sorted(sudokuGridLookupDict[prKey].keys())
	
		for pcKey in plateColPools:
			PCSubElement = ET.SubElement(PRSubElement, 'PC')
			PCSubElement.set('plateCol', pcKey)
		
			sudokuPlate = sudokuGridLookupDict[prKey][pcKey]
			ExportSudokuPlateAsXMLSubElement(PCSubElement, sudokuPlate, rowPools=rowPools, \
			colPools=colPools, useSudokuColonyPurifiedWells=useSudokuColonyPurifiedWells)
						
		
	indent(sudokuGridTreeRoot)
	sudokuGridTree = ET.ElementTree(sudokuGridTreeRoot)
	sudokuGridTree.write(fileName)
	
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ExportSudokuPlateAsXMLSubElement(parent, plate, rowPools=[], colPools=[], \
useSudokuColonyPurifiedWells=False):
	
	import xml.etree.ElementTree as ET
	
	plateName = plate.plateName
		
	plateSubElement = ET.SubElement(parent, 'plate')
	plateSubElement.set('plateName', str(plateName))
	
	plateSubElement.set('plateRow', str(plate.plateRow))
	plateSubElement.set('plateCol', str(plate.plateCol))
	
	wellGrid = plate.wellGrid
	
	if len(rowPools) == 0:
		rows = sorted(wellGrid.keys())
	else:
		rows = rowPools
	
	if len(colPools) == 0:
		cols = wellGrid[0].keys()
	else:
		cols = colPools
	
	plateSubElement.set('rows', ExportListForXML(rows))
	plateSubElement.set('cols', ExportListForXML(cols))
	
	
	for row in rows:
		rowSubElement = ET.SubElement(plateSubElement, 'row')
		rowSubElement.set('row', str(row))
			
		for col in cols:
			colSubElement = ET.SubElement(rowSubElement, 'col')
			colSubElement.set('col', str(col))
			
			well = plate.wellGrid[row][col]
			
			if useSudokuColonyPurifiedWells == False:
				ExportSudokuWellAsXMLSubElement(colSubElement, well)
			else:
				ExportSudokuColonyPurifiedWellAsXMLSubElement(colSubElement, well)
	
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ExportSudokuColonyPurifiedWellAsXMLSubElement(parent, well):
	
	import xml.etree.ElementTree as ET
	
	plateName = well.plateName
	row = well.row
	col = well.col
	libraryAddress = well.libraryAddress
	OD = well.OD
	readAlignmentCoords = well.readAlignmentCoords
	plateRow = well.plateRow
	plateCol = well.plateCol
	addressDict = well.addressDict
	addressDictKeys = addressDict.keys()
	
	# Variables specific to SudokuColonyPurifiedWell
	
	hasPredictionForContents = well.hasPredictionForContents
	predictionsForContents = well.predictionsForContents
	predictionCorrect = well.predictionCorrect
	hopedForPresent = well.hopedForPresent
	hopedForCoord = well.hopedForCoord
	simplifiedReadAlignmentCoords = well.simplifiedReadAlignmentCoords
	simplifiedLikelyReadAlignmentCoords = well.simplifiedLikelyReadAlignmentCoords
	likelyReadAlignmentCoords = well.likelyReadAlignmentCoords
	progenitorContents = well.progenitorContents
	progenitorLocatabilities = well.progenitorLocatabilities
	condensationType = well.condensationType
	

	wellSubElement = ET.SubElement(parent, 'well')
	
	wellSubElement.set('plateName', str(plateName))
	wellSubElement.set('row', str(row))
	wellSubElement.set('col', str(col))
	wellSubElement.set('libraryAddress', str(libraryAddress))
	wellSubElement.set('od', str(OD))
	wellSubElement.set('plateRow', str(plateRow))
	wellSubElement.set('plateCol', str(plateCol))
	
	# Sub elements specific to SudokuColonyPurifiedWell
	wellSubElement.set('condensationType', str(condensationType))
	wellSubElement.set('hopedForCoord', str(hopedForCoord))
	wellSubElement.set('hopedForPresent', str(hopedForPresent))
	
	
	for key in addressDictKeys:
		addressSubElement = ET.SubElement(wellSubElement, 'addressSystem')
		addressSubElement.set('name', str(key))
		coordSubElement = ET.SubElement(addressSubElement, 'coords')
		for coordKey in addressDict[key].keys():
			coordSubElement.set(coordKey, str(addressDict[key][coordKey]))

	for gCoord in readAlignmentCoords:
		ExportSudokuCoordAsXMLSubElement(wellSubElement, gCoord)	
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ExportSudokuWellAsXMLSubElement(parent, well):
	
	import xml.etree.ElementTree as ET
	
	plateName = well.plateName
	row = well.row
	col = well.col
	libraryAddress = well.libraryAddress
	OD = well.OD
	readAlignmentCoords = well.readAlignmentCoords
	plateRow = well.plateRow
	plateCol = well.plateCol
	addressDict = well.addressDict
	addressDictKeys = addressDict.keys()
	
	wellSubElement = ET.SubElement(parent, 'well')
	
	wellSubElement.set('plateName', str(plateName))
	wellSubElement.set('row', str(row))
	wellSubElement.set('col', str(col))
	wellSubElement.set('libraryAddress', str(libraryAddress))
	wellSubElement.set('od', str(OD))
	wellSubElement.set('plateRow', str(plateRow))
	wellSubElement.set('plateCol', str(plateCol))
	
	for key in addressDictKeys:
		addressSubElement = ET.SubElement(wellSubElement, 'addressSystem')
		addressSubElement.set('name', str(key))
		coordSubElement = ET.SubElement(addressSubElement, 'coords')
		for coordKey in addressDict[key].keys():
			coordSubElement.set(coordKey, str(addressDict[key][coordKey]))

	for gCoord in readAlignmentCoords:
		ExportSudokuCoordAsXMLSubElement(wellSubElement, gCoord)	
	
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ExportSudokuCoordAsXMLSubElement(parent, gCoord):
	
	import xml.etree.ElementTree as ET
	
	gCoordSubElement = ET.SubElement(parent, 'genomicCoord')
		
	gCoordSubElement.set('coord', str(gCoord.coord))
	gCoordSubElement.set('locatability', str(gCoord.locatability))
	gCoordSubElement.set('locatabilityScore', str(gCoord.locatabilityScore))
	gCoordSubElement.set('readCount', str(gCoord.readCount))
	
	gCoordSubElement.set('featureName', ExportListForXML(gCoord.featureName, delimeter=','))
	
	gCoordSubElement.set('distanceFromFeatureTranslationStart', \
	ExportListForXML(gCoord.distanceFromFeatureTranslationStart, delimeter=','))
	
	gCoordSubElement.set('fracDistanceFromFeatureTranslationStart', \
	ExportListForXML(gCoord.fracDistanceFromFeatureTranslationStart, delimeter=','))
	
	gCoordSubElement.set('locusTag', ExportListForXML(gCoord.locusTag, delimeter=','))

	
	return
# ------------------------------------------------------------------------------------------------ #












# ------------------------------------------------------------------------------------------------ #
def ExportFeatureTagDictAsXMLSubElement(parent, tagDict):
	
	import xml.etree.ElementTree as ET
	
	tagDictKeys  = tagDict.keys()

	tagDictSubElement = ET.SubElement(parent, 'tagDict')

	for key in tagDictKeys:
		tagSubElement = ET.SubElement(tagDictSubElement, key)
		tagSubElement.text = ExportListForXML(tagDict[key], delimeter=',')
			
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ExportGeneFeatureArrayToXML(geneFeatureArray, fileName):
	
	import xml.etree.ElementTree as ET
	
	geneFeatureTreeRoot = ET.Element('genes')
	
	i = 0
	
	while i < len(geneFeatureArray):
		
		feature = geneFeatureArray[i]
		
		coordinates = feature.coordinates
		featureType = feature.featureType
		featureName = feature.featureName
		
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		startTranslation = feature.startTranslation
		endTranslation = feature.endTranslation
		
		tagDict = feature.tagDict
		

		featureSubElement = ET.SubElement(geneFeatureTreeRoot, 'gene')
		
		featureSubElement.set('coordinates', coordinates)
		featureSubElement.set('featureType', featureType)
		featureSubElement.set('featureName', featureName)
		featureSubElement.set('startCoord', str(startCoord))
		featureSubElement.set('endCoord', str(endCoord))
		featureSubElement.set('startTranslation', str(startTranslation))
		featureSubElement.set('endTranslation', str(endTranslation))
		
		for well in feature.sudokuGridEntries:
			ExportSudokuWellAsXMLSubElement(featureSubElement, well)
			
		for well in feature.rearrayedGridEntries:
			ExportSudokuRearrayedWellAsXMLSubElement(featureSubElement, well)
			
		ExportFeatureTagDictAsXMLSubElement(featureSubElement, feature.tagDict)
	
		i += 1
	
	indent(geneFeatureTreeRoot)
	geneFeatureTree = ET.ElementTree(geneFeatureTreeRoot)
	geneFeatureTree.write(fileName)
		
	
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportTagDictFromXML(xmlTagDicts):
	
	import pdb
	
	tagDict = {}
	
	for dict in xmlTagDicts:
		for child in dict:
			tagDictKeys = tagDict.keys()
			if child.tag not in tagDictKeys:
				tagDict[child.tag] = [child.text]
			elif child.tag in tagDictKeys:
				tagDict[child.tag].append(child.text)
	
# 	pdb.set_trace()	

	return tagDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportFeatureArrayFromXML(fileName):
	
	import xml.etree.ElementTree as ET
	
	tree = ET.parse(fileName)
	root = tree.getroot()
	
	geneFeatureArray = []
	
	features = root.findall('gene')
	
	for feature in features:
		
		coordinates = feature.attrib['coordinates']
		featureType = feature.attrib['featureType']

		geneFeature = Feature3(featureType, coordinates)
		
		geneFeature.featureName = feature.attrib['featureName']
		geneFeature.startCoord = int(feature.attrib['startCoord'])
		geneFeature.endCoord = int(feature.attrib['endCoord'])
		geneFeature.startTranslation = int(feature.attrib['startTranslation'])
		geneFeature.endTranslation = int(feature.attrib['endTranslation'])
		
		sudokuWells = feature.findall('well')
		
		for well in sudokuWells:
			sudokuWell = ImportSudokuWellFromXML(well)
			geneFeature.sudokuGridEntries.append(sudokuWell)
			
		sudokuRearrayedWells = feature.findall('rearrayedWell')
		
		for well in sudokuRearrayedWells:
			sudokuRearrayedWell = ImportSudokuRearrayedWellFromXML(well)
			geneFeature.rearrayedGridEntries.append(sudokuRearrayedWell)
		
		
		tagDicts = feature.findall('tagDict')
		
		tagDict = ImportTagDictFromXML(tagDicts)
		
		geneFeature.tagDict = tagDict
		
		
		geneFeatureArray.append(geneFeature)
		
	return geneFeatureArray
# ------------------------------------------------------------------------------------------------ #