# ---------------------------------------------------------------------------- #
def generateOutputMatrix(vectorList, delimeter='\t'):
	
	outputMatrix = []
	outputString = ''
	i = 0
	
	
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def generateOutputMatrixWithHeaders(vectorList, headers, delimeter='\t'):
	
	if len(headers) != len(vectorList):
		print("Header array length != vector list length")
		return
	
	outputMatrix = []
	
	k = 0
	headerString = ''
	while k < len(headers):
		headerString += str(headers[k]) + delimeter
		k += 1
	
	headerString += '\n'
	outputMatrix.append(headerString)
	
	outputString = ''
	i = 0
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def writeOutputMatrix(fileName, outputMatrix):
	outputFile = open(fileName, 'w')
	i = 0
	while i < len(outputMatrix):
		outputFile.write(outputMatrix[i])
		i+=1
	return	
# ---------------------------------------------------------------------------- #
