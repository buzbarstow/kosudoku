
# ------------------------------------------------------------------------------------------------ #
def CalculatePlotAndSavePoolPresenceTableTotalReadCountHistogram(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, cumulativePlotFileName, cumulativePlotFileNameNC, \
threshold):
# Generate, plot and save the pool presence table total read count per line histogram

	import numpy
	import pdb

	[totalReadCountPerLineArray, totalReadCountPerLineNoControlsArray, \
	totalRowPoolCountArray, totalColPoolCountArray, totalPRPoolCountArray, totalPCPoolCountArray, \
	rowPoolCountArray, colPoolCountArray, prPoolCountArray, pcPoolCountArray] = \
	CalculatePoolPresenceTableReadCountHistograms(\
	poolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools, threshold)
	

	
	# Make and Plot histograms
	
	MakeAndPlotHistogramOfReadCounts(totalReadCountPerLineArray, plotCumulative=True, \
	plotHistogram=True, \
	title="Distribution of Total Read Count", \
	xlabel="Total Read Count Per Line", ylabel="Lines", \
	cumTitle="Cumulative Distribution of Total Read Count", \
	cumXlabel="Total Read Count for Line", \
	cumYlabel="Fraction of All Non-Zero Pool Presence Table Entries")
	
	
	MakeAndPlotHistogramOfReadCounts(totalReadCountPerLineNoControlsArray, plotCumulative=True, \
	plotHistogram=True, \
	title="Distribution of Total Non Control Read Count for All Pool Presence Table Entries", \
	xlabel="Total Read Count Per Line", ylabel="Lines", \
	cumTitle="Cumulative Distribution of Total Non Control Read Count", \
	cumXlabel="Total Read Count for Line", \
	cumYlabel="Fraction of All Non-Zero Pool Presence Table Entries")
	
	MakeAndPlotHistogramOfReadCounts(totalRowPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Total Row Pool Read Count", \
	xlabel="Total Row Pool Read Count Per Line", ylabel="Lines")

	MakeAndPlotHistogramOfReadCounts(totalColPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Total Col Pool Read Count", \
	xlabel="Total Col Pool Read Count Per Line", ylabel="Lines")

	MakeAndPlotHistogramOfReadCounts(totalPRPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Total PR Pool Read Count", \
	xlabel="Total PR Pool Read Count Per Line", ylabel="Lines")


	MakeAndPlotHistogramOfReadCounts(totalPCPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Total PC Pool Read Count", \
	xlabel="Total PC Pool Read Count Per Line", ylabel="Lines")
	
	
	MakeAndPlotHistogramOfReadCounts(rowPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Row Pool Read Count", \
	xlabel="Row Pool Read Count Per Line", ylabel="Lines")
	
	MakeAndPlotHistogramOfReadCounts(colPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of Col Pool Read Count", \
	xlabel="Col Pool Read Count Per Line", ylabel="Lines")
	
	MakeAndPlotHistogramOfReadCounts(prPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of PR Pool Read Count", \
	xlabel="PR Pool Read Count Per Line", ylabel="Lines")
	
	MakeAndPlotHistogramOfReadCounts(pcPoolCountArray, plotCumulative=False, \
	plotHistogram=True, \
	title="Distribution of PC Pool Read Count", \
	xlabel="PC Pool Read Count Per Line", ylabel="Lines")
	
	
	
	
	# Cumulative density plots 

# 	bins = numpy.arange(0, max(totalReadCountPerLineArray), 1)
# 	binsNC = numpy.arange(0, max(totalReadCountPerLineNoControlsArray), 1)
# 
# 	values, base = numpy.histogram(totalReadCountPerLineArray, bins=bins)
# 	valuesNC, baseNC = numpy.histogram(totalReadCountPerLineNoControlsArray, bins=binsNC)
# 
# 	cumulative = numpy.cumsum(values)
# 	cumulativeNC = numpy.cumsum(valuesNC)
# 	
# 	PlotCumulativeDensityOfTotalReadCounts(base, cumulative, baseNC, cumulativeNC)
# 	
# 	# Save cumulative density plot data
# 	
# 	vectorList = [base[:-1], cumulative, cumulative/max(cumulative)]
# 	vectorListNC = [baseNC[:-1], cumulativeNC, cumulativeNC/max(cumulativeNC)]
# 	
# 	headers = ["Total Reads", "Cum. Lines", "Frac Cum. Lines"]
# 	headersNC = ["Total Reads NC", "Cum. Lines NC", "Frac Cum. Lines NC"]
# 	
# # 	pdb.set_trace()
# 
# 	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
# 	oMatrixNC = generateOutputMatrixWithHeaders(vectorListNC, headersNC)
# 	
# 	writeOutputMatrix(cumulativePlotFileName, oMatrix)
# 	writeOutputMatrix(cumulativePlotFileNameNC, oMatrixNC)
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateCoordNumbersForAmbiguousAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):

	import pdb
	
	# Initialize the number arrays
	
	nRowCoordsArray = []
	nColCoordsArray = []
	nPRCoordsArray = []
	nPCCoordsArray = []
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		# Figure out if the line is ambiguous or not
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
		
		# Here's the test for ambiguous mapping to multiple library addresses
		
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			
			[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
			CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)
			
			nRowCoordsArray.append(nRowCoords)
			nColCoordsArray.append(nColCoords)
			nPRCoordsArray.append(nPRCoords)
			nPCCoordsArray.append(nPCCoords)
		
		i += 1

	return [nRowCoordsArray, nColCoordsArray, nPRCoordsArray, nPCCoordsArray]

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, threshold):

	import pdb
	
	# Initialize the number arrays
	
	nPossibleAddressesArray = []
	nTotalCoordsArray = []
	maxNumberCoordsArray = []
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		# Figure out if the line is ambiguous or not
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
		
		# Here's the test for ambiguous mapping to multiple library addresses
		
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			
			nPossibleAddressesArray.append(lenPossibleAddresses)
			
			[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
			CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)
			
			nTotalCoordsArray.append(nRowCoords + nColCoords + nPRCoords + nPCCoords)
			
			maxNumberCoordsArray.append(max([nRowCoords, nColCoords, nPRCoords, nPCCoords]))
		
		i += 1

	return [nPossibleAddressesArray, nTotalCoordsArray, maxNumberCoordsArray]

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateAndPlotPossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, threshold):

	import matplotlib.pyplot as pyplot
	from numpy import arange, histogram, cumsum

	[nPossibleAddressesArray, nTotalCoordsArray, maxNumberCoordsArray] = \
	CalculatePossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
	colPools, prPools, pcPools, controlPools, threshold)
	
	nPossibleAddressesBins = arange(0,max(nPossibleAddressesArray),1)
	nTotalCoordsBins = arange(0,max(nTotalCoordsArray),1)
	maxNumberCoordsBins = arange(0,max(maxNumberCoordsArray),1)


	pyplot.figure()
	pyplot.hist(nPossibleAddressesArray, bins=nPossibleAddressesBins)
	pyplot.xlabel("Number Possible Addresses")
	pyplot.ylabel("Number of Ambiguous Lines")

	pyplot.figure()
	pyplot.hist(nTotalCoordsArray, bins=nTotalCoordsBins)
	pyplot.xlabel("Total Number of Coords")
	pyplot.ylabel("Number of Ambiguous Lines")
	
	
	pyplot.figure()
	pyplot.hist(maxNumberCoordsArray, bins=maxNumberCoordsBins)
	pyplot.xlabel("Maximum Number of Coords in a Pool Axis")
	pyplot.ylabel("Number of Ambiguous Lines")

	pyplot.show()
	
	
	binsPossAdd = arange(0, max(nPossibleAddressesArray), 1)
	valuesPossAdd, basePossAdd = histogram(nPossibleAddressesArray, bins=binsPossAdd)
	cumulativePossAdd = cumsum(valuesPossAdd)
	
	binsTotalCoords = arange(0, max(nTotalCoordsArray), 1)
	valuesTotalCoords, baseTotalCoords = histogram(nTotalCoordsArray, bins=binsTotalCoords)
	cumulativeTotalCoords = cumsum(valuesTotalCoords)
	
	binsMaxCoords = arange(0, max(maxNumberCoordsArray), 1)
	valuesMaxCoords, baseMaxCoords = histogram(maxNumberCoordsArray, bins=binsMaxCoords)
	cumulativeMaxCoords = cumsum(valuesMaxCoords)

	pyplot.figure()
	pyplot.plot(basePossAdd[:-1], cumulativePossAdd/max(cumulativePossAdd))
	pyplot.title("Cumulative Distribution of Possible Addresses for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Possible Addresses")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	
	pyplot.figure()
	pyplot.plot(baseTotalCoords[:-1], cumulativeTotalCoords/max(cumulativeTotalCoords))
	pyplot.title("Cumulative Distribution of Total Coords for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Total Coordinates")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	
	pyplot.figure()
	pyplot.plot(baseMaxCoords[:-1], cumulativeMaxCoords/max(cumulativeMaxCoords))
	pyplot.title("Cumulative Distribution of Maximum Pool Coordinates for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Maximum Number of Coordinates in a Pool")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	

	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateReadCountHistograms(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):
	
	import pdb
	
	# Total number of reads per line, with and without controls
	totalReadCountPerLineArray = []
	totalReadCountPerLineNoControlsArray = []
	
	# Number of counts in each row, col, pr and pc coord pool per line
	rowPoolCountArray = []
	colPoolCountArray = []
	prPoolCountArray = []
	pcPoolCountArray = []
	
	# Total number of reads in row, col, pr and pc pools per line
	totalRowPoolCountArray = []
	totalColPoolCountArray = []
	totalPRPoolCountArray = []
	totalPCPoolCountArray = []
	
	
	i = 0
	while i < len(poolPresenceTable):
	
		lineEntries = []

		poolPresenceTableLine = poolPresenceTable[i]
		
		rowEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, rowPools, threshold=threshold)
		colEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, colPools, threshold=threshold)
		prEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, prPools, threshold=threshold)
		pcEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, pcPools, threshold=threshold)
		
# 		pdb.set_trace()
		
		if controlPools != None:
			controlEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, controlPools)
		else:
			controlEntries = []
		
		
		totalRowCounts = 0
		totalColCounts = 0
		totalPRCounts = 0
		totalPCCounts = 0
		
		totalLineCounts = 0
		totalLineCountsNC = 0
		
		for entry in rowEntries:
			rowPoolCountArray.append(entry)
			totalRowCounts += entry
			totalLineCounts += entry
			totalLineCountsNC += entry
		totalRowPoolCountArray.append(totalRowCounts)
		
		for entry in colEntries:
			colPoolCountArray.append(entry)
			totalColCounts += entry
			totalLineCounts += entry
			totalLineCountsNC += entry
		totalColPoolCountArray.append(totalColCounts)
		
		for entry in prEntries:
			prPoolCountArray.append(entry)
			totalPRCounts += entry
			totalLineCounts += entry
			totalLineCountsNC += entry
		totalPRPoolCountArray.append(totalPRCounts)
		
		for entry in pcEntries:
			pcPoolCountArray.append(entry)
			totalPCCounts += entry
			totalLineCounts += entry
			totalLineCountsNC += entry
		totalPCPoolCountArray.append(totalPCCounts)
		
		for entry in controlEntries:
			totalLineCounts += entry
		
		
		totalReadCountPerLineArray.append(totalLineCounts)
		totalReadCountPerLineNoControlsArray.append(totalLineCountsNC)
		
		i += 1

	return [totalReadCountPerLineArray, totalReadCountPerLineNoControlsArray, \
	totalRowPoolCountArray, totalColPoolCountArray, totalPRPoolCountArray, totalPCPoolCountArray, \
	rowPoolCountArray, colPoolCountArray, prPoolCountArray, pcPoolCountArray]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindEntriesInPoolPresenceLine(poolPresenceTableLine, poolsList, threshold=1):
	
	entries = []
	
	j = 0
	while j < len(poolsList):
		if 	poolPresenceTableLine[poolsList[j]] >= threshold:
			entries.append(poolPresenceTableLine[poolsList[j]])
		j += 1

	return entries
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def gaussian(x, sigma, mu):
	from numpy import pi,sqrt, exp
	
	p = (1./(sigma*sqrt(2.*pi)))*exp(-1.*(x-mu)**2/(2.*sigma**2))
	
	return p
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Functions for performing Gaussian fit to data
def peval(x, p): 
	fit = p[0]*(exp( - ((x-p[1])**2) / p[2])) \
		+ p[3]*(exp( - ((x-p[4])**2) / p[5])) \
		+ p[6]*(exp( - ((x-p[7])**2) / p[8])) \
		+ p[9]*(exp( - ((x-p[10])**2) / p[11]))  
	return fit 

def residuals(p, y, x): 
	err = y-peval(x,p) 
	return err

def fitgaussian(lowerLimit, upperLimit, data, datalabel, g):
	import scipy
	from scipy.optimize import leastsq
	
	y = data[lowerLimit:upperLimit,1]
	x = data[lowerLimit:upperLimit,0]
	
	g.plot(data)
	g.set_range('yrange',(0,100))
	g.replot()
	
	A1_0=90
	A2_0=525
	A3_0=20
	A4_0=75
	A5_0=530
	A6_0=15
	A7_0=40
	A8_0=545
	A9_0=25
	A10_0=40
	A11_0=560
	A12_0=100
	
	p0 = array([A1_0, A2_0, A3_0, A4_0, A5_0, A6_0,A7_0,A8_0,A9_0,A10_0,A11_0,A12_0])
	plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000)
	
	fit = []
	
	for i in  x:
		fit.append((i,peval(i,plsq[0])))
	
	g.replot(fit)
	
	print("Final parameters for " + datalabel)
	i = 0
	while i < len(p0):
		print("%s = %.4f " % (i, plsq[0][i]))
		i = i+1
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def MakeAndPlotReadCountHistogram(readCountArray, plotCumulative=False, plotHistogram=False, \
title='', xlabel='', ylabel='', cumTitle='', cumXlabel='', cumYlabel='', binWidth=1):
	
	import numpy
	import matplotlib.pyplot as plt
	
	bins = numpy.arange(0, max(readCountArray)+binWidth, binWidth)
	values, base = numpy.histogram(readCountArray, bins=bins)
	cumulative = numpy.cumsum(values)
	
	binCenters = ConvertBinEdgesToCenters(base)	
	
	if plotCumulative:
		plt.figure()
		plt.semilogx(base[:-1],cumulative/max(cumulative))
		plt.title(cumTitle)
		plt.xlabel(cumXlabel)
		plt.ylabel(cumYlabel) 
		plt.show()
	
	if plotHistogram:
		plt.figure()
		plt.plot(binCenters, values)
		plt.title(title)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel) 
		plt.show()

	
	return [bins, values, base, cumulative]

# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateScoreHistogramForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict,\
threshold, cumulativeScoreDistributionFileName):
	
	from numpy import histogram, log, arange, cumsum
	from matplotlib.pyplot import hist
	import matplotlib.pyplot as plt
	import pdb
	
	i = 0
	scoreArray = []
	
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
		
		possibleAddressesAndScores = \
		CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, rowPools, \
		colPools, prPools, pcPools, logReadNumberRatioHistogramFitDict, \
		logReadNumberRatioHistogramIntegralDict, threshold)
		
		# Calculate if the line is single occupancy
		lenPossibleAddresses = len(possibleAddressesAndScores)
	
		if lenPossibleAddresses == 1:	
			scoreArray.append(possibleAddressesAndScores[0][2])

		i += 1
		
	logScoreArray = log(scoreArray)
# 	pdb.set_trace()
	
	logBins = arange(min(logScoreArray),max(scoreArray),0.1)
	bins = arange(0,max(scoreArray),0.001)
	
	values, base = histogram(logScoreArray, bins=logBins)
	cumulative = cumsum(values)
	
	plt.figure()
	hist(scoreArray, bins)
	plt.xlabel("Score")
	plt.ylabel("Number of Single Address Lines")
	
	plt.figure()
	plt.plot(base[:-1],cumulative/max(cumulative))
	plt.title("Cumulative Score Distribution of Single Address Pool Presence Table Entries")
	plt.xlabel("Natural Log of Score")
	plt.ylabel("Fraction of Single Address Pool Presence Table Entries") 
	plt.show()
	

	return [scoreArray, cumulative, base[:-1]]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def FitReadNumberRatiosAndPlot(logNr2ncArray, logBins, xAxisLabel, plsq0, integrate=False):

	import numpy
	import matplotlib.pyplot as plt
	import pdb
	from scipy.integrate import simps
	from kosudoku.pool import voigtFit

	valueslogNr2nc, baselogNr2nc = numpy.histogram(logNr2ncArray, bins=logBins)
	binCentersLogNr2nc = ConvertBinEdgesToCenters(baselogNr2nc)
	
	

	[plsqlogNr2nc, voigtlogNr2nc] = voigtFit(binCentersLogNr2nc, valueslogNr2nc, plsq0)
	
	plt.figure()
	
	
	plt.plot(binCentersLogNr2nc, voigtlogNr2nc, linewidth=1)
	plt.scatter(binCentersLogNr2nc, valueslogNr2nc)
	
	plt.xlabel(xAxisLabel)
	plt.ylabel("Number of Pool Presence Table Lines")
	
	plt.grid()
	
	plt.xlim(min(logBins),max(logBins))
	
# 	pdb.set_trace()
	
	plt.ylim(-100,2500)

	plt.show()
	
	if integrate == True:
		integral = simps(voigtlogNr2nc, binCentersLogNr2nc)
		return [valueslogNr2nc, baselogNr2nc, integral, plsqlogNr2nc, \
		voigtlogNr2nc, binCentersLogNr2nc]
	else:
		return [valueslogNr2nc, baselogNr2nc] 
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculateRatioOfReadNumbersForSingleAddressLines(poolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, threshold, logReadNumberRatioHistogramsFileName, \
fitFileName):
	
	from numpy import array, arange
	from kosudoku.output import generateOutputMatrixWithHeaders, writeOutputMatrix
	
	import pdb
	
	[nr2ncArray, nr2nprArray, nr2npcArray, nc2nprArray, nc2npcArray, npr2npcArray, \
	logNr2ncArray, logNr2nprArray, logNr2npcArray, logNc2nprArray, logNc2npcArray, \
	logNpr2npcArray] = \
	CalculateReadNumberRatiosForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
	pcPools, controlPools, threshold)

	a0 = 4000
	c0 = -0.5
	delta0 = 0.35
	sigma0 = 0.5

	plsq0 = array([a0, c0, delta0, sigma0])

	logBins = arange(-10,10.2,0.2)
	
	
	[valueslogNr2nc, baselogNr2nc, integralNr2nc, plsqNr2nc, voigtlogNr2nc, binCentersLogNr2nc] = \
	FitReadNumberRatiosAndPlot(logNr2ncArray, logBins, "log(n_row/n_col)", plsq0, integrate=True)

	[valueslogNr2npr, baselogNr2npr, integralNr2npr, plsqNr2npr, voigtlogNr2npr, \
	binCentersLogNr2npr] = \
	FitReadNumberRatiosAndPlot(logNr2nprArray, logBins, "log(n_row/n_pr)", plsq0, integrate=True)

	[valueslogNr2npc, baselogNr2npc, integralNr2npc, plsqNr2npc, voigtlogNr2npc, \
	binCentersLogNr2npc] = \
	FitReadNumberRatiosAndPlot(logNr2npcArray, logBins, "log(n_row/n_pc)", plsq0, integrate=True)
 
	[valueslogNc2npr, baselogNc2npr, integralNc2npr, plsqNc2npr, voigtlogNc2npr, \
	binCentersLogNc2npr] = \
	FitReadNumberRatiosAndPlot(logNc2nprArray, logBins, "log(n_col/n_pr)", plsq0, integrate=True)
 
	[valueslogNc2npc, baselogNc2npc, integralNc2npc, plsqNc2npc, voigtlogNc2npc, \
	binCentersLogNc2npc] = \
	FitReadNumberRatiosAndPlot(logNc2npcArray, logBins, "log(n_col/n_pc)", plsq0, integrate=True)

	[valueslogNpr2npc, baselogNpr2npc, integralNpr2npc, plsqNpr2npc, voigtlogNpr2npc, \
	binCentersLogNpr2npc] = \
	FitReadNumberRatiosAndPlot(logNpr2npcArray, logBins, "log(n_pr/n_pc)", plsq0, integrate=True)

	
	headers = [\
	'baselogNr2nc', 'valueslogNr2nc', \
	'baselogNr2npr', 'valueslogNr2npr', \
	'baselogNr2npc', 'valueslogNr2npc', \
	'baselogNc2npr', 'valueslogNc2npr', \
	'baselogNc2npc', 'valueslogNc2npc', \
	'baselogNpr2npc', 'valueslogNpr2npc' ]
	
	fitHeaders = [\
	'binCenterNr2nc', 'fitValueslogNr2nc', \
	'binCenterNr2npr', 'fitValueslogNr2npr', \
	'binCenterNr2npc', 'fitValueslogNr2npc', \
	'binCenterNc2npr', 'fitValueslogNc2npr', \
	'binCenterNc2npc', 'fitValueslogNc2npc', \
	'binCenterNpr2npc', 'fitValueslogNpr2npc' ]


	vectorList = [\
	binCentersLogNr2nc, valueslogNr2nc, \
	binCentersLogNr2npr, valueslogNr2npr, \
	binCentersLogNr2npc, valueslogNr2npc, \
	binCentersLogNc2npr, valueslogNc2npr, \
	binCentersLogNc2npc, valueslogNc2npc, \
	binCentersLogNpr2npc, valueslogNpr2npc]
	
	fitVectorList = [\
	binCentersLogNr2nc, voigtlogNr2nc, \
	binCentersLogNr2npr, voigtlogNr2npr, \
	binCentersLogNr2npc, voigtlogNr2npc, \
	binCentersLogNc2npr, voigtlogNc2npr, \
	binCentersLogNc2npc, voigtlogNc2npc, \
	binCentersLogNpr2npc, voigtlogNpr2npc ]
	
	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
	writeOutputMatrix(logReadNumberRatioHistogramsFileName, oMatrix)
	
	oMatrixFit = generateOutputMatrixWithHeaders(fitVectorList, fitHeaders, delimeter=',')
	writeOutputMatrix(fitFileName, oMatrixFit)
	
	
	logReadNumberRatioHistogramFitDict = {\
	'nr2nc':plsqNr2nc, 'nr2npr':plsqNr2npr, \
	'nr2npc':plsqNr2npc, 'nc2npr':plsqNc2npr, \
	'nc2npc':plsqNc2npc, 'npr2npc':plsqNpr2npc }
	
	logReadNumberRatioHistogramIntegralDict = {\
	'nr2nc':integralNr2nc, 'nr2npr':integralNr2npr, \
	'nr2npc':integralNr2npc, 'nc2npr':integralNc2npr, \
	'nc2npc':integralNc2npc, 'npr2npc':integralNpr2npc }
	
	return [logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PlotEffectOfReadThresholdOnPoolPresenceTableTaxonomy(poolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools, \
startThreshold, maxThreshold, poolPresenceTaxonomyTableFileName):

	import matplotlib.pyplot as plt
	from pdb import set_trace
	from kosudoku.pool import CalculatePoolPresenceTableTaxonomyDict
	from kosudoku.output import generateOutputMatrixWithHeaders, writeOutputMatrix

	# Count the number of pool presence table lines with only control wells present

	thresholds = []
	linesThatMapToLibraryAddresses = []
	linesThatMapToSingleLibraryAddresses = []
	linesThatMapToMultipleLibraryAddresses = []
	linesThatMapToUnambiguousLibraryAddresses = []
	linesThatMapToAmbiguousLibraryAddresses = []
	linesThatDoNotMapToLibraryAddresses = []
	linesThatHaveNoReadsAboveThresholdInAnyPool = []
	linesThatMapToControlIndexesOnly = []
	linesThatHaveCoordinatesInNoPoolAxis = []
	linesThatHaveCoordinatesInOnlyOnePoolAxis = []
	linesThatHaveCoordinatesInOnlyTwoPoolAxes = []
	linesThatHaveCoordinatesInOnlyThreePoolAxes = []

	threshold = startThreshold
	
	while threshold <= maxThreshold:
		
		poolPresenceTaxonomyDict = CalculatePoolPresenceTableTaxonomyDict(poolPresenceTable, \
		rowPools, colPools, prPools, pcPools, controlPools, threshold)
	
		thresholds.append(threshold)
		linesThatMapToLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToLibraryAddresses'])
		
		linesThatMapToSingleLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToSingleLibraryAddresses'])
		
		linesThatMapToMultipleLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToMultipleLibraryAddresses'])
		
		linesThatMapToUnambiguousLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToUnambiguousLibraryAddresses'])
		
		linesThatMapToAmbiguousLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToAmbiguousLibraryAddresses'])
		
		linesThatDoNotMapToLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatDoNotMapToLibraryAddresses'])
		
		linesThatHaveNoReadsAboveThresholdInAnyPool.append(\
		poolPresenceTaxonomyDict['linesThatHaveNoReadsAboveThresholdInAnyPool'])
		
		linesThatMapToControlIndexesOnly.append(\
		poolPresenceTaxonomyDict['linesThatMapToControlIndexesOnly'])
		
		linesThatHaveCoordinatesInNoPoolAxis.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInNoPoolAxis'])
		
		linesThatHaveCoordinatesInOnlyOnePoolAxis.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyOnePoolAxis'])
		
		linesThatHaveCoordinatesInOnlyTwoPoolAxes.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyTwoPoolAxes'])
		
		linesThatHaveCoordinatesInOnlyThreePoolAxes.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyThreePoolAxes'])
	
	
		threshold += 1


	headers = ['threshold', 
	'linesThatMapToLibraryAddresses', \
	'linesThatMapToSingleLibraryAddresses', \
	'linesThatMapToMultipleLibraryAddresses', \
	'linesThatMapToUnambiguousLibraryAddresses', \
	'linesThatMapToAmbiguousLibraryAddresses', \
	'linesThatDoNotMapToLibraryAddresses', \
	'linesThatHaveNoReadsAboveThresholdInAnyPool', \
	'linesThatMapToControlIndexesOnly', \
	'linesThatHaveCoordinatesInNoPoolAxis', \
	'linesThatHaveCoordinatesInOnlyOnePoolAxis', \
	'linesThatHaveCoordinatesInOnlyTwoPoolAxes', \
	'linesThatHaveCoordinatesInOnlyThreePoolAxes']

	vectorList = [\
	thresholds, \
	linesThatMapToLibraryAddresses, \
	linesThatMapToSingleLibraryAddresses, \
	linesThatMapToMultipleLibraryAddresses, \
	linesThatMapToUnambiguousLibraryAddresses, \
	linesThatMapToAmbiguousLibraryAddresses, \
	linesThatDoNotMapToLibraryAddresses, \
	linesThatHaveNoReadsAboveThresholdInAnyPool, \
	linesThatMapToControlIndexesOnly, \
	linesThatHaveCoordinatesInNoPoolAxis, \
	linesThatHaveCoordinatesInOnlyOnePoolAxis, \
	linesThatHaveCoordinatesInOnlyTwoPoolAxes, \
	linesThatHaveCoordinatesInOnlyThreePoolAxes]


	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
	writeOutputMatrix(poolPresenceTaxonomyTableFileName, oMatrix)
	
	plt.figure()
	plt.plot(thresholds, linesThatMapToLibraryAddresses, \
	label='Lines That Map To Collection Addresses')
	
	plt.plot(thresholds, linesThatMapToSingleLibraryAddresses, \
	label='Lines That Map To Single Collection Addresses')
	
	plt.plot(thresholds, linesThatMapToMultipleLibraryAddresses, \
	label='Lines That Map To Multiple Collection Addresses')
	
	plt.plot(thresholds, linesThatMapToUnambiguousLibraryAddresses, \
	label='Lines That Map Unambiguously to Collection Addresses')
	
	plt.plot(thresholds, linesThatMapToAmbiguousLibraryAddresses, \
	label='Lines That Map To Ambiguously to Collection Addresses')
	
	plt.plot(thresholds, linesThatDoNotMapToLibraryAddresses, \
	label='Lines That Do Not Map To Collection Addresses')
	
	plt.plot(thresholds, linesThatHaveNoReadsAboveThresholdInAnyPool, \
	label='Lines That Have No Reads Above Threshold In Any Pool')
	
	plt.plot(thresholds, linesThatMapToControlIndexesOnly, \
	label='Lines That Map To Control Indexes Only')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInNoPoolAxis, \
	label='Lines That Have Coordinates In No Pool Axis')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyOnePoolAxis, \
	label='Lines That Have Coordinates In Only One Pool Axis')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyTwoPoolAxes, \
	label='Lines That Have Coordinates In Only Two Pool Axes')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyThreePoolAxes, \
	label='Lines That Have Coordinates In Only Three Pool Axes')
	
	plt.grid()
	plt.legend()
	plt.xlabel('Read Count Threshold')
	plt.title('Effect of Read Count Threshold on Taxonomy of Pool Presence Table Solution')
	plt.show()
	
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
rowPools, colPools, prPools, pcPools):
	
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import figure, plot
	
	[readAlignmentCoords, totalReads] = \
	GenerateHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
	rowPools, colPools, prPools, pcPools)
	
	figure()
	
	plot(readAlignmentCoords, totalReads, marker='o')
	
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ConvertBinEdgesToCenters(binEdges):
	
	from numpy import zeros
	
	binCenters = zeros(len(binEdges)-1)
	
	i = 0
	
	while i < len(binEdges) - 1:
		
		binCenter = binEdges[i] + (binEdges[i + 1] - binEdges[i])
		
		binCenters[i] = binCenter
		i += 1
		
	
	return binCenters
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateReadNumberRatiosForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):
	
	import pdb
	from numpy import log
	
	# Initialize the ratio arrays
	
	nr2ncArray = []
	nr2nprArray = []
	nr2npcArray = []
	nc2nprArray = []
	nc2npcArray = []
	npr2npcArray = []
	
	logNr2ncArray = []
	logNr2nprArray = []
	logNr2npcArray = []
	logNc2nprArray = []
	logNc2npcArray = []
	logNpr2npcArray = []
	
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		# Calculate if the line is single occupancy
		
		lenPossibleAddresses = len(possibleAddresses)
	
		if lenPossibleAddresses == 1:
	
			[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
			CalculatePoolCoordsForLine(\
			poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, \
			threshold=threshold)
			
			
			nr = addresses_r[1][addresses_r[0][0]]
			nc = addresses_c[1][addresses_c[0][0]]
			npr = addresses_pr[1][addresses_pr[0][0]]
			npc = addresses_pc[1][addresses_pc[0][0]]
					
			nr2nc = nr/nc
			nr2npr = nr/npr
			nr2npc = nr/npc
			nc2npr = nc/npr
			nc2npc = nc/npc
			npr2npc = npr/npc

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
			
			logNr2ncArray.append(log(nr2nc))
			logNr2nprArray.append(log(nr2npr)) 
			logNr2npcArray.append(log(nr2npc))
			logNc2nprArray.append(log(nc2npr))
			logNc2npcArray.append(log(nc2npc))
			logNpr2npcArray.append(log(npr2npc))

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
		
		i += 1

	return [nr2ncArray, nr2nprArray, nr2npcArray, nc2nprArray, nc2npcArray, npr2npcArray, \
	logNr2ncArray, logNr2nprArray, logNr2npcArray, logNc2nprArray, logNc2npcArray, logNpr2npcArray]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculateReadNumberRatiosForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):
	
	import pdb
	from numpy import log
	from kosudoku.pool import CalculateLibraryAddressesForPoolPresenceTableLine, \
	CalculatePoolCoordsForLine
	
	# Initialize the ratio arrays
	
	nr2ncArray = []
	nr2nprArray = []
	nr2npcArray = []
	nc2nprArray = []
	nc2npcArray = []
	npr2npcArray = []
	
	logNr2ncArray = []
	logNr2nprArray = []
	logNr2npcArray = []
	logNc2nprArray = []
	logNc2npcArray = []
	logNpr2npcArray = []
	
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		# Calculate if the line is single occupancy
		
		lenPossibleAddresses = len(possibleAddresses)
	
		if lenPossibleAddresses == 1:
	
			[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
			CalculatePoolCoordsForLine(\
			poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, \
			threshold=threshold)
			
			
			nr = addresses_r[1][addresses_r[0][0]]
			nc = addresses_c[1][addresses_c[0][0]]
			npr = addresses_pr[1][addresses_pr[0][0]]
			npc = addresses_pc[1][addresses_pc[0][0]]
					
			nr2nc = nr/nc
			nr2npr = nr/npr
			nr2npc = nr/npc
			nc2npr = nc/npr
			nc2npc = nc/npc
			npr2npc = npr/npc

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
			
			logNr2ncArray.append(log(nr2nc))
			logNr2nprArray.append(log(nr2npr)) 
			logNr2npcArray.append(log(nr2npc))
			logNc2nprArray.append(log(nc2npr))
			logNc2npcArray.append(log(nc2npc))
			logNpr2npcArray.append(log(npr2npc))

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
		
		i += 1

	return [nr2ncArray, nr2nprArray, nr2npcArray, nc2nprArray, nc2npcArray, npr2npcArray, \
	logNr2ncArray, logNr2nprArray, logNr2npcArray, logNc2nprArray, logNc2npcArray, logNpr2npcArray]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculateReadCountHistogramsForSingleAddressLines(poolPresenceTable, rowPools, colPools, \
prPools, pcPools, threshold):
	
	import pdb
	from kosudoku.pool import CalculateLibraryAddressesForPoolPresenceTableLine
	
	# Total number of reads per line, with and without controls
	totalReadCountPerLineArray = []
	totalReadCountPerLineNoControlsArray = []
	
	# Number of counts in each row, col, pr and pc coord pool per line
	rowPoolCountArray = []
	colPoolCountArray = []
	prPoolCountArray = []
	pcPoolCountArray = []
	
	# Total number of reads in row, col, pr and pc pools per line
	totalRowPoolCountArray = []
	totalColPoolCountArray = []
	totalPRPoolCountArray = []
	totalPCPoolCountArray = []
	
	
	i = 0
	while i < len(poolPresenceTable):
	
		poolPresenceTableLine = poolPresenceTable[i]
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTableLine, rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		# Calculate if the line is single occupancy
		
		lenPossibleAddresses = len(possibleAddresses)
	
		if lenPossibleAddresses == 1:
		
			rowEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, rowPools, \
			threshold=threshold)
			colEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, colPools, \
			threshold=threshold)
			prEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, prPools, \
			threshold=threshold)
			pcEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, pcPools, \
			threshold=threshold)
		
					
			for entry in rowEntries:
				rowPoolCountArray.append(entry)
						
			for entry in colEntries:
				colPoolCountArray.append(entry)
				
			for entry in prEntries:
				prPoolCountArray.append(entry)
				
			for entry in pcEntries:
				pcPoolCountArray.append(entry)
				
		
		i += 1

	return [rowPoolCountArray, colPoolCountArray, prPoolCountArray, pcPoolCountArray]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WriteBayesianInferenceParameters(fitDict, integralDict, bayesianParameterFileName):
# Code for writing out read count ratio fit parameters for use in Bayesian inference	
	import pdb
	import numpy
	from kosudoku.output import generateOutputMatrixWithHeaders, writeOutputMatrix
	
	fileHandle = open(bayesianParameterFileName, 'w')
	
	fitKeys = fitDict.keys()
	
	for key in fitKeys:
		
		fitArray = fitDict[key]
		
		outputStr = key + ','
		
		i = 0
		fitStr = ''
		while i < len(fitArray):
			fitStr += str(fitArray[i]) + ','
			i += 1
		
		outputStr = key + ',' + fitStr + str(integralDict[key]) + '\n'
		fileHandle.write(outputStr)
		

	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ReadBayesianInferenceParameters(fileName):
# Code for reading in Bayesian inference parameters. Assumes that key is in the first column of each
# line, that fit parameters are in the first 4 columns and that normalization constant is in the
# 5th column. 
	
	import pdb
	from numpy import array, float
	
	fileHandle = open(fileName, 'r')
	lines = fileHandle.readlines()
	
	integralDict = {}
	fitDict = {}
	
	i = 0
	while i < len(lines):
		lineData = lines[i].split(',')
		key = lineData[0]
		
		j = 1
		fitArray = []
		while j < len(lineData) - 1:
			fitArray.append(float(lineData[j]))
			j += 1
		
		integralDict[key] = float(lineData[-1])
		fitDict[key] = array(fitArray, float)
		
		i += 1
	
	return [fitDict, integralDict]
# ------------------------------------------------------------------------------------------------ #
