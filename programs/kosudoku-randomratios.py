#!/sw/bin/python3.5

# Random ratios code


# ------------------------------------------------------------------------------------------------ #
# Functions for calculating a Voigt function to a read ratio histogram
def voigtFunction(x, p):
	
	from scipy.special import erfc
	from numpy import exp
	from numpy import sqrt
	from numpy import pi, float64

	a = p[0]
	c = p[1]
	delta = p[2]
	sigma = p[3]
	
	firstArg = ((-1.j)*(-c + x) + delta)/(sqrt(2)*sigma)
	secondArg = ((1.j)*(-c + x) + delta)/(sqrt(2)*sigma)

	voigtEquation = a\
	*(exp(firstArg**2)*erfc(firstArg) \
	+ exp(secondArg**2)*erfc(secondArg) ) \
	/ (2*sqrt(2*pi)*sigma)
	
	voigtEquation = float64(voigtEquation)
	
	return voigtEquation

	
def voigtResiduals(p, y, x): 
	err = y - voigtFunction(x,p) 
	return err

def voigtFit(x,y, p0):
	import scipy
	from scipy.optimize import leastsq
	
	plsq = leastsq(voigtResiduals, p0, args=(y, x), maxfev=2000)
	
	return [plsq[0], voigtFunction(x, plsq[0])]

# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
# Functions for performing Gaussian fit to data
def gaussianFunction(x, p): 
	fit = p[0]*(exp( - ((x-p[1])**2) / p[2]))
	return fit

def gaussianResiduals(p, y, x): 
	err = y - gaussianFunction(x,p) 
	return err

def gaussianFit(x, y, p0):
	import scipy
	from scipy.optimize import leastsq
	
	plsq = leastsq(gaussianResiduals, p0, args=(y, x), maxfev=2000)
	
	return [plsq[0], gaussianFunction(x, plsq[0])]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ConvertBinEdgesToCenters(binEdges):
	
	from numpy import zeros
	
	binCenters = zeros(len(binEdges)-1)
	
	i = 0
	
	while i < len(binEdges) - 1:
		
		binCenter = binEdges[i] + (binEdges[i + 1] - binEdges[i])/2
		
		binCenters[i] = binCenter
		i += 1
		
	
	return binCenters
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotAndFitRatio(ratioC, binSize=0.1, title='', xlabel='', ylabel='', outputFileName='', \
outputData=False): 

	import matplotlib.pyplot as plt 
	import pdb
	from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix
	

	binSize = 0.1
	bins = numpy.arange(min(ratioC)-binSize, max(ratioC)+binSize, binSize)
	values, base = numpy.histogram(ratioC, bins=bins)
	binCenters = ConvertBinEdgesToCenters(base)
	integral = simps(values, binCenters)
	normValues = values/integral

	[plsqG0, gaussFit] = gaussianFit(binCenters, normValues, [1,1,1])
	[plsqV0, vFit] = voigtFit(binCenters, normValues, array([0.1, -0.5, 0.35, 0.5]))


	plt.figure()
	plt.scatter(binCenters, normValues, marker='o')
	plt.plot(binCenters, gaussFit, color='r')
	plt.plot(binCenters, vFit, color='g')
	plt.grid()
	plt.xlabel(xlabel)
	plt.title(title)
	plt.ylabel(ylabel)
	
	
	if outputData == True and len(outputFileName) > 0:
# 		pdb.set_trace()
		headers = [xlabel, ylabel, ylabel + '_gf', ylabel + '_vf']
		vectorList = [binCenters, normValues, gaussFit, vFit]
		oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
		writeOutputMatrix(outputFileName, oMatrix)
	
	
	return
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def PlotRowsVsColumsOnScatterPlot(rowReads, colReads, xlabel='', ylabel='', outputFileName='', \
outputData=False):

	import matplotlib.pyplot as plt 
	import pdb
	from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix
	


	plt.figure()
	plt.scatter(log(rowReads), log(colReads))
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	
	if outputData == True and len(outputFileName) > 0:
		vectorList = [log(rowReads), log(colReads)]
		headers = [xlabel, ylabel]
		oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
		writeOutputMatrix(outputFileName, oMatrix)
		
	return

# ------------------------------------------------------------------------------------------------ #




from numpy.random import lognormal
from numpy import sum
import numpy
from scipy.integrate import simps
from scipy.linalg import cholesky, eigh
from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix


sampleSize=1000000


# Generate read count ratios when completely uncorrelated
rowReads = lognormal(mean=6., sigma=1, size=sampleSize)
colReads = lognormal(mean=6, sigma=1, size=sampleSize)
ratio = log(colReads/rowReads)

r = numpy.array([[1.0, 0.0], [0.0, 12/8.]])

# evals, evecs = eigh(r)
# c = np.dot(evecs, numpy.diag(np.sqrt(evals)))

c = cholesky(r, lower=True)

# Convert the data to correlated random variables. 
[rowReadsC, colReadsC] = numpy.dot(c, [rowReads, colReads])
ratioC = log(rowReadsC/colReadsC)


r2 = numpy.array([[1.0, 0.0], [0.7, 1.0]])
c2 = cholesky(r2, lower=True)

[rowReadsC2, colReadsC2] = numpy.dot(c2, [log(rowReadsC), log(colReadsC)])

rowReadsC2 = exp(rowReadsC2)
colReadsC2 = exp(colReadsC2)
ratioC2 = log(rowReadsC2/colReadsC2)


i = 0
threshRowReadsC2 = []
threshColReadsC2 = []
threshold = 5
while i < len(rowReadsC2):	
	if rowReadsC2[i] >= threshold and colReadsC2[i] >= threshold:
		threshRowReadsC2.append(rowReadsC2[i])
		threshColReadsC2.append(colReadsC2[i])
	i += 1

threshColReadsC2 = array(threshColReadsC2)
threshRowReadsC2 = array(threshRowReadsC2)
ratioThreshC2 = log(threshRowReadsC2/threshColReadsC2)



outputData = True

PlotAndFitRatio(ratio, binSize=0.1, \
title='Ratio for non-correlated Row and Cols', xlabel='log(col/row)', ylabel='pdf', \
outputFileName='no_correlation.csv', outputData=outputData)

PlotAndFitRatio(ratioC, binSize=0.1, \
title='Ratio for weakly correlated Row and Cols', xlabel='log(colC/rowC)', ylabel='pdfC', \
outputFileName='mean_correlation.csv', outputData=outputData)

PlotAndFitRatio(ratioC2, binSize=0.1, \
title='Ratio for more strongly correlated Row and Cols', xlabel='log(colC2/rowC2)', ylabel='pdfC2', \
outputFileName='tight_correlation.csv', outputData=outputData)

PlotAndFitRatio(ratioThreshC2, binSize=0.1, \
title='Ratio for more strongly correlated and thresholded Row and Cols', \
xlabel='log(threshColC2/threshRowC2)', ylabel='pdfC2T', \
outputFileName='threshold_tight_correlation.csv', outputData=outputData)


PlotRowsVsColumsOnScatterPlot(rowReads, colReads, xlabel='row counts', ylabel='col counts', \
outputFileName='row_vs_col.csv', outputData=outputData)

PlotRowsVsColumsOnScatterPlot(rowReadsC, colReadsC, xlabel='row counts c', ylabel='col counts c', \
outputFileName='row_vs_col_c.csv', outputData=outputData)

PlotRowsVsColumsOnScatterPlot(rowReadsC2, colReadsC2, xlabel='row counts c2', \
ylabel='col counts c2', outputFileName='row_vs_col_c2.csv', outputData=outputData)

PlotRowsVsColumsOnScatterPlot(threshRowReadsC2, threshColReadsC2, xlabel='row counts c2t', \
ylabel='col counts c2t', outputFileName='row_vs_col_c2t.csv', outputData=outputData)


show()