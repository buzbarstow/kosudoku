#!/sw/bin/python3.6

# ----------------------------------------------------------------------------------------------- #
# kosudoku-collectionmc
# Created by Buz Barstow 2016-05-06
# Last modified by Buz Barstow 2018-11-09

# Monte Carlo code that simulates picking of a random progenitor library and determines gene
# knockout coverage.
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
from utilities.montecarlo import SimulateMultiplePickings, SimulatePicking, PoissonEstimateOfGenesHit
from utilities.output import generateOutputMatrixWithHeaders, writeOutputMatrix
from utilities.input import get_input
from utilities.utils import ensure_dir

import sys
from numpy import ones, float
from matplotlib.pyplot import ion, ioff, show, figure, plot, grid, xlabel, ylabel, legend

# ------------------------------------------------------------------------------------------------ #


# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file
## Import the input parameter file
inputParameters = ['numberOfTrials', 'maxMutants', 'outputFileName',
'transposonCoordToFeatureIndexFileName']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/pre-pool-files/'
file_name = file_intro + 'collectionmc.inp'
inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
numberOfTrials = int(inputParameterValues['numberOfTrials'])
maxMutants = int(inputParameterValues['maxMutants'])
transposonCoordToFeatureIndexFileName = file_intro + inputParameterValues['transposonCoordToFeatureIndexFileName']
outputFileName = file_intro + inputParameterValues['outputFileName']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
ensure_dir(outputFileName)
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Do the simulation for the maximum non essential set

[iAxis, averageFeatureHitCount, sdFeatureHitCount, featureHitCountUpperBound, \
 featureHitCountLowerBound, noUniqHittableFeatures] = \
    SimulateMultiplePickings(transposonCoordToFeatureIndexFileName, numberOfTrials, maxMutants)

noUniqHittableFeaturesArray = ones(len(iAxis), float) * noUniqHittableFeatures
noUniqHittableFeaturesArray99 = ones(len(iAxis), float) * noUniqHittableFeatures * 0.99
noUniqHittableFeaturesArray95 = ones(len(iAxis), float) * noUniqHittableFeatures * 0.95

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Do a poisson simulation as well
poissonUniqueGenesHit = PoissonEstimateOfGenesHit(iAxis, noUniqHittableFeatures)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
headers = ['Mutants', 'AvUniqFeatures', 'SDUniq', 'FeaturesUpper', 'FeaturesLower']
vectorList = [iAxis, averageFeatureHitCount, sdFeatureHitCount, \
              featureHitCountUpperBound, featureHitCountLowerBound]
oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
writeOutputMatrix(outputFileName, oMatrix)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output simulation data
figure()
plot(iAxis, averageFeatureHitCount, c='purple', label='Monte Carlo estimate')
plot(iAxis, featureHitCountUpperBound, c='blue')
plot(iAxis, featureHitCountLowerBound, c='blue')

plot(iAxis, noUniqHittableFeaturesArray, c='black')
plot(iAxis, noUniqHittableFeaturesArray99, c='black')
plot(iAxis, noUniqHittableFeaturesArray95, c='black')

plot(iAxis, poissonUniqueGenesHit, c='red', label='Poisson estimate')
grid()
xlabel('Mutants Picked')
ylabel('Non-essential Genes Hit')
legend(loc=4)
show()
# ------------------------------------------------------------------------------------------------ #