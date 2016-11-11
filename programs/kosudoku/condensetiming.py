# ------------------------------------------------------------------------------------------------ #
def CalculateRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, sealAlTime, \
aerasealTime, numberOfHumans, numberOfRobots, avMutRearrayPerPlate, avMutCPPerPlate, \
avPicksPerCPMutant, mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94):
# This function calculates the rate of rearraying progenitor collection plates

	import pdb
	
	if numberOfRobots == 0:
		rearrayAndStreakOutRate, rearrayCondPlateGenerationRate, CPCondPlateGenerationRate, \
		coloniesToPickGenerationRate = \
		CalculateManualRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, \
		sealAlTime, \
		aerasealTime, numberOfHumans, avMutRearrayPerPlate, avMutCPPerPlate, \
		avPicksPerCPMutant, mutantsPerRearrayCondPlate=mutantsPerRearrayCondPlate, \
		mutantsPerCPCondPlate=mutantsPerCPCondPlate)
	
	elif numberOfRobots > 0:
		
		rearrayAndStreakOutRate, rearrayCondPlateGenerationRate, \
		CPCondPlateGenerationRate, coloniesToPickGenerationRate, bestRobottingHumans = \
		CalculateRoboticallyAssistedRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, \
		sealAlTime, \
		aerasealTime, numberOfHumans, numberOfRobots, avMutRearrayPerPlate, avMutCPPerPlate, \
		avPicksPerCPMutant, mutantsPerRearrayCondPlate=mutantsPerRearrayCondPlate, \
		mutantsPerCPCondPlate=mutantsPerCPCondPlate)
	
	
	return rearrayAndStreakOutRate, rearrayCondPlateGenerationRate, \
	CPCondPlateGenerationRate, coloniesToPickGenerationRate
# ------------------------------------------------------------------------------------------------ #











# ------------------------------------------------------------------------------------------------ #
def CalculateManualRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, sealAlTime, \
aerasealTime, numberOfHumans, avMutRearrayPerPlate, avMutCPPerPlate, \
avPicksPerCPMutant, mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94):
# This function calculates the rate of rearraying progenitor collection plates

	import pdb
	
	# This part calculates how long it takes to deal with a single progenitor collection plate
	rearrayAndStreakOutTimePerProgPlate = unsealAlTime \
	+ avMutRearrayPerPlate/rearrayRate \
	+ avMutCPPerPlate*streakOutTime \
	+ sealAlTime
	
	# pdb.set_trace()
	
	# Figure out how many simply rearrayed condensed plates come from one progenitor plate
	# (on average)
	noRearrayCondPlatesPerProgPlate = avMutRearrayPerPlate/mutantsPerRearrayCondPlate
	
	# Figure out how many colony purified condensed plates come from one progenitor plate
	# (on average)
	noCPCondPlatesPerProgPlate = avMutCPPerPlate*avPicksPerCPMutant/mutantsPerCPCondPlate
	noColoniesToPickPerProgPlate = avMutCPPerPlate*avPicksPerCPMutant
	
	# Calculate the fractional time needed to seal up the simple rearray plate. Remember,
	# for every rearray condensed plate that you seal up, you will have unsealed and sealed 
	# a lot of progenitor plates. We can divide this time between these progenitor plates, and then
	# add this fractional time to the time needed to process a progenitor plate
	condPlateSealTimePerProgPlate = aerasealTime/noRearrayCondPlatesPerProgPlate

	# Calculate the total average amount of time it takes to unseal a progenitor plate, 
	# rearray and streak out from it, seal it back up, and the fractional time needed to 
	# seal up the simple rearray plate. 
	
	totalRearrayAndStreakOutTimePerProgPlate = rearrayAndStreakOutTimePerProgPlate \
	+ condPlateSealTimePerProgPlate
	
	# Calculate the rate of rearraying and streaking out from progenitor collection plates
	# Calculate this in plates per second
	manualRearrayAndStreakOutRate = numberOfHumans/totalRearrayAndStreakOutTimePerProgPlate

	# Then, calculate the rate of generation of rearrayed condensed plates, and eventually
	# of colony purified condensed plates per second
	manualCPCondPlateGenerationRate = manualRearrayAndStreakOutRate*noCPCondPlatesPerProgPlate
	manualRearrayCondPlateGenerationRate = \
	manualRearrayAndStreakOutRate*noRearrayCondPlatesPerProgPlate
	
	coloniesToPickGenerationRate = manualRearrayAndStreakOutRate*noColoniesToPickPerProgPlate
	

	return manualRearrayAndStreakOutRate, manualRearrayCondPlateGenerationRate, \
	manualCPCondPlateGenerationRate, coloniesToPickGenerationRate
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateRoboticallyAssistedRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, \
sealAlTime, \
aerasealTime, noHumans, noRobots, avMutRearrayPerPlate, avMutCPPerPlate, \
avPicksPerCPMutant, mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94):
# This function calculates the rate of rearraying progenitor collection plates

	import pdb
	
	robottingHumans = 1
	rearrayStreakOutRateArray = []
	robottingHumansArray = []
	
	while robottingHumans <= min(noRobots, noHumans - 1):
		
		manualHumans = noHumans-robottingHumans
		
		# This part calculates how long it takes to deal with a single progenitor collection plate
		rearrayAndStreakOutTimePerProgPlate = (unsealAlTime + sealAlTime)/manualHumans \
		+ avMutRearrayPerPlate/(rearrayRate*robottingHumans) \
		+ avMutCPPerPlate*streakOutTime/manualHumans \
	
		# pdb.set_trace()
	
		# Figure out how many simply rearrayed condensed plates come from one progenitor plate
		# (on average)
		noRearrayCondPlatesPerProgPlate = avMutRearrayPerPlate/mutantsPerRearrayCondPlate
	
		# Figure out how many colony purified condensed plates come from one progenitor plate
		# (on average)
		noCPCondPlatesPerProgPlate = avMutCPPerPlate*avPicksPerCPMutant/mutantsPerCPCondPlate
		noColoniesToPickPerProgPlate = avMutCPPerPlate*avPicksPerCPMutant
	
		# Calculate the fractional time needed to seal up the simple rearray plate. Remember,
		# for every rearray condensed plate that you seal up, you will have unsealed and sealed 
		# a lot of progenitor plates. We can divide this time between these progenitor plates, and then
		# add this fractional time to the time needed to process a progenitor plate
		condPlateSealTimePerProgPlate = aerasealTime/\
		(noRearrayCondPlatesPerProgPlate*(noHumans-robottingHumans))

		# Calculate the total average amount of time it takes to unseal a progenitor plate, 
		# rearray and streak out from it, seal it back up, and the fractional time needed to 
		# seal up the simple rearray plate. 
	
		totalRearrayAndStreakOutTimePerProgPlate = rearrayAndStreakOutTimePerProgPlate \
		+ condPlateSealTimePerProgPlate
				
		# Calculate the rate of rearraying and streaking out from progenitor collection plates
		# Calculate this in plates per second
		rearrayAndStreakOutRate = 1.0/totalRearrayAndStreakOutTimePerProgPlate
	
		robottingHumansArray.append(robottingHumans)
		rearrayStreakOutRateArray.append(rearrayAndStreakOutRate)
		
		robottingHumans += 1

	try:
		bestRearrayStreakOutRate = max(rearrayStreakOutRateArray)
	except:
		pdb.set_trace()
	
	bestRearrayStreakOutRateIndex = \
	rearrayStreakOutRateArray.index(bestRearrayStreakOutRate)
	
	bestRobottingHumans = robottingHumansArray[bestRearrayStreakOutRateIndex]
	
	
	# Then, calculate the rate of generation of rearrayed condensed plates, and eventually
	# of colony purified condensed plates per second
	CPCondPlateGenerationRate = bestRearrayStreakOutRate*noCPCondPlatesPerProgPlate
	
	rearrayCondPlateGenerationRate = \
	bestRearrayStreakOutRate*noRearrayCondPlatesPerProgPlate

	coloniesToPickGenerationRate = bestRearrayStreakOutRate*noColoniesToPickPerProgPlate
	

	return bestRearrayStreakOutRate, rearrayCondPlateGenerationRate, \
	CPCondPlateGenerationRate, coloniesToPickGenerationRate, bestRobottingHumans
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def CalcNoProgPlatesThatCanBeManuallyRearrayedAndStruckOut(rearrayRate, streakOutTime, \
unsealAlTime, sealAlTime, aerasealTime, hoursAvailable, \
noHumans, avMutRearrayPerPlate, avMutCPPerPlate, avPicksPerCPMutant, \
mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94):
# Calculate the number of progenitor plates that can be rearrayed and struck out in a given
# amount of time, and how many rearrayed condensed plates this will generate, how many colonies
# will need to be picked, and how colony purified condensed plates it will generate.

	import pdb

	manualRearrayAndStreakOutRate, manualRearrayCondPlateGenerationRate, \
	manualCPCondPlateGenerationRate, coloniesToPickGenerationRate \
	= CalculateManualRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, sealAlTime, \
	aerasealTime, noHumans, avMutRearrayPerPlate, \
	avMutCPPerPlate, avPicksPerCPMutant, mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94)
	
# 	pdb.set_trace()
	
	noProgPlates = hoursAvailable*3600*manualRearrayAndStreakOutRate
	
	noRerrayCondPlates = hoursAvailable*3600*manualRearrayCondPlateGenerationRate
	noCPCondPlates = hoursAvailable*3600*manualCPCondPlateGenerationRate
	noColoniesToPick = hoursAvailable*3600*coloniesToPickGenerationRate
	
	return noProgPlates, noRerrayCondPlates, noCPCondPlates, noColoniesToPick
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalcNoProgPlatesThatCanBeRearrayedAndStruckOut(rearrayRate, streakOutTime, \
unsealAlTime, sealAlTime, aerasealTime, hoursAvailable, \
noHumans, noRobots, avMutRearrayPerPlate, avMutCPPerPlate, avPicksPerCPMutant, \
mutantsPerRearrayCondPlate=94, mutantsPerCPCondPlate=94):
# Calculate the number of progenitor plates that can be rearrayed and struck out in a given
# amount of time, and how many rearrayed condensed plates this will generate, how many colonies
# will need to be picked, and how colony purified condensed plates it will generate.

	import pdb

	rearrayAndStreakOutRate, rearrayCondPlateGenerationRate, \
	CPCondPlateGenerationRate, coloniesToPickGenerationRate = \
	CalculateRearrayAndStreakOutRate(rearrayRate, streakOutTime, unsealAlTime, sealAlTime, \
	aerasealTime, noHumans, noRobots, avMutRearrayPerPlate, avMutCPPerPlate, \
	avPicksPerCPMutant, mutantsPerRearrayCondPlate=mutantsPerRearrayCondPlate, \
	mutantsPerCPCondPlate=mutantsPerCPCondPlate)	

# 	pdb.set_trace()
	
	noProgPlates = hoursAvailable*3600*rearrayAndStreakOutRate
	
	noRerrayCondPlates = hoursAvailable*3600*rearrayCondPlateGenerationRate
	noCPCondPlates = hoursAvailable*3600*CPCondPlateGenerationRate
	noColoniesToPick = hoursAvailable*3600*coloniesToPickGenerationRate
	
	return noProgPlates, noRerrayCondPlates, noCPCondPlates, noColoniesToPick
# ------------------------------------------------------------------------------------------------ #









# ------------------------------------------------------------------------------------------------ #
def CalculateCryoPreservationTime(noPipettorOperators, noSealUnsealOperators, \
pipettingTimePerPlate, aerasealUnsealTime, sealAlTime, noPlates):

	pipettingRate = noPipettorOperators/pipettingTimePerPlate
	totalPipettingTime = noPlates/pipettingRate
		
	unsealingAndSealingRate = noSealUnsealOperators/(aerasealUnsealTime + sealAlTime)
	totalUnsealingAndSealingTime = noPlates/unsealingAndSealingRate
		
	longestCryoPreservationTaskTime = max(totalPipettingTime, totalUnsealingAndSealingTime)
	
	return longestCryoPreservationTaskTime, totalPipettingTime, totalUnsealingAndSealingTime
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateMinimizedCryoPreservationTime(noHumans, noPipettors, pipettingTimePerPlate, \
aerasealUnsealTime, sealAlTime, noPlates):
# This is in no way an ideal optimizer, but it should help a little
	
	import pdb

	if noHumans == 1:
		longestCryoPreservationTaskTime, totalPipettingTime, totalUnsealingAndSealingTime = \
		CalculateCryoPreservationTime(1, 1, pipettingTimePerPlate, aerasealUnsealTime, sealAlTime, \
		noPlates)
		
		bestCryopreservationTime = totalPipettingTime + totalUnsealingAndSealingTime
		bestPipettingHumans = 1
		bestSealingHumans = 1
		
		# pdb.set_trace()
		
	elif noHumans == 2:
		# Assume if we have 2 humans that one is on pipetting duty and one on unsealing and 
		# sealing duty. The total time will be the longer of the unsealing and sealing time
		# and the pipetting. For now, I'm choosing to not do anything fancy, where say after
		# pipetting, another person switches over to sealing. 

		longestCryoPreservationTaskTime, totalPipettingTime, totalUnsealingAndSealingTime = \
		CalculateCryoPreservationTime(1, 1, pipettingTimePerPlate, aerasealUnsealTime, sealAlTime, \
		noPlates)
		
		bestCryopreservationTime = longestCryoPreservationTaskTime
		bestPipettingHumans = 1
		bestSealingHumans = noHumans - 1
		
	elif noHumans > 2:
		if noPipettors == 1:
			# Assume that we put surplus people on sealing and unsealing duty
			totalCryopreservationTime, totalPipettingTime, totalUnsealingAndSealingTime = \
			CalculateCryoPreservationTime(1, noHumans-1, pipettingTimePerPlate, \
			aerasealUnsealTime, sealAlTime, noPlates)
			
			bestCryopreservationTime = totalCryopreservationTime
			bestPipettingHumans = 1
			bestSealingHumans = noHumans - 1

		elif noPipettors > 1:
			# If there is more than 1 pipettor available, try to find the best balance
			# between parallelizing the pipetting and the unsealing-sealing.
			pipettingHumans = 1
			totalCryopreservationTimeArray = []
			pipettingHumansArray = []
			while pipettingHumans < min(noPipettors, noHumans - 1):
				
				totalCryopreservationTime, totalPipettingTime, totalUnsealingAndSealingTime = \
				CalculateCryoPreservationTime(pipettingHumans, noHumans-pipettingHumans, \
				pipettingTimePerPlate, aerasealUnsealTime, sealAlTime, noPlates)

				totalCryopreservationTimeArray.append(totalCryopreservationTime)
				pipettingHumansArray.append(pipettingHumans)
				pipettingHumans += 1
			
			# Find the lowest total cryopreservation time
			bestCryopreservationTime = min(totalCryopreservationTimeArray)
			bestCryopreservationTimeIndex = \
			totalCryopreservationTimeArray.index(bestCryopreservationTime)
			bestPipettingHumans = pipettingHumansArray[bestCryopreservationTimeIndex]
			bestSealingHumans = noHumans - bestPipettingHumans
			
	return bestCryopreservationTime, bestPipettingHumans, bestSealingHumans
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculatePickingTime(noHumans, noColonies, pickingTime):
	
	totalPickingTime = noColonies*pickingTime/noHumans
	
	return totalPickingTime
# ------------------------------------------------------------------------------------------------ #
