# ------------------------------------------------------------------------------------------------ #
# batch module
# Created by Buz Barstow 2016-05-06
# ------------------------------------------------------------------------------------------------ #

# ---------------------------------------------------------------------------- #
def CalculateNPick(tLiquid, tSolid, rPick=1500, tIncubate=18):

	# Inputs to program
	# rPick - Colony picking rate per hour of picking device
	# tIncubate - time needed for colonies to reach saturation in liquid media
	# tLiquid - survival time of picked colonies after reaching saturation in 
	# liquid media in hours
	# tSolid - survival time of plated colonies after appearing on solid media
	# in hours
	
	# Picking time per day and how many days to we pick for?
	
	if tIncubate > 24:
		tPickPerDay = 15
		pickDays = 1
	else:
		tPickPerDay = 13 + 1/3
		
		if tSolid < tLiquid:
			pickDays = tSolid/24
		else:
			pickDays = tLiquid/24
	
	
	# Maximum number of colonies that you can pick in a single day
	rPickPerDay = rPick*tPickPerDay
	
	# Total number of colonies that you can pick 
	Npick = pickDays*rPickPerDay
	
	# Print out some results
	
	
	return Npick, pickDays, tPickPerDay
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def CalculateNPoolDays(TLiquidShelfLife, restDays, pickDays, tIncubate):
	
	from numpy import ceil
	
	# TLiquidShelfLife - shelf life of organism in liquid media in days
	# restDays - number of days we take to rest between picking and pooling
	# pickDays - number of days spent picking
	# tIncubate - number of hours needed to reach saturation
	
	NPoolDays = TLiquidShelfLife - restDays - pickDays + ceil(tIncubate/24)
	
	return NPoolDays
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def CalculateNPool(rpool, tPoolPerDay, NPoolDays, NcoloniesPerPlate):

	# rpool - pooling rate in plates per hour
	# tPoolPerDay - amount of time spent pooling per day in hours
	# NPoolDays - amount of days spent pooling in days
	# NcoloniesPerPlate - number of colonies per plate
	
	NPool = rpool*NcoloniesPerPlate*tPoolPerDay*NPoolDays
	
	return NPool
# ---------------------------------------------------------------------------- #




