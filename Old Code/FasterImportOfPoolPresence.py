
#!/sw/bin/python3.5

# ------------------------------------------------------------------------------------------------ #
def ImportPoolPresenceTableFaster(poolPresenceTableFileName, poolColumns):
	
	import numpy
	import pdb
	import csv
	
	from numpy import genfromtxt
	
	
	
	
	
	fHandle = open(poolPresenceTableFileName, 'r')
	header = fHandle.readline()
	headerCols = header.strip().split(',')
	fHandle.close()
	
	names = []
	colsToUse = []
	i = 0
	while i < len(headerCols):
		col = headerCols[i]
		if col in poolColumns:
			names.append(col)
			colsToUse.append(i)
	
	
	poolPresenceTable = np.genfromtxt(poolPresenceTableFileName, delimiter=',', useCols=colsToUse, \
	names=names)

		
	return poolPresenceTable

# ------------------------------------------------------------------------------------------------ #



# ----------------------------------------------------------------------------------------------- #
# FasterImportOfPoolPresenceTable.py
# ----------------------------------------------------------------------------------------------- #


from sudokuutils13.input import get_input

from sudokuutils13.grid import ImportSudokuGridLayout

from sudokuutils13.pool import ImportPoolPresenceTable





# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

argv = sys.argv

inputParameters = [\
'rowPools', 'colPools', 'controlPools', \
'poolPresenceTableFileName', 'sudokuGridLayout']

inputParameterFileName = \
'/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Files/Version 13/' \
+ 'SudokuAnalyze13_PoolPresenceAnalyze.inp'

inputParameterValues = get_input(['', inputParameterFileName], inputParameters)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data

# Some input parameters that we don't need right now, but might be useful


controlPools = inputParameterValues['controlPools'].split(',')
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')

poolPresenceTableFileName = inputParameterValues['poolPresenceTableFileName']
sudokuGridLayout = inputParameterValues['sudokuGridLayout']

# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Get the sudoku grid plate layout and the plate row and plate column pools from the sudoku
# grid definition file
[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
colPools)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Calculate the sudoku pool columns

if controlPools[0].lower() == 'none':
	sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools
	controlPools = None
else:
	sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools \
	+ controlPools


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table
poolPresenceTable = ImportPoolPresenceTableFaster(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #






