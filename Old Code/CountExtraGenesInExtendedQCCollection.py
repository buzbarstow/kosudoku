from scipy import unique

fileHandle = open('/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-05-10 - Sudoku Sequencing of Shewanella Sudoku Library/Input Data/TransposonsInExtendedQCCollection.py', 'r')

data = fileHandle.readlines()

fileHandle.close()

transposons = []

i = 0
while i < len(data):
	transposons.append(int(data[i].strip()))
	i += 1
	
uniqTransposons = unique(transposons)

lenUniquTransposons = len(uniqTransposons)

print('Unique extra transposons: ' + str(lenUniquTransposons))