# from itertools import product
# from collections import deque


# def needleman_wunsch(x, y):
#     """Run the Needleman-Wunsch algorithm on two sequences.

#     x, y -- sequences.

#     Code based on pseudocode in Section 3 of:

#     Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
#     "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
#     https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
#     """
#     N, M = len(x), len(y)
#     s = lambda a, b: int(a == b)

#     DIAG = -1, -1
#     LEFT = -1, 0
#     UP = 0, -1

#     # Create tables F and Ptr
#     F = {}
#     Ptr = {}

#     F[-1, -1] = 0
#     for i in range(N):
#         F[i, -1] = -i
#     for j in range(M):
#         F[-1, j] = -j

#     option_Ptr = DIAG, LEFT, UP
#     for i, j in product(range(N), range(M)):
#         option_F = (
#             F[i - 1, j - 1] + s(x[i], y[j]),
#             F[i - 1, j] - 1,
#             F[i, j - 1] - 1,
#         )
#         F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

#     # Work backwards from (N - 1, M - 1) to (0, 0)
#     # to find the best alignment.
#     alignment = deque()
#     i, j = N - 1, M - 1
#     while i >= 0 and j >= 0:
#         direction = Ptr[i, j]
#         if direction == DIAG:
#             element = i, j
#         elif direction == LEFT:
#             element = i, None
#         elif direction == UP:
#             element = None, j
#         alignment.appendleft(element)
#         di, dj = direction
#         i, j = i + di, j + dj
#     while i >= 0:
#         alignment.appendleft((i, None))
#         i -= 1
#     while j >= 0:
#         alignment.appendleft((None, j))
#         j -= 1

#     return list(alignment)




class ScoreParams:
	'''
	Define scores for each parameter
	'''
	def __init__(self,gap,match,mismatch):
		self.gap = gap
		self.match = match
		self.mismatch = mismatch

	def misMatchChar(self,x,y):
		if x != y:
			return self.mismatch
		else:
			return self.match

def getMatrix(sizeX,sizeY,gap):
	'''
	Create an initial matrix of zeros, such that its len(x) x len(y)
	'''
	matrix = []
	for i in range(len(sizeX)+1):
		subMatrix = []
		for j in range(len(sizeY)+1):
			subMatrix.append(0)
		matrix.append(subMatrix)

	# Initializing the first row and first column with the gap values
	for j in range(1,len(sizeY)+1):
		matrix[0][j] = j*gap
	for i in range(1,len(sizeX)+1):
		matrix[i][0] = i*gap
	return matrix

def getTraceBackMatrix(sizeX,sizeY):
	'''
	Create an initial matrix of zeros, such that its len(x) x len(y)
	'''
	matrix = []
	for i in range(len(sizeX)+1):
		subMatrix = []
		for j in range(len(sizeY)+1):
			subMatrix.append('0')
		matrix.append(subMatrix)

	# Initializing the first row and first column with the up or left values
	for j in range(1,len(sizeY)+1):
		matrix[0][j] = 'left'
	for i in range(1,len(sizeX)+1):
		matrix[i][0] = 'up'
	matrix[0][0] = 'done'
	return matrix


def globalAlign(x,y,score):
	'''
	Fill in the matrix with alignment scores
	'''
	matrix = getMatrix(x,y,score.gap)
	traceBack = getTraceBackMatrix(x,y)

	for i in range(1,len(x)+1):
		for j in range(1,len(y)+1):
			left = matrix[i][j-1] + score.gap
			up = matrix[i-1][j] + score.gap
			diag = matrix[i-1][j-1] + score.misMatchChar(x[i-1],y[j-1])
			matrix[i][j] = max(left,up,diag)
			if matrix[i][j] == left:
				traceBack[i][j] = 'left'
			elif matrix[i][j] == up:
				traceBack[i][j] = 'up'
			else:
				traceBack[i][j] = 'diag'
	return matrix,traceBack

def getAlignedSequences(x,y,matrix,traceBack):
	'''
	Obtain x and y globally aligned sequence arrays using the bottom-up approach
	'''
	xSeq = []
	ySeq = []
	i = len(x)
	j = len(y)
	while(i > 0 or j > 0):
		if traceBack[i][j] == 'diag':
			# Diag is scored when x[i-1] == y[j-1]
			xSeq.append(x[i-1])
			ySeq.append(y[j-1])
			i = i-1
			j = j-1
		elif traceBack[i][j] == 'left':
			# Left holds true when '-' is added from x string and y[j-1] from y string
			xSeq.append('-')
			ySeq.append(y[j-1])
			j = j-1
		elif traceBack[i][j] == 'up':
			# Up holds true when '-' is added from y string and x[j-1] from x string
			xSeq.append(x[i-1])
			ySeq.append('-')
			i = i-1
		elif traceBack[i][j] == 'done':
			# Break condition when we reach the [0,0] cell of traceback matrix
			break
	return xSeq,ySeq

def printMatrix(matrix):
	'''
	Create a custom function to print the matrix
	'''
	for i in range(len(matrix)):
		print(matrix[i])
	print()

'''
Driver Code:
'''
x = 'GAATTC'
y = 'GAATTA'
print('Input sequences are: ')
print(x)
print(y)
print()
score = ScoreParams(-2,2,-1)
matrix,traceBack = globalAlign(x,y,score)

print('Printing the score matrix:')
printMatrix(matrix)

print('Printing the trace back matrix:')
printMatrix(traceBack)

xSeq,ySeq = getAlignedSequences(x,y,matrix,traceBack)

print('The globally aligned sequences are:')
print(*xSeq[::-1])
print(*ySeq[::-1])