# some functions for mesh grading in 1D

def getGradedMesh( x1, x2, n, ratio):	#	define function that grades mesh
	#	x1: first node coord, 
	#	x2: last node coord
	#	n: number of elements
	#	ratio: ratio of last element/first element
	
	floatn = float(n)	# we use a float version of n for some parts
	coordList = []		#	this is the list we store the coords in
	if(n==1):
		coordList = [x1, x2];
		return coordList
	
	if(ratio == 1):
		alpha = 1.0
		factor = 1.0 / floatn
	else:
		texp = 1.0 / ( floatn - 1.0)
		alpha = ratio**texp
		factor = (1.0 - alpha) / (1.0-alpha**floatn)
	
	deltax = ( x2 - x1 ) * factor
	coordList.append(x1)
	
	for i in range(1,n+1):
		coordi = coordList[i-1] + deltax
		coordList.append(coordi)
		deltax = deltax * alpha
	
	return coordList
		
# 	.................. end of getGradedMesh ..............................