#!/usr/bin/python

from numpy import *
from scipy.linalg import solve
import pylab
import forcingTerm, shapeFunctions, meshGrading, gaussNodes


# 	This is a 1D FE code to solve Poisson's equation for semiconductor modelling
#	We'll start with a simple 1D semiconductor where a potential is prescribed at
#	both ends

#	-------------------------- Device paramters ----------------------

q = 1.60219e-19						#	electronic charge (C)
eps_0 = 8.85419e-14					# 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 11.9						# 	relative permittivity of silicon dioxide
k = 1.38062e-23						# 	Boltzmann constant
T = 300.0							# 	temperature (K)
ni = 1.45e10						# 	intrinsic concentration (cm^-3)
#Na = 5.0e18						# 	concentration of acceptor atoms (cm^-3)
kT = k * T							#	simply k * T
epsilon = eps_0*eps_sr				#	the permittivity of silicon

def main(dopingConc, appliedVoltage, plotting):			#	the main function that runs our code
											#	pass in: doping Conc, applied V (LHS) and plotting flag
	Na = dopingConc
	phi_n = kT/q * math.log( Na / ni )   	#	quasi fermi levels
	phi_p = phi_n
	
	L = 100.0e-7						#	length of semiconductor (cm)

	#	-------------------------- Mesh parameters -------------------------

	M = 100						# 	number of elements
	x0 = 0.0						#	position of first node
	xL = x0 + L						# 	position of last node
	mesh_ratio = 1.0					# 	mesh ratio = size last element / size of first element


	#	-------------------------- Boundary conditions ---------------------

	boundCond = array([[0, appliedVoltage], [M, 0.0]])	#	the nodal location of boundary conditions

	#	-------------------------- Mesh generation ------------------------

	nodalCoords = meshGrading.getGradedMesh(x0, xL, M, mesh_ratio)
	nodArray = array(nodalCoords)				#	create an array (useful functions)

	#	-------------------------- Node connectivity -------------------------

	#	we need to know which nodes are connected to which elements

	connectivity = zeros((M,2))
	for element in range(M):
		connectivity[element,0] = element
		connectivity[element,1] = element + 1

	#   -------------------------- psi interpolation --------------------------

	def getPsi_n(element, localXi):
		globalNode1 = connectivity[element,0]		#	global node num of local node 1
		globalNode2 = connectivity[element,1]		#	global node num of local node 2
		psi = shapeFunctions.N1(localXi) * psi_n[globalNode1] \
				+ shapeFunctions.N2(localXi) * psi_n[globalNode2]
		return psi
	
	#	-------------------------- NR + matrix setup -------------------------

	#	The Newton-Raphson procedure requires an initial guess to allow the 
	#	increment (Delta psi) to be found. We simply make this vector zero while
	#	making sure to satisfy the boundary conditions

	iteration = 0						# 	iteration number is zero at start
	residualRef = 0.0					#	the reference residual (put psi_n = {0}) in function
	residual = 0.0						#	the residual for the present iteration
	residualRatio = 1.0					#	simply the ratio of residual / residualRef
	residualPlotValues = zeros((20, M+1))
	tolerance = 1e-6					#	the NR tolerance
	max_iterations = 1000

	psi_n = zeros((size(nodArray), 1),dtype=float)
	psi_n = psi_n + 0.0

	for n in boundCond:
		psi_n[n[0]] = n[1]				#	set the initial guess to satisfy BCs
	
	#	NR ITERATION LOOP
	while residualRatio > tolerance:
	
		if iteration > max_iterations:
			print "max no. of iterations reached"
			break

		T = zeros((M+1, M+1))					# 	zero Tangent matrix
		Mass = zeros((M+1, M+1))
		K = zeros((M+1, M+1))					# 	zero K matrix
		Fb = zeros((M+1,1))					# 	zero body force vector	

		#	-------------------------- Matrix computation -------------------------
	
		numGPs = 2					# 	define num. Gauss points
		gpt,gwt = gaussNodes.gaussNodes(numGPs) 	# 	gauss points and weights

		for element in range(M):
			nodes = connectivity[element]		#	get nodes for current element
			x1 = nodArray[nodes[0]]			#	the node coords for element
			x2 = nodArray[nodes[1]]
			length = abs(x2 - x1)			# 	length of element
		
			#	STIFFNESS MATRIX COMPUTATION (analytical expression)
			K_el = 1.0 / length * array([[1,-1],[-1,1]]) * epsilon
	
			#	NUMERICAL INTEGRATION (mass matrix, body force vector)
			m11 = m12 = 0.0				# 	initialise mass matrix terms
			fb1 = fb2 = 0.0				# 	initialise body force terms
			for gp in range(numGPs):
				gxi = gpt[gp]
				gw =  gwt[gp]
			
				psi = getPsi_n(element, gxi)
				E = forcingTerm.E(psi, q, epsilon, phi_n, phi_p, kT, ni, Na)
				Ederiv = forcingTerm.Ederiv(psi, q, epsilon, phi_n, phi_p, kT, ni, Na)
				#print E, Ederiv
				N1 = shapeFunctions.N1(gxi)
				N2 = shapeFunctions.N2(gxi) 
				detJacob = length / 2.0
			
				#	MASS MATRIX COMPUTATION
				m11 = m11 + N1 * Ederiv * N1 * gw * detJacob
				m12 = m12 + N1 * Ederiv * N2 * gw * detJacob
			
				#	BODY FORCE COMPUTATION
				fb1 = fb1 + N1 * E * gw * detJacob
				fb2 = fb2 + N2 * E * gw * detJacob
				
				#gp = 2
				# .............. end loop over Gauss points
				
			M_el = array([[m11, m12], [m12, m11]])
			Fb_el = array([[fb1],[fb2]])
			#print K_el, M_el, Fb_el
		
			# 	now put the element matrices in the global matrices
			for i in range(size(nodes)):
				globali = connectivity[element,i]
				Fb[globali] = Fb[globali] + Fb_el[i]
				for j in range(size(nodes)):
					globalj = connectivity[element,j]
					T[globali,globalj] = T[globali,globalj] + K_el[i,j] + M_el[i,j]
					K[globali,globalj] = K[globali,globalj] + K_el[i,j]
					Mass[globali,globalj] = Mass[globali,globalj] + M_el[i,j]
				
		# .............. end loop over elements
	
		#	-------------------------- Solution ----------------------------------

		RHS = -dot(K, psi_n) - Fb

		if iteration == 0: 
			for term in RHS:			#	reference residual 
				residualRef = residualRef + term**2		
			residualRef = residualRef**0.5

		for i in boundCond[:,0]:			# 	apply boundary conditions
			for j in range(size(nodArray)):
				T[i,j] = 0.0
			T[i, i] = 1.0
			RHS[i] = 0.0
		#if iteration < 20: RHS[0] = 0.1
		
		solution = solve(T,RHS)				#	solve the system of equations
	
		for term in range(size(solution)):		#	get the new solution using the increment
			#if term in boundCond[:,0]:
				#continue
			psi_n[term] =  psi_n[term] + solution[term]
	
		#	-------------------------- Residual ----------------------------------
	
		residual = 0.0					#	calculate the residual and if on first iteration, the 
		for term in RHS:				#	reference residual
			residual = residual + term**2		
		residual = residual**0.5
		residualRatio = residual / residualRef
		
		if iteration < 20: 
			for term in range(size(RHS)):
				residualPlotValues[iteration, term] = RHS[term]
		
		print iteration, residualRatio
		
		#pylab.plot(nodArray, RHS)
		#pylab.savefig('Resiudal_plot.eps')
		#pylab.show()

		iteration = iteration + 1

		# .............. end loop over NR iteration

	psiList = zeros(size(psi_n))
	for term in range(size(psi_n)):
		psiList[term] = psi_n[term,0]

#	-------------------------- Plotting ---------------------------------

	if(plotting):
		y = zeros(size(nodArray,0))		#	array of zeros for plotting
		#pylab.plot(nodArray, y, 'ko-')		#	plot the mesh
		#pylab.plot(nodArray, psiList, 'ko-')
		#pylab.show()
	
		
		for count in range(10):			#	residual plotting
			pylab.plot(nodArray[0:10], residualPlotValues[count,0:10], label='iteration %s' % count)
			pylab.xlabel('x')
			pylab.ylabel('residual')
			pylab.legend()
		pylab.savefig('residual_6GP_potential_2V.eps')	
		pylab.show()
		
			
	# .............. main functionhttp://matplotlib.sourceforge.net/
	
	return nodArray, psiList

#	let's run the main function

#main(5.0e18, 2.0, True)



	

	










