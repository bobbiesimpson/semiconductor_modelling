#!/usr/bin/python

from numpy import array, zeros, dot
from scipy.linalg import solve
import pylab
from math import exp, log
import meshGrading

# This is a simple bit of code which computes the Poisson equation for semicondcutor modelling
# using finite differences

#	-------------------------- Device paramters ----------------------

q = 1.60219e-19						#	electronic charge (C)
eps_0 = 8.85419e-14					# 	absolute permittivity (C V^-1 cm^-1)
eps_sr = 11.9						# 	relative permittivity of silicon dioxide
k = 1.38062e-23						# 	Boltzmann constant
T = 300.0							# 	temperature (K)
ni = 1.45e10						# 	intrinsic concentration (cm^-3)
Na = 5.0e18						# 	concentration of acceptor atoms (cm^-3)
kT = k * T							#	simply k * T
epsilon = eps_0*eps_sr				#	the permittivity of silicon
phi_n = kT/q * log( Na / ni )   	#	quasi fermi levels
phi_p = phi_n	
ft = k* T / q

L = 100.0e-7						#	length of semiconductor (cm)

#	-------------------------- Forcing functions ----------------------

def getForcingTerm(psi):
	n = ni * exp( (  ( psi - phi_n ) ) / ft )
	p = ni * exp( (  ( phi_p - psi ) ) / ft ) 
	E =  q  * (- n - Na + p)
	return E
		
def getForcingTermDeriv(psi): 
	Ederiv = - q  * q / kT * ( ni * exp( q * ( psi - phi_n ) / kT ) \
					+ ni * exp( q * ( phi_p - psi ) / kT ) )
	return Ederiv		
	
#	-------------------------- Mesh parameters -------------------------

M = 4							# 	number of elements
Nnodes = M + 1					#	number of nodes
x0 = 0.0						#	position of first node
xL = x0 + L						# 	position of last node
mesh_ratio = 1.0				# 	mesh ratio = size last element / size of first element
deltaX = L / M					#	assume uniform nodal spacing

#	-------------------------- Boundary conditions ---------------------

boundCondLoc = [0, M]			#	the nodal location of boundary conditions
boundCondValue = [2.0, 0.0]		#	and the values of the boundary conditions

#	-------------------------- Mesh generation ------------------------

nodalCoords = meshGrading.getGradedMesh(x0, xL, M, mesh_ratio)
nodArray = array(nodalCoords)				#	create an array (useful functions)

#	-------------------------- NR procedure ------------------------

tolerance = 1.0e-8
error = 1.0
iterationCount = 0

K = zeros((Nnodes, Nnodes))			#	create the K matrix
RHS = zeros((Nnodes, 1))			# 	create the RHS vector

psi_n = zeros((Nnodes, 1), dtype=float)	#	create initial guess for psi
psi_n = psi_n + 0.0

for n in boundCondLoc:
	psi_n[n] = boundCondValue[boundCondLoc.index(n)]		#	set the initial guess to satisfy BCs

while error > tolerance:

	if iterationCount > 0: break
	
	for i in range(0, Nnodes):			# 	clear the matrices
		RHS[i] = 0.0
		for j in range(0, Nnodes):
			K[i][j] = 0.0
	
	for point in range(0,Nnodes):	# 	loop over all the points in the mesh
	
		if point in boundCondLoc:	# 	if we are at a Dirichlet boundary condition
			BCvalue = boundCondValue[boundCondLoc.index(point)]
			K[point, point] = 1.0
			RHS[point] = BCvalue
		else:		
			K[point, point - 1] += epsilon / deltaX
			K[point, point] += - 2.0 * epsilon / deltaX
			K[point, point + 1] += epsilon / deltaX
			
			potential_n = psi_n[point]			#	get the potential for the previous iteration
			Ederiv = getForcingTermDeriv(potential_n)
			K[point, point] += Ederiv * deltaX
		
			RHS[point] += deltaX * ( - getForcingTerm(potential_n) + Ederiv * potential_n )
		
	print K
	print RHS
		
	newPsi = solve(K, RHS)

	error = 0.0
	for term in range(0, Nnodes):
		error = error + abs(newPsi[term] - psi_n[term]) / ( 1.0 + psi_n[term])
		
	psi_n = newPsi
	iterationCount += 1
	
print iterationCount, " iterations"
pylab.plot(nodArray, psi_n, 'ko-')
pylab.show()


