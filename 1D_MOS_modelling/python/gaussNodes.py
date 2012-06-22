# filename: gaussNodes.py
#
# x,A = gaussNodes(m,tol=10e-9)
# Returns nodal abscissas {x} and weights {A} of 
# Gauss-Legendre m-point quadrature
#
# Sundar, ITACM, Cardiff University
# 2010
#----------------------------------------------------------------------
from math import cos,pi
from numpy import *

def gaussNodes(m,tol=10e-9):

	def legendre(t,m):
		p0 = 1.0; p1 = t
		for k in range(1,m):
			p = ((2.0*k+1.0)*t*p1 - k*p0)/(1.0+k)
			p0 = p1; p1 = p
		dp = m*(p0-t*p1)/(1.0 - t**2)
		return p,dp

	A = zeros( (m),dtype=float)
	x = zeros( (m),dtype=float)
	nRoots = (m+1)/2
	for i in range(nRoots):
		t = cos(pi*(i+0.75)/(m+0.5))
		for j in range(30):
			p,dp = legendre(t,m)
			dt = -p/dp; t = t+dt
			if abs(dt) < tol:
				x[i] = t; x[m-i-1] =-t
				A[i] = 2.0/(1.0 - t**2)/(dp**2)
				A[m-i-1] = A[i]
				break
	return x,A
