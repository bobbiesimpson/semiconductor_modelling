# 	the forcing term (ie the RHS of Poisson's equation) is required for calculation
#	of various matrix terms. 

from math import exp, cos, sin

def E(psi, q, epsilon, phi_n, phi_p, kT, ni, Na):
	ft = kT / q	
	n = ni * exp( (  ( psi - phi_n ) ) / ft )
	p = ni * exp( (  ( phi_p - psi ) ) / ft ) 
	E =  - q  * (- n - Na + p)	
	return  E
						  
def Ederiv(psi, q, epsilon, phi_n, phi_p, kT, ni, Na):
	ft = kT / q	
	Ederiv =  q  * q / kT * ( ni * exp( ( psi - phi_n ) / ft ) \
					+ ni * exp( ( phi_p - psi ) / ft ) ) 
	#n = ni * exp( (  ( psi - phi_n ) ))# / ft )
	#p = ni * exp( (  ( phi_p - psi ) ))# / ft ) 
	#Ederiv =   q  * ( n  + p)							
	return Ederiv

