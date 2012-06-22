#!/usr/bin/env python 

#/software/devmod/python/2.6.1/bin/python

import sys
import os
import math
from matplotlib import *
from pylab import *
from point import Point

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})



#sys.path[:0]=['/usr/lib64/python2.6/lib-tk/']
#print sys.path

from Tkinter import *

class Poisson:
    def __init__(self, parent):
    
        self.bgcolour='#EEF'
    
        self.mpParent=parent
	self.myContainer1 = Frame(parent, background=self.bgcolour) 
	self.myContainer1.pack(expand=YES, fill=BOTH)
	
	#======================================================================#
	
	self.vg_scale = Scale(self.myContainer1, background=self.bgcolour,
	                      highlightbackground=self.bgcolour)
	self.vg_scale['from']=2.0
	self.vg_scale['to']=0.0
	self.vg_scale['resolution']=0.1
	#self.vg_scale['tickinterval']=0.5
	self.vg_scale['label']='Vg'
	self.vg_scale.pack()
	
	#======================================================================#
	
	self.oxideFrame = Frame(self.myContainer1, background=self.bgcolour)
	self.oxideFrame.pack()
	Label(self.oxideFrame, text="Tox [nm]", justify=LEFT, 
	      background=self.bgcolour).pack(side=LEFT)
	self.oxide = StringVar()
	self.oxide.set('1.0')
	self.oxideEntry = Entry(self.oxideFrame, textvariable=self.oxide,
	                         width=6, justify=CENTER, background=self.bgcolour,
	                         highlightbackground=self.bgcolour)
	self.oxideEntry.focus_force()
	self.oxideEntry.pack(side=RIGHT, expand=YES, anchor=CENTER)
	
	#======================================================================#
	
	self.dopingFrame = Frame(self.myContainer1, background=self.bgcolour)
	self.dopingFrame.pack()
	Label(self.dopingFrame, text="Substrate\nDoping", justify=LEFT, 
	      background=self.bgcolour).pack(side=LEFT)
	self.doping = StringVar()
	self.doping.set('5e18')
	self.dopingEntry = Entry(self.dopingFrame, textvariable=self.doping,
	                         width=6, justify=CENTER, background=self.bgcolour,
	                         highlightbackground=self.bgcolour)
	self.dopingEntry.focus_force()
	self.dopingEntry.pack(side=RIGHT, expand=YES, anchor=CENTER)
	
	#======================================================================#
	
	self.gracePlot = IntVar()
	self.gracePlot.set(0)
	self.graph_plot = Checkbutton(self.myContainer1, text='Plot graph?',
	                              variable=self.gracePlot, 
				      background=self.bgcolour,
	                              highlightbackground=self.bgcolour)
	self.graph_plot.pack()
	
	#======================================================================#
	
        self.button_run = Button(self.myContainer1, text="Run!")
	self.button_run["background"]='green'
	self.button_run["activebackground"]='yellow'
	self.button_run.pack()
	self.button_run.bind("<Button-1>", self.run)   
	self.button_run.bind("<Return>", self.run) 
	
    
    def run(self,event):
	self.button_run["background"]='red'
	self.button_run["activebackground"]='red'
        self.vg=self.vg_scale.get()
        self.solve()
	self.button_run["background"]='green'
	self.button_run["activebackground"]='yellow'
    
    
    def solve(self):
	# physical constants
	qel=1.60219e-19		#electronic charge
	bkb=1.38062e-23		#Boltamann constant
	am0=0.910956e-30	#electron mass
	hbr=1.05459e-32		#h bar
	e00=8.85419e-14		#permitivity of free space
	esr=11.9		#rel. permitivity of Si
	eor=3.9			#rel. permitivity of oxide
	afs=4.05		#electron affinity of Si
	afo=0.95		#electron affinity of oxide
	tab=300			#temperature
	ft=bkb*tab/qel		#kT/q
	cni=1.45e10		#intrinsic concentration

	#device parameters
	dx=1e-8
	tdev=100e-7					# ie the length of the device is 100 nm or 100e-7cm
	tox=float(self.oxide.get())*1e-7 # get the thickness of the oxide layer
	acceptors=float(self.doping.get())	# the doping concentrating of donor atoms
	print 'Acceptor doping:',acceptors

	tol=1e-8				# tolerance we set for newton iteration scheme
	niter=0
	vg=self.vg				
	error=1.0

	nox=int(tox/dx)	      #grid line for interface	
	ntot=int((tox+tdev)/dx)	# the last grid line on the model

	print 'nox =',nox
	print 'ntot =',ntot

	fn=ft*math.log(acceptors/cni)      #quasi fermi level

	mesh=[]			# list of points that we create

	for i in range(ntot+1):
	    mesh.append(Point())	# create a point and put it in the mesh
	    mesh[i].x=i*dx			# set the xcoord of the point
	    mesh[i].p=0.0			# set the potential of the point = 0
	    if i<nox:				# if we are in the oxide
	       mesh[i].eps=e00*eor  # set the permitivity in oxide
	    else:
	       mesh[i].eps=e00*esr	# otherwise set permitivity in Si
	       mesh[i].ca=acceptors # set the acceptor doping concentration in Si

	    #print i,mesh[i].x,mesh[i].eps,mesh[i].ca

	# calculate discretisation
	for i in range(ntot):					# loop over all elements
	    mesh[i].disc_r += mesh[i].eps/dx	# 
	    mesh[i].disc_m += -mesh[i].eps/dx
	    if i>=nox: mesh[i].disc_rhs += dx/2.0

	    mesh[i+1].disc_l += mesh[i].eps/dx
	    mesh[i+1].disc_m += -mesh[i].eps/dx
	    if i>=nox: mesh[i+1].disc_rhs += dx/2.0

	# boundary conditions
	mesh[0].fix_point(vg)
	mesh[ntot].fix_point(0.0)
	print vg

	n=0

	# Gaussian elimination solver
	def gaus3d(a,b,c,d,r):
	    k=len(a)	   
	    d[0]=d[0]/b[0]  		     
	    for i in range(1,k): 
        	# print 'loop 1:',i		     
		b[i]=b[i]-a[i]*c[i-1]/b[i-1]      
		d[i]=(d[i]-a[i]*d[i-1])/b[i] 			
	    r[k-1]=d[k-1]			
	    for i in range(k-2,-1,-1):	
        	# print 'loop 2:',i		  	
		r[i]=d[i]-c[i]*r[i+1]/b[i]	

	# Do Newton iterations
	while error>tol:
	    niter+=1
	    a=[mp.disc_l for mp in mesh]
	    b=[]
	    c=[mp.disc_r for mp in mesh]
	    d=[]
	    v=[mp.p for mp in mesh]
	    # Go through all mesh points to generate final discretisation arrays
	    for mp in mesh:
        	if mp.fixed:
        	    b.append(mp.disc_m)
        	    d.append(mp.disc_rhs)
		else:
		    mp.cn=cni*math.exp((mp.p-fn)/ft)
		    mp.cp=cni*math.exp((fn-mp.p)/ft)
		    # calculate RHS charge
        	    bb=qel*(-mp.ca-mp.cn+mp.cp)*mp.disc_rhs
        	    # calculate Jacobian
        	    bj=-(qel/ft)*(mp.cn+mp.cp)*mp.disc_rhs
        	    b.append(mp.disc_m+bj)
        	    d.append(-bb+bj*mp.p)

	    gaus3d(a,b,c,d,v)

	    #for i in range(1000):
	    #    for i in range(1,ntot):
	    #        pnew=-(a[i]*v[i-1]+c[i]*v[i+1])/b[i]
	    #	    v[i]+=1.9*(pnew-v[i])

	    error=0.0
	    for i,mp in enumerate(mesh):
        	error=error+abs(v[i]-mp.p)/(1.0+mp.p)

	    print 'Iteration:',niter,' -- Error:',error

	    for i,mp in enumerate(mesh):
        	mp.p=v[i]
		if i<nox:
		    mp.cn=1.0
        	    mp.cp=1.0
		else:
		    mp.cn=cni*math.exp((mp.p-fn)/ft)
        	    mp.cp=cni*math.exp((fn-mp.p)/ft)
		#print i,mp.p,mp.cn,mp.cp
		

	#with open('out.dat', 'w') as f:
	f=open('out.dat', 'w')
	x=[] 
	y=[]
	for mp in mesh:
            #f.write('{0:6.2f} {1:12.8f} {2:18.10e} {3:18.10e}\n'.format((mp.x-mesh[nox].x)*1e7,mp.p,mp.cn,mp.cp))
            f.write('%6.2f %12.8f %18.10e %18.10e \n' % ((mp.x-mesh[nox].x)*1e7,mp.p,mp.cn,mp.cp))
            x.append((mp.x-mesh[nox].x)*1e7)
            y.append(mp.p)
    	f.close()
	
	if self.gracePlot.get(): os.system('xmgrace out.dat&')
	
	# plot 
	plot(x, y)
	xlabel(r'distance ($nm$)')
	ylabel(r'electron concentration $(cm^{-3})$')
	savefig('electronConcentration.eps')
	show()
	
# instatiate the application class so it will run
root = Tk()
poisson=Poisson(root)
root.mainloop()
