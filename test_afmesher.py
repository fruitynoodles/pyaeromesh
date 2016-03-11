#!/usr/bin/env python
import os
from math import *
#import numpy
from PyFoam.Applications.Runner import Runner
from PyFoam.Applications.PlotRunner import PlotRunner
from PyFoam.Applications.CloneCase import CloneCase
from aerofoilmesher import *


# clone wing_template case
#CloneCase(args=["wing_template",case])
# generate eliptical wing planform
threed = False # set this to False to generate a two-dimensional aerofoil
blunt_te = True # set this to True to make trailing edge blunt instead of radius
num_c = 5 # number of chord sections for the wing
span = 3.6 # wingspan
c_0 = 0.5 # root chord
th = 0.1 # max thickness-to-chord ratio
alpha = 0.0 # angle of attack
camber = 10.0 # camber angle
stack = 0.25 # stacking line as a fration of chord
sweep = 40.0
span_c = []
chords = []
dihedral = []
if(threed):
	for i in range(num_c):
		span_c.append(span/2.0*i/(num_c-1))
		chords.append(max(c_0*sqrt(1.0-(2.0*span_c[i]/span)**2),0.1))
		if(span_c[i]<span/6):
			dihedral.append(-10.0)
		else:
			dihedral.append(7.0)
	# area beyond wingtip
	span_c.append(span/2.0*2)
	chords.append(chords[-1])
	dihedral.append(0.0)
else:
	num_c = 2
	chords = [c_0,c_0]
	span_c = [0.0,c_0/10]
	dihedral = [0.0,0.0]

for i in range(0,1):
	alpha = float(i)
	# generate root aerofoil
	#[cline,p,s] = genNACA65(c_0,th,camber,alpha)
	[p,s] = read_ps_profile("clarky.csv",c_0,alpha)
	farfield = 5.0 #number of chord lengths to farfield boundary
	case = "wing2d_a"+str(i) #case directory name
	CloneCase(args=["wing2d",case])
	# create blockmesh file for wing
	write_blockmesh(case,threed,blunt_te,span_c,chords,stack,sweep,dihedral,p,s,grid_w = 60)
	# run blockMesh on the case
	Runner(args=["--clear","blockMesh", "-case",case])
	# ask the user if they want to examine the mesh in paraview
	#viewit = raw_input("Do you want to view the mesh in ParaView? (y/N): ")
	#if(viewit in ["y","Y","yes","Yes","YES"]):
	#	Runner(args=["paraFoam","-case",case])
	# run the simulation
	PlotRunner(args=["--progress","simpleFoam","-case",case])
	# calculate y+ wall values
	Runner(args=["yPlusRAS","-case",case])
	# ask the user if they want to examine the mesh in paraview
	#viewit = raw_input("Do you want to view the results in ParaView? (y/N): ")
	#if(viewit in ["y","Y","yes","Yes","YES"]):
	#	Runner(args=["paraFoam","-case",case])
