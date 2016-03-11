#!/usr/bin/env python
import os
from math import *
#import numpy
from PyFoam.Applications.Runner import Runner
from PyFoam.Applications.PlotRunner import PlotRunner
from PyFoam.Applications.CloneCase import CloneCase

#from aerofoil_mesher import read_ps_profile

#generate a set of vertices for a NACA-65 aerofoil with a circular arc camberline

# chord fraction
C_frac = (0.000, 0.0050, 0.0075, 0.0125, 0.0254, 0.050, 0.075,
	0.100, 0.150, 0.200, 0.250, 0.300,
	0.350, 0.400, 0.450, 0.500, 0.550,
	0.600, 0.650, 0.700, 0.750, 0.800,
	0.850,0.900, 0.950, 0.9800, 0.9877,
	0.9941, 1.0000)

# thickness  as a fraction of max thickness
Th_frac = ( 0.0000, 0.0772, 0.0932, 0.1169, 0.1574, 0.2177, 0.2647,
	0.3040, 0.3666, 0.4143,	0.4503, 0.4760,
	0.4924, 0.4996, 0.4963, 0.4812, 0.4530,
	0.4146, 0.3732, 0.3318, 0.2904, 0.2490,
	0.2076, 0.1662, 0.1248, 0.1000, 0.0924,
	0.0707, 0.0000)

"""	generate a set of coordinates for a NACA-65 aerofoil profile with a circular-arc camberline
	Arguments:
	c:      chord-length (m)
	t_c:    maximum thickness as a fraction of chord length
	theta:  camber angle (degrees)
	gamma:  stagger angle (degrees)
	returns: [camberline[x][y], pressure surface[x][y], suction surface[x][y]]
"""
def genNACA65(c,t_c,theta,gamma):
	# convert angles to radians
	theta_rad = theta*pi/180.0
	gamma_rad = gamma*pi/180.0
	# initialise lists for camber line, pressure and suction surface coordinates 
	camberline = [[],[]]
	suction = [[],[]]
	pressure = [[],[]]
	# generate geometery with chord line parallel with x-axis
	r = c/2.0/sin(theta_rad/2.0) # camber line arc radius
	ycentre = c/2.0/tan(theta_rad/2.0) #distance between chord and camber center
	for i in range(len(C_frac)):
		x = c*C_frac[i]
        	# calculate camber line coordinates
		camberline[0].append(x-c/2.0)
		camberline[1].append(sqrt(r**2-camberline[0][i]**2)- ycentre)
		# calculate suction and pressure surface coordinates
		suction[0].append(camberline[0][i]+Th_frac[i]*t_c*c*camberline[0][i]/r)
		suction[1].append(camberline[1][i]+Th_frac[i]*t_c*c*(camberline[1][i]+
			ycentre)/r)
		pressure[0].append(camberline[0][i]-Th_frac[i]*t_c*c*camberline[0][i]/r)
		pressure[1].append(camberline[1][i]-Th_frac[i]*t_c*c*(camberline[1][i]+
			ycentre)/r)
		# rotate geometery through stagger angle
		# rotate camber line
		cx = camberline[0][i]*cos(gamma_rad)-camberline[1][i]*sin(gamma_rad)
		cy = camberline[0][i]*sin(gamma_rad)+camberline[1][i]*cos(gamma_rad)
		camberline[0][i] = cx
		camberline[1][i] = cy
		# rotate pressure surface
		px = pressure[0][i]*cos(gamma_rad)-pressure[1][i]*sin(gamma_rad)
		py = pressure[0][i]*sin(gamma_rad)+pressure[1][i]*cos(gamma_rad)
		pressure[0][i] = px
		pressure[1][i] = py
		# rotate suction surface
		sx = suction[0][i]*cos(gamma_rad)-suction[1][i]*sin(gamma_rad)
		sy = suction[0][i]*sin(gamma_rad)+suction[1][i]*cos(gamma_rad)
		suction[0][i] = sx
		suction[1][i] = sy    
	return [camberline,pressure,suction]

# Write a blockMesh file to the constant/polymesh directory of specified case
# usage: write_blockmesh(case,threed,blunt_te,rad,chords,stack
#	farfield=5.0,grid_r=60,grid_c=40, grid_s=8,grid_le=50,grid_w=60,
#	exp_r=1000,exp_lete=4.0,exp_wall=100.0,foam_version="2.0"):
# case: name of case directory
# Wing geometry:
# threed: a boolean value denoting whether the case is a 3D wing (True) or a 2D aerofoil profile (False)
# blunt_te: a boolean value denoting whether the trailing edge is blunt (True) or an arc (False)
# rad: a list of the spanwise distance from the wing symmetry line to each chord 
# chords: a list of the all the chord lengths
# stack: the wing stacking line as a fraction of chord
# p: list of pressure-surface coordinates; p[0] contains x-coords and p[1] contains y-coords
# s: list of suction-surface coordinates; s[0] contains x-coords and s[1] contains y-coords
#    Note: The first and last coordinates of p and s must represent the point of intersection
#          of the camber line and leading/trailing edge radius, and the second/second last point
#          is where the leading/trailing edge arc ends
# Number of cells along principle edges:
# grid_r: grid thickness around wing
# grid_c: number of cells along the wing chord
# grid_s: number of cells along the blade span
# grid_le: number of cells along the leading and trailing edge arcs
# grid_w: number of cells downstream of trailing edge (along wake)
# expansion ratios
# exp_r: expansion ratio from blade to boundary
# exp_lete: expansion ratio from mid-chord to leading and trailing edge
# exp_wall: expansion ratio from tip to edge of domain
# foam_version: String containing the OpenFOAM version number

def write_blockmesh(case,threed,blunt_te,n_blades,rad,chords,camber,stagger,t_c,stack,
	farfield=5.0,grid_r=60,grid_c=40, grid_s=8,grid_le=50,grid_w=60,
	exp_r=100,exp_lete=1.0,exp_wall=10.0,foam_version="2.0"):
	c_0 = chords[0]
	passage_angle = 2*pi/n_blades

	# x and y offsets of aerofoil
	def get_xoff(chord_no,stack_frac,farfield):
		return c_0*(farfield+stack_frac*cos(stagger[0]*pi/180.0))

	def get_yoff(chord_no):
		return rad[chord_no]*passage_angle/2.0

	#vertex layout scheme:
	#	6------7-----8---9
	#	|\  B  |  C  | D |
	#	| 0----1-----2---10
	#	|A(LE wing TE) E |
	#	| 5----4-----3---11
	#	|/  H  |  G  | F |
	#	15-----14----13--12

	def generate_geom(i):
		xoff = get_xoff(i,stack,farfield)
		yoff = get_yoff(i)
		[cline,p,s] = genNACA65(chords[i],t_c,camber[i],stagger[i])
		#leading edge gradient
		legrad = ((s[1][1]+p[1][1])/2.0-s[1][0])/((s[0][1]+p[0][1])/2.0-s[0][0])
		#trailing edge gradient
		tegrad = (s[1][-1]-(s[1][-2]+p[1][-2])/2.0)/(s[0][-1]-(s[0][-2]+p[0][-1])/2.0)

		geom = [(s[0][1],s[1][1]),			#0
			(s[0][len(s[0])/2],s[1][len(s[0])/2]),	#1
			(s[0][-2],s[1][-2]),			#2
			(p[0][-2],p[1][-2]),			#3
			(p[0][len(p[0])/2],p[1][len(p[0])/2]),	#4
			(p[0][1],p[1][1]),			#5
			#(s[0][0]-farfield*chords[i],yoff+s[1][0]-farfield*chords[i]*legrad),	#6
			(s[0][0]-farfield*chords[i],yoff+s[1][0]-farfield*chords[i]*legrad*0.5),	#6
			#(s[0][len(s[0])/2],yoff*0.8),					#7
			(s[0][len(s[0])/2],yoff+cline[1][len(cline[1])/2-1]),					#7
			(s[0][-1],yoff+s[1][-1]),					#8
			(s[0][-1]+farfield*chords[i],yoff+s[1][-1]+farfield*chords[i]*tegrad*0.5),	#9
			(s[0][-1]+farfield*chords[i],yoff*0.1+s[1][-1]+farfield*chords[i]*tegrad*0.5),	#10
			(s[0][-1]+farfield*chords[i],-yoff*0.1+s[1][-1]+farfield*chords[i]*tegrad*0.5),	#11
			(s[0][-1]+farfield*chords[i],-yoff+s[1][-1]+farfield*chords[i]*tegrad*0.5),	#12
			(s[0][-1],-yoff+s[1][-1]),					#13
			#(s[0][len(s[0])/2],-yoff*1.2),					#14
			(s[0][len(s[0])/2],-yoff+cline[1][len(cline[1])/2-1]),					#14
			#(s[0][0]-farfield*chords[i],-yoff+s[1][0]-farfield*chords[i]*legrad)]	#15
			(s[0][0]-farfield*chords[i],-yoff+s[1][0]-farfield*chords[i]*legrad*0.5)]	#15
		return geom
	
	#return an interpolated point for arc edges
	def arc_interp(pt1,pt2,i):
		x = (pt1[0]+pt2[0])/2.0+xoff
		theta = (pt1[1]+pt2[1])/2.0/rad[i]
		y = rad[i]*sin(theta)
		z = rad[i]*cos(theta)
		return [x,y,z]

	xoff = get_xoff(0,stack,farfield)
	yoff = get_yoff(0)
	# set values for 2d case
	if(not threed):
		num_r = 2
		chords = [c_0,c_0]
		rad = [rad[0],rad[0]]
		stagger = [stagger[0],stagger[0]]
		camber = [camber[0],camber[0]]
		grid_s = 1

	#start writing blockMeshDict file
	f = file(os.path.join(case,"constant","polyMesh","blockMeshDict"),"w")
	f.write("""FoamFile
{
	version	"""+foam_version+""";
	format	ascii;
	class	dictionary;
	object	blockMeshDict;
}
convertToMeters	1.0;

vertices
(
""")

	#vertex layout scheme:
	#	6------7-----8---9
	#	|\  B  |  C  | D |
	#	| 0----1-----2---10
	#	|A(LE wing TE) E |
	#	| 5----4-----3---11
	#	|/  H  |  G  | F |
	#	15-----14----13--12

	block_A = (0,6,15,5)
	block_B = (1,7,6,0)
	block_C = (2,8,7,1)
	block_D = (10,9,8,2)
	block_E = (3,11,10,2)
	block_F = (3,13,12,11)
	block_G = (4,14,13,3)
	block_H = (5,15,14,4)

	swirl = farfield*c_0*tan(stagger[len(stagger)/2]*pi/180)/rad[len(rad)/2]
	swirl_mid = farfield*c_0*tan((1-stack)/farfield*stagger[len(stagger)/2]*pi/180)/rad[len(rad)/2]
	for i in range(len(chords)):
		# generate aerofoil
		xoff = get_xoff(i,stack,farfield)
		yoff = get_yoff(i)
		[cline,p,s] = genNACA65(chords[i],t_c,camber[i],stagger[i])
		points = generate_geom(i)
		"""
		for j in range(6):
			f.write("\t("+str(points[j][0]*chords[i]/c_0+xoff)+" "+str(points[j][1]*chords[i]/c_0)+
				" "+str(rad[i])+")\n")
		for j in range(6,16):
			f.write("\t("+str(points[j][0]+xoff)+" "+str(points[j][1])+" "+str(rad[i])+")\n")
		"""
		for j in range(16):
			px = points[j][0]+xoff
			if(threed):
				theta = points[j][1]/rad[i]
				py = rad[i]*sin(theta)
				pz = rad[i]*cos(theta)
			else:
				py = points[j][1]
				pz = i*chords[i]/10.0
			f.write("\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")

	# tip mesh points
	#    0--------1------2
	#   / \   R   |  S  / \
	#  /   A0____A1___A2   \
	# ( Q  |  W   |  X | T  )
	#  \   A5----A4---A3   /
	#   \ /   V   |  U  \ /
	#    5--------4------3
	"""
	if(threed):
		xoff = get_xoff(i-1,stack,farfield)
		f.write("\t("+str(xoff+0.98*(0.7*points[0][0]+0.3*points[5][0])*chords[i-1]/c_0)+" "+
			str((0.7*points[0][1]+0.3*points[5][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A0
		f.write("\t("+str(xoff+(0.7*points[1][0]+0.3*points[4][0])*chords[i-1]/c_0)+" "+
			str((0.7*points[1][1]+0.3*points[4][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A1
		f.write("\t("+str(xoff+0.98*(0.7*points[2][0]+0.3*points[3][0])*chords[i-1]/c_0)+" "+
			str((0.7*points[2][1]+0.3*points[3][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A2
		f.write("\t("+str(xoff+0.98*(0.3*points[2][0]+0.7*points[3][0])*chords[i-1]/c_0)+" "+
			str((0.3*points[2][1]+0.7*points[3][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A3
		f.write("\t("+str(xoff+(0.3*points[1][0]+0.7*points[4][0])*chords[i-1]/c_0)+" "+
			str((0.3*points[1][1]+0.7*points[4][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A4
		f.write("\t("+str(xoff+0.98*(0.3*points[0][0]+0.7*points[5][0])*chords[i-1]/c_0)+" "+
			str((0.3*points[0][1]+0.7*points[5][1])*chords[i-1]/c_0)+" "+
			str(rad[i-1])+")\n")#A5

		xoff = get_xoff(i,stack,farfield)
		f.write("\t("+str(xoff+0.98*(0.7*points[0][0]+0.3*points[5][0])*chords[i]/c_0)+" "+
			str((0.7*points[0][1]+0.3*points[5][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A6
		f.write("\t("+str(xoff+(0.7*points[1][0]+0.3*points[4][0])*chords[i]/c_0)+" "+
			str((0.7*points[1][1]+0.3*points[4][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A7
		f.write("\t("+str(xoff+0.98*(0.7*points[2][0]+0.3*points[3][0])*chords[i]/c_0)+" "+
			str((0.7*points[2][1]+0.3*points[3][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A8
		f.write("\t("+str(xoff+0.98*(0.3*points[2][0]+0.7*points[3][0])*chords[i]/c_0)+" "+
			str((0.3*points[2][1]+0.7*points[3][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A9
		f.write("\t("+str(xoff+(0.3*points[1][0]+0.7*points[4][0])*chords[i]/c_0)+" "+
			str((0.3*points[1][1]+0.7*points[4][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A10
		f.write("\t("+str(xoff+0.98*(0.3*points[0][0]+0.7*points[5][0])*chords[i]/c_0)+" "+
			str((0.3*points[0][1]+0.7*points[5][1])*chords[i]/c_0)+" "+
			str(rad[i])+")\n")#A11
"""
	f.write(");\n")

	#write edges
	f.write("""
edges
(
""")

	for i in range(len(chords)):
		xoff = get_xoff(i,stack,farfield)
		yoff = get_yoff(i)
		[cline,p,s] = genNACA65(chords[i],t_c,camber[i],stagger[i])
		points = generate_geom(i)

		# Annulus arcs
		if(threed):
			#inlet (pts 6-15)
			[ix,iy,iz] = arc_interp(points[6],points[15],i)
			f.write("\tarc "+str(6+i*16)+" "+str(15+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 0-6
			[ix,iy,iz] = arc_interp(points[0],points[6],i)
			f.write("\tarc "+str(0+i*16)+" "+str(6+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 5-15
			[ix,iy,iz] = arc_interp(points[5],points[15],i)
			f.write("\tarc "+str(5+i*16)+" "+str(15+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 1-7
			[ix,iy,iz] = arc_interp(points[1],points[7],i)
			f.write("\tarc "+str(1+i*16)+" "+str(7+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 4-14
			[ix,iy,iz] = arc_interp(points[4],points[14],i)
			f.write("\tarc "+str(4+i*16)+" "+str(14+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 2-8
			[ix,iy,iz] = arc_interp(points[2],points[8],i)
			f.write("\tarc "+str(2+i*16)+" "+str(8+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			# pts 3-13
			[ix,iy,iz] = arc_interp(points[3],points[13],i)
			f.write("\tarc "+str(3+i*16)+" "+str(13+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			#outlet (pts 10-11)
			[ix,iy,iz] = arc_interp(points[10],points[11],i)
			f.write("\tarc "+str(10+i*16)+" "+str(11+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			#outlet (pts 9-10)
			[ix,iy,iz] = arc_interp(points[9],points[10],i)
			f.write("\tarc "+str(9+i*16)+" "+str(10+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")
			#outlet (pts 11-12)
			[ix,iy,iz] = arc_interp(points[11],points[12],i)
			f.write("\tarc "+str(11+i*16)+" "+str(12+i*16)+" ("+str(ix)+" "+str(iy)+" "+str(iz)+")\n")

		# Leading edge arcs
		if(threed):
			f.write("\tarc "+str(0+i*16)+" "+str(5+i*16)+" ("+str(xoff+s[0][0])+" "+
				str(rad[i]*sin(s[1][0]/rad[i]))+" "+str(rad[i]*cos(s[1][0]/rad[i]))+")\n")
		else:
			f.write("\tarc "+str(0+i*16)+" "+str(5+i*16)+" ("+str(xoff+s[0][0])+" "+
				str(s[1][0])+" "+str(i*chords[i]/10.0)+")\n")
		# !!!!FIXME!!!! Trailing edge arc
		if(not blunt_te):
			f.write("\tarc "+str(2+i*16)+" "+str(3+i*16)+" ("+str(xoff+s[0][-1]*chords[i]/c_0)+" "+
				str(s[1][-1]*chords[i]/c_0)+" "+str(rad[i])+")\n")

	#pressure and suction surface splines
	for i in range(len(chords)):
		[cline,p,s] = genNACA65(chords[i],t_c,camber[i],stagger[i])
		points = generate_geom(i)
		xoff = get_xoff(i,stack,farfield)
		yoff = get_yoff(i)
		# suction surfaces
		f.write("\tspline "+str(0+i*16)+" "+str(1+i*16)+" (\n")
		for n in range(1,len(s[0])/2):
			#f.write("\t\t("+str(xoff+s[0][n]*chords[i]/c_0)+" "+str(s[1][n]*chords[i]/c_0)+" "+str(rad[i])+")\n")
			px = xoff+s[0][n]
			if(threed):
				theta = s[1][n]/rad[i]
				py = rad[i]*sin(theta)
				pz = rad[i]*cos(theta)
			else:
				py = s[1][n]
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")
		f.write("\tspline "+str(1+i*16)+" "+str(2+i*16)+" (\n")
		for n in range(len(s[0])/2,len(s[0])-1):
			#f.write("\t\t("+str(xoff+s[0][n]*chords[i]/c_0)+" "+str(s[1][n]*chords[i]/c_0)+" "+str(rad[i])+")\n")
			px = xoff+s[0][n]
			if(threed):
				theta = s[1][n]/rad[i]
				py = rad[i]*sin(theta)
				pz = rad[i]*cos(theta)
			else:
				py = s[1][n]
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		# pressure surfaces
		f.write("\tspline "+str(5+i*16)+" "+str(4+i*16)+" (\n")
		for n in range(1,len(p[0])/2):
			#f.write("\t\t("+str(xoff+p[0][n]*chords[i]/c_0)+" "+str(p[1][n]*chords[i]/c_0)+" "+str(rad[i])+")\n")
			px = xoff+p[0][n]
			if(threed):
				theta = p[1][n]/rad[i]
				py = rad[i]*sin(theta)
				pz = rad[i]*cos(theta)
			else:
				py = p[1][n]
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")
		f.write("\tspline "+str(4+i*16)+" "+str(3+i*16)+" (\n")
		for n in range(len(p[0])/2,len(p[0])-1):
			#f.write("\t\t("+str(xoff+p[0][n]*chords[i]/c_0)+" "+str(p[1][n]*chords[i]/c_0)+" "+str(rad[i])+")\n")
			px = xoff+p[0][n]
			if(threed):
				theta = p[1][n]/rad[i]
				py = rad[i]*sin(theta)
				pz = rad[i]*cos(theta)
			else:
				py = p[1][n]
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		# cyclic surfaces
		f.write("\tspline "+str(6+i*16)+" "+str(7+i*16)+" (\n")
                # first point before LE
		px = (points[6][0]+points[7][0])/2.0+xoff
		if(threed):
			theta = (2*points[6][1]+points[7][1])/3.0/rad[i]
			py = rad[i]*sin(theta)
			pz = rad[i]*cos(theta)
		else:
			py = (2*points[6][1]+points[7][1])/3.0
			pz = i*chords[i]/10.0
		f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		for n in range(1,len(cline[0])/2):
			px = xoff+cline[0][n]
			if(threed):
				theta = cline[1][n]/rad[i]
				py = rad[i]*sin(theta)+yoff
				pz = rad[i]*cos(theta)
			else:
				py = cline[1][n]+yoff
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		f.write("\tspline "+str(15+i*16)+" "+str(14+i*16)+" (\n")
                # first point before LE
		px = (points[15][0]+points[14][0])/2.0+xoff
		if(threed):
			theta = (2*points[15][1]+points[14][1])/3.0/rad[i]
			py = rad[i]*sin(theta)
			pz = rad[i]*cos(theta)
		else:
			py = (2*points[15][1]+points[14][1])/3.0
			pz = i*chords[i]/10.0
		f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		for n in range(1,len(cline[0])/2):
			px = xoff+cline[0][n]
			if(threed):
				theta = cline[1][n]/rad[i]
				py = rad[i]*sin(theta)-yoff
				pz = rad[i]*cos(theta)
			else:
				py = cline[1][n]-yoff
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		f.write("\tspline "+str(15+i*16)+" "+str(14+i*16)+" (\n")
		for n in range(1,len(cline[0])/2):
			px = xoff+cline[0][n]
			if(threed):
				theta = cline[1][n]/rad[i]
				py = rad[i]*sin(theta)-yoff
				pz = rad[i]*cos(theta)
			else:
				py = cline[1][n]-yoff
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		f.write("\tspline "+str(7+i*16)+" "+str(8+i*16)+" (\n")
		for n in range(len(cline[0])/2,len(cline[0])-1):
			px = xoff+cline[0][n]
			if(threed):
				#theta = p[1][n]/rad[i]
				theta = cline[1][n]/rad[i]
				py = rad[i]*sin(theta)+yoff
				pz = rad[i]*cos(theta)
			else:
				py = cline[1][n]+yoff
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")

		f.write("\tspline "+str(14+i*16)+" "+str(13+i*16)+" (\n")
		for n in range(len(cline[0])/2,len(cline[0])-1):
			px = xoff+cline[0][n]
			if(threed):
				#theta = p[1][n]/rad[i]
				theta = cline[1][n]/rad[i]
				py = rad[i]*sin(theta)-yoff
				pz = rad[i]*cos(theta)
			else:
				py = cline[1][n]-yoff
				pz = i*chords[i]/10.0
			f.write("\t\t("+str(px)+" "+str(py)+" "+str(pz)+")\n")
		f.write("\t)\n")
                
	f.write(");\n")

	# write blocks
	f.write("""
blocks
(
""")

	for i in range(len(chords)-1):
		#block A
		f.write("\thex (")
		for v in block_A:
			f.write(" "+str(v+i*16)+" ")
		for v in block_A:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_le)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" 1 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" 1 "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" 1 "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" 1 "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" 1 1)\n")

		#block B
		f.write("\thex (")
		for v in block_B:
			f.write(" "+str(v+i*16)+" ")
		for v in block_B:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" 1)\n")

		#block C
		f.write("\thex (")
		for v in block_C:
			f.write(" "+str(v+i*16)+" ")
		for v in block_C:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")

		#block D
		f.write("\thex (")
		for v in block_D:
			f.write(" "+str(v+i*16)+" ")
		for v in block_D:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_w)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" 1)\n")

		#block E
		f.write("\thex (")
		for v in block_E:
			f.write(" "+str(v+i*16)+" ")
		for v in block_E:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_w)+" "+str(grid_le)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_lete)+" 1 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_lete)+" 1 "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_lete)+" 1 "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_lete)+" 1 "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_lete)+" 1 1)\n")

		#block F
		f.write("\thex (")
		for v in block_F:
			f.write(" "+str(v+i*16)+" ")
		for v in block_F:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_w)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")

		#block G
		f.write("\thex (")
		for v in block_G:
			f.write(" "+str(v+i*16)+" ")
		for v in block_G:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" 1)\n")

		#block H
		f.write("\thex (")
		for v in block_H:
			f.write(" "+str(v+i*16)+" ")
		for v in block_H:
			f.write(" "+str(v+(i+1)*16)+" ")
		f.write(")\n")
		f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		"""
		if((not threed)or(i<len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		"""
		if(threed and (i==0)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_wall)+")\n")
		elif(threed and (i==len(chords)-2)):
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(1.0/exp_wall)+")\n")
		else:
			f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")

	#wingtip blocks
	"""
	if(threed):
		#block Q
		f.write("\thex (")
		f.write(str(5+i*16)+" "+str(5+(i+2)*16)+" "+str(0+(i+2)*16)+" "+str(0+i*16)+" "+
			str(5+(i+1)*16)+" "+str(11+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(0+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_le)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 1 "+str(exp_wall)+")\n")

		#block R
		f.write("\thex (")
		f.write(str(0+i*16)+" "+str(0+(i+2)*16)+" "+str(1+(i+2)*16)+" "+str(1+i*16)+" "+
			str(0+(i+1)*16)+" "+str(6+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(1+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_wall)+")\n")

		#block S
		f.write("\thex (")
		f.write(str(1+i*16)+" "+str(1+(i+2)*16)+" "+str(2+(i+2)*16)+" "+str(2+i*16)+" "+
			str(1+(i+1)*16)+" "+str(7+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(2+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")

		#block T
		f.write("\thex (")
		f.write(str(2+i*16)+" "+str(2+(i+2)*16)+" "+str(3+(i+2)*16)+" "+str(3+i*16)+" "+
			str(2+(i+1)*16)+" "+str(8+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(3+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_le)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 1 "+str(exp_wall)+")\n")

		#block U
		f.write("\thex (")
		f.write(str(3+i*16)+" "+str(3+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(4+i*16)+" "+
			str(3+(i+1)*16)+" "+str(9+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(4+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_wall)+")\n")

		#block V
		f.write("\thex (")
		f.write(str(4+i*16)+" "+str(4+(i+2)*16)+" "+str(5+(i+2)*16)+" "+str(5+i*16)+" "+
			str(4+(i+1)*16)+" "+str(10+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(5+(i+1)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")

		#block W
		f.write("\thex (")
		f.write(str(0+(i+2)*16)+" "+str(5+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(1+(i+2)*16)+" "+
			str(6+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(7+(i+2)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_wall)+")\n")

		#block X
		f.write("\thex (")
		f.write(str(1+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(3+(i+2)*16)+" "+str(2+(i+2)*16)+" "+
			str(7+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(8+(i+2)*16))
		f.write(")\n")
		f.write("\t("+str(grid_le)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
		f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_wall)+")\n")
		"""

	f.write(");\n")

	#patches
	f.write("""
boundary
(
""")

	f.write("""
	blade
	{
		type wall;
		faces
		(
""")
	if(not threed):
		i = 0
		f.write("\t\t\t("+str(0+i*16)+" "+str(0+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(5+i*16)+")\n")
		f.write("\t\t\t("+str(1+i*16)+" "+str(1+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(0+i*16)+")\n")
		f.write("\t\t\t("+str(2+i*16)+" "+str(2+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(1+i*16)+")\n")
		f.write("\t\t\t("+str(3+i*16)+" "+str(3+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(2+i*16)+")\n")
		f.write("\t\t\t("+str(4+i*16)+" "+str(4+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(3+i*16)+")\n")
		f.write("\t\t\t("+str(5+i*16)+" "+str(5+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(4+i*16)+")\n")

	else:
		#for i in range(len(chords)-2): #put this back for tip-gap!
		for i in range(len(chords)-1):
			f.write("\t\t\t("+str(0+i*16)+" "+str(0+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(5+i*16)+")\n")
			f.write("\t\t\t("+str(1+i*16)+" "+str(1+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(0+i*16)+")\n")
			f.write("\t\t\t("+str(2+i*16)+" "+str(2+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(1+i*16)+")\n")
			f.write("\t\t\t("+str(3+i*16)+" "+str(3+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(2+i*16)+")\n")
			f.write("\t\t\t("+str(4+i*16)+" "+str(4+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(3+i*16)+")\n")
			f.write("\t\t\t("+str(5+i*16)+" "+str(5+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(4+i*16)+")\n")

	"""
		# blade tip
		f.write("\t\t\t("+str(0+(i+3)*16)+" "+str(0+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(5+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(1+(i+3)*16)+" "+str(1+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(0+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(2+(i+3)*16)+" "+str(2+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(1+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(3+(i+3)*16)+" "+str(3+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(2+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(4+(i+3)*16)+" "+str(4+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(3+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(5+(i+3)*16)+" "+str(5+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(4+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(0+(i+3)*16)+" "+str(5+(i+3)*16)+" "+str(4+(i+3)*16)+" "+str(1+(i+3)*16)+")\n")
		f.write("\t\t\t("+str(1+(i+3)*16)+" "+str(4+(i+3)*16)+" "+str(3+(i+3)*16)+" "+str(2+(i+3)*16)+")\n")
"""

	f.write("""
		);
	}
""")

	f.write("""
	inlet
	{
		type inlet;
		faces
		(
""")
	for i in range(len(chords)-1):
		f.write("\t\t\t("+str(6+(i+1)*16)+" "+str(6+i*16)+" "+str(15+i*16)+" "+str(15+(i+1)*16)+")\n")
	f.write("""
		);
	}
""")
	# write createPatchDict file for cyclic boundaries
	fc = file(os.path.join(case,"system","createPatchDict"),"w")
	fc.write("""FoamFile
{
	version	"""+foam_version+""";
	format	ascii;
	class	dictionary;
	object	createPatchDict;
}
pointSync false;

patches
(
	{
		name right_periodic;
		patchInfo
		{
			type cyclic;
			neighbourPatch left_periodic;
""")
	if(threed):
		fc.write("""
		transform	rotational;
		rotationAxis	(1 0 0);
		rotationCentre	(0 0 0);
""")
	else:
		fc.write("""
		transform       translational;
		separationVector (0 """+str(2*yoff)+""" 0);
""")
	fc.write("""
			matchTolerance 0.1;
		}
		constructFrom patches;
		patches ( r_periodic );
	}

	{
		name left_periodic;
		patchInfo
		{
			type cyclic;
			neighbourPatch right_periodic;
""")
	if(threed):
		fc.write("""
			transform	rotational;
			rotationAxis	(1 0 0);
			rotationCentre	(0 0 0);
""")
	else:
		fc.write("""
			transform       translational;
			separationVector (0 """+str(2*yoff)+""" 0);
""")
	fc.write("""
			matchTolerance 0.1;
		}
		constructFrom patches;
		patches ( l_periodic );
	}
);
""");
	fc.close() 
	# periodic boundaries
	f.write("""
	r_periodic
	{
		type patch;
		""");
	f.write("""
		faces
		(
""")
	for i in range(len(chords)-1):
		f.write("\t\t\t("+str(7+(i+1)*16)+" "+str(7+i*16)+" "+str(6+i*16)+" "+str(6+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(8+(i+1)*16)+" "+str(8+i*16)+" "+str(7+i*16)+" "+str(7+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(9+(i+1)*16)+" "+str(9+i*16)+" "+str(8+i*16)+" "+str(8+(i+1)*16)+")\n")
	f.write("""
		);
	}
""")
	f.write("""
	l_periodic
	{
		type patch;
		faces
		(
""")

	for i in range(len(chords)-1):
		f.write("\t\t\t("+str(14+(i+1)*16)+" "+str(14+i*16)+" "+str(15+i*16)+" "+str(15+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(13+(i+1)*16)+" "+str(13+i*16)+" "+str(14+i*16)+" "+str(14+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(12+(i+1)*16)+" "+str(12+i*16)+" "+str(13+i*16)+" "+str(13+(i+1)*16)+")\n")

	f.write("""
		);
	}
""")
	# outside of domain
	f.write("""
	casing	
	{
""")
	if(threed):
		f.write("""
	type wall;
""")
	else:
		f.write("""
	type empty;
""")
	f.write("""
		faces
		(
""")

	if(threed):
		f.write("\t\t\t("+str(6+(i+1)*16)+" "+str(15+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(0+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(7+(i+1)*16)+" "+str(6+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(1+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(8+(i+1)*16)+" "+str(7+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(2+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(9+(i+1)*16)+" "+str(8+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(10+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(11+(i+1)*16)+" "+str(10+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(13+(i+1)*16)+" "+str(12+(i+1)*16)+" "+str(11+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(14+(i+1)*16)+" "+str(13+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(4+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(15+(i+1)*16)+" "+str(14+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(5+(i+1)*16)+")\n")

		"""
		# "blade tip" outside of domain
		f.write("\t\t\t("+str(0+(i+1)*16)+" "+str(6+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(5+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(1+(i+1)*16)+" "+str(7+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(0+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(2+(i+1)*16)+" "+str(8+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(1+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(3+(i+1)*16)+" "+str(9+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(2+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(4+(i+1)*16)+" "+str(10+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(3+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(5+(i+1)*16)+" "+str(11+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(4+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(11+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(10+(i+2)*16)+")\n")
		f.write("\t\t\t("+str(10+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(9+(i+2)*16)+")\n")
"""

	f.write("""
		);
	}
""")
	f.write("""
	outlet
	{
		type outlet;
		faces
		(
""")
	for i in range(len(chords)-1):
		f.write("\t\t\t("+str(10+(i+1)*16)+" "+str(10+i*16)+" "+str(9+i*16)+" "+str(9+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(11+(i+1)*16)+" "+str(11+i*16)+" "+str(10+i*16)+" "+str(10+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(12+(i+1)*16)+" "+str(12+i*16)+" "+str(11+i*16)+" "+str(11+(i+1)*16)+")\n")
	f.write("""
		);
	}
""")
	f.write("""
	hub
	{
""")
	if(threed):
		f.write("""
	type wall;
""")
	else:
		f.write("""
	type empty;
""")
	f.write("""
		faces
		(
			(0 5 15 6)
			(1 0 6 7)
			(2 1 7 8)
			(10 2 8 9)
			(3 11 10 2)
			(3 11 12 13)
			(4 3 13 14)
			(5 4 14 15)
""")
	if(not threed):
		f.write("\t\t\t("+str(6+(i+1)*16)+" "+str(15+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(0+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(7+(i+1)*16)+" "+str(6+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(1+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(8+(i+1)*16)+" "+str(7+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(2+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(9+(i+1)*16)+" "+str(8+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(10+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(11+(i+1)*16)+" "+str(10+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(13+(i+1)*16)+" "+str(12+(i+1)*16)+" "+str(11+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(14+(i+1)*16)+" "+str(13+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(4+(i+1)*16)+")\n")
		f.write("\t\t\t("+str(15+(i+1)*16)+" "+str(14+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(5+(i+1)*16)+")\n")
	f.write("""
		);
	}
""")
	f.write(");\n")
	f.close()
	return

#MAIN
threed = False # set this to False to generate a two-dimensional aerofoil
blunt_te = True # set this to True to make trailing edge blunt instead of radius
n_blades = 43  # number of blades
th = 0.1 # max thickness-to-chord ratio
stack = 0.25 # stacking line as a fration of chord
"""
rad = [0.150,0.165,0.180,0.195,0.21]
chords = [0.03,0.03,0.03,0.03,0.03]
camber = [31.04,23.48,17.93,13.85,10.90]
stagger = [38.0,45.0,49.40,53.00,56.10]
"""
rad = [0.180]
chords = [0.03]
camber = [17.93]
stagger = [30.00]

farfield = 0.3 #number of chord lengths to farfield boundary
case = "cascade_test" #case directory name
# create blockmesh file for wing
write_blockmesh(case,threed,blunt_te,n_blades,rad,chords,camber,stagger,th,stack,farfield=farfield,grid_le=50,grid_c = 100,grid_w = 80,exp_lete=10,grid_s=50)
# run blockMesh on the case
Runner(args=["--clear","blockMesh", "-case",case])
# run createPatch to sort out cyclic boundaries
Runner(args=["--clear","createPatch", "-overwrite", "-case",case])
# ask the user if they want to examine the mesh in paraview
viewit = raw_input("Do you want to view the mesh in ParaView? (y/N): ")
if(viewit in ["y","Y","yes","Yes","YES"]):
	Runner(args=["paraFoam","-case",case])
# run the simulation
PlotRunner(args=["--progress","simpleFoam","-case",case])
# calculate y+ wall values
Runner(args=["yPlusRAS","-case",case])
# ask the user if they want to examine the mesh in paraview
viewit = raw_input("Do you want to view the results in ParaView? (y/N): ")
if(viewit in ["y","Y","yes","Yes","YES"]):
	Runner(args=["paraFoam","-case",case])
