#!/usr/bin/env python
# This script generates a text stl of a single axial flow turbomachine (compressor or turbine) blades.
# It is intended for use with cut-cell meshers (eg. cfMesh with OpenFOAM). You may need correct the facet normals.
# Meshlab is useful for this purpose

import os
from math import *
#import numpy

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
    0.3040, 0.3666, 0.4143,    0.4503, 0.4760,
    0.4924, 0.4996, 0.4963, 0.4812, 0.4530,
    0.4146, 0.3732, 0.3318, 0.2904, 0.2490,
    0.2076, 0.1662, 0.1248, 0.1000, 0.0924,
    0.0707, 0.0000)

"""    generate a set of coordinates for a NACA-65 aerofoil profile with a circular-arc camberline
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

#MAIN
threed = True # set this to False to generate a two-dimensional aerofoil
n_blades = 43  # number of blades
th = 0.1 # max thickness-to-chord ratio
stack = 0.25 # stacking line as a fration of chord
tip_gap = True # set to False if there is no tip gap
hub_gap = False # set to False if there is no hub gap

rad = [0.148,0.165,0.180,0.195,0.2095]
#rad = [0.150,0.165,0.180,0.195,0.21]
#chords = [0.03,0.03,0.03,0.03,0.03]
chords = [0.1,0.1,0.08,0.06,0.03]
camber = [31.04,23.48,17.93,13.85,10.90]
stagger = [38.0,45.0,49.40,53.00,56.10]

farfield = 0.0 #number of chord lengths to farfield boundary
xoff = farfield
yoff = 0
xpts = []
ypts = []
zpts = []

#pressure and suction surface points
for i in range(len(chords)):
    [cline,p,s] = genNACA65(chords[i],th,camber[i],stagger[i])
    xpts.append([])
    ypts.append([])
    zpts.append([])
    # suction surfaces
    for n in range(0,len(s[0])):
        theta = s[1][n]/rad[i]
        px = rad[i]*sin(theta)
        py = rad[i]*cos(theta)
        pz = xoff+s[0][n]
        xpts[-1].append(px)
        ypts[-1].append(py)
        zpts[-1].append(pz)

    # pressure surfaces
    for n in range(len(p[0])-1,-1,-1):
        theta = p[1][n]/rad[i]
        px = rad[i]*sin(theta)
        py = rad[i]*cos(theta)
        pz = xoff+p[0][n]
        xpts[-1].append(px)
        ypts[-1].append(py)
        zpts[-1].append(pz)

filename="rotorblade.stl"
def writeSTLfile(filename,xpts,ypts,zpts,n_blades,hub_gap=False,tip_gap=True,clockwise=True):
    cw = 1.0
    if(not clockwise):
        cw = -1.0
    #open stl file for output
    f = file(filename,"w")
    f.write("solid rotorblade\n")
    # loop over blades in blade row
    for m in range(n_blades):
        # calculate angle between blades
        theta=m*2*pi/n_blades
        #loop over spanwise profiles
        for c in range(1,len(xpts)):
            #loop over profile points
            for n in range(len(xpts[0])):
                # optionally rotate points around z-axis (turbomachine axis of symmetry)
                xs1=cw*(xpts[c-1][n-1]*cos(theta)-ypts[c-1][n-1]*sin(theta))
                ys1=xpts[c-1][n-1]*sin(theta)+ypts[c-1][n-1]*cos(theta)
                xs2=cw*(xpts[c-1][n]*cos(theta)-ypts[c-1][n]*sin(theta))
                ys2=xpts[c-1][n]*sin(theta)+ypts[c-1][n]*cos(theta)
                xn1=cw*(xpts[c][n-1]*cos(theta)-ypts[c][n-1]*sin(theta))
                yn1=xpts[c][n-1]*sin(theta)+ypts[c][n-1]*cos(theta)
                xn2=cw*(xpts[c][n]*cos(theta)-ypts[c][n]*sin(theta))
                yn2=xpts[c][n]*sin(theta)+ypts[c][n]*cos(theta)

                dx1 = xs2-xs1
                dx2 = xn1-xs1
                dy1 = ys2-ys1
                dy2 = yn1-ys1
                dz1 = zpts[c-1][n]-zpts[c-1][n-1]
                dz2 = zpts[c][n-1]-zpts[c-1][n-1]
                # calculate cross-products of two edges of facet
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                # use cross-product to calculate vector normals
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xs1)+" "+str(ys1)+" "+str(zpts[c-1][n-1])+"\n")
                f.write("        vertex "+str(xs2)+" "+str(ys2)+" "+str(zpts[c-1][n])+"\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[c][n-1])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")

                dx1 = xs2-xn2
                dx2 = xn1-xn2
                dy1 = ys2-yn2
                dy2 = yn1-yn2
                dz1 = zpts[c][n]-zpts[c-1][n]
                dz2 = zpts[c][n]-zpts[c][n-1]
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xn2)+" "+str(yn2)+" "+str(zpts[c][n])+"\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[c][n-1])+"\n")
                f.write("        vertex "+str(xs2)+" "+str(ys2)+" "+str(zpts[c-1][n])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")
        # optionally write facets on blade hub
        if(hub_gap):
            for n in range(2,len(xpts[-1])/2):
                xn1=cw*(xpts[0][n-1]*cos(theta)-ypts[0][n-1]*sin(theta))
                yn1=xpts[0][n-1]*sin(theta)+ypts[0][n-1]*cos(theta)
                xn2=cw*(xpts[0][n]*cos(theta)-ypts[0][n]*sin(theta))
                yn2=xpts[0][n]*sin(theta)+ypts[0][n]*cos(theta)
                xs1=cw*(xpts[0][len(xpts[0])-n]*cos(theta)-ypts[0][len(xpts[0])-n]*sin(theta))
                ys1=xpts[0][len(xpts[0])-n]*sin(theta)+ypts[0][len(xpts[0])-n]*cos(theta)
                xs2=cw*(xpts[0][len(xpts[0])-n+1]*cos(theta)-ypts[0][len(xpts[0])-n+1]*sin(theta))
                ys2=xpts[0][len(xpts[0])-n+1]*sin(theta)+ypts[0][len(xpts[0])-n+1]*cos(theta)

                dx1 = xs2-xs1
                dx2 = xn1-xs1
                dy1 = ys2-ys1
                dy2 = yn1-ys1
                dz1 = zpts[0][n]-zpts[0][n-1]
                dz2 = zpts[0][n]-zpts[0][n-1]
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xs1)+" "+str(ys1)+" "+str(zpts[0][len(xpts[0])-n])+"\n")
                f.write("        vertex "+str(xs2)+" "+str(ys2)+" "+str(zpts[0][len(xpts[0])-n+1])+"\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[0][n-1])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")

                dx1 = xn2-xs2
                dx2 = xn2-xn1
                dy1 = yn2-ys2
                dy2 = yn2-yn1
                dz1 = zpts[0][n]-zpts[0][n]
                dz2 = zpts[0][n]-zpts[0][n]
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[0][n-1])+"\n")
                f.write("        vertex "+str(xn2)+" "+str(yn2)+" "+str(zpts[0][n])+"\n")
                f.write("        vertex "+str(xs1)+" "+str(ys1)+" "+str(zpts[0][len(xpts[0])-n])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")

        # optionally write facets on blade tip
        if(tip_gap):
            for n in range(2,len(xpts[-1])/2):
                xn1=cw*(xpts[-1][n-1]*cos(theta)-ypts[-1][n-1]*sin(theta))
                yn1=xpts[-1][n-1]*sin(theta)+ypts[-1][n-1]*cos(theta)
                xn2=cw*(xpts[-1][n]*cos(theta)-ypts[-1][n]*sin(theta))
                yn2=xpts[-1][n]*sin(theta)+ypts[-1][n]*cos(theta)
                xs1=cw*(xpts[-1][len(xpts[-1])-n]*cos(theta)-ypts[-1][len(xpts[-1])-n]*sin(theta))
                ys1=xpts[-1][len(xpts[-1])-n]*sin(theta)+ypts[-1][len(xpts[-1])-n]*cos(theta)
                xs2=cw*(xpts[-1][len(xpts[-1])-n+1]*cos(theta)-ypts[-1][len(xpts[-1])-n+1]*sin(theta))
                ys2=xpts[-1][len(xpts[-1])-n+1]*sin(theta)+ypts[-1][len(xpts[-1])-n+1]*cos(theta)

                dx1 = xs2-xs1
                dx2 = xn1-xs1
                dy1 = ys2-ys1
                dy2 = yn1-ys1
                dz1 = zpts[-1][n]-zpts[-1][n-1]
                dz2 = zpts[-1][n]-zpts[-1][n-1]
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xs1)+" "+str(ys1)+" "+str(zpts[-1][len(xpts[-1])-n])+"\n")
                f.write("        vertex "+str(xs2)+" "+str(ys2)+" "+str(zpts[-1][len(xpts[-1])-n+1])+"\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[-1][n-1])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")

                dx1 = xn2-xs2
                dx2 = xn2-xn1
                dy1 = yn2-ys2
                dy2 = yn2-yn1
                dz1 = zpts[-1][n]-zpts[-1][n]
                dz2 = zpts[-1][n]-zpts[-1][n]
                nx = dy1*dz2-dz1*dy2
                ny = dz2*dx1-dx2*dz1
                nz = dx1*dy2-dx2*dy1
                vek = sqrt(nx**2+ny**2+nz**2)
                if(vek!=0):
                    nx/=vek
                    ny/=vek
                    nz/=vek
                f.write("facet normal "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
                f.write("    outer loop\n")
                f.write("        vertex "+str(xn1)+" "+str(yn1)+" "+str(zpts[-1][n-1])+"\n")
                f.write("        vertex "+str(xn2)+" "+str(yn2)+" "+str(zpts[-1][n])+"\n")
                f.write("        vertex "+str(xs1)+" "+str(ys1)+" "+str(zpts[-1][len(xpts[-1])-n])+"\n")
                f.write("    endloop\n")
                f.write("endfacet\n")
    f.write("endsolid rotorblade\n")
    f.close()

writeSTLfile(filename,xpts,ypts,zpts,n_blades,tip_gap=True,hub_gap=True,clockwise=False)
