#!/usr/bin/env python
import os
from math import *
#import numpy
from string import split

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
    alpha:  angle of incidence (degrees)
    returns: [camberline[x][y], pressure surface[x][y], suction surface[x][y]]
"""
def genNACA65(c,t_c,theta,alpha):
    # convert angles to radians
    theta_rad = theta*pi/180.0
    alpha_rad = -alpha*pi/180.0
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
        cx = camberline[0][i]*cos(alpha_rad)-camberline[1][i]*sin(alpha_rad)
        cy = camberline[0][i]*sin(alpha_rad)+camberline[1][i]*cos(alpha_rad)
        camberline[0][i] = cx
        camberline[1][i] = cy
        # rotate pressure surface
        px = pressure[0][i]*cos(alpha_rad)-pressure[1][i]*sin(alpha_rad)
        py = pressure[0][i]*sin(alpha_rad)+pressure[1][i]*cos(alpha_rad)
        pressure[0][i] = px
        pressure[1][i] = py
        # rotate suction surface
        sx = suction[0][i]*cos(alpha_rad)-suction[1][i]*sin(alpha_rad)
        sy = suction[0][i]*sin(alpha_rad)+suction[1][i]*cos(alpha_rad)
        suction[0][i] = sx
        suction[1][i] = sy    
    return [camberline,pressure,suction]

"""    Reads in the pressure and suction surface coordinates from a text file containing
    columns representing chord fraction, pressure and suction surface coordinates.
    It then scales it based on chord-length, and rotates to specified angle of incidence
    Arguments:
    datafilename: the name of the data file
    c:      chord-length (m)
    alpha:  angle of incidence (degrees)
    delimiter: text file delimiter (defaults to a comma)
    returns: [pressure surface[x][y], suction surface[x][y]]
"""
def read_ps_profile(datafilename,c,alpha,delimiter=","):
    alpha_rad = -alpha*pi/180.0
    i = -1
    pressure = [[],[]]
    suction = [[],[]]
    datafile = file(datafilename,"r")
    for line in datafile:
        data = split(line,delimiter)
        if(len(data)>0):
            i+=1
            pressure[0].append((float(data[0])-0.5)*c)
            pressure[1].append(float(data[1])*c)
            suction[0].append((float(data[0])-0.5)*c)
            suction[1].append(float(data[2])*c)
            # rotate pressure surface
            px = pressure[0][i]*cos(alpha_rad)-pressure[1][i]*sin(alpha_rad)
            py = pressure[0][i]*sin(alpha_rad)+pressure[1][i]*cos(alpha_rad)
            pressure[0][i] = px
            pressure[1][i] = py
            # rotate suction surface
            sx = suction[0][i]*cos(alpha_rad)-suction[1][i]*sin(alpha_rad)
            sy = suction[0][i]*sin(alpha_rad)+suction[1][i]*cos(alpha_rad)
            suction[0][i] = sx
            suction[1][i] = sy 
    datafile.close()
    return [pressure, suction]

# Write a blockMesh file to the constant/polymesh directory of specified case
# usage: write_blockmesh(case,threed,blunt_te,span_c,chords,stack,sweep,dihedral,p,s,
#    farfield=5.0,grid_r=60,grid_c=40, grid_s=8,grid_le=50,grid_w=60,
#    exp_r=1000,exp_lete=4.0,exp_tip=100.0,foam_version="2.0"):
# case: name of case directory
# Wing geometry:
# threed: a boolean value denoting whether the case is a 3D wing (True) or a 2D aerofoil profile (False)
# blunt_te: a boolean value denoting whether the trailing edge is blunt (True) or an arc (False)
# span_c: a list of the spanwise distance from the wing symmetry line to each chord 
# chords: a list of the all the chord lengths
# stack: the wing stacking line as a fraction of chord
# sweep: the wing sweep in degrees
# dihedral: a list of the wing dihedrals outboard of each chord position in degrees
# p: list of pressure-surface coordinates; p[0] contains x-coords and p[1] contains y-coords
# s: list of suction-surface coordinates; s[0] contains x-coords and s[1] contains y-coords
#    Note: The first and last coordinates of p and s must represent the point of intersection
#          of the camber line and leading/trailing edge radius, and the second/second last points
#          are where the leading/trailing edge arcs end
# Number of cells along principle edges:
# grid_r: grid thickness around wing
# grid_c: number of cells along the wing chord
# grid_s: number of cells along the blade span
# grid_le: number of cells along the leading and trailing edge arcs
# grid_w: number of cells downstream of trailing edge (along wake)
# expansion ratios
# exp_r: expansion ratio from blade to boundary
# exp_lete: expansion ratio from mid-chord to leading and trailing edge
# exp_tip: expansion ratio from tip to edge of domain
# foam_version: String containing the OpenFOAM version number

def write_blockmesh(case,threed,blunt_te,span_c,chords,stack,sweep,dihedral,p,s,
    farfield=5.0,grid_r=60,grid_c=40,grid_s=8,grid_le=50,grid_w=60,
    exp_r=1000,exp_lete=4.0,exp_tip=100.0,foam_version="2.0"):

    # x and y offsets of aerofoil
    def get_xoff(chord_no,stack_frac,sweep_angle,farfield):
        if(type(sweep_angle)in(list,tuple)):
            x_off = chords[0]*farfield
            for i in range(1,chord_no+1):
                x_off+=stack_frac*(chords[i]-chords[i-1])+tan(sweep_angle[i-1]*pi/180.0)*(span_c[i]-span_c[i-1])
            return x_off
        else:
            return chords[0]*farfield+stack_frac*(chords[chord_no]-chords[0])+tan(sweep_angle*pi/180.0)*span_c[chord_no]

    def get_yoff(chord_no,dihedral,farfield):
        if(type(dihedral)in(list,tuple)):
            y_off = chords[0]*farfield
            for i in range(1,chord_no+1):
                y_off+=(span_c[i]-span_c[i-1])*tan(dihedral[i-1]*pi/180.0)
            return y_off
        else:
            return chords[0]*farfield+span_c[chord_no]*tan(dihedral*pi/180.0)

    xoff = get_xoff(0,stack,sweep,farfield)
    yoff = get_yoff(0,dihedral,farfield)
    # set values for 2d case
    if(not threed):
        num_c = 2
        chords = [chords[0],chords[0]]
        span_c = [0.0,chords[0]/10]
        dihedral = [0.0,0.0]
        grid_s = 1

    #start writing blockMeshDict file
    f = file(os.path.join(case,"constant","polyMesh","blockMeshDict"),"w")
    f.write("""FoamFile
{
    version    """+foam_version+""";
    format    ascii;
    class    dictionary;
    object    blockMeshDict;
}
convertToMeters    1.0;

vertices
(
""")

    #vertex layout scheme:
    #    6------7-----8---9
    #    |\  B  |  C  | D |
    #    | 0----1-----2---10
    #    |A(LE wing TE) E |
    #    | 5----4-----3---11
    #    |/  H  |  G  | F |
    #    15-----14----13--12

    block_A = (0,6,15,5)
    block_B = (1,7,6,0)
    block_C = (2,8,7,1)
    block_D = (10,9,8,2)
    block_E = (3,11,10,2)
    block_F = (3,13,12,11)
    block_G = (4,14,13,3)
    block_H = (5,15,14,4)

    points = [(s[0][1],s[1][1]),            #0
        (s[0][len(s[0])/2],s[1][len(s[0])/2]),    #1
        (s[0][-2],s[1][-2]),            #2
        (p[0][-2],p[1][-2]),            #3
        (p[0][len(p[0])/2],p[1][len(p[0])/2]),    #4
        (p[0][1],p[1][1]),            #5
        (s[0][0]-farfield/sqrt(2)*chords[0],farfield/sqrt(2)*chords[0]),    #6
        (s[0][len(s[0])/2],farfield*chords[0]),                    #7
        (s[0][-2]+farfield/2/sqrt(2)*chords[0],farfield*chords[0]),        #8
        (s[0][-1]+farfield*chords[0],farfield*chords[0]),            #9
        (s[0][-1]+farfield*chords[0],chords[0]*0.25),                #10
        (p[0][-1]+farfield*chords[0],chords[0]*-0.25),                #11
        (p[0][-1]+farfield*chords[0],-farfield*chords[0]),            #12
        (p[0][-2]+farfield/2/sqrt(2)*chords[0],-farfield*chords[0]),        #13
        (p[0][len(p[0])/2],-farfield*chords[0]),                #14
        (p[0][0]-farfield/sqrt(2)*chords[0],-farfield/sqrt(2)*chords[0])]    #15

    for i in range(len(chords)):
        xoff = get_xoff(i,stack,sweep,farfield)
        yoff = get_yoff(i,dihedral,farfield)
        for j in range(6):
            f.write("\t("+str(points[j][0]*chords[i]/chords[0]+xoff)+" "+str(points[j][1]*chords[i]/chords[0]+yoff)+" "+str(span_c[i])+")\n")
        for j in range(6,16):
            f.write("\t("+str(points[j][0]+xoff)+" "+str(points[j][1]+yoff)+" "+str(span_c[i])+")\n")

    # tip mesh points
    #    0--------1------2
    #   / \   R   |  S  / \
    #  /   A0____A1___A2   \
    # ( Q  |  W   |  X | T  )
    #  \   A5----A4---A3   /
    #   \ /   V   |  U  \ /
    #    5--------4------3
    if(threed):
        xoff = get_xoff(i-1,stack,sweep,farfield)
        yoff = get_yoff(i-1,dihedral,farfield)
        f.write("\t("+str(xoff+0.98*(0.7*points[0][0]+0.3*points[5][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.7*points[0][1]+0.3*points[5][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A0
        f.write("\t("+str(xoff+(0.7*points[1][0]+0.3*points[4][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.7*points[1][1]+0.3*points[4][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A1
        f.write("\t("+str(xoff+0.98*(0.7*points[2][0]+0.3*points[3][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.7*points[2][1]+0.3*points[3][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A2
        f.write("\t("+str(xoff+0.98*(0.3*points[2][0]+0.7*points[3][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.3*points[2][1]+0.7*points[3][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A3
        f.write("\t("+str(xoff+(0.3*points[1][0]+0.7*points[4][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.3*points[1][1]+0.7*points[4][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A4
        f.write("\t("+str(xoff+0.98*(0.3*points[0][0]+0.7*points[5][0])*chords[i-1]/chords[0])+" "+
            str(yoff+(0.3*points[0][1]+0.7*points[5][1])*chords[i-1]/chords[0])+" "+
            str(span_c[i-1])+")\n")#A5

        xoff = get_xoff(i,stack,sweep,farfield)
        yoff = get_yoff(i,dihedral,farfield)
        f.write("\t("+str(xoff+0.98*(0.7*points[0][0]+0.3*points[5][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.7*points[0][1]+0.3*points[5][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A6
        f.write("\t("+str(xoff+(0.7*points[1][0]+0.3*points[4][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.7*points[1][1]+0.3*points[4][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A7
        f.write("\t("+str(xoff+0.98*(0.7*points[2][0]+0.3*points[3][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.7*points[2][1]+0.3*points[3][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A8
        f.write("\t("+str(xoff+0.98*(0.3*points[2][0]+0.7*points[3][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.3*points[2][1]+0.7*points[3][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A9
        f.write("\t("+str(xoff+(0.3*points[1][0]+0.7*points[4][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.3*points[1][1]+0.7*points[4][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A10
        f.write("\t("+str(xoff+0.98*(0.3*points[0][0]+0.7*points[5][0])*chords[i]/chords[0])+" "+
            str(yoff+(0.3*points[0][1]+0.7*points[5][1])*chords[i]/chords[0])+" "+
            str(span_c[i])+")\n")#A11
    f.write(");\n")

    #write edges
    f.write("""
edges
(
""")

    for i in range(len(chords)):
        xoff = get_xoff(i,stack,sweep,farfield)
        yoff = get_yoff(i,dihedral,farfield)
        # Inlet arcs
        f.write("\tarc "+str(6+i*16)+" "+str(15+i*16)+" ("+str(-farfield*chords[0]+xoff)+" "+str(yoff+s[1][0])+" "+str(span_c[i])+")\n")
        # Leading edge arcs
        f.write("\tarc "+str(0+i*16)+" "+str(5+i*16)+" ("+str(xoff+s[0][0]*chords[i]/chords[0])+" "+
            str(yoff+s[1][0]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")
        # Trailing edge arc
        if(not blunt_te):
            f.write("\tarc "+str(2+i*16)+" "+str(3+i*16)+" ("+str(xoff+s[0][-1]*chords[i]/chords[0])+" "+
                str(yoff+s[1][-1]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")

    #pressure and suction surface splines
    for i in range(len(chords)):
        xoff = get_xoff(i,stack,sweep,farfield)
        yoff = get_yoff(i,dihedral,farfield)
        # suction surfaces
        f.write("\tspline "+str(0+i*16)+" "+str(1+i*16)+" (\n")
        for n in range(1,len(s[0])/2):
            f.write("\t\t("+str(xoff+s[0][n]*chords[i]/chords[0])+" "+str(yoff+s[1][n]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")
        f.write("\t)\n")
        f.write("\tspline "+str(1+i*16)+" "+str(2+i*16)+" (\n")
        for n in range(len(s[0])/2,len(s[0])-1):
            f.write("\t\t("+str(xoff+s[0][n]*chords[i]/chords[0])+" "+str(yoff+s[1][n]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")
        f.write("\t)\n")

        # pressure surfaces
        f.write("\tspline "+str(5+i*16)+" "+str(4+i*16)+" (\n")
        for n in range(1,len(p[0])/2):
            f.write("\t\t("+str(xoff+p[0][n]*chords[i]/chords[0])+" "+str(yoff+p[1][n]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")
        f.write("\t)\n")
        f.write("\tspline "+str(4+i*16)+" "+str(3+i*16)+" (\n")
        for n in range(len(p[0])/2,len(p[0])-1):
            f.write("\t\t("+str(xoff+p[0][n]*chords[i]/chords[0])+" "+str(yoff+p[1][n]*chords[i]/chords[0])+" "+str(span_c[i])+")\n")
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
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" 1 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" 1 "+str(exp_tip)+")\n")

        #block B
        f.write("\thex (")
        for v in block_B:
            f.write(" "+str(v+i*16)+" ")
        for v in block_B:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+"  "+str(1.0/exp_lete)+" "+str(exp_tip)+")\n")

        #block C
        f.write("\thex (")
        for v in block_C:
            f.write(" "+str(v+i*16)+" ")
        for v in block_C:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_tip)+")\n")

        #block D
        f.write("\thex (")
        for v in block_D:
            f.write(" "+str(v+i*16)+" ")
        for v in block_D:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_w)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/(exp_lete))+" "+str(exp_tip)+")\n")

        #block E
        f.write("\thex (")
        for v in block_E:
            f.write(" "+str(v+i*16)+" ")
        for v in block_E:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_w)+" "+str(grid_le)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_lete)+" 1 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_lete)+" 1 "+str(exp_tip)+")\n")

        #block F
        f.write("\thex (")
        for v in block_F:
            f.write(" "+str(v+i*16)+" ")
        for v in block_F:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_w)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_tip)+")\n")

        #block G
        f.write("\thex (")
        for v in block_G:
            f.write(" "+str(v+i*16)+" ")
        for v in block_G:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(1.0/exp_lete)+" "+str(exp_tip)+")\n")

        #block H
        f.write("\thex (")
        for v in block_H:
            f.write(" "+str(v+i*16)+" ")
        for v in block_H:
            f.write(" "+str(v+(i+1)*16)+" ")
        f.write(")\n")
        f.write("\t("+str(grid_r)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        if((not threed)or(i<len(chords)-2)):
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" 1)\n")
        else:
            f.write("\tsimpleGrading ("+str(exp_r)+" "+str(exp_lete)+" "+str(exp_tip)+")\n")

    #wingtip blocks
    if(threed):
        #block Q
        f.write("\thex (")
        f.write(str(5+i*16)+" "+str(5+(i+2)*16)+" "+str(0+(i+2)*16)+" "+str(0+i*16)+" "+
            str(5+(i+1)*16)+" "+str(11+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(0+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_le)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 1 "+str(exp_tip)+")\n")

        #block R
        f.write("\thex (")
        f.write(str(0+i*16)+" "+str(0+(i+2)*16)+" "+str(1+(i+2)*16)+" "+str(1+i*16)+" "+
            str(0+(i+1)*16)+" "+str(6+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(1+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_tip)+")\n")

        #block S
        f.write("\thex (")
        f.write(str(1+i*16)+" "+str(1+(i+2)*16)+" "+str(2+(i+2)*16)+" "+str(2+i*16)+" "+
            str(1+(i+1)*16)+" "+str(7+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(2+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_tip)+")\n")

        #block T
        f.write("\thex (")
        f.write(str(2+i*16)+" "+str(2+(i+2)*16)+" "+str(3+(i+2)*16)+" "+str(3+i*16)+" "+
            str(2+(i+1)*16)+" "+str(8+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(3+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_le)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 1 "+str(exp_tip)+")\n")

        #block U
        f.write("\thex (")
        f.write(str(3+i*16)+" "+str(3+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(4+i*16)+" "+
            str(3+(i+1)*16)+" "+str(9+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(4+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_tip)+")\n")

        #block V
        f.write("\thex (")
        f.write(str(4+i*16)+" "+str(4+(i+2)*16)+" "+str(5+(i+2)*16)+" "+str(5+i*16)+" "+
            str(4+(i+1)*16)+" "+str(10+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(5+(i+1)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le/3)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_tip)+")\n")

        #block W
        f.write("\thex (")
        f.write(str(0+(i+2)*16)+" "+str(5+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(1+(i+2)*16)+" "+
            str(6+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(7+(i+2)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(exp_lete)+" "+str(exp_tip)+")\n")

        #block X
        f.write("\thex (")
        f.write(str(1+(i+2)*16)+" "+str(4+(i+2)*16)+" "+str(3+(i+2)*16)+" "+str(2+(i+2)*16)+" "+
            str(7+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(8+(i+2)*16))
        f.write(")\n")
        f.write("\t("+str(grid_le)+" "+str(grid_c/2)+" "+str(grid_s)+")\n")
        f.write("\tsimpleGrading (1 "+str(1.0/exp_lete)+" "+str(exp_tip)+")\n")

    f.write(");\n")

    #patches
    f.write("""
boundary
(
""")

    f.write("""
    wall
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
        for i in range(len(chords)-2):
            f.write("\t\t\t("+str(0+i*16)+" "+str(0+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(5+i*16)+")\n")
            f.write("\t\t\t("+str(1+i*16)+" "+str(1+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(0+i*16)+")\n")
            f.write("\t\t\t("+str(2+i*16)+" "+str(2+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(1+i*16)+")\n")
            f.write("\t\t\t("+str(3+i*16)+" "+str(3+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(2+i*16)+")\n")
            f.write("\t\t\t("+str(4+i*16)+" "+str(4+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(3+i*16)+")\n")
            f.write("\t\t\t("+str(5+i*16)+" "+str(5+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(4+i*16)+")\n")

        # blade tip
        f.write("\t\t\t("+str(0+(i+3)*16)+" "+str(0+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(5+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(1+(i+3)*16)+" "+str(1+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(0+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(2+(i+3)*16)+" "+str(2+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(1+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(3+(i+3)*16)+" "+str(3+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(2+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(4+(i+3)*16)+" "+str(4+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(3+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(5+(i+3)*16)+" "+str(5+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(4+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(0+(i+3)*16)+" "+str(5+(i+3)*16)+" "+str(4+(i+3)*16)+" "+str(1+(i+3)*16)+")\n")
        f.write("\t\t\t("+str(1+(i+3)*16)+" "+str(4+(i+3)*16)+" "+str(3+(i+3)*16)+" "+str(2+(i+3)*16)+")\n")

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
        f.write("\t\t\t("+str(7+(i+1)*16)+" "+str(7+i*16)+" "+str(6+i*16)+" "+str(6+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(8+(i+1)*16)+" "+str(8+i*16)+" "+str(7+i*16)+" "+str(7+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(9+(i+1)*16)+" "+str(9+i*16)+" "+str(8+i*16)+" "+str(8+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(13+(i+1)*16)+" "+str(13+i*16)+" "+str(12+i*16)+" "+str(12+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(14+(i+1)*16)+" "+str(14+i*16)+" "+str(13+i*16)+" "+str(13+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(15+(i+1)*16)+" "+str(15+i*16)+" "+str(14+i*16)+" "+str(14+(i+1)*16)+")\n")

    # outside of domain
    if(threed):
        f.write("\t\t\t("+str(6+(i+1)*16)+" "+str(15+(i+1)*16)+" "+str(5+(i+1)*16)+" "+str(0+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(7+(i+1)*16)+" "+str(6+(i+1)*16)+" "+str(0+(i+1)*16)+" "+str(1+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(8+(i+1)*16)+" "+str(7+(i+1)*16)+" "+str(1+(i+1)*16)+" "+str(2+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(9+(i+1)*16)+" "+str(8+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(10+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(11+(i+1)*16)+" "+str(10+(i+1)*16)+" "+str(2+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(13+(i+1)*16)+" "+str(12+(i+1)*16)+" "+str(11+(i+1)*16)+" "+str(3+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(14+(i+1)*16)+" "+str(13+(i+1)*16)+" "+str(3+(i+1)*16)+" "+str(4+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(15+(i+1)*16)+" "+str(14+(i+1)*16)+" "+str(4+(i+1)*16)+" "+str(5+(i+1)*16)+")\n")

        # "blade tip" outside of domain
        f.write("\t\t\t("+str(0+(i+1)*16)+" "+str(6+(i+2)*16)+" "+str(11+(i+2)*16)+" "+str(5+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(1+(i+1)*16)+" "+str(7+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(0+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(2+(i+1)*16)+" "+str(8+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(1+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(3+(i+1)*16)+" "+str(9+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(2+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(4+(i+1)*16)+" "+str(10+(i+2)*16)+" "+str(9+(i+2)*16)+" "+str(3+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(5+(i+1)*16)+" "+str(11+(i+2)*16)+" "+str(10+(i+2)*16)+" "+str(4+(i+1)*16)+")\n")
        f.write("\t\t\t("+str(11+(i+2)*16)+" "+str(6+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(10+(i+2)*16)+")\n")
        f.write("\t\t\t("+str(10+(i+2)*16)+" "+str(7+(i+2)*16)+" "+str(8+(i+2)*16)+" "+str(9+(i+2)*16)+")\n")


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
    root
    {
""")
    if(threed):
        f.write("""
    type symmetryPlane;
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
