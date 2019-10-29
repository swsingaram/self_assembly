#!/usr/bin/env python

import os, sys, re, time
import math, itertools, random
from config_dict import config_dict
from hoomd import data
import gsd.hoomd

#this is modified from Oren's gen_multi.py
# TODO: use numpy arrays instead of respinning
Vdif = lambda a,b: (a[0]-b[0],a[1]-b[1],a[2]-b[2])
Vsum = lambda a,b: (a[0]+b[0],a[1]+b[1],a[2]+b[2])
Vlsq = lambda d:   d[0]*d[0] + d[1]*d[1] + d[2]*d[2]
VDistSq = lambda a,b: Vlsq(Vdif(a,b))


#Fill a grid with capsomers and polymer
def LatticeFill(N,minLateral,minZ,box,polLength,polSigma,straightPolLength):

    #for linear polyelectrolyte
	
    lowestZ = 0.5*minZ+.001 - box[-1]/2                 #last lattice position for a capsomer 
    numZ = int(box[-1]/minZ)				#number lattice positions in one dimension; equivalently, number of layers
    print "Packing, Number of rungs in z direction  %d, length of z-dimension %f, min spacing %f" %(numZ, box[-1], minZ)
    numPerLayer = int(N/numZ)+2				#Number of particles per layer
    numPerLine = math.ceil(math.sqrt(numPerLayer))	#Number of particles in a row per layer, assuming square grid at each z rung
    print "Number of particles/Layer %d, Number of particles/Line %d" %(numPerLayer,numPerLine)

    #this spreads out capsids in box a little more
    #what is the purpose of this if statement?
    if (polLength>0):    
	tempZ =  (box[2]-(minZ*2))/(int((N/((numPerLine*numPerLine)-1))))
	if (tempZ>minZ):
		minZ=tempZ
		print "New rescaled z-row separation %f" %minZ

    polymerSpot = int(numPerLine/2)                   # ?
    minLateralUse = float(box[0])/float(numPerLine+1) #how much space each particle takes up laterally
    if (polSigma>minLateral):
	minLateral=polSigma				#If the polymer particles have a larger radii, use that instead
    if minLateralUse < minLateral:
        sys.exit("particles are too close, separation %f, min separation %f" %(minLateralUse,minLateral) )
	
    #positions for capsomers
    pos = []
    for zi in xrange(numZ):
        z = zi*minZ + lowestZ          #position of z; has the form z(n) = n*del_z+ z_0
        staggerDist = 0.25*minLateralUse*(zi%2) #add stagger to every other rung. don't know why
        xi=0
        while xi < numPerLine and len(pos)<N:
            x = (xi+0.5)*minLateralUse+staggerDist     #position of x; has the form x(n) = n*del_x+x_0; x_0=0
            yi=0
            while yi < numPerLine and len(pos)<N:
                y = (yi+0.5)*minLateralUse+staggerDist  #position of y; has the form y(n) = n*del_x+x_0; x_0=0
		if ((yi != polymerSpot) or (xi != polymerSpot)):
                	pos.append([x-0.5*box[0],y-0.5*box[1],z])
                #print "xyz", xi,yi,zi
                yi+=1
            xi+=1
    if len(pos) != N:
        sys.exit("LatticeFill failed, Positions assigned %d, Total # of capsomers %d" %(len(pos),N))

    #now do polymer positions
    if (straightPolLength < box[-1]):
        for zi in xrange(polLength):
	    z = zi*polSigma + lowestZ
            #z = zi*math.sqrt(0.3333334)*polSigma + lowestZ
	    staggerDist = (zi%2)*polSigma	
    	    #staggerDist = (zi%2)*math.sqrt(0.3333334)*polSigma		#simplest staggering, using body distance
            xi=polymerSpot
            yi=polymerSpot
            x = (xi+0.5)*minLateralUse+staggerDist
            y = (yi+0.5)*minLateralUse+staggerDist
            pos.append([x-0.5*box[0],y-0.5*box[1],z])
    elif ((polSigma*2)<minLateralUse):
	xmax = ((polymerSpot+1.5)*minLateralUse)-(2*minLateral)
	xmin = ((polymerSpot-0.5)*minLateralUse)+(2*minLateral)
	ymax = ((polymerSpot+1.5)*minLateralUse)-(2*minLateral)
	ymin = ((polymerSpot-0.5)*minLateralUse)+(2*minLateral)
	zmin = lowestZ+(2*minLateral)
	zmax = box[-1]-(2*minLateral)
	xi=xmin
	yi=ymin
	zi=zmin
	xDirection=1
	yDirection=1
	scale=1.0
	pos.append([xi-0.5*box[0],yi-0.5*box[1],zi])
	print "Xmin= %s, Xmax= %s, Ymin= %s, Ymax= %s, Zmin= %s, Zmax= %s, Bond=%s" % \
	(xmin, xmax, ymin, ymax, zmin, zmax, polSigma)
	#check to see if the number of positions is equal to the 
 	#number of particles
	for ppp in xrange(polLength-1):
		xi=xi+(polSigma*scale*xDirection)
		if ((xi>xmax) or (xi<xmin)):
			xi=xi-(polSigma*scale*xDirection)
			yi=yi+(polSigma*scale*yDirection)
			xDirection=xDirection*-1
			if ((yi>ymax) or (yi<ymin)):
				yi=yi-(polSigma*scale*yDirection)
				yDirection=yDirection*-1
				zi=zi+(polSigma*scale)			
				if (zi>zmax):
					sys.exit("Grid filled up after %s units" % (ppp))
		#print "%s %s %s" % (xi, yi, zi)
		pos.append([xi-0.5*box[0],yi-0.5*box[1],zi])
    else:
	sys.exit("Polymer too large to fit in box!")

    #print "Length of pos array is %d" % len(pos) 
    return pos


class struct:
    def __init__ (self, **argw):
		self.__dict__.update (argw)
	
def GenMulti(config,seed):

	num_multi1 = int(config['Number of MultiSpheres1'])
	num_multi2 = int(config['Number of MultiSpheres2'])
	num_multi = num_multi1 + num_multi2

	atomPositionsFile1 = config['Atom Positions File1']
	atomPositionsFile2 = config['Atom Positions File2']

	num_atoms1 = int(config['Number of atoms in Species1'])
	num_atoms2 = int(config['Number of atoms in Species2'])

	box = map(float,config['Dimensions of Box'].split())
	#box = [b,b,b]
	#number of monomers
        polLength = int(float(config['Pol length']))
	#size of monomer
        polSigma = float(config['Pol sigma'])
	#charge of monomer
	polCharge = float(config['Pol charge'])

	straightPolLength = polSigma*polLength       	
 	
	minD = 6.0   #vertical spacing; where does this number come from?    
	minLateral = 2.1   #lateral spacing; where does this number comes from?
	#make the grid to populate with CPs and polymer
	multi_pos = LatticeFill(num_multi,minLateral,minD,box,polLength,polSigma,straightPolLength)

	#multi_pos = LatticeFill(num_multi,minD,box,polLength,polSigma)
	
	#list of particles
	#'C' is the central particle of the CP
	#'P' is the monomer
	particles = []
	#particle ids, 'C0'=0,'C1'=1,'P'=2
	partid = []
	#particle posititions
	partpos =[]
	#particle charges
	partcharges = []
	#bonds of polymer
	bondList=[]
	#angles of polymer
	angleList=[]
	#list of principle moments of inertia
	moment_inertia = []
	#calculated from combo.dat
	moment_diamond_bc = [93.9,73.67,30.53]
	moment_diamond_aa = [104.54,91.18,23.67]
	moment_pointpart = [0.025,0.025,0.025]

	
        for i in xrange(0, num_multi+polLength):
		#Central particles C0 = 0
		if i < num_multi1:
			particles.append('C0')
			partid.append(0)
			partpos.append(multi_pos[i])
			partcharges.append(0)
			moment_inertia.append(moment_diamond_bc)
		#Central particles C1 = 1
		elif num_multi1 <= i < num_multi:
			particles.append('C1')
			partid.append(1)
			partpos.append(multi_pos[i])
			partcharges.append(0)
			moment_inertia.append(moment_diamond_aa)
		#Polymer positions P = 2
		else:
			particles.append('P')
			partid.append(2)
			partpos.append(multi_pos[i])
			partcharges.append(polCharge)
			moment_inertia.append(moment_pointpart)
		
			if i != num_multi:
				firstIndex = len(particles)-2
               			bondList.append([firstIndex,firstIndex+1])
        		if ((i != num_multi) and (i != (num_multi+1))):
                		firstIndex = len(particles)-3
                		secondIndex = len(particles)-2
                		angleList.append([firstIndex,firstIndex+1,firstIndex+2])
	
	#Create the gsd
	s = gsd.hoomd.Snapshot()
	s.particles.N=len(particles)
	s.particles.types=['C0','C1','P','TC','TB','A1B','A2C','A2B','A3C','A5B','A5C','A3B','A4C','A6B','A6C','A4B','TA','A1A','A2A','A5A','A3A','A4A','A6A','X1','X2','Q'] #distinct type of particles, I use 'C' to denote central particle
	s.particles.typeid=partid #CPs are labelled 0 and polymer units labelled 1
	s.particles.position = multi_pos #center of CPs and polymer units
#	s.particles.charge = partcharges
	s.particles.moment_inertia = moment_inertia
	s.configuration.box=[box[0],box[1],box[2],0,0,0]
	#bonds
	s.bonds.N = len(bondList) #polLength-1
	s.bonds.types=['bondPP']
	s.bonds.group = bondList
	s.bonds.typeid = [0]*len(bondList) #(polLength-1)
	#angles
	s.angles.N = len(angleList) #polLenth-2
	s.angles.types = ['anglePPP']
	s.angles.group = angleList
	s.angles.typeid = [0]*len(angleList) #polLength-2
	#print
	gsd.hoomd.create(name='sd%s.gsd' % seed, snapshot=s)


def mol2generator(s,seed):
	p=s.particles
	molfile=open('sd%s.mol2'%seed,'w')
        molfile.write('@<TRIPOS>MOLECULE\n')
        molfile.write('GENERATED BY MOl2GeNERATOR\n')
        molfile.write('%s 1\n'%len(p))
        molfile.write('NO_CHARGES\n')
        molfile.write('@<TRIPOS>ATOM\n')
        for eachP in p:
                molfile.write('%s %s %s %s %s %s\n'%(eachP.tag+1,eachP.type,str(eachP.position[0]),str(eachP.position[1]),str(eachP.position[2]),eachP.type))
        molfile.write('@<TRIPOS>BOND\n')
        #for now
	molfile.write('1 1 2 1\n')	 
	#it should be updated by bonds information if any
	st = s.take_snapshot(bonds=True)
        bonds=st.bonds.group
        bondN=1
        for eachbond in bonds:
        	molfile.write('%s %s %s 1\n'%(bondN,eachbond[0],eachbond[1]))
        	bondN+=1

if __name__ == '__main__':
	try:
		seed = int(sys.argv[2])
		random.seed(seed)
		config = config_dict(sys.argv[1])
		#b = int(sys.argv[3])
		#box = [b,b,b]
	except:
		sys.exit("Usage: %s <input_file> <seed>" % sys.argv[0])
	
	GenMulti(config, seed)
