#!/usr/bin/env python
import fileinput
import math
import sys
import csv
import itertools
import networkx as nx
import numpy as np
BOX = 172 #cubic box with length 40
TRAJFILE = './temp.mol2'
FRAMEA = 0
FRAMEB = 899 #899 frames in total
NUMCP = 450 #158
ATOMSINCP = 36

#variables for dumping
EBOND = float(sys.argv[1])
POLLENGTH = float(sys.argv[2]) #800
DEBYE = float(sys.argv[3])
ITR = int(sys.argv[4])
numboundfile = 'numbound_%s_%s_%s_%s.dat' % (EBOND, POLLENGTH,DEBYE,ITR)
clusterfile = 'cluster_%s_%s_%s_%s.dat' % (EBOND, POLLENGTH,DEBYE,ITR)
######
CUTOFFCPRNA = 3.0
CUTOFFCPCP  = 4*0.1
#Beads, 2A-2C and 4A-4C are too close
#to other diamonds so overcounting results
#this is a temp fix.  clean up adjacency matrix
#after
CUTOFFCPCP2 = CUTOFFCPCP #0.9
CUTOFFCPCP3 = CUTOFFCPCP #0.87

#initialize system
lineA = FRAMEA*(10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP) - 1
lineB = FRAMEB*(10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP) - 1
if FRAMEA == 0:
   lineA = 1
class Bead:
    def __init__(self):
        self.indices = []
        self.posns = []
    def add_index(self,index):
        self.indices.append(index)
    def add_posn(self,posn):
        self.posns.append(posn)
def Vect(v):
    return map(float,v)
def Vlen(v):
    return math.sqrt(v[0]**2+v[1]**2+v[2]**2)
def Vmod(a):
    if (a/BOX - math.floor(a/BOX)) < 0.5:
        return (a - BOX*(math.floor(a/BOX)))
    else:
        return (a - BOX*(math.ceil(a/BOX)))
def Vdef(u,v):
    return ([Vmod(u[0]-v[0]),Vmod(u[1]-v[1]),Vmod(u[2]-v[2])])
def QPBound(filename,linenumA,QPcutoff):
    Qposns = []
    Pposns = []
    BoundQs = []
    lineend = linenumA + 10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP - 1
    with open(filename) as fin:
        lines = itertools.islice(fin,linenumA,lineend)
        for line in lines:
           line = line.strip()
           col = line.split()
           if (len(col) > 2) and (col[1] == 'C0'):
               Qposns.append(Vect([col[2],col[3],col[4]]))
	   if (len(col) > 2) and (col[1] == 'C1'):
	       Qposns.append(Vect([col[2],col[3],col[4]]))
           elif (len(col) > 2) and (col[1] == 'P'):
               Pposns.append(Vect([col[2],col[3],col[4]]))
    for i in Qposns:
        for j in Pposns:
            if Vlen(Vdef(i,j)) <= QPcutoff:
                BoundQs.append(i) #keep record of list for analysis
    #get rid of duplicates
    BoundQs.sort()
    BoundQs = list(BoundQs for BoundQs,_ in itertools.groupby(BoundQs))
    return (len(BoundQs))
def CPBound(filename,linenumA,CPCPcutoff):
    A1A = Bead() 
    A1B = Bead()
    A2A = Bead()
    A2B = Bead()
    A2C = Bead()
    A3A = Bead()
    A3B = Bead()
    A3C = Bead()
    A4A = Bead()
    A4B = Bead()
    A4C = Bead()
    A5A = Bead()
    A5B = Bead()
    A5C = Bead()
    A6A = Bead()
    A6B = Bead()
    A6C = Bead()
    lineend = linenumA + 10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP - 1
    with open(filename) as fin:
        lines = itertools.islice(fin,linenumA,lineend)
        cpindex = -1 #cpindex will start at 0;same numbering as cpmat
        counter = 0
        for line in lines:
            line = line.strip()
            col = line.split()
            if (len(col) > 2) and (col[1] == 'TA' or col[1] == 'TB' or col[1] == 'TC'):
                counter += 1
                if counter == 2:
                    cpindex += 1
                    counter = 0
            if (len(col) > 2) and (col[1] == 'A1A'):
                A1A.add_posn(Vect([col[2],col[3],col[4]]))
                A1A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A1B'):
                A1B.add_posn(Vect([col[2],col[3],col[4]]))
                A1B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A2A'):
                A2A.add_posn(Vect([col[2],col[3],col[4]]))
                A2A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A2B'):
                A2B.add_posn(Vect([col[2],col[3],col[4]]))
                A2B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A2C'):
                A2C.add_posn(Vect([col[2],col[3],col[4]]))
                A2C.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A3A'):
                A3A.add_posn(Vect([col[2],col[3],col[4]]))
                A3A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A3B'):
                A3B.add_posn(Vect([col[2],col[3],col[4]]))
                A3B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A3C'):
                A3C.add_posn(Vect([col[2],col[3],col[4]]))
                A3C.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A4A'):
                A4A.add_posn(Vect([col[2],col[3],col[4]]))
                A4A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A4B'):
                A4B.add_posn(Vect([col[2],col[3],col[4]]))
                A4B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A4C'):
                A4C.add_posn(Vect([col[2],col[3],col[4]]))
                A4C.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A5A'):
                A5A.add_posn(Vect([col[2],col[3],col[4]]))
                A5A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A5B'):
                A5B.add_posn(Vect([col[2],col[3],col[4]]))
                A5B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A5C'):
                A5C.add_posn(Vect([col[2],col[3],col[4]]))
                A5C.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A6A'):
                A6A.add_posn(Vect([col[2],col[3],col[4]]))
                A6A.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A6B'):
                A6B.add_posn(Vect([col[2],col[3],col[4]]))
                A6B.add_index(cpindex)
            elif (len(col) > 2) and (col[1] == 'A6C'):
                A6C.add_posn(Vect([col[2],col[3],col[4]]))
                A6C.add_index(cpindex)
    #A1-A1; A2-A2; A3-A5; A4-A6
    #Verified in VMD.  Optimal distances are:
    #     A1A-A1B --> 0.4
    #     A1B-A1B --> 0.2
    #     A3A-A5C --> 0.4
    #     A3B-A5B --> 0.4
    #     A4A-A6C --> 0.4
    #     A4B-A6B --> 0.4
    ctmat = np.zeros( (NUMCP,NUMCP), dtype=np.int16 )
   # ctmat = [[0 for i in range(NUMCP)] for j in range(NUMCP)]
    #loop over A1-A1 interactins
    for i in range(len(A1A.posns)):
        for j in range(len(A1B.posns)):
            if not(A1A.indices[i] == A1B.indices[j]): #not on same cp!
                if Vlen(Vdef(A1A.posns[i],A1B.posns[j])) <= CPCPcutoff:
                    ctmat[A1A.indices[i]][A1B.indices[j]] = 1
                    ctmat[A1B.indices[j]][A1A.indices[i]] = 1
#    for i in range(len(A1B.posns)):
 #       for j in range(len(A1B.posns)):
  #          if not(A1B.indices[i] == A1B.indices[j]):
   #             if Vlen(Vdef(A1B.posns[i],A1B.posns[j])) <= 0*CPCPcutoff:
    #                ctmat[A1B.indices[i]][A1B.indices[j]] = 1
     #               ctmat[A1B.indices[j]][A1B.indices[i]] = 1
#    for i in range(len(A2A.posns)):
 #       for j in range(len(A2C.posns)):
 #           if not(A2A.indices[i] == A2C.indices[j]):
  #              if Vlen(Vdef(A2A.posns[i],A2C.posns[j])) <= 0*CPCPcutoff: #0.87fix this
   #                 ctmat[A2A.indices[i]][A2C.indices[j]] = 1
    #                ctmat[A2C.indices[j]][A2A.indices[i]] = 1
  #  for i in range(len(A2B.posns)):
   #     for j in range(len(A2B.posns)):
    #        if not(A2B.indices[i] == A2B.indices[j]):
     #           if Vlen(Vdef(A2B.posns[i],A2B.posns[j])) <= 0*CPCPcutoff:  #0.87
      #              ctmat[A2B.indices[i]][A2B.indices[j]] = 1
       #             ctmat[A2B.indices[j]][A2B.indices[i]] = 1
    for i in range(len(A3A.posns)):
        for j in range(len(A5C.posns)):
            if not(A3A.indices[i] == A5C.indices[j]):
                if Vlen(Vdef(A3A.posns[i],A5C.posns[j])) <= CPCPcutoff:
                    ctmat[A3A.indices[i]][A5C.indices[j]] = 1
                    ctmat[A5C.indices[j]][A3A.indices[i]] = 1
    for i in range(len(A3B.posns)):
        for j in range(len(A5B.posns)):
            if not(A3B.indices[i] == A5B.indices[j]):
                if Vlen(Vdef(A3B.posns[i],A5B.posns[j])) <= CPCPcutoff:
                    ctmat[A3B.indices[i]][A5B.indices[j]] = 1
                    ctmat[A5B.indices[j]][A3B.indices[i]] = 1
#    for i in range(len(A4A.posns)):
 #       for j in range(len(A6C.posns)):
  #          if not(A4A.indices[i] == A6C.indices[j]):
   #             if Vlen(Vdef(A4A.posns[i],A6C.posns[j])) <= 0*CPCPcutoff: # 0.9fix this 
    #                ctmat[A4A.indices[i]][A6C.indices[j]] = 1
     #               ctmat[A6C.indices[j]][A4A.indices[i]] = 1
#    for i in range(len(A4B.posns)):
 #       for j in range(len(A6B.posns)):
  #          if not(A4B.indices[i] == A6B.indices[j]):
  #              if Vlen(Vdef(A4B.posns[i],A6B.posns[j])) <= 0*CPCPcutoff:
   #                 ctmat[A4B.indices[i]][A6B.indices[j]] = 1
    #                ctmat[A6B.indices[j]][A4B.indices[i]] = 1
    return(ctmat)
#main program
#loop over QPBound to compute the number of bound
#CPs to the polymer from frameA through frameB.
frame=FRAMEA
data = []
clstrdata = []
line = lineA
if line == 1:
    numbound = QPBound(TRAJFILE,line,CUTOFFCPRNA)
    adjmat = CPBound(TRAJFILE,line,CUTOFFCPCP)
    graph = nx.from_numpy_matrix(adjmat)
    clusters = [len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
    data.append([frame,numbound,NUMCP,BOX**3])
    clstrdata.append(clusters)
    frame += 1
    line = line + 10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP - 2
while line <= lineB:
    numbound = QPBound(TRAJFILE,line,CUTOFFCPRNA)
    adjmat = CPBound(TRAJFILE,line,CUTOFFCPCP)
    graph = nx.from_numpy_matrix(adjmat)
    clusters = [len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
    data.append([frame,numbound,NUMCP,BOX**3])
    clstrdata.append(clusters)
    frame += 1
    line = line + 10 + 2 + NUMCP + POLLENGTH + NUMCP*ATOMSINCP
    
with open(numboundfile,'w') as csvfile:
    writer = csv.writer(csvfile,delimiter=' ')
    for row in data:
        writer.writerow(row)
with open(clusterfile,'w') as csvfile:
    writer = csv.writer(csvfile,delimiter=' ')
    for row in clstrdata:
        writer.writerow(row)

