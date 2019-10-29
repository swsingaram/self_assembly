#!/share/apps/md/bin/pyt2.7
import os, sys, re, time, math
from config_dict import config_dict
from hoomd import *
from hoomd import md
from hoomd import data
from hoomd import deprecated
context.initialize()

try:
	seed = int(sys.argv[2])
	config = config_dict(sys.argv[1])
except:
	print "Usage: %s <input_file> <seed>" % sys.argv[0]
	raise
start=0
#load gsd file to initialize system
if not os.path.exists('sd%s.gsd' % seed):
	from gen_wales_diamonds_poly import GenMulti
	GenMulti(config, seed)
	start=1
	system = init.read_gsd(filename='sd%s.gsd' % seed, restart='sd%s_a.gsd')
else:
	print(os.getcwd() + "\n")
	with open('sd%s.log'%seed) as f:
		timeStep = int(f.readlines()[-1].split('\t')[0])+1
	print('restarting from %d'%timeStep)
	system = init.read_gsd(filename='sd%s_a.gsd' % seed, restart='sd%s_b.gsd',time_step=timeStep)	

#load atomic positions
atomPositionsFile1 = config['Atom Positions File1']
atom_pos1=[]
#with open(atomPositionsFile,'r') as fPos:
fPos =  open(atomPositionsFile1,'r')
for line in fPos:
	data = line.strip().split()
	atom_pos1.append( map(float,data) )
#list of atoms in CP B-C
nameList1= ['TC']+['TB'] + ['X1'] + ['X1'] + ['A1B']+['A1B'] + ['A2C']+['A2B'] + ['A3C']+['A5B']+['A5C']+['A3B'] + ['A4C']+['A6B']+['A6C']+['A4B']+['X1']*4+['X2']*4+['Q']*12
beadsizeList1 = [1.0]*4 + [0.1]*12 + [1.0]*4 + [0.8]*4 + [0.5]*12
#load atomic positions
atomPositionsFile2 = config['Atom Positions File2']
atom_pos2=[]
#with open(atomPositionsFile,'r') as fPos:
fPos =  open(atomPositionsFile2,'r')
for line in fPos:
	data = line.strip().split()
	atom_pos2.append( map(float,data) )
#list of atoms in CP A-A
nameList2 = ['TA']*2 + ['X1']*2 + ['A1A']+['A1A'] + ['A2A']+['A2A'] + ['A5A']+['A3A']+['A3A']+['A5A'] + ['A6A']+['A4A']+['A4A']+['A6A']+['X1']*4+['X2']*4+['Q']*12
beadsizeList2 = [1.0]*4 + [0.1]*12 + [1.0]*4 + [0.8]*4 + [0.5]*12

#set up twp rigid bodies
#rigid body 1
rigid = md.constrain.rigid()
rigid.set_param('C0', positions=atom_pos1, types=nameList1, diameters=beadsizeList1)
rigid.set_param('C1', positions=atom_pos2, types=nameList2, diameters=beadsizeList2)
rigid.create_bodies(create=True)

#dump mol2 file
if start==1:
	from gen_wales_diamonds_poly import mol2generator
	mol2generator(system,seed)

ebondCC = float(config['MultiSphere Bond Strength'])

ebondCCvert = ebondCC/2.0


###Neighbor list
nl=md.nlist.tree() #bonds are exlcuded from short range interaction by default
nl.reset_exclusions(exclusions = ['bond']) 
#parameters for morse potential

#define groups
rigid = group.rigid()
centralrigid = group.rigid_center() #CPs
nonrigid = group.nonrigid() #RNA
integrableparts = group.union(name='int-particles',a=centralrigid, b=nonrigid)

r0 = 0.2
rho = float(config['rho'])
morsecut=2.0
alpha = rho/r0
morse = md.pair.morse(r_cut=morsecut, nlist=nl)
morse.set_params(mode="shift")	#This line may have been left out previously!!!
morse.pair_coeff.set(['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'],['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'],D0=0,alpha=alpha,r0=r0,r_cut=0.0)



#LJ Interactions:
sigmaTop=2.25	#completely arbitary, to prevent overlap
pol_sigma = float(config['Pol sigma'])
polX1Cut=0.5*(pol_sigma+1.0)
polX2Cut=0.5*(pol_sigma+0.8)
topRepulsion=(ebondCC/4)
topCut=sigmaTop
lj = md.pair.lj(r_cut=topCut, nlist=nl)
lj.set_params(mode="shift")
lj.pair_coeff.set(['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'],['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'],alpha=0,epsilon=0,r_cut=0,sigma=0)
lj.pair_coeff.set('TB','TB',alpha=0,epsilon=topRepulsion,r_cut=sigmaTop,sigma=sigmaTop)
lj.pair_coeff.set('TA','TC',alpha=0,epsilon=topRepulsion,r_cut=sigmaTop,sigma=sigmaTop)
#lj.pair_coeff.set(['P'],['BA','BB','BC'],alpha=1,epsilon=ebondCP,r_cut=particleSigma*3.0,sigma=particleSigma)

#even though A-A don't bind, TA-TA repulsion may be necessary to prevent wierd overlapping thing.
lj.pair_coeff.set('TA','TA',alpha=0,epsilon=topRepulsion,r_cut=1.75,sigma=1.75)
lj.pair_coeff.set('TC','TC',alpha=0,epsilon=topRepulsion,r_cut=1.75,sigma=1.75)
lj.pair_coeff.set('TB','TC',alpha=0,epsilon=topRepulsion,r_cut=1.75,sigma=1.75)
lj.pair_coeff.set('TB','TA',alpha=0,epsilon=topRepulsion,r_cut=1.75,sigma=1.75)
#lj.pair_coeff.set('P','BA',alpha=0,epsilon=1,r_cut=polQCut,sigma=polQCut)
#lj.pair_coeff.set('P','BB',alpha=0,epsilon=1,r_cut=polQCut,sigma=polQCut)
#lj.pair_coeff.set('P','BC',alpha=0,epsilon=1,r_cut=polQCut,sigma=polQCut)
lj.pair_coeff.set('P','Q',alpha=0,epsilon=1,r_cut=pol_sigma, sigma=pol_sigma)
lj.pair_coeff.set('P','P',alpha=0,epsilon=1,r_cut=pol_sigma, sigma=pol_sigma)
lj.pair_coeff.set('P',['C0','C1','X1'],alpha=0,epsilon=1,r_cut=polX1Cut,sigma=polX1Cut)
lj.pair_coeff.set('P','X2',alpha=0,epsilon=1,r_cut=polX2Cut,sigma=polX2Cut)
lj.pair_coeff.set('X1','X1',alpha=0,epsilon=1,r_cut=1.0, sigma=1.0)
lj.pair_coeff.set('X2','X2',alpha=0,epsilon=1,r_cut=0.8,sigma=0.8)
lj.pair_coeff.set('X1','X2',alpha=0,epsilon=1,r_cut=0.9,sigma=0.9)



#parameters for polymer
unitLength=1 
pol_charge = float(config['Pol charge'])
polBend=0

#polymer
harmonic = md.bond.harmonic()
harmonic.bond_coeff.set('bondPP', k=330.0*unitLength**2, r0=pol_sigma)
harmonic = md.angle.harmonic()
harmonic.angle_coeff.set('anglePPP', k=polBend, t0=3.1415)

#yukawa parameters
chargeQ = float(config['Bead charge'])
lDebye = float(config['Debye'])
kappa = 1.0/lDebye
bjerrum= 0.7
yukawaPrefactor = bjerrum * lDebye * math.exp(pol_sigma/lDebye)/(lDebye+pol_sigma)
rc = 3.0*lDebye
ron = 2.0*lDebye
print "For Debye length of %f rb, Yukawa Kappa = %f and Yukawa Epsilon = %f" % (lDebye, kappa, yukawaPrefactor)

#polymer-CP interactions
yukawa = md.pair.yukawa(r_cut=rc, nlist=nl)
yukawa.set_params(mode="xplor")
yukawa.pair_coeff.set(['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'],['C0','C1','P','TA','TB','TC','A1A','A1B','A2A','A2B','A3A','A3B','A4A','A4B','A5A','A5B','A6A','A6B','A1C','A2C','A3C','A4C','A5C','A6C','X1','X2','Q'], epsilon=0, kappa=0, r_cut=0, r_on=0)

yukawa.pair_coeff.set('P', 'P', epsilon=yukawaPrefactor, kappa=kappa, r_cut=rc, r_on=ron)


#### NOW RUN POLYMER + CPS WITHOUT ANY ATTRACTIVE INTERACTIONS
#### TO GET THE SYSTEM INTO A MORE REASONABLE
#### CONFIGURATION
# Set the integrator
ts = float(config['Time Step Length'])
if start == 1:
        firepoly = md.integrate.mode_minimize_fire(dt=ts*0.1, ftol=1e-2, Etol=1e-7)
	nve = md.integrate.nve(group=nonrigid)
        #while not(firepoly.has_converged()):
        #    run(100)
        run(5000)
	nve.disable()
	del nve
        imode = md.integrate.mode_standard(dt=ts,aniso=True)
        langevin = md.integrate.langevin(group=integrableparts,  kT=1.0, seed=seed, dscale=1.0)
        anim_period = 100000.0/ts
        dcd = dump.dcd(filename="sd%s.dcd" % seed, period = anim_period*0.025)
        run(50000000)
        langevin.disable()
        del langevin #hoomd complains if you simply disable the integrator


###TURN ON INTERACTIONS

morse.pair_coeff.set('A1','A1',D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set('A2','A2',D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set('A3','A5',D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set('A4','A6',D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)

#Pentameric
morse.pair_coeff.set(['A1B'],['A1B'],D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A2B'],['A2B'],D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A3B'],['A5B'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A4B'],['A6B'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)

#Hexameric
#there is no A1C...
morse.pair_coeff.set(['A1A'],['A1B'],D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A2A'],['A2C'],D0=ebondCCvert,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A3A'],['A5C'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A4A'],['A6C'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A5A'],['A3C'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)
morse.pair_coeff.set(['A6A'],['A4C'],D0=ebondCC,alpha=alpha,r0=r0,r_cut=2.0)


#yukawa
yukawa.pair_coeff.set('Q', 'Q', epsilon=yukawaPrefactor*chargeQ**2, kappa=kappa, r_cut=rc, r_on=ron)
#yukawa.pair_coeff.set('P', 'P', epsilon=yukawaPrefactor, kappa=kappa, r_cut=rc, r_on=ron)
yukawa.pair_coeff.set('P', 'Q', epsilon=-yukawaPrefactor*chargeQ, kappa=kappa, r_cut=rc, r_on=ron)


imode = md.integrate.mode_standard(dt=ts,aniso=True)
langevin = md.integrate.langevin(group=integrableparts,  kT=1.0, seed=seed, dscale=1.0)
if start != 1:
        anim_period = 100000.0/ts
        dcd = dump.dcd(filename="sd%s.dcd" % seed, period = anim_period*0.025)
logger = analyze.log(filename="sd%s.log" % seed, quantities=['temperature', 'potential_energy'], period=anim_period)
if start==1:
        gsd_restart = dump.gsd(filename='sd%s_a.gsd' % seed, truncate=True, group=group.all(),period = anim_period , phase=0)
else:
        gsd_restart1 = dump.gsd(filename='sd%s_a.gsd' % seed, truncate=True,  group=group.all(), period = anim_period , phase=0)
        gsd_restart2 = dump.gsd(filename='sd%s_b.gsd' % seed, truncate=True, group=group.all(),period = anim_period , phase=0)
runlen = anim_period*20
#runlen = anim_period*10
run(runlen, limit_hours=59.9, limit_multiple=10)
#run(runlen)
