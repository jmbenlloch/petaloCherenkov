import random as rnd
from Scintillator import *
from Geometry import *
from Util import *
from Centella.histoManager import *
from Centella.treeManager import *
from Centella.messenger import *
from system_of_units import *
from NumpyAlgebra import *
import array 
import sys


"""
Generates photons in a point located at coordinates (0,0,0) in a system
of reference in which box1 and box2 are located as in the description below.
The box coordinates are defined by 8 vertices.
Each vertex is an (x,y,z) point. 
Define the box with the following convention:
        
        v1 = (xl,yb,zb), where xl = lefmost x coordinate,
                               yb = bottom y coordinate
                               zb = back z coordinate
        x-------x (xr,yt)
        |       |
        |       |
        x-------x            zb ----- zf
        (xl,yb)

        v2 = (xl,yt,zb), where xl = lefmost x coordinate,
                               yt = top y coordinate
                          

        v3 = (xr,yb,zb), where xr = rightmost x coordinate
        v4 = (xr,yt,zb)
        v5 = (xl,yb,zf), where zf = front z coordinate
        v6 = (xl,yt,zf)
        v7 = (xl,yt,zf)
        v8 = (xr,yb,zf)                           


"""
BOX1 =[
 [-12.8,-12.8,-1000.],[-12.8,12.8,-1000.],[12.8,-12.8,-1000.],[12.8,12.8,-1000.],
 [-12.8,-12.8,-1050.],[-12.8,12.8,-1050.],[12.8,-12.8,-1050.],[12.8,12.8,-1050.]
 ]

BOX2 =[
 [-12.8,-12.8,1000.],[-12.8,12.8,1000.],[12.8,-12.8,1000.],[12.8,12.8,1000.],
 [-12.8,-12.8,1050.],[-12.8,12.8,1050.],[12.8,-12.8,1050.],[12.8,12.8,1050.]
 ]

WBOX1 =[
 [-13.8,-13.8,-100.],[-13.8,13.8,-100.],[13.8,-13.8,-100.],[13.8,13.8,-100.],
 [-13.8,-13.8,-150.],[-13.8,13.8,-150.],[13.8,-13.8,-150.],[13.8,13.8,-150.]
 ]

WBOX2 =[
 [-13.8,-13.8,100.],[-13.8,13.8,100.],[13.8,-13.8,100.],[13.8,13.8,100.],
 [-13.8,-13.8,150.],[-13.8,13.8,150.],[13.8,-13.8,150.],[13.8,13.8,150.]
 ]


EPHOT = 511*keV
STEP = 0.2*mm

class PhotonGenerator:
	def __init__(self,box1,box2,scint='LYSO',level=0, debug=0):
		"""
		box1, box2: instances of boxes
		level = debug level
		 
		"""
		self.m = Messenger(level)
		self.box1 = box1
		self.box2 = box2 
		if scint == "LYSO":
			self.scint = LYSO()
		elif scint == "LXe":
			self.scint = LXe()
		else:
			print "scintillator %s not implemented yet"%scint
			sys.exit()

		self.debug = debug

		# distance in z from the interaction point to box1 or box2
		z = self.box2.zmin #assumes that boxes are at the same distance from IP
		x = self.box2.x/2.
		y = self.box2.y/2.

		self.txmax = x/z
		self.tymax = y/z

		self.m.log(1, "Box1 --", self.box1)
		self.m.log(1, "Box2 ---", self.box2)

		self.m.log(1, "x = %7.2f y = %7.2f z = %7.2f txmax = %7.2f tymax = %7.2f"%(
			x,y,z,self.txmax,self.tymax))


		self.ProbStep =self.scint.PhotoelectricEfficiencyAt511KeV(STEP)
		

		self.m.log(1, "Photoelectric fraction for 511 keV = %7.2f"%(
			self.scint.PhotoelectricFractionAt511KeV()))
		self.m.log(1, "Probability interaction per mm at 511 keV = %7.2g"%(
			self.scint.EfficiencyAt511KeV(STEP)))
		self.m.log(1, "Probability PE per mm at 511 keV = %7.2g"%(
			self.ProbStep))

		if self.debug > 0:
			wait()

	def GenerateEvent(self):
		"""
		Generates:
		gamma1 and gamma2 with normalized momentum P1 and P2
		by convention, P1 has negative Pz and P2 positive Pz
		"""
		P2 = self.generateMomentum(lvl=3)
		P1 = -1.*P2

		self.m.log(3," P1 =%s, P2 =%s"%(P1,P2))
			
		self.m.log(3," |P2| =%7.2f, |P1| =%7.2f, P1*P2 =%7.2f "%(
			P2.Norm(),P1.Norm(),P1*P2))

		#extrapolate P1 to box1 (z <0)
		z1 = self.box1.zmin
		x1,y1 = extrapToZ(z1,(0,0,0),(P1[0],P1[1],P1[2]))

		#extrapolate P2 to box2 (z >0)
		z2 = self.box2.zmin
		x2,y2 = extrapToZ(z2,(0,0,0),(P2[0],P2[1],P2[2]))

		
		self.m.log(3, 'x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm '%(
					x1/mm,y1/mm,z1/mm))

		self.m.log(3, 'x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm '%(
					x2/mm,y2/mm,z2/mm))

		#safety check

		if self.box1.Active((x1,y1,z1)) == False:
			self.m.log(1, '**generated point in box 1 out of box 1*** ')
			self.m.log(1, ' x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm '%(
					x1/mm,y1/mm,z1/mm))
			return False

		if self.box2.Active((x2,y2,z2)) == False:
			self.m.log(1, '**generated point in box 2 out of box 2*** ')
			self.m.log(1, ' x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm '%(
					x2/mm,y2/mm,z2/mm))
			return False

		#compute gammas time
		d1 = distance((x1,y1,z1),(0,0,0))
		d2 = distance((x2,y2,z2),(0,0,0))

		t1 = d1/c_light
		t2 = d2/c_light

		#compute path in box1 and in box2

		path1 = pathInBox((x1,y1,z1),(P1[0],P1[1],P1[2]),
			self.box1)

		path2 = pathInBox((x2,y2,z2),(P2[0],P2[1],P2[2]),
			self.box2)
				
		
		
		self.m.log(3,'Prob of int for g1 (511 keV), path =%7.2f mm) = %7.2f'%(
				path1/mm,
				self.scint.EfficiencyAt511KeV(path1)*
				self.scint.PhotoelectricFractionAt511KeV()))

		self.m.log(3,'Prob of int for g2 (511 keV), path =%7.2f mm) = %7.2f'%(
				path1/mm,
				self.scint.EfficiencyAt511KeV(path2)*
				self.scint.PhotoelectricFractionAt511KeV()))

		
		#loop until photons interact

		step1=0
		step2=0
		ntr=0
		while True: 
			ntr+=1
			#propagate photon 1:
			inter =0
			istep =0
			for step in drange(0, path1, STEP):
				self.m.log(5,'step =%7.2f mm, prob step =%7.2g '%(
					step/mm,self.ProbStep))
				#interacts?
				if rnd.uniform(0.,1.) <= self.ProbStep:
					inter =1
					istep = step
					break

			if inter == 0:
				self.m.log(4,' Photon 1 did not interact, start over ')
				continue  #fail in box1 no need to try in box2


			self.m.log(4,' Photon 1 interacts in step =%7.2f '%(step))
			step1 = step

			#propagate photon 2:
			inter =0
			istep =0
		
			for step in drange(0, path2, STEP):
				self.m.log(5,'step =%7.2f mm '%(step/mm))
				#interacts?
				if rnd.uniform(0.,1.) <= self.ProbStep:
					inter =1
					istep = step
					break

			if inter == 0:
				self.m.log(4,' Photon 2 did not interact, start over ')
				continue  

			self.m.log(4,' Photon 2 interacts in step =%7.2f  '%(step))
			step2 = step
			break

		self.m.log(2,' *Interaction* step1 = %7.2f step2 =%7.2f tried =%d time'%(
			step1,step2,ntr))

		#finally, find the interaction point in each box.
		xi1,yi1,zi1 = propagateInBox((x1,y1,z1), 
			(P1[0],P1[1],P1[2]), step1)

		xi2,yi2,zi2 = propagateInBox((x2,y2,z2), 
			(P2[0],P2[1],P2[2]), step2) 

		self.m.log(2, 'xi1 =%7.2f mm, yi1 =%7.2f mm, zi1 =%7.2f mm '%(
					xi1/mm,yi1/mm,zi1/mm))

		self.m.log(2, 'xi2 =%7.2f mm, yi2 =%7.2f mm, zi2 =%7.2f mm '%(
					xi2/mm,yi2/mm,zi2/mm))

		#gamma times

		d1 = distance((xi1,yi1,zi1),(x1,y1,z1))
		d2 = distance((xi2,yi2,zi2),(x2,y2,z2))

		ti1 = t1 + d1/c_light
		ti2 = t2 + d2/c_light

		self.m.log(2, 't1 =%7.2f ps, t2 =%7.2f ps '%(
					t1/ps,t2/ps))

		self.X1 = [x1,y1,z1,t1]
		self.X2 = [x2,y2,z2,t2]
		self.Xi1 = [xi1,yi1,zi1,ti1]
		self.Xi2 = [xi2,yi2,zi2,ti2]

		if self.debug > 0:
			wait()
		return True
		

		
	def generateMomentum(self,lvl=0):
		"""
		generate (px,py,pz) for photon 1
		px and py are generated inside box solid angle
		tx = tan(tx) = (x/2)/(z/2) = x/z  = px/pz
		where x is the length of the box and z the distance
		between the two boxes.
		ty = tan(ty) = (y/2)/(z/2) = y/z  = py/pz
		
		pz = e/sqrt(1+tx**2+ty**2)
		px = pz * tx
		py = pz * ty

		"""

		# theta = rnd.uniform(0.,2*pi)
		# phi = rnd.uniform(0.,pi)
		tx = rnd.uniform(0.,self.txmax)
		if rnd.uniform(0.,1.) > 0.5: 
			tx = -tx

		ty = rnd.uniform(0.,self.tymax)
		if rnd.uniform(0.,1.) > 0.5: 
			ty = -ty

		pz = EPHOT/sqrt(1 + tx**2 + ty**2)
		px = pz *tx
		py = pz *ty

		self.m.log(lvl, " in kev: px = %7.2g py = %7.2g pz = %7.2g phot = %7.2g "%(
			px,py,pz,EPHOT))

		self.m.log(lvl, " px**2 + py**2 + pz**2 -EPHOT = %7.2g"%(
			px**2+py**2+pz**2-EPHOT**2))

		#pz = sqrt(EPHOT**2 - px**2 - py**2)

		self.m.log(lvl, " normalized: px = %7.2g py = %7.2g pz = %7.2g "%(
			px/EPHOT,py/EPHOT,pz/EPHOT))

		P = NPVector([px/EPHOT,py/EPHOT,pz/EPHOT])

		if self.debug > 0:
			wait() 
		return P

	def __str__(self):

		s= """
		Photon Generator:
		  box 1  = %s
		  box 2 = %s
        
			"""%(self.box1,self.box2)

		return s

def Histograms(hman,box1,box2):

	# hman.h3("XYZBox1", "XYZBox1", 
	# 	10, box1.xmin, box1.xmax,
	# 	10, box1.ymin, box1.ymax,
	# 	10, box1.zmin, box1.zmax)
	# hman.fetch("XYZBox1").GetXaxis().SetTitle(
	# 	"XYZ interaction box 1 (mm)")

	hman.h2("XYBox1", "XYBox1", 
		25, box1.xmin, box1.xmax,
		25, box1.ymin, box1.ymax)
	hman.fetch("XYBox1").GetXaxis().SetTitle(
		"XY interaction box 1 (mm)")

	hman.h1("XBox1", "XBox1", 
		50, box1.xmin, box1.xmax)
	hman.fetch("XBox1").GetXaxis().SetTitle(
		"X interaction box 1 (mm)")

	hman.h1("YBox1", "YBox1", 
		50, box1.ymin, box1.ymax)
	hman.fetch("YBox1").GetXaxis().SetTitle(
		"Y interaction box 1 (mm)")

	hman.h1("ZBox1", "ZBox1", 
		50, 0, 50)
	hman.fetch("ZBox1").GetXaxis().SetTitle(
		"Z interaction box 1 (mm)")

	hman.h1("TBox1", "TBox1", 
		50, 0., 500.)
	hman.fetch("TBox1").GetXaxis().SetTitle(
		"Time interaction box 1 (ps)")

	# hman.h3("XYZBox2", "XYZBox2", 
	# 	10, box2.xmin, box2.xmax,
	# 	10, box2.ymin, box2.ymax,
	# 	10, box2.zmin, box2.zmax)
	# hman.fetch("XYZBox2").GetXaxis().SetTitle(
	# 	"XYZ interaction box 2 (mm)")

	hman.h2("XYBox2", "XYBox2", 
		25, box2.xmin, box2.xmax,
		25, box2.ymin, box2.ymax)
	hman.fetch("XYBox2").GetXaxis().SetTitle(
		"XY interaction box 2 (mm)")

	hman.h1("XBox2", "XBox2", 
		50, box2.xmin, box2.xmax)
	hman.fetch("XBox2").GetXaxis().SetTitle(
		"X interaction box 1 (mm)")

	hman.h1("YBox2", "YBox2", 
		50, box2.ymin, box2.ymax)
	hman.fetch("YBox2").GetXaxis().SetTitle(
		"Y interaction box 1 (mm)")

	hman.h1("ZBox2", "ZBox2", 
		50, 0, 50)
	hman.fetch("ZBox2").GetXaxis().SetTitle(
		"Z interaction box 1 (mm)")

	hman.h1("TBox2", "TBox2", 
		50, 0, 500.)
	hman.fetch("TBox2").GetXaxis().SetTitle(
		"Time interaction box 2 (ps)")

	# hman.h2("TBox12", "TBox12", 
	# 	25, 300, 900,
	# 	25, 300, 900)
	# hman.fetch("TBox12").GetXaxis().SetTitle(
	# 	"Time interaction box1-box2 (ps)")

	hman.h1("TDiffBox12", "TDiffBox12", 
		100, -300, 300.)
	hman.fetch("TDiffBox12").GetXaxis().SetTitle(
		"Time Diff interaction box 1-2 (ps)")

	hman.h1("XBoxG2", "XBoxG2", 
		100, box2.xmin, box2.xmax)
	hman.fetch("XBoxG2").GetXaxis().SetTitle(
		"X interaction box 1 (mm)")

	hman.h1("YBoxG2", "YBoxG2", 
		100, box2.ymin, box2.ymax)
	hman.fetch("YBoxG2").GetXaxis().SetTitle(
		"Y interaction box 1 (mm)")



if __name__ == '__main__':

	nevents = 10000
	nprint = 100
	Debug = 0
	lvl = 1
	scint = "LXe"

	m = Messenger(0)


	hman =HistoManager() 
	tman =TreeManager() 

	box1 = Box(BOX1)
	box2 = Box(BOX2)
	wbox1 = Box(WBOX1)
	wbox2 = Box(WBOX2)
	
	Histograms(hman,wbox1,wbox2)

	px1 = array.array('f',[0.])
	py1 = array.array('f',[0.])
	pz1 = array.array('f',[0.])
	pt1 = array.array('f',[0.])
	px2 = array.array('f',[0.])
	py2 = array.array('f',[0.])
	pz2 = array.array('f',[0.])
	pt2 = array.array('f',[0.])
	

	tman.book('tpg',"photon generator tree")
	tman.addBranch('tpg','px1',px1,dim=1)
	tman.addBranch('tpg','py1',py1,dim=1)
	tman.addBranch('tpg','pz1',pz1,dim=1)
	tman.addBranch('tpg','pt1',pt1,dim=1)
	tman.addBranch('tpg','px2',px2,dim=1)
	tman.addBranch('tpg','py2',py2,dim=1)
	tman.addBranch('tpg','pz2',pz2,dim=1)
	tman.addBranch('tpg','pt2',pt2,dim=1)
	
		
	pg = PhotonGenerator(box1,box2,scint=scint,level=lvl, debug=Debug)
	print pg

	nfail = 0
	nOK = 0
	
	for event in range(0,nevents):
		if event%nprint == 0:
			print 'event ', event
		evt =pg.GenerateEvent()
		if evt == False:
			nfail+=1
			continue
		nOK+=1
		x1,y1,z1,t1 = pg.Xi1
		x2,y2,z2,t2 = pg.Xi2

		xg1,yg1,zg1,tg1 = pg.X1
		xg2,yg2,zg2,tg2 = pg.X2

		m.log(1, 'event =%d x1 =%7.2f mm, y1 =%7.2f mm, z1 =%7.2f mm'%(
					event,x1/mm,y1/mm,z1/mm))

		m.log(1, 'event =%d xg1 =%7.2f mm, yg1 =%7.2f mm, zg1 =%7.2f mm'%(
					event,xg1/mm,yg1/mm,zg1/mm))

		m.log(1, 'event =%d x2 =%7.2f mm, y2 =%7.2f mm, z2 =%7.2f mm'%(
					event,x2/mm,y2/mm,z2/mm))

		m.log(1, 'event =%d xg2 =%7.2f mm, yg2 =%7.2f mm, zg2 =%7.2f mm'%(
					event,xg2/mm,yg2/mm,zg2/mm))

		m.log(1, 'event =%d t1 =%7.2f ns, t2 =%7.2f ns dt = %7.2g ps'%(
			event,t1/ps,t2/ps, abs(t1-t2)/ps))

		m.log(1, 'event =%d z1 =%7.2f mm, zg1=%7.2f mm dz = %7.2g mm'%(
			event,z1/mm,zg1/mm, abs(z1-zg1)/mm))


		if Debug > 0:
			wait() 
		

		#hman.fill("XYZBox1", x1/mm,y1/mm,z1/mm)
		hman.fill("XYBox1", x1/mm,y1/mm)
		hman.fill("XBox1", x1/mm)
		hman.fill("YBox1", y1/mm)
		hman.fill("ZBox1", abs(z1 - zg1)/mm)
		hman.fill("TBox1", (t1-tg1)/ps)

		#hman.fill("XYZBox2", x2/mm,y2/mm,z2/mm)
		hman.fill("XYBox2", x2/mm,y2/mm)
		hman.fill("XBox2", x2/mm)
		hman.fill("YBox2", y2/mm)
		hman.fill("ZBox2", abs(z2 -zg2)/mm)
		hman.fill("TBox2", (t2 - tg2)/ps)
		hman.fill("XBoxG2", xg2/mm)
		hman.fill("YBoxG2", yg2/mm)

		#hman.fill("TBox12", t1/ps,t2/ps)
		hman.fill("TDiffBox12", (t1-t2)/ps)

		px1[0]=x1/mm
		py1[0]=y1/mm
		pz1[0]=abs(z1 - zg1)/mm
		pt1[0]=(t1-tg1)/ps

		px2[0]=x2
		py2[0]=y2
		pz2[0]=abs(z2 -zg2)/mm
		pt2[0]=(t2 - tg2)/ps

		tman.fill('tpg')

	m.log(0,"nevents =%d, nfail = %d nOK =%d"%(nevents,nfail,nOK))
	pathFile = '/Users/jjgomezcadenas/Documents/Development/PETALO/histo/'
	fileName = 'vertices_'+scint+'_128_128_50.root'
	hfile = pathFile+fileName
	pathFile = '/Users/jjgomezcadenas/Documents/Development/PETALO/tree/'
	fileName = 'vertices_'+scint+'_128_128_50.root'
	tfile = pathFile+fileName
	
	hman.save(file_name=hfile)
	tman.save(file_name=tfile)


	
	

