from Centella.AAlgo import AAlgo

from Particles import *
from Util import *
from Geometry import *
from TOF import *


"""
This algorithm defines the fiducial volume for the CRT analysis. 
"""


class AFiducial(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		AFiducial  Algorithm
		"""
		
		self.name = 'AFiducial'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    

	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')

		### PARAMETERS
		
		self.debug = self.ints["Debug"]  #used to stop program at key break points

		self.xres = self.doubles["XRES"]  #Resolution in the X coordinate in FHWM
		self.yres = self.doubles["YRES"]  #Resolution in the Y coordinate in FWHM
		self.zres = self.doubles["ZRES"]  #Resolution in the Z coordinate in FWHM

		#load coordinates of box and fiducial box

		boxCoord1 =self.loadCoord("Box1V")
  		boxCoord2 =self.loadCoord("Box2V")
  		fboxCoord1 =self.loadCoord("FBox1V")
  		fboxCoord2 =self.loadCoord("FBox2V")

		self.box1 = Box(boxCoord1)
		self.box2 = Box(boxCoord2)
		self.fbox1 = Box(fboxCoord1)
		self.fbox2 = Box(fboxCoord2)

		self.HBOX = self.ints["HBOX"] # if 0 no histos, if 1 only box 1 if 2 also box2
		
		self.m.log(1, "Box1 --", self.box1)
		self.m.log(1, "Box2 ---", self.box2)

		self.m.log(1, "Fiducial Box1 --", self.box1)
		self.m.log(1, "Fiducial Box2 ---", self.box2)
	
		if self.debug == 1:
			wait()
		
		### Defining histos
		
		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

		#Fiducial cut. There must be one photon in the fiducial volue
		#(defined by parameters) in each box. return the interaction vertex
		#of the photon

		fiducial = self.Fiducial(event)
		self.hman.fill(self.Fiducial_histo_name,fiducial)
		
		if self.debug == 1:
			wait()

		if fiducial != 2:
			return False

		self.numOutputEvents += 1
		
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input Events: ', self.numInputEvents)
		self.m.log(1, 'Pass Fiducial cuts Events: ', self.numOutputEvents)
		
		return

############################################################
	def Fiducial(self,event):
		"""
		This method checks that there is one photon in the fiducial volue
		(defined by parameters) in each box. It returns 
		fid = 0 if no photon is found in the fiduical volume
		fid = 1 if one photon found in one box
		fid = 2 if one photon found in each box
		vertexBox1, vertexBox2 the vertices of the photons
		"""

		primaryParticles = PrimaryParticles(event)
		self.m.log(2, ' number of primary Particles =%d '%(len(primaryParticles)))
		
		vertexBox1=Point4D()
		vertexBox2=Point4D()
		fid = 0
		for pparticle in primaryParticles:
			self.m.log(2, '\n+++primary particle+++\n')
			ei,ef = particleKineticEnergy(pparticle)

			self.m.log(2,'name = %s t =%7.2f ps E = %7.2f keV'%(
			particleName(pparticle), particleTime(pparticle)/ps,ei/keV))

			x0,y0,z0 = particleInitialVtx(pparticle)
			x,y,z = particleFinalVtx(pparticle)
			t = particleTime(pparticle)
			
			self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm, t = %7.2f ps '%(
			x0/mm,y0/mm,z0/mm, t/ps))
			self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm, t = %7.2f ps '%(
			x/mm,y/mm,z/mm, t/ps))

			if self.fbox1.Active((x,y,z)) == True:
				self.m.log(2,'gamma found in box1')
				
				if self.HBOX ==1 :
					self.hman.fill(self.XYZBox1_histo_name, 
					x/mm,y/mm,z/mm)
					self.hman.fill(self.XYBox1_histo_name, 
					x/mm,y/mm)
					self.hman.fill(self.XBox1_histo_name, 
					x/mm)
					self.hman.fill(self.YBox1_histo_name, 
					y/mm)
					self.hman.fill(self.ZBox1_histo_name, 
					z/mm)
					self.hman.fill(self.TBox1_histo_name, 
					t/ps)

                                xSmear, ySmear, zSmear = SmearVertex((x,y,z,), self.xres, self.yres, self.zres)
                                #print "x: %s, xS: %s, y: %s, yS: %s, z: %s, zS: %s" %(x,xSmear,y,ySmear,z,zSmear)
				#vertexBox1.x = x 
				#vertexBox1.y = y
				#vertexBox1.z = z
				vertexBox1.x = xSmear
				vertexBox1.y = ySmear
				vertexBox1.z = zSmear
				vertexBox1.t = t

				fid+=1
			
			elif self.fbox2.Active((x,y,z)) == True:
				self.m.log(2,'gamma found in box2')
				
				if self.HBOX ==2 :
					self.hman.fill(self.XYZBox2_histo_name, 
					x/mm,y/mm,z/mm)
					self.hman.fill(self.XYBox2_histo_name, 
					x/mm,y/mm)
					self.hman.fill(self.ZBox2_histo_name, 
					z/mm)
					self.hman.fill(self.XBox2_histo_name, 
					x/mm)
					self.hman.fill(self.YBox2_histo_name, 
					y/mm)
					self.hman.fill(self.TBox2_histo_name, 
					t/ps)

                                xSmear, ySmear, zSmear = SmearVertex((x,y,z,), self.xres, self.yres, self.zres)
                                #print "x: %s, xS: %s, y: %s, yS: %s, z: %s, zS: %s" %(x,xSmear,y,ySmear,z,zSmear)
				#vertexBox2.x = x
				#vertexBox2.y = y
				#vertexBox2.z = z
				vertexBox2.x = xSmear
				vertexBox2.y = ySmear
				vertexBox2.z = zSmear
				vertexBox2.t = t

				fid+=1
			else:
				self.m.log(2,'gamma not found in box1 or in box2')
				self.m.log(2, ' x0 =%7.2f mm, y0 =%7.2f mm, z0 =%7.2f mm '%(
			x0/mm,y0/mm,z0/mm))
				self.m.log(2, ' xf =%7.2f mm, yf =%7.2f mm, zf =%7.2f mm '%(
			x/mm,y/mm,z/mm))
				self.m.log(2, "fBox1 --", self.fbox1)
				self.m.log(2, "fBox2 ---", self.fbox2)
		
		self.logman["USER"].gparam["vertexBox"]=[vertexBox1,vertexBox2]
		
		self.m.log(2,' fid =%7.2f'%(fid))
		
		return fid


	############################################################
	def loadCoord(self, blb):
		boxCoord =[]
		for i in range(1,9):
			lbl = blb+str(i)
			self.m.log(3, "loading parameter", lbl)
			coord  = self.vdoubles[lbl]
			try:
				boxCoord.append(coord);
			except KeyError:
				self.m.log(1, "WARNING!! Parameter %s not defined."%(lbl))
				exit(0)
		return boxCoord


	############################################################		
	def Histos(self):
		"""
		book the histograms for the algo 
		"""
		Fiducial_histo_desc = "Fiducial"
		self.Fiducial_histo_name = self.alabel(Fiducial_histo_desc)
		self.hman.h1(self.Fiducial_histo_name, Fiducial_histo_desc, 
			10, 0, 5)
		self.hman.fetch(
			self.Fiducial_histo_name).GetXaxis().SetTitle(
			"Fiducial interactions")

		if self.HBOX ==1 :
			XYZBox1_histo_desc = "XYZBox1"
			self.XYZBox1_histo_name = self.alabel(XYZBox1_histo_desc)
			self.hman.h3(self.XYZBox1_histo_name, XYZBox1_histo_desc, 
				10, self.box1.xmin, self.box1.xmax,
           		10, self.box1.ymin, self.box1.ymax,
           		10, self.box1.zmin, self.box1.zmax)
			self.hman.fetch(
				self.XYZBox1_histo_name).GetXaxis().SetTitle(
				"XYZ interaction box 1 (mm)")

			XYBox1_histo_desc = "XYBox1"
			self.XYBox1_histo_name = self.alabel(XYBox1_histo_desc)
			self.hman.h2(self.XYBox1_histo_name, XYBox1_histo_desc, 
			25, self.box1.xmin, self.box1.xmax,
           	25, self.box1.ymin, self.box1.ymax)  	
			self.hman.fetch(
			self.XYBox1_histo_name).GetXaxis().SetTitle(
			"XY interaction box 1 (mm)")

			XBox1_histo_desc = "XBox1"
			self.XBox1_histo_name = self.alabel(XBox1_histo_desc)
			self.hman.h1(self.XBox1_histo_name, XBox1_histo_desc, 
			25, self.box1.xmin, self.box1.xmax)
			self.hman.fetch(
			self.XBox1_histo_name).GetXaxis().SetTitle(
			"X interaction box 1 (mm)")

			YBox1_histo_desc = "YBox1"
			self.YBox1_histo_name = self.alabel(YBox1_histo_desc)
			self.hman.h1(self.YBox1_histo_name, YBox1_histo_desc, 
			25,self.box1.ymin, self.box1.ymax)
			self.hman.fetch(
			self.YBox1_histo_name).GetXaxis().SetTitle(
			"Y interaction box 1 (mm)")

			ZBox1_histo_desc = "ZBox1"
			self.ZBox1_histo_name = self.alabel(ZBox1_histo_desc)
			self.hman.h1(self.ZBox1_histo_name, ZBox1_histo_desc, 
			25, self.box1.zmin, self.box1.zmax)
			self.hman.fetch(
			self.ZBox1_histo_name).GetXaxis().SetTitle(
			"Z interaction box 1 (mm)")

			TBox1_histo_desc = "TBox1"
			self.TBox1_histo_name = self.alabel(TBox1_histo_desc)
			self.hman.h1(self.TBox1_histo_name, TBox1_histo_desc, 
			50, 300, 600)
			self.hman.fetch(
			self.TBox1_histo_name).GetXaxis().SetTitle(
			"Time of gamma in  box 1 (ps)")
		
		if self.HBOX ==2 :
			XYZBox2_histo_desc = "XYZBox2"
			self.XYZBox2_histo_name = self.alabel(XYZBox2_histo_desc)
			self.hman.h3(self.XYZBox2_histo_name, XYZBox2_histo_desc, 
			10, self.box2.xmin, self.box2.xmax,
           	10, self.box2.ymin, self.box2.ymax,
           	10, self.box2.zmin, self.box2.zmax)
			self.hman.fetch(
			self.XYZBox2_histo_name).GetXaxis().SetTitle(
			"XYZ interaction box 2 (mm)")

			XYBox2_histo_desc = "XYBox2"
			self.XYBox2_histo_name = self.alabel(XYBox2_histo_desc)
			self.hman.h2(self.XYBox2_histo_name, XYBox2_histo_desc, 
			25, self.box2.xmin, self.box2.xmax,
           	25, self.box2.ymin, self.box2.ymax)
			self.hman.fetch(
			self.XYBox2_histo_name).GetXaxis().SetTitle(
			"XY interaction box 2 (mm)")

			XBox2_histo_desc = "XBox2"
			self.XBox2_histo_name = self.alabel(XBox2_histo_desc)
			self.hman.h1(self.XBox2_histo_name, XBox2_histo_desc, 
			25, self.box2.xmin, self.box2.xmax)
			self.hman.fetch(
			self.XBox2_histo_name).GetXaxis().SetTitle(
			"X interaction box 2 (mm)")

			YBox2_histo_desc = "YBox2"
			self.YBox2_histo_name = self.alabel(YBox2_histo_desc)
			self.hman.h1(self.YBox2_histo_name, YBox2_histo_desc, 
			25,self.box2.ymin, self.box2.ymax)
			self.hman.fetch(
			self.YBox2_histo_name).GetXaxis().SetTitle(
			"Y interaction box 2 (mm)")

			ZBox2_histo_desc = "ZBox2"
			self.ZBox2_histo_name = self.alabel(ZBox2_histo_desc)
			self.hman.h1(self.ZBox2_histo_name, ZBox2_histo_desc, 
			25, self.box1.zmin, self.box1.zmax)
			self.hman.fetch(
			self.ZBox2_histo_name).GetXaxis().SetTitle(
			"Z interaction box 2 (mm)")

			TBox2_histo_desc = "TBox2"
			self.TBox2_histo_name = self.alabel(TBox2_histo_desc)
			self.hman.h1(self.TBox2_histo_name, TBox2_histo_desc, 
			50, 300, 600)
			self.hman.fetch(
			self.TBox2_histo_name).GetXaxis().SetTitle(
			"Time of gamma in  box 2 (ps)")

		
	
	
