from Centella.AAlgo import AAlgo

from Util import *
from Geometry import *
from TOF import *
from Scintillator import *

"""
This algorithm computes the Coincidence Resolution Time. 
"""


class DTOFWL(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		DTOFWL Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'DTOFWL'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS

		self.VEL = self.doubles["VEL"]  #used to stop program at key break points
		
		self.debug = self.ints["Debug"]  #used to stop program at key break points
		self.SCINT = self.strings["SCINTILLATOR"]
		self.INTER = self.strings["INTER"] # cher or scint
  		self.NINDEX = self.strings["NINDEX"] # FIX or VAR refraction index

		self.csv = self.strings["CSV"]
                self.csvFile = open(self.csv,'w')
                dataEvt = 'dt,dtg,dbox1,dbox2,wl1,wl2\n'
                self.csvFile.write(dataEvt)


		self.profileROOT = self.strings["VELHIST"] # ROOT file with vel hists
                if self.profileROOT != 'None':
		    rootFile = ROOT.TFile.Open(self.profileROOT, "read")
		    rootFile.PhVelTime.Rebin2D(40,40)
		    self.profileVel = rootFile.PhVelTime.ProfileX()
		    # Needed to avoid   'PyROOT_NoneType' object error
		    self.profileVel.SetDirectory(0)
                else:
                    self.profileVel = False
		#self.vel = photonVelocity(self.profileVel, self.SCINT, self.NINDEX, self.INTER)
		self.vel = self.VEL

		if self.debug == 1:
			wait()


	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos
		# Event energy histogram

		self.Histos()

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')	

		#get time map from previous algo
		self.timeMap = self.logman["USER"].gparam["timeMap"]
		self.vertexBox = self.logman["USER"].gparam["vertexBox"]
                self.sensors = self.logman["USER"].gparam["sensors"]
		
		#Compute DT
                print "Event {0}".format(event.GetEventID())
		self.ComputeDTOFWL() 
		
		if self.debug == 1:
			wait()	
		
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)

		self.logman["USER"].ints[self.alabel("InputEvents")] = self.numInputEvents
		self.logman["USER"].ints[self.alabel("OutputEvents")] = self.numOutputEvents

		return


###########################################################
	def ComputeDTOFWL(self):
		"""
		Compute DT for the first PE using directly gamma information
		DTOF = (1/2)*(dt - dtg - dpg)
		where: DT is the difference in TOF
		dt = t2 -t1 difference of time stamp of first PE 
		dtg = tg2 -tg1 difference of time stamp of gamma 1 and tGamma2 
		dpg = difference of path1 and path2, where 
		path = (n/c)*d and distance is the distance between gamma and SiPM vertex. 
		"""

		self.m.log(2," ---ComputeDTOFWL---")

		vertexBox1 =self.vertexBox[0]
		vertexBox2 =self.vertexBox[1]

		self.hman.fill(self.vertex1_histo,vertexBox1.t*1000)
		self.hman.fill(self.vertex2_histo,vertexBox2.t*1000)

                hitsBox1 = self.timeMap[0]
                hitsBox2 = self.timeMap[1]

		dt  = hitsBox2[0].time() - hitsBox1[0].time()
#                print hitsBox1[0].time()
#                print hitsBox2[0].time()
#                print dt

		dtg = vertexBox2.t - vertexBox1.t
#                print vertexBox1.t
#                print vertexBox2.t
#                print dtg

                sipm1 = self.sensors[hitsBox1[0].sensorID]
                sipm2 = self.sensors[hitsBox2[0].sensorID]
                print "sipm1ID: {0}".format(sipm1.id)
                print "sipm2ID: {0}".format(sipm2.id)

                dbox1 = distance(sipm1.XYZ(), vertexBox1.XYZ())
                dbox2 = distance(sipm2.XYZ(), vertexBox2.XYZ())
#                print dbox1
#                print dbox2
                
                tpath1 = dbox1 / self.vel
                tpath2 = dbox2 / self.vel
#                print tpath1
#                print tpath2

		dpg = tpath2 - tpath1
    
#                print "dbox1: {}".format(dbox1)
#                print "dbox2: {}".format(dbox2)

		dtof = 0.5*(dt - dtg - dpg)
#                print dtof

		self.hman.fill(self.dtof_histo_name,dtof/ps)
		self.hman.fill(self.dtof2_histo_name,dtof/ps)
		self.hman.fill(self.dtof3_histo_name,dtof/ps)

		self.hman.fill(self.hist_dt,dt*1000)
		self.hman.fill(self.hist_dtg,dtg*1000)
		self.hman.fill(self.hist_dpg,dpg*1000)
                print "dtof: %s" % dtof
                print "dt: %s" % dt
                print "dtg: %s" % dtg
                print "dpg: %s" % dpg
                print "dbox1: %s" % dbox1
                print "dbox2: %s" % dbox2
                print "v1: {0}".format(vertexBox1.XYZ())
                print "v2: {0}".format(vertexBox2.XYZ())
                print "sipm1: {0}".format(sipm1.XYZ())
                print "sipm2: {0}".format(sipm2.XYZ())

                dataEvt = "{},{},{},{},{},{}\n".format(dt,dtg,dbox1,dbox2,hitsBox1[0].wl,hitsBox2[0].wl)
                print "csv: {}".format(dataEvt)
                self.csvFile.write(dataEvt)


	############################################################		
	def Histos(self):
		"""
		book the histograms for the algo 
		"""
		self.dtof_histo_name = self.defineHisto("DTOF",
										 " DTOF (ps)",
										 1,[50],[-500],[500])

		self.dtof2_histo_name = self.defineHisto("DTOF2",
										 " DTOF2 (ps)",
										 1,[80],[-200],[200])

		self.dtof3_histo_name = self.defineHisto("DTOF3",
										 " DTOF3 (ps)",
										 1,[80],[-100],[100])

		self.vertex1_histo = self.defineHisto("Vertex1", "Vertex1 (ps)",1,[80],[0],[1000])
		self.vertex2_histo = self.defineHisto("Vertex2", "Vertex1 (ps)",1,[80],[0],[1000])

		self.hist_dt = self.defineHisto("dt", "dt (ps)",1,[80],[0],[1000])
		self.hist_dpg = self.defineHisto("dpg", "dpg (ps)",1,[80],[0],[1000])
		self.hist_dtg = self.defineHisto("dtg", "dtg (ps)",1,[80],[0],[1000])
		

	############################################################
	def defineHisto(self,histoName,histoTitle,histoType,nbin,xmin,xmax):

		histo_desc = histoName
		histo_name = self.alabel(histo_desc)
		if histoType == 1:
			self.hman.h1(histo_name, histo_desc,nbin[0],xmin[0],xmax[0])
		elif histoType == 2:
			self.hman.h2(histo_name, histo_desc,
				nbin[0],xmin[0],xmax[0],
				nbin[1],xmin[1],xmax[1])
		elif histoType == 3:
			self.hman.h3(histo_name, histo_desc,
				nbin[0],xmin[0],xmax[0],
				nbin[1],xmin[1],xmax[1],
				nbin[2],xmin[2],xmax[2])
		else:
			print "not implemented"
			sys.exit()
			
		self.hman.fetch(histo_name).GetXaxis().SetTitle(histoTitle)
		return histo_name
			
	
	
