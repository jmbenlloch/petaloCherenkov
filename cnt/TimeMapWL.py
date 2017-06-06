from Centella.AAlgo import AAlgo

from Particles import *
from Util import *
from Geometry import *
from TOF import *


"""
This algorithm defines the fiducial volume for the CRT analysis. 
"""


class TimeMapWL(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		TimeMapWL  Algorithm
		"""
		
		self.name = 'TimeMapWL'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

  		self.QE = self.doubles["QE"]  #quantum efficiency
  		self.SPTR = self.doubles["SPTR"]*ps  #single photon time resolution
  		self.ASIC= self.doubles["ASIC"]*ps  #ASIC contribution

  		self.minWL= self.doubles["minWL"]  #Min wavelength
  		self.maxWL= self.doubles["maxWL"]  #Max wavelength

		self.csv = self.strings["CSV"]
                self.csvFile = open(self.csv,'w')
                dataEvt = 'evtID,vtx1X,vtx1Y,vtx1Z,vtx1T,vtx2X,vtx2Y,vtx2Z,vtx2T\n'
                self.csvFile.write(dataEvt)

    

	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')

		### PARAMETERS
		
    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0

                self.histos()
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

                self.sensors = self.logman["USER"].gparam["sensors"]

                self.timeMap = [[],[]]

                self.ComputeTimeMap(event)

                if len(self.timeMap[0]) == 0 or len(self.timeMap[1]) == 0:
                    return False

                #Sort them according to time
                map(lambda tm: tm.sort(key=lambda s: s.time()), self.timeMap)
                #print "Event: {0}, box1: {1}, box2: {2}".format(event.GetEventID(), len(self.timeMap[0]),len(self.timeMap[1]))

                map(lambda t: self.hman.fill(self.wf1, t.time()*1000), self.timeMap[0])
                map(lambda t: self.hman.fill(self.wf2, t.time()*1000), self.timeMap[1])
                #print map(lambda t: t.time(), self.timeMap[0])
                #print map(lambda t: t.time(), self.timeMap[1])

#                print map(lambda t: (t.sensorID,t.wl,t.time()), self.timeMap[1])

                #Print file
		self.vertexBox = self.logman["USER"].gparam["vertexBox"]
		vtx1 =self.vertexBox[0]
		vtx2 =self.vertexBox[1]
                dataEvt = '{},'.format(event.GetEventID())
                dataEvt += '{},{},{},{},'.format(vtx1.x,vtx1.y,vtx1.z,vtx1.t)
                dataEvt += '{},{},{},{},'.format(vtx2.x,vtx2.y,vtx2.z,vtx2.t)
                dataEvt += 'box1,'
                for ph in self.timeMap[0]:
                    dataEvt += '{},{},{},'.format(ph.time(),ph.wl,ph.sensorID)
                dataEvt += 'box2,'
                for ph in self.timeMap[1]:
                    dataEvt += '{},{},{},'.format(ph.time(),ph.wl,ph.sensorID)
                dataEvt += '\n'
#                print "csv: {}".format(dataEvt)
                self.csvFile.write(dataEvt)


		self.numOutputEvents += 1

		self.logman["USER"].gparam["timeMap"] = self.timeMap
		
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)
		
		return


############################################################
	def ComputeTimeMap(self,event):

            b1Times = []
            b2Times = []
            ng1 = 0
            ng2 = 0
            for s in event.names_dvstore():

                #Apply PDE
                if rnd.uniform(0.,1.) <= self.QE:
                    sensorID = int(s.split('_')[0])
                    ph = event.fetch_dvstore(s)
                    
                    realTime = ph[0]
                    timeASIC = self.ASICTime(ph[0])
                    timeJitter = SmearTime(timeASIC,self.SPTR,self.ASIC)
                    wl = ph[1]
                    
                    #print "sensorID: {0}\ttime: {1}".format(sensorID,timeJitter)

#                    print realTime, timeASIC, timeJitter
                    photon = Photon(sensorID,wl,realTime,timeASIC,timeJitter)
                    
                    histWL = self.name + '.' + self.histWL
                    self.hman.fill(histWL, wl)

                    if self.minWL <= wl <= self.maxWL:
                        histFiltered = self.name + '.' + self.histWLFiltered
                        self.hman.fill(histFiltered, wl)
                
                        box = self.sensors[sensorID].box - 1
                        self.timeMap[box].append(photon)

                        if box == 0:
                            ng1 += 1
                        if box == 1:
                            ng2 += 1

                        if box == 0:
                            b1Times.append((sensorID, timeJitter))
                        if box == 1:
                            b2Times.append((sensorID, timeJitter))

                        
            self.hman.fill(self.ngammas1, ng1)
            self.hman.fill(self.ngammas2, ng2)
            b1Times.sort(key=lambda a: a[1])
            b2Times.sort(key=lambda a: a[1])
            #print "b1Times"
            #print b1Times
            #print "b2Times"
            #print b2Times
            return


############################################################
        def histos(self):
            self.histWL = "Wavelength"
            self.defineHisto(self.histWL,"Wavelength (nm)",1,[200],[0],[1000])

            self.histWLFiltered = "WavelengthFiltered"
            self.defineHisto(self.histWLFiltered,"Wavelength (nm)",1,[200],[0],[1000])

            self.ngammas1 = self.defineHisto("NGBOX1", "ngammas",1,[200],[0],[200])
            self.ngammas2 = self.defineHisto("NGBOX2", "ngammas",1,[200],[0],[200])

            self.wf1 = self.defineHisto("wf1", "wf1",1,[200],[0],[2000])
            self.wf2 = self.defineHisto("wf2", "wf2",1,[200],[0],[2000])

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
			
############################################################
        #TODO Add resolution instead of only 5ps
        #TODO Check correctness
        def ASICTime(self, t):
            t0 = (t // 10) * 10
            t1 = t0 + 5
#            t2 = t0 + 10
            if t0 <= t < t1:
                return t0/1000
            else:
                return t1/1000
                

