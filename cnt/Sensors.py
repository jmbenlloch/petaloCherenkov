from Centella.AAlgo import AAlgo
from Util import *
from Geometry import *
from TOF import *
from Scintillator import *
import random as rnd


"""
This algorithm get all sensors and their positions
"""


class Sensors(AAlgo):

	############################################################
	def __init__(self, param=False, level=1, label="", **kargs):

		"""
		Sensors Algorithm
		"""
		#self.m.log(1, 'Constructor()')

		### GENERAL STUFF
		self.name = 'Sensors'
		
		AAlgo.__init__(self, param, level, self.name, 0, label, kargs)

    ### PARAMETERS
		self.debug = self.ints["Debug"]  #used to stop program at key break points

  		self.box1ID=self.vints["Box1Id"]  #ids of SiPMs in box1
  		self.box2ID=self.vints["Box2Id"] 

                #DataFrame
		self.csv = self.strings["CSVSensors"]
                self.csvFile = open(self.csv,'w')
                dataEvt = 'sensorID,x,y,z,box\n'
                self.csvFile.write(dataEvt)

		if self.debug == 1:
			wait()


	############################################################		
	def initialize(self):

		self.m.log(1, 'Initialize()')
		
		### Defining histos

    ### Counters:
		self.numInputEvents = 0
		self.numOutputEvents = 0
		
		return


	############################################################
	def execute(self, event=""):

		self.m.log(2, 'Execute()')		
	
		self.numInputEvents += 1   

                self.sensors = {}

		self.GetSensors(event)

		if self.debug == 1:
			wait()	
		
		self.numOutputEvents += 1
		self.logman["USER"].gparam["sensors"]=self.sensors
		return True

############################################################
	def finalize(self):

		self.m.log(1, 'Finalize()')

		self.m.log(1, 'Input  Events: ', self.numInputEvents)
		self.m.log(1, 'Output Events: ', self.numOutputEvents)

		return



############################################################
	def GetSensors(self,event):
		"""
		Get sensors
		"""

		sensorhits =  event.GetMCSensHits()
		self.m.log(3, " Get sensors: event has %d sensor hits = "%(
			len(sensorhits)))

		for hit in sensorhits:
			hid = hit.GetSensorID()
			x = hit.GetPosition().x()
			y = hit.GetPosition().y()
			z = hit.GetPosition().z()
                        box = self.boxId(hid)
                        
                        sensor = SiPM(hid,x,y,z,box)
                        self.sensors[hid] = sensor

                        dataEvt = "{},{},{},{},{}\n".format(hid,x,y,z,box)
                        #print "csv: {}".format(dataEvt)
                        self.csvFile.write(dataEvt)


	############################################################
	def boxId(self, hid):
		b1f1min = self.box1ID[0]
		b1f1max = self.box1ID[1]
		b1f2min = self.box1ID[2]
		b1f2max = self.box1ID[3]

		b2f1min = self.box2ID[0]
		b2f1max = self.box2ID[1]
		b2f2min = self.box2ID[2]
		b2f2max = self.box2ID[3]

		if inRange(hid,b1f1min,b1f1max) or inRange(hid,b1f2min,b1f2max):
			return 1

		if inRange(hid,b2f1min,b2f1max) or inRange(hid,b2f2min,b2f2max):
			return 2

		return 0

