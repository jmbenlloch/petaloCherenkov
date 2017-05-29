import pandas as pd
import numpy as np
import random as rnd
from math import sqrt

from sensor import *
from wfProcess import *

class PetaloAnalysis(object):
    def __init__(self,sensorsFile,eventsFile):
        self.sensors = pd.read_csv(sensorsFile)
        eventsF = open(eventsFile,'r')

        self.sensorsIndex = {}
        for index,s in self.sensors.iterrows():
            self.sensorsIndex[int(s['sensorID'])] = index
            #access s.iloc[[255]]

        self._events = pd.DataFrame()

        eventIDs = []
        vtx1X = []
        vtx1Y = []
        vtx1Z = []
        vtx1T = []
        vtx2X = []
        vtx2Y = []
        vtx2Z = []
        vtx2T = []
        wf1 = []
        wf2 = []

        eventsF.readline()
        for line in eventsF:
            data = line.split(',')
            eventIDs.append(int(data[0]))
            vtx1X.append(float(data[1]))
            vtx1Y.append(float(data[2]))
            vtx1Z.append(float(data[3]))
            vtx1T.append(float(data[4]))
            vtx2X.append(float(data[5]))
            vtx2Y.append(float(data[6]))
            vtx2Z.append(float(data[7]))
            vtx2T.append(float(data[8]))
            if data[9] != 'box1':
                print 'Error'
            count = 10
            times1 = []
            while data[count] != 'box2':
                times1.append((float(data[count]),float(data[count+1]),int(data[count+2])))
                count += 3
            wf1.append(times1)
            times1 = []
            
            count += 1
            times2 = []
            while count+2 < len(data)-1:
                times2.append((float(data[count]),float(data[count+1]),int(data[count+2])))
                count += 3
            wf2.append(times2)
            times2 = []
            
        self._events['eventIDs'] = eventIDs
        self._events['vtx1X'] = vtx1X
        self._events['vtx1Y'] = vtx1Y
        self._events['vtx1Z'] = vtx1Z
        self._events['vtx1T'] = vtx1T
        self._events['vtx2X'] = vtx2X
        self._events['vtx2Y'] = vtx2Y
        self._events['vtx2Z'] = vtx2Z
        self._events['vtx2T'] = vtx2T
        self._events['wf1'] = wf1
        self._events['wf2'] = wf2

    def processWaveform(self, minWL, maxWL, pde, sptr, asic):
        self.minWL = minWL
        self.maxWL = maxWL
        self.pde = pde
        self.sptr = sptr
        self.asic = asic
        self.eventsProcessed = self._events.copy()

        self.eventsProcessed['wf1'] = self.eventsProcessed['wf1'].apply(lambda photons: wfProcess3(minWL,maxWL,pde,sptr,asic,photons))
        self.eventsProcessed['wf2'] = self.eventsProcessed['wf2'].apply(lambda photons: wfProcess3(minWL,maxWL,pde,sptr,asic,photons))

    def computeDTOF(self, velocity):
        self.eventsProcessed['dtg'] = self.eventsProcessed['vtx2T'] - self.eventsProcessed['vtx1T']

        self.eventsProcessed['d1'] = self.eventsProcessed.apply(distances(self.sensors,self.sensorsIndex,1), axis=1)
        self.eventsProcessed['d2'] = self.eventsProcessed.apply(distances(self.sensors,self.sensorsIndex,2), axis=1)
        self.eventsProcessed['dpg'] = (self.eventsProcessed['d2'] - self.eventsProcessed['d1']) / (velocity*1.0)

        self.eventsProcessed['t1'] = self.eventsProcessed['wf1'].apply(firstTime)
        self.eventsProcessed['t2'] = self.eventsProcessed['wf2'].apply(firstTime)
        self.eventsProcessed['wl1'] = self.eventsProcessed['wf1'].apply(firstWL)
        self.eventsProcessed['wl2'] = self.eventsProcessed['wf2'].apply(firstWL)

        self.eventsProcessed['dt'] = self.eventsProcessed['t2'] - self.eventsProcessed['t1']

        self.eventsProcessed['dtof'] = 0.5 * (self.eventsProcessed['dt'] - self.eventsProcessed['dtg'] - self.eventsProcessed['dpg']) * 1000 #in ps

        self.result = self.eventsProcessed[self.eventsProcessed['dtof'].notnull()][['wl1','wl2','dtof']]
        return self.result


    def computeAVG(self, velocity, n):
        self.eventsProcessed['avg'] = self.eventsProcessed.apply(lambda row: avg(self.sensors, self.sensorsIndex, row, velocity, n), axis=1)

        self.result = self.eventsProcessed[self.eventsProcessed['avg'].notnull()]['avg']
        return self.result


    def buildHist(self):
        hhist, hedges = np.histogram(self.eventsProcessed['dtof'], bins=80,range=(-100,100))
        result = dict()
        result['bottom'] = 0
        result['left'] = hedges[:-1].tolist()
        result['right'] = hedges[1:].tolist()
        result['top'] = hhist.tolist()

        return result


    def buildScatter(self):
        dtof = self.eventsProcessed['dtof'].tolist()
        dtof.extend(dtof)
        wl1 = self.eventsProcessed['wl1'].tolist()
        wl2 = self.eventsProcessed['wl2'].tolist()
        wl1.extend(wl2)
        result = dict()
        result['x'] = wl1
        result['y'] = dtof
        
        return result

    def buildBox(self):
        dtof = self.eventsProcessed['dtof'].tolist()
        dtof.extend(dtof)
        wl1 = self.eventsProcessed['wl1'].tolist()
        wl2 = self.eventsProcessed['wl2'].tolist()
        wl1.extend(wl2)
        wl1 = [str(int(w//100)*100) for w in wl1]
    
        catData = pd.DataFrame(dict(wl=wl1,dtof=dtof))
    #    cats = map(str,range(100,1300,100))
        start = max(100,minwl)
        end = min(1300,maxwl)
        cats = map(str,range(100*int(start/100),end,100))
    
        # find the quartiles and IQR for each category
        groups = catData.groupby('wl')
        q1 = groups.quantile(q=0.25)
        q2 = groups.quantile(q=0.5)
        q3 = groups.quantile(q=0.75)
        iqr = q3 - q1
        upper = q3 + 1.5*iqr
        lower = q1 - 1.5*iqr
        
        # find the outliers for each category
        def outliers(group):
            cat = group.name
            return group[(group.dtof > upper.loc[cat][0]) | (group.dtof < lower.loc[cat][0])]['dtof']
        out = groups.apply(outliers).dropna()
        
        # prepare outlier data for plotting, we need coordinates for every outlier.
        outx = []
        outy = []
        for cat in cats:
            # only add outliers if they exist
            if not out.loc[cat].empty:
                for value in out[cat]:
                    outx.append(cat)
                    outy.append(value)
        
        qmin = groups.quantile(q=0.00)
        qmax = groups.quantile(q=1.00)
        upper.dtof = [min([x,y]) for (x,y) in zip(list(qmax.iloc[:,0]),upper.dtof) ]
        lower.dtof = [max([x,y]) for (x,y) in zip(list(qmin.iloc[:,0]),lower.dtof) ]

        result = dict()
        result['q1'] = q1.dtof.tolist()
        result['q2'] = q2.dtof.tolist()
        result['q3'] = q3.dtof.tolist()
        result['upper'] = upper.dtof.tolist()
        result['lower'] = lower.dtof.tolist()
        result['outx'] = outx
        result['outy'] = outy
        result['cats'] = cats

        return result

    def buildColor(self):
        #-100 is red, 100 is green, 0 is black
        def color(dtof):
            color = [0,0,0]
            if dtof < 0:
                color[0] = int(-2.55 * dtof)
            else:
                color[1] = int(2.55 * dtof)
            return "#%02x%02x%02x" % (color[0], color[1], color[2])

        x = self.eventsProcessed['wl1'].as_matrix()
        y = self.eventsProcessed['wl2'].as_matrix()

        radii = map(lambda dtof: abs(dtof)*0.12+1.5,self.eventsProcessed['dtof'].as_matrix())
        colors = [color(dtof) for dtof in self.eventsProcessed['dtof'].as_matrix()]

        result = dict()
        result['radii'] = radii
        result['colors'] = colors

        return result
