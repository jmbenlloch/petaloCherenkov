import pandas as pd
import numpy as np
import random as rnd
from math import sqrt
from wfProcess  import distance

def avg(sensors, sensorsIndex, row, velocity, n):
    np1 = len(row['wf1'])
    np2 = len(row['wf2'])
    result = np.nan
    if np1 > 0 and np2 > 0:
        lastNPE = min(np1,np2)
        lastNPE = min(n,lastNPE)

        t1 = 0
        t2 = 0
        result = []
        for i in xrange(lastNPE):
            ph1 = row['wf1'][i]
            ph2 = row['wf2'][i]
            sensor1 = sensors.iloc[[ sensorsIndex[ph1[2]] ]]
            sensor2 = sensors.iloc[[ sensorsIndex[ph2[2]] ]]

            d1 = distance(row['vtx1X'],row['vtx1Y'],row['vtx1Z'],sensor1['x'],sensor1['y'],sensor1['z'])
            d2 = distance(row['vtx2X'],row['vtx2Y'],row['vtx2Z'],sensor2['x'],sensor2['y'],sensor2['z'])

            t1 += row['wf1'][i][0] - d1/(velocity*1.0)
            t2 += row['wf2'][i][0] - d2/(velocity*1.0)

            t1Avg = t1/(i+1) - row['vtx1T']
            t2Avg = t2/(i+1) - row['vtx2T']
            result.append(500 * (t2Avg - t1Avg)) #Expect ns
    return result

def dtof(sensors, sensorsIndex, row, velocity):

    if len(row['wf1']) > 0 and len(row['wf2']) > 0:
        dtg = row['vtx2T'] - row['vtx1T']
    
        sensor1 = sensors.iloc[[ sensorsIndex[row['wf1'][0][2]] ]]
        sensor2 = sensors.iloc[[ sensorsIndex[row['wf2'][0][2]] ]]

        d1 = distance(row['vtx1X'],row['vtx1Y'],row['vtx1Z'],sensor1['x'],sensor1['y'],sensor1['z'])
        d2 = distance(row['vtx2X'],row['vtx2Y'],row['vtx2Z'],sensor2['x'],sensor2['y'],sensor2['z'])
        dpg = (d2-d1) / (velocity*1.0)

        t1 = row['wf1'][0][0]
        t2 = row['wf2'][0][0]
        dt = t2-t1
        return 0.5 * (dt - dtg - dpg) * 1000

    else:
        return np.nan
    
def distances(sensors,indexes,box):
    def distance1(row):
        try:
            sensor = sensors.iloc[[ indexes[row['wf1'][0][2]] ]]
            return distance(row['vtx1X'],row['vtx1Y'],row['vtx1Z'],sensor['x'],sensor['y'],sensor['z'])
        except:
            return np.nan
    
    def distance2(row):
        try:
            sensor = sensors.iloc[[ indexes[row['wf2'][0][2]] ]]
            return distance(row['vtx2X'],row['vtx2Y'],row['vtx2Z'],sensor['x'],sensor['y'],sensor['z'])
        except:
            return np.nan

    if box == 1:
        return distance1
    elif box == 2:
        return distance2
    else:
        raise ValueError

def firstTime(photons):
    try:
        return photons[0][0]
    except IndexError:
        return np.nan
    
def firstWL(photons):
    try:
        return photons[0][1]
    except IndexError:
        return np.nan
