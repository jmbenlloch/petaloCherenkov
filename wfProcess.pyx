#import random as rnd
import numpy as np
cimport numpy as np
#from math import sqrt

cdef extern from "math.h":
    double sqrt(double x)

def wfProcess2(int minWL, int maxWL, double pde, int sptr, int asic, photons):
    result = []
    for photon in photons:
        if minWL < photon[1] < maxWL:
            if np.random.uniform(0.,1.) <= pde:
                time = photon[0]
                if sptr > 0:
                    time += np.random.normal(0, sptr)
                if asic > 0:
                    time += np.random.normal(0, asic)
                result.append((time, photon[1],photon[2])) 
    return sorted(result, key=lambda p:p[0])



cpdef double processPhoton(int minWL, int maxWL, double pde, double sptr, double asic, double time, double wl):
    cdef double result = -142857
    if minWL < wl < maxWL:
        if np.random.uniform(0.,1.) <= pde:
            result = time
            if sptr > 0:
                result += np.random.normal(0, sptr)
            if asic > 0:
                result += np.random.normal(0, asic)
    return result

def wfProcess3(int minWL, int maxWL, double pde, double sptr, double asic, photons):
    result = []
    for ph in photons:
        newtime = processPhoton(minWL, maxWL, pde/100., sptr/1000., asic/1000., ph[0], ph[1])
        if newtime != -142857:
            result.append((newtime, ph[1],ph[2])) 
    return sorted(result, key=lambda p:p[0])


cpdef double  distance(double x1, double y1, double z1, double x2, double y2, double z2):
    return sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


