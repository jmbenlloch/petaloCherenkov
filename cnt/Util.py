"""
Utilities
"""
from __future__ import print_function
from math import *
import sys, os
from Centella.physical_constants import *

import matplotlib.pyplot as plt
import numpy as np

ps = picosecond
nm = nanometer 
mol = mole
micron = micrometer

def wait():
	raw_input("Press a key...")


def drange(start, stop, step):
	"""
	a range of doubles
	"""
	r = start
	while r < stop:
		yield r
		r += step

def lin(x,x0,y0,x1,y1):
	"""
	lineal extrapolation
	"""
	
	return y0 + (x-x0)*(y1-y0)/(x1-x0)

def inRange(x,xmin,xmax):
	if x >= xmin and x <= xmax:
		return True
	else:
		return False

