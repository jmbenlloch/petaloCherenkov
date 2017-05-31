#import petaloAnalysis
from ROOT import *
import numpy as np
import pandas as pd
from array import array

from pymongo import MongoClient

#map(processData, [(210,pde,0,1500,0,0) for pde in xrange(10,105,5)])
#connect to db
client = MongoClient('127.0.0.1', 27017)
db = client['petalo']
collection = db['cherenkov']

#old version -> compute here
#petalo = petaloAnalysis.PetaloAnalysis('sensors.csv','events.csv')
#petalo.processWaveform(0,1500,1,0,0)

c1 = TCanvas("c1","multipads",900,700);

pdes = []
sigmas = []
errors = []

for pde in xrange(10,105,5):
#    result = petalo.computeDTOF(vel)
    result = collection.find_one({'vel': 210, 'pde': pde, 'minwl': 150, 'maxwl': 1500, 'asic': 0, 'sptr': 0, "avg": {"$exists": False}})

    gStyle.SetOptStat(1110)
    gStyle.SetOptFit(1)
    h1 = TH1F("h1", "PETALO", 80,-100,100)
    for n in result['dtof']:
        h1.Fill(n)

    gauF1 = TF1("gauF1","gaus",-200,200)
    h1.Fit("gauF1","","e",-10,10);
    h1.Draw()
    c1.Print('plots/pde/' + str(pde) + '.png')

    sigma = h1.GetFunction("gauF1").GetParameter(2)
    error = h1.GetFunction("gauF1").GetParError(2)
    sigmas.append(sigma)
    errors.append(error)
    pdes.append(pde)


pdes = array('f', pdes)
sigmas = array('f', sigmas)
errors = array('f', errors)
zeros = array('f', np.zeros(len(pdes)))

    
c2 = TCanvas("c2","multipads",900,700);
pdePlot = TGraphErrors(len(pdes), pdes, sigmas,zeros, errors)
pdePlot.SetTitle("");
pdePlot.GetXaxis().SetTitle("PDE (%)");
pdePlot.GetXaxis().SetLimits(0,100);
pdePlot.GetYaxis().SetTitle("CRT (ps)");
pdePlot.SetMinimum(0.);
pdePlot.SetMaximum(10.);



pdePlot.Draw()
c2.Print('plots/pde/pde.png')
