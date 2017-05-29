import petaloAnalysis
from ROOT import *
import numpy as np
import pandas as pd
from array import array

from pymongo import MongoClient

#map(processData, [(210,pde,0,1500,0,sptr) for pde in xrange(10,105,5) for sptr in xrange(0,90,10)])
#connect to db
client = MongoClient('127.0.0.1', 27017)
db = client['petalo']
collection = db['cherenkov']

#old version -> compute here
#petalo = petaloAnalysis.PetaloAnalysis('sensors.csv','events.csv')
#petalo.processWaveform(0,1500,1,0,0)

c1 = TCanvas("c1","multipads",900,700);

sigmas = {}
errors = {}

for sptr in xrange(0,90,10):
    sigmas[sptr] = []
    errors[sptr] = []
    for pde in xrange(10,105,5):
#    result = petalo.computeDTOF(vel)
        result = collection.find_one({'vel': 210, 'pde': pde, 'minwl': 150, 'maxwl': 1500, 'asic': 0, 'sptr': sptr})
        
        gStyle.SetOptStat(1110)
        gStyle.SetOptFit(1)
        
        if sptr > 40:
            h1 = TH1F("h1", "PETALO", 80,-200,200)
        else:
            h1 = TH1F("h1", "PETALO", 80,-100,100)
        
        for n in result['dtof']:
            h1.Fill(n)
        
        gauF1 = TF1("gauF1","gaus",-200,200)
        if sptr > 40:
            h1.Fit("gauF1","","e",-70,70);
        elif sptr > 0:
            h1.Fit("gauF1","","e",-40,40);
        else:
            h1.Fit("gauF1","","e",-15,15);
        h1.Draw()
        c1.Print('plots/jitter/sptr_' + str(sptr) + '_pde_' + str(pde) + '.png')
    
        sigma = h1.GetFunction("gauF1").GetParameter(2)
        error = h1.GetFunction("gauF1").GetParError(2)
        sigmas[sptr].append(sigma)
        errors[sptr].append(error)
    

#Plot
c2 = TCanvas("c2","multipads",900,700);
colors = [ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack]

plots = []

leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)

for index,sptr in enumerate(xrange(0,90,10)):
    pdes = array('f', range(10,105,5))
    sigma = array('f', sigmas[sptr])
    error = array('f', errors[sptr])
    zeros = array('f', np.zeros(len(pdes)))
    
    pdePlot = TGraphErrors(len(pdes), pdes, sigma,zeros, error)

    pdePlot.SetLineColor(colors[index%20]);
    pdePlot.SetLineWidth(2);
    pdePlot.SetLineStyle(index%10);

    plots.append(pdePlot)
    leg.AddEntry(plots[index],"Jitter {} ps".format(sptr), "lp")

    if index == 0:
        plots[index].SetTitle("");
        plots[index].GetXaxis().SetTitle("PDE (%)");
        plots[index].GetXaxis().SetLimits(0,100);
        plots[index].GetYaxis().SetTitle("CRT (ps)");
        plots[index].SetMinimum(0.);
        plots[index].SetMaximum(100.);
        plots[index].Draw()
    else:
        plots[index].Draw('same')


leg.Draw('same')
c2.Print('plots/jitter/jitter.png')
