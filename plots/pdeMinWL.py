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

leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)

for minwl in xrange(150,450,50):
    sigmas[minwl] = []
    errors[minwl] = []
    for pde in xrange(10,105,5):
#    result = petalo.computeDTOF(vel)
        result = collection.find_one({'vel': 210, 'pde': pde, 'minwl': minwl, 'maxwl': 1500, 'asic': 0, 'sptr': 0})
        
        gStyle.SetOptStat(1110)
        gStyle.SetOptFit(1)
        
        h1 = TH1F("h1", "PETALO", 80,-100,100)
        
        for n in result['dtof']:
            h1.Fill(n)
        
        gauF1 = TF1("gauF1","gaus",-200,200)
        h1.Fit("gauF1","","e",-10,10);
        h1.Draw()
        c1.Print('plots/minwl/minwl_' + str(minwl) + '_pde_' + str(pde) + '.png')
    
        sigma = h1.GetFunction("gauF1").GetParameter(2)
        error = h1.GetFunction("gauF1").GetParError(2)
        sigmas[minwl].append(sigma)
        errors[minwl].append(error)
    

#Plot
c2 = TCanvas("c2","multipads",900,700);
colors = [ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack]

plots = []

for index,minwl in enumerate(xrange(150,450,50)):
    pdes = array('f', range(10,105,5))
    sigma = array('f', sigmas[minwl])
    error = array('f', errors[minwl])
    zeros = array('f', np.zeros(len(pdes)))
    
    pdePlot = TGraphErrors(len(pdes), pdes, sigma,zeros, error)

    pdePlot.SetLineColor(colors[index%20]);
    pdePlot.SetLineWidth(2);
    pdePlot.SetLineStyle(index%10);

    plots.append(pdePlot)
    leg.AddEntry(plots[index],"minwl {} nm".format(minwl), "lp")

    if index == 0:
        plots[index].SetTitle("");
        plots[index].GetXaxis().SetTitle("PDE (%)");
        plots[index].GetXaxis().SetLimits(0,100);
        plots[index].GetYaxis().SetTitle("CRT (ps)");
        plots[index].SetMinimum(0.);
        plots[index].SetMaximum(10.);
        plots[index].Draw()
    else:
        plots[index].Draw('same')

leg.Draw('same')

c2.Print('plots/minwl/pde_minwl.png')
