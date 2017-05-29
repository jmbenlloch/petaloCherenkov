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


def getDTOF(times, npe):
    result = []
    for t in times:
        try:
            result.append(t[npe-1])
        except:
            pass
    return result


pdes = [30,50,70,100]
for pde in pdes:
    sigmas[pde] = []
    errors[pde] = []
    for npe in xrange(1,11):
        result = collection.find_one({'vel': 210, 'pde': pde, 'minwl': 150, 'maxwl': 1500, 'asic': 0, 'sptr': 0, 'npe': 10})
        
        gStyle.SetOptStat(1110)
        gStyle.SetOptFit(1)
        
        h1 = TH1F("h1", "PETALO", 80,-100,100)
        
        for n in getDTOF(result['avg'],npe):
            h1.Fill(n)
        
        gauF1 = TF1("gauF1","gaus",-200,200)
        h1.Fit("gauF1","","e",-10,10);
        h1.Draw()
        c1.Print('plots/avg/npe_' + str(npe) + '_pde_' + str(pde) + '.png')
    
        sigma = h1.GetFunction("gauF1").GetParameter(2)
        error = h1.GetFunction("gauF1").GetParError(2)
        sigmas[pde].append(sigma)
        errors[pde].append(error)
    

#Plot
c2 = TCanvas("c2","multipads",900,700);
colors = [ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack, ROOT.kGreen+3, ROOT.kBlue, ROOT.kGreen, ROOT. kCyan+3, ROOT.kRed, ROOT.kBlack]

plots = []

for index,pde in enumerate(pdes):
    npes = array('f', range(1,10))
    sigma = array('f', sigmas[pde])
    error = array('f', errors[pde])
    zeros = array('f', np.zeros(len(npes)))
    
    pdePlot = TGraphErrors(len(npes), npes, sigma,zeros, error)

    pdePlot.SetLineColor(colors[index%20]);
    pdePlot.SetLineWidth(2);
    pdePlot.SetLineStyle(index%10);

    plots.append(pdePlot)
    leg.AddEntry(plots[index],"PDE {}%".format(pde), "lp")

    if index == 0:
        plots[index].SetTitle("");
        plots[index].GetXaxis().SetTitle("#PE");
        plots[index].GetXaxis().SetLimits(0,10);
        plots[index].GetYaxis().SetTitle("CRT (ps)");
        plots[index].SetMinimum(0.);
        plots[index].SetMaximum(50.);
        plots[index].Draw()
    else:
        plots[index].Draw('same')

leg.Draw('same')

c2.Print('plots/avg/pde_avg.png')
