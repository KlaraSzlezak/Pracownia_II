#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "histos.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

c1 = TCanvas('cHisto','cHisto',800,600)

c1.SetLeftMargin(0.14)

histo = gROOT.FindObject('histo')
histo.SetStats(0)
histo.SetLineColor( 1  )
histo.SetLineWidth( 2  )
histo.SetFillColor( 19 )
histo.Draw()
c1.Print("histo.jpg")
input('press enter to exit')
