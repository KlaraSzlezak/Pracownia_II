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

c1 = TCanvas('cHisto','cHisto',600,600)
histo = gROOT.FindObject('histo')
histo.Draw()
c1.Print("histo.pdf")
input('press enter to exit')
