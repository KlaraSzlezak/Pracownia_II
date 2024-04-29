#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *

def DrawHisto ( histo, color):
    histo.SetStats(0)
    histo.SetLineColor( 1  )
    histo.SetLineWidth( 2  )
    histo.SetFillColor( color )

print ("Hello ROOT")
fileName = "jobs/histo.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

fG = TF1('fG' , '[4] * TMath::Gaus(x, [0], [1]) + [2]*x + [3] ' , 2.9 , 3.3 )
fG.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fG.SetParNames ( 'Mass' , 'Sigma' , 'Slope' , 'Offset', 'Normalization' )
fG.SetLineColor ( 4 )
fG.SetLineWidth ( 2 )

fBW = TF1('fG' , '[4] * TMath::BreitWigner(x, [0], [1]) + [2]*x + [3] ' , 2.9 , 3.3 )
fBW.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fBW.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'Normalization' )
fBW.SetLineColor ( 3 )
fBW.SetLineWidth ( 2 )

# Slimmed Muons - invariant mass
cInvM = TCanvas( 'cInvM' , 'cInvM' , 1000 , 600 )

cInvM.SetLeftMargin(0.14)
hInvMass = gROOT.FindObject('histo')

DrawHisto(hInvMass, 19)

cInvM.cd()
hInvMass.Draw('hist')
cInvM.Draw()

cInvM.Print("Plots/SMinv2704.png")

# J/Psi Fit
cZoom = TCanvas("cZoom", "cZoom", 800, 600)
hJPsi = hInvMass.Clone('hJPsi')
hJPsi.SetTitle( 'J/#Psi' )
hJPsi.GetXaxis().SetRangeUser( 2.85 , 3.35 )
hJPsi.Fit( fG , 'R' ,"", 2.8 , 3.25 )
hJPsi.Fit( fBW , 'R' ,"", 2.8 , 3.2 )

lZoom = TLegend( 0.6 , 0.7 , 0.9 , 0.9 )
lZoom.SetTextSize   (0.035)
lZoom.AddEntry( hJPsi    , 'Data' , 'f' )
lZoom.AddEntry( fBW    , 'Breit-Wigner'  , 'l' )
lZoom.AddEntry( fG , 'Gauss'         , 'l' )

cZoom.cd()
hJPsi.Draw('hist')
fG.DrawCopy( 'same' )
fBW.DrawCopy( 'same' )
lZoom.Draw()
cZoom.Draw()

# Particle IDs
cID = TCanvas('cID','cID',1000,600)
cID.SetLeftMargin(0.14)

cID.cd()

hID = gROOT.FindObject('hID')
DrawHisto(hID, 19)
hID.GetXaxis().SetRangeUser(-5.5, hID.GetXaxis().GetXmax())

hID.Draw('hist')

cID.Draw()

#### Slimmed Displaced Muons
cDisp = TCanvas('cDisp','cDisp',1000,600)
cDisp.SetLeftMargin(0.14)

cDisp.cd()

hDisp = gROOT.FindObject('hDisp')
DrawHisto(hDisp, 19)

hDisp.Draw('hist')
#fG.DrawCopy( 'same' )
#fBW.DrawCopy( 'same' )
cDisp.Draw()


#### Comparison of invariant mass histogram for different collections
cComp = TCanvas('cComp','cComp',1000,600)
cComp.cd()

hSM = hInvMass.Clone('hSM')
hSDM = hDisp.Clone('hSDM')


hSM.SetLineColor( 0  )
hSM.SetLineWidth( 0  )
hSM.SetFillColorAlpha( 38 , 0.4 )


hSDM.SetLineColor( 0  )
hSDM.SetLineWidth( 0  )
hSDM.SetFillColorAlpha( 90 , 0.571)

lComp = TLegend( 0.55 , 0.75 , 0.9 , 0.9 )
lComp.SetTextSize   (0.035)
lComp.AddEntry( hSM   , 'Slimmed Muons' , 'f' )
lComp.AddEntry( hSDM    , 'Slimmed Displaced Muons'  , 'f' )


hSM.DrawCopy()
hSDM.DrawCopy('same')
lComp.Draw()
cComp.Draw()

#### Vertex 
cVertex = TCanvas( 'cVertex' , 'cVertex' , 1000 , 600 )
cVertex.SetSupportGL(True)
cVertex.SetLeftMargin(0.14)

hVz = gROOT.FindObject('hVz')
DrawHisto(hVz, 19)

cVertex.cd()
hVz.Draw('hist')
cVertex.Draw()

input('press enter to exit')
