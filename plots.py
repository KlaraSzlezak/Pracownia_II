#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
import numpy as np
from ROOT import *
from ROOT import TMath


def DrawHisto ( histo, color, line=1, width=2):
    histo.SetStats(0)
    histo.SetLineColor( line  )
    histo.SetLineWidth( width  )
    histo.SetFillColor( color )



def Background( a , b , fun , color ):
    
    fun.SetParameters( a , b)
    fun.SetLineColor( color )  
    fun.SetLineStyle(7)  
    fun.SetLineWidth(2)  
    

print ("Hello ROOT")
fileName = "jobs/part21h.root"
print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

xMin = 2.97
xMax = 3.23

# FUNCTIONS
fG = TF1('fG' , '[4] * TMath::Gaus(x, [0], [1]) + [2]*x + [3] ' , xMin , xMax )
fG.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fG.SetParNames ( 'Mass' , 'Sigma' , 'Slope' , 'Offset', 'Normalization' )
fG.SetLineColor ( 4 )
fG.SetLineWidth ( 2 )

fBW = TF1('fBW' , '[4] * TMath::BreitWigner(x, [0], [1]) + [2]*x + [3] ' , xMin , xMax )
fBW.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fBW.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'Normalization' )
fBW.SetLineColor ( 3 )
fBW.SetLineWidth ( 2 )

sigmaToGamma = 1/(2 * np.sqrt( 2 * np.log(2) ) )
fBWG = TF1('fBWG' , '[4] *( TMath::BreitWigner(x, [0], [1]) *  TMath::Gaus(x, [0], [1] * [5] ) )+ [2]*x + [3] ' , xMin , xMax )
fBWG.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fBWG.FixParameter(5, sigmaToGamma)
fBWG.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'Normalization' )
fBWG.SetLineColor ( 798 )
fBWG.SetLineWidth ( 2 )


####### Slimmed Muons - invariant mass
cInvM = TCanvas( 'cInvM' , 'cInvM' , 1000 , 600 )

cInvM.SetLeftMargin(0.14)
hInvMass = gROOT.FindObject('hInv')

DrawHisto(hInvMass, 19)

cInvM.cd()
hInvMass.Draw('hist')
cInvM.Draw()

#cInvM.Print("Plots/SMinv.png")

####### J/Psi Fit
cZoom = TCanvas("cZoom", "cZoom", 1200, 800)
hJPsi = hInvMass.Clone('hJPsi')
hJPsi.SetTitle( 'J/#Psi' )
hJPsi.GetXaxis().SetRangeUser( 2.85 , 3.35 )
hJPsi.Fit( fG , 'R' ,"", 2.8 , 3.25 )
hJPsi.Fit( fBW , 'R' ,"", 2.8 , 3.2 )
hJPsi.Fit( fBWG , 'R' ,"", 2.8 , 3.2 )

lZoom = TLegend( 0.6 , 0.7 , 0.9 , 0.9 )
lZoom.SetTextSize   (0.035)
lZoom.AddEntry( hJPsi    , 'Data' , 'f' )
lZoom.AddEntry( fBW    , 'Breit-Wigner'  , 'l' )
lZoom.AddEntry( fG     , 'Gauss'         , 'l' )
lZoom.AddEntry( fBWG   , 'Gauss+Breit-Wigner' , 'l' )

fBWGbckg = TF1('fBWGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
fGbckg = TF1('fGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
fBWbckg = TF1('fBWbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)

Background(fBWG.GetParameter(2), fBWG.GetParameter(3), fBWGbckg , fBWG.GetLineColor() )
Background(fBW.GetParameter(2), fBW.GetParameter(3), fBWbckg , fBW.GetLineColor() )
Background(fG.GetParameter(2), fG.GetParameter(3), fGbckg , fG.GetLineColor() )

cZoom.cd()
hJPsi.Draw('hist')
fG.DrawCopy( 'same' )
fBW.DrawCopy( 'same' )
fBWGbckg.DrawCopy( 'same' )
fBWbckg.DrawCopy( 'same' )
fGbckg.DrawCopy( 'same' )
fBWG.DrawCopy( 'same' )
lZoom.Draw()
cZoom.Draw()


#######  Particle IDs from 25 events
cID = TCanvas('cID','cID',1000,600)
cID.SetLeftMargin(0.14)

cID.cd()

hID = gROOT.FindObject('hID')
DrawHisto(hID, 19)
hID.GetXaxis().SetRangeUser(-5.5, hID.GetXaxis().GetXmax())

hID.Draw('hist')

cID.Draw()


#######  Slimmed Displaced Muonsh
hDisp = gROOT.FindObject('hDisp')
'''
cDisp = TCanvas('cDisp','cDisp',1000,600)
cDisp.SetLeftMargin(0.14)

cDisp.cd()


DrawHisto(hDisp, 19)

hDisp.Draw('hist')
#fG.DrawCopy( 'same' )
#fBW.DrawCopy( 'same' )
cDisp.Draw()
'''

#######  Comparison of invariant mass histogram for different collections
cComp = TCanvas('cComp','cComp',1200,600)
cComp.cd()

hCan = gROOT.FindObject('hCan')
hInvG = gROOT.FindObject('hInvG')

hSM = hInvMass.Clone('hSM')
hSDM = hDisp.Clone('hSDM')

DrawHisto(hCan, 0, 3, 1)
DrawHisto(hSM, 0, 4, 1)
#hSM.SetFillColorAlpha( 0 , 0.4 )
DrawHisto(hSDM, 0, 798 , 2)
#hSDM.SetFillColorAlpha( 0, 0.571)
DrawHisto(hInvG, 0, 880 , 2)

lComp = TLegend( 0.53 , 0.7 , 0.9 , 0.9 )
lComp.SetTextSize   (0.035)
lComp.AddEntry( hSM   , 'Slimmed Muons' , 'l' )
lComp.AddEntry( hSDM    , 'Slimmed Displaced Muons'  , 'l' )
lComp.AddEntry( hCan    , 'Packed Candidates (Muons)'  , 'l' )
lComp.AddEntry( hInvG    , 'Global muons (Slimmed Muons)'  , 'l' )


hSM.DrawCopy()
hSDM.DrawCopy('same')
hCan.DrawCopy('same')
hInvG.DrawCopy('same')
lComp.Draw()
cComp.Draw()

####### Vertex (distance for muons)
cVertex = TCanvas( 'cVertex' , 'cVertex' , 1000 , 600 )
cVertex.SetSupportGL(True)
cVertex.SetLeftMargin(0.14)

hVz = gROOT.FindObject('hVz')
DrawHisto(hVz, 19)

cVertex.cd()
hVz.Draw('hist')
cVertex.Draw()



cProb = TCanvas( 'cProb' , 'cProb' , 1000 , 600 )

cProb.SetLeftMargin(0.14)
hProb = gROOT.FindObject('hProb')

DrawHisto(hProb, 19)

cProb.cd()
hProb.Draw('hist')
cProb.Draw()

cJpsiProb = TCanvas( 'cJpsiProb' , 'cJpsiProb' , 1000 , 600 )

cJpsiProb.SetLeftMargin(0.14)
hJpsiProb = gROOT.FindObject('hJpsiProb')

DrawHisto(hJpsiProb, 19)

cJpsiProb.cd()
hJpsiProb.Draw('hist')
cJpsiProb.Draw()

cVzPC = TCanvas( 'cVzPC' , 'cVzPC' , 1000 , 600 )

cVzPC.SetLeftMargin(0.14)
hVzPC = gROOT.FindObject('hVzPC')
hVzPC.GetEntries()

DrawHisto(hVzPC, 19)

cVzPC.cd()
hVzPC.Draw('hist')
cVzPC.Draw()

####### Vertices distance for muons from J/Psi decay for different probability cuts
cVz0x = TCanvas( 'cVz0x' , 'cVz0x' , 1000 , 600 )

cVzPC.SetLeftMargin(0.14)
hVz01 = gROOT.FindObject('hVz01')
hVz08 = gROOT.FindObject('hVz08')

print("AVERAGE DISCTANCE BETWEEN VERTICES FOR P>= 0.8, jpsi: ", hVz08.GetMean() )
print("AVERAGE DISCTANCE BETWEEN VERTICES FOR P>= 0.1, jpsi: ", hVz01.GetMean() )

DrawHisto(hVz01, 0, 4)

DrawHisto(hVz08, 0, 3)

lVz0x = TLegend( 0.55 , 0.8 , 0.9 , 0.9 )
lVz0x.SetTextSize   (0.035)
lVz0x.AddEntry( hVz01   , 'P > 0.1' , 'l' )
lVz0x.AddEntry( hVz08    , 'SP > 0.8'  , 'l' )

hVz01.GetXaxis().SetRangeUser( -0.02 ,0.1)
cVz0x.cd()
hVz01.Draw('hist')
hVz08.Draw('same')
lVz0x.Draw()
cVz0x.Draw()


####### Invariant mass of K+ K-
cInvK = TCanvas( 'cInvK' , 'cInvK' , 1000 , 600 )

cInvK.SetLeftMargin(0.14)
hInvK = gROOT.FindObject('hInvK')

DrawHisto(hInvK, 19)

cInvK.cd()
hInvK.Draw('hist')
cInvK.Draw()

####### Probability of a common vertex for a pair of kaons
ckkProb = TCanvas( 'ckkProb' , 'ckkProb' , 1000 , 600 )

ckkProb.SetLeftMargin(0.14)
hkkProb = gROOT.FindObject('hkkProb')

DrawHisto(hkkProb, 19)

ckkProb.cd()
hkkProb.Draw('hist')
ckkProb.Draw()

print("ENTRIES, hInvK: ", hInvK.GetEntries() )

print("ENTRIES, hkkProb: ", hkkProb.GetEntries() )

'''####### Invariant mass of K+ K- mu+ mu-
cMuMuKK = TCanvas( 'cMuMuKK' , 'cMuMuKK' , 1000 , 600 )

cMuMuKK.SetLeftMargin(0.14)
hMuMuKK = gROOT.FindObject('hMuMuKK')

DrawHisto(hMuMuKK, 19, 1, 5)

cMuMuKK.cd()
hMuMuKK.Draw('hist')
cMuMuKK.Draw()

print("ENTRIES, hMuMuKK: ", hMuMuKK.GetEntries() )
'''
input('press enter to exit')
