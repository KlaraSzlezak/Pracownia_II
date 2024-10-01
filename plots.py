#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
import numpy as np
from ROOT import *
from ROOT import gStyle
from ROOT import TPaveText


pdgJPsi = 3.096900 
pdgJPsiErr = 0.000006

pdgPsi = 3.68610
pdgPsiErr = 0.00006

pdgUps1 = 9.46040 
pdgUps1Err = 0.0001

pdgUps2 = 10.0234 
pdgUps2Err = 0.0005

pdgUps3 = 10.3551
pdgUps3Err = 0.0005

pdgBp = 5.27934 
pdgBpErr = 0.00012

pdgX = 4.474 
pdgXErr = 0.004

def DrawHisto ( histo, color, line=1, width=2):
    histo.SetStats(0)
    histo.SetLineColor( line  )
    histo.SetLineWidth( width  )
    histo.SetFillColor( color )

def Background( a , b , fun , color ):
    
    fun.SetParameters( a , b)
    fun.SetLineColor( color )  
    fun.SetLineStyle(9)  
    fun.SetLineWidth(2)  
    

print ("Hello ROOT")
fileName = "jobs/histo1.root"
print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

xMin = 2.95
xMax = 3.25

# FUNCTIONS
fG = TF1('fG' , '[4] * TMath::Gaus(x, [0], [1]) + [2]*x + [3] ' , xMin , xMax )
fG.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fG.SetParNames ( 'Mass' , 'Sigma' , 'Slope' , 'Offset', 'Normalization' )
fG.SetLineColor ( 2 )
fG.SetLineWidth ( 3 )

fBW = TF1('fBW' , '[4] * TMath::BreitWigner(x, [0], [1]) + [2]*x + [3] ' , xMin , 3.25 )
fBW.SetParameters(3.1 , 1. , -4. , 150. , 30000.)
fBW.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'Normalization' )
fBW.SetLineColor ( 880 )
fBW.SetLineWidth ( 3 )

sigmaToGamma = 1/(2 * np.sqrt( 2 * np.log(2) ) )
fBWG = TF1('fBWG' , '[4] * TMath::BreitWigner(x, [0], [1]) * TMath::Gaus(x, [0], [5] ) + [2]*x + [3] ' , xMin , xMax )
fBWG.SetParameters(3.1 , 1. , -4. , 150. , 30000., 0.5 )
#fBWG.FixParameter(5, sigmaToGamma)
fBWG.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'NormalizationBW', 'Width2'  )
fBWG.SetLineColor ( 4 )
fBWG.SetLineWidth ( 2 )

fGplusG = TF1('fGplusG' , '[5] * TMath::Gaus(x, [0], [1]) + [6] *  TMath::Gaus(x, [0], [2] ) + [3]*x + [4] ' , xMin , xMax )
fGplusG.SetParameters(3.1 , 1. , 1. , -4. , 150. , 3000., 3000.)

fGplusG.SetParNames ( '#mu' , '#sigma_{1}' , '#sigma_{2}', 'Slope' , 'Offset', 'N_{1}', 'N_{2}' )
fGplusG.SetLineColor ( 2 )
fGplusG.SetLineWidth ( 3 )

fGG = TF1('fGG' , '[4] * TMath::Gaus(x, [0], [1]) + [5] *  TMath::Gaus(x, [7], [6] ) + [2]*x + [3] ' , xMin , xMax )
fGG.SetParameters(3.1 , 1. , -4. , 150. , 3000., 3000., 1.)
fGG.FixParameter(7, 3.0968)
fGG.SetParNames ( 'Mass' , 'Width' , 'Slope' , 'Offset', 'Normalization1', 'Normalization2', 'Width2',  'Mass2')
fGG.SetLineColor ( 7 )
fGG.SetLineWidth ( 2 )

####### Slimmed Muons - invariant mass
cInvM = TCanvas( 'cInvM' , 'cInvM' , 1000 , 600 )

cInvM.SetLeftMargin(0.14)
hInvMass = gROOT.FindObject('hInv')
hInvMass.GetXaxis().SetTitle("M_{inv} [GeV]")
DrawHisto(hInvMass, 19, 1, 2)

cInvM.cd()
hInvMass.Draw('hist')
cInvM.Draw()

cInvM.Print("Plots/MuMuInv25_05.png")

####### J/Psi Fit
cZoom = TCanvas("cZoom", "cZoom", 1200, 800)
hJPsi = hInvMass.Clone('hJPsi')
hJPsi.SetTitle( 'J/#Psi' )
hJPsi.GetXaxis().SetRangeUser( 2.9 , 3.3 )


print('----------------------- Gauss')
hJPsi.Fit( fG , 'R' ,"", 2.9 , 3.25 )

#print('----------------------- Breit-Wigner')
#hJPsi.Fit( fBW , 'R' ,"", 2.95 , 3.3 )
print('----------------------- Gauss * Breit-Wigner (same mass, different widths)')
hJPsi.Fit( fBWG , 'R' ,"", 2.93 , 3.23 )
print('----------------------- Gauss + Gauss (same mass, different widths)')
hJPsi.Fit( fGplusG , 'R' ,"", 2.8 , 3.25 )


sigmaJPsiSMG = np.abs( fGplusG.GetParameter(0) - pdgJPsi ) / np.sqrt( pdgJPsiErr**2 + fGplusG.GetParError(0)**2 )
relDiffJPsiSMG= 100* (pdgJPsi-fGplusG.GetParameter(0))/pdgJPsi

sigmaJPsiSM = np.abs( fBWG.GetParameter(0) - pdgJPsi ) / np.sqrt( pdgJPsiErr**2 + fBWG.GetParError(0)**2 )
relDiffJPsiSM = 100* (pdgJPsi-fBWG.GetParameter(0))/pdgJPsi

lZoom = TLegend( 0.6 , 0.65 , 0.9 , 0.9 )
lZoom.SetTextSize   (0.03)
lZoom.AddEntry( hJPsi    , 'Data' , 'f' )
#lZoom.AddEntry( fBW    , 'Breit-Wigner'  , 'l' )
lZoom.AddEntry( fG     , 'Gaussian'   , 'l' )
lZoom.AddEntry( fBWG   , 'Convolution', 'l' )
lZoom.AddEntry( hJPsi  , 'Gaussian and Breit-Wigner', '' )
#lZoom.AddEntry( fGplusG     , 'Sum of Gaussians'         , 'l' )


fBWGbckg = TF1('fBWGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
fGbckg = TF1('fGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
#fBWbckg = TF1('fBWbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
fGplusGbckg = TF1('fGplusGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)


Background(fBWG.GetParameter(2), fBWG.GetParameter(3), fBWGbckg , fBWG.GetLineColor() )
Background(fG.GetParameter(2), fG.GetParameter(3), fGbckg , fG.GetLineColor() )
Background(fGplusG.GetParameter(3), fGplusG.GetParameter(4), fGplusGbckg, fGplusG.GetLineColor() )

lZoom2 = TLegend( 0.57, 0.75, 0.9, 0.9 )
lZoom2.SetTextSize   (0.03)
lZoom2.AddEntry( hJPsi    , 'Data'              , 'f' )
lZoom2.AddEntry( fBWG     , 'Convolution'       , 'l' )
lZoom2.AddEntry( fBWGbckg , 'Linear background' , 'l' )

aveText = TPaveText(0.57, 0.55, 0.9, 0.75, "NDC")
aveText.SetFillColor(0)
aveText.SetTextSize(0.032)
aveText.SetBorderSize(1)

aveText.AddText(" ")
aveText.AddText("Fit Parameters:")
aveText.AddText("#mu :     (%.6f #pm %.6f) GeV" % (fBWG.GetParameter(0), fBWG.GetParError(0)))
aveText.AddText("#sigma :     (%.5f #pm %.5f) GeV" % (fBWG.GetParameter(5), fBWG.GetParError(5)))
aveText.AddText("#Gamma :     (%.5f #pm %.5f) GeV" % (fBWG.GetParameter(1), fBWG.GetParError(1)))
#aveText.AddText("Slope: %.3f #pm %.3f" % (fGplusG.GetParameter(3), fGplusG.GetParError(3)))
aveText.SetTextAlign(10)
aveText.SetTextFont(42) 

cZoom.cd()
hJPsi.Draw('hist')
fG.DrawCopy( 'same' )
fGplusG.DrawCopy( 'same' )
#fGbckg.DrawCopy( 'same' )
#fBWGGbckg.DrawCopy('same')
fBWG.DrawCopy('same')
fBWGbckg.DrawCopy('same')
lZoom.Draw()
cZoom.Draw()
cZoom.Print("Plots/JPsi_SM_25_05.png")

cZoom.Clear()

hJPsi.Draw('hist')
fBWGbckg.DrawCopy('same')
fBWG.DrawCopy('same')
aveText.Draw()
lZoom2.Draw()

cZoom.Print("Plots/JPsi_SM_25_05_fit_parameters.png")

#######  Particle IDs from 25 events
cID = TCanvas('cID','cID',1000,600)
cID.SetLeftMargin(0.14)
hID = gROOT.FindObject('hID')
hID.SetTitle('Particle ID (Sum of 25 Events)')
cID.cd()
DrawHisto(hID, 19)
hID.GetXaxis().SetRangeUser(-5.5, hID.GetXaxis().GetXmax())

hID.Draw('hist')

cID.Draw()
cID.Print("Plots/Particle_ID_26_05.png")

#######  Slimmed Displaced Muonsh
hDisp = gROOT.FindObject('hDisp')
hDisp.SetTitle( 'Invariant mass of #mu^{-}#mu^{+} pairs (Slimmed Displaced Muons)' )
cDisp = TCanvas('cDisp','cDisp',1000,600)
cDisp.SetLeftMargin(0.14)
hDisp.GetXaxis().SetTitle("M_{inv} [GeV]")
cDisp.cd()


DrawHisto(hDisp, 19)

hDisp.Draw('hist')
#fG.DrawCopy( 'same' )
#fBW.DrawCopy( 'same' )
cDisp.Draw()
cDisp.Print("Plots/MuMuInv_SDM_26_05.png")

####### Vertex (distance for muons)
cVertex = TCanvas( 'cVertex' , 'cVertex' , 1000 , 600 )
cVertex.SetSupportGL(True)
cVertex.SetLeftMargin(0.14)

hVz = gROOT.FindObject('hVz')
hVzPC = gROOT.FindObject('hVzPC')
DrawHisto(hVz, 0, 4)
DrawHisto(hVzPC, 0, 3)

lVz = TLegend( 0.53 , 0.7 , 0.9 , 0.9 )
lVz.SetTextSize   (0.035)
lVz.AddEntry( hVz  , 'Slimmed Muons' , 'l' )
lVz.AddEntry( hVzPC  , 'Packed Candidates'  , 'l' )

cVertex.cd()
cVertex.SetLogy()
hVz.Draw('hist')
hVzPC.Draw('same')
lVz.Draw()
cVertex.Draw()
cVertex.Print("Plots/vertex_dist_PC_SM_26_05.png")

cProb = TCanvas( 'cProb' , 'cProb' , 1000 , 600 )

cProb.SetLeftMargin(0.14)
hProb = gROOT.FindObject('hProb')

DrawHisto(hProb, 19)

cProb.cd()
hProb.Draw('hist')
cProb.Draw()
cProb.Print("Plots/MuMu_SM_probCommonVx_26_05.png")

cJpsiProb = TCanvas( 'cJpsiProb' , 'cJpsiProb' , 1000 , 600 )

cJpsiProb.SetLeftMargin(0.14)
hJpsiProb = gROOT.FindObject('hJpsiProb')

DrawHisto(hJpsiProb, 19)

cJpsiProb.cd()
hJpsiProb.SetTitle('The probability of a common vertex for #mu^{+}#mu^{-} pairs, M_{inv} = m_{J/#Psi}')
hJpsiProb.Draw('hist')
cJpsiProb.Draw()
cJpsiProb.Print("Plots/JPsi_SM_probCommonVx_26_05.png")

cProbComp = TCanvas( 'cProbComp' , 'cProbComp' , 1000 , 600 )
hPall = hProb.Clone('hPall')
hPJpsi = hJpsiProb.Clone('hPJpsi')
cProbComp.SetLeftMargin(0.14)
#hProb = gROOT.FindObject('hProb')
#hJpsiProb = gROOT.FindObject('hJpsiProb')
DrawHisto(hPall, 0, 4 )
DrawHisto(hPJpsi,0,  3)

lProb = TLegend( 0.56 , 0.7 , 0.9 , 0.9 )
lProb.SetTextSize   (0.045)
lProb.AddEntry( hPall  , 'All #mu^{+}#mu^{-} pairs' , 'l' )
lProb.AddEntry( hPJpsi  , 'M_{inv} = m_{J/#Psi}'  , 'l' )

cProbComp.cd()
cProbComp.SetLogy()
hPall.GetYaxis().SetRangeUser(10**4, 2*10**7)
hPall.Draw('hist')
hPJpsi.Draw('same')
lProb.Draw()
cProbComp.Draw()
cProbComp.Print("Plots/Prob_commonV_JPsiVsAll_26_05.png")

#cJpsiProb = TCanvas( 'cJpsiProb' , 'cJpsiProb' , 1000 , 600 )

####### Vertices distance for muons from J/Psi decay for different probability cuts
cVz0x = TCanvas( 'cVz0x' , 'cVz0x' , 1000 , 600 )

hVz01 = gROOT.FindObject('hVz01')
hVz08 = gROOT.FindObject('hVz08')

print("AVERAGE DISCTANCE BETWEEN VERTICES FOR P>= 0.8, jpsi: ", hVz08.GetMean() )
print("AVERAGE DISCTANCE BETWEEN VERTICES FOR P>= 0.1, jpsi: ", hVz01.GetMean() )

DrawHisto(hVz01, 0, 4)

DrawHisto(hVz08, 0, 3)

lVz0x = TLegend( 0.7 , 0.75 , 0.9 , 0.9 )
lVz0x.AddEntry(hVzPC, "Probability cuts", "")

lVz0x.SetTextSize   (0.035)
lVz0x.AddEntry( hVz01   , 'P > 0.1' , 'l' )
lVz0x.AddEntry( hVz08    , 'P > 0.8'  , 'l' )

hVz01.GetXaxis().SetRangeUser( -0.02 ,0.1)
cVz0x.cd()
hVz01.SetTitle('The distance along the Z-axis between vertices for #mu^{+}#mu^{-} pairs, M_{inv} = m_{J/#Psi}')
hVz01.Draw('hist')
hVz08.Draw('same')
lVz0x.Draw()
cVz0x.Draw()
cVz0x.Print("Plots/Vertex_dist_prob_26_05.png")

hCan = gROOT.FindObject('hCan')
print ("Hello ROOT")
fileName5 = "jobs2/histos.root"
print ('Read data from: ', fileName5)
gROOT.Reset()
f5 = TFile(fileName5);
#f5.ls();

cComp = TCanvas('cComp','cComp',1200,600)
cComp.cd()


hInvG = gROOT.FindObject('hInvG')
hInvGC = hInvG.Clone('hInvGC')
hSM = hInvMass.Clone('hSM')
hSDM = hDisp.Clone('hSDM')

DrawHisto(hCan, 0, 3)
DrawHisto(hSM, 0, 4)
#hSM.SetFillColorAlpha( 0 , 0.4 )
DrawHisto(hSDM, 0, 798 )
#hSDM.SetFillColorAlpha( 0, 0.571)
DrawHisto(hInvGC, 0, 880 )

lComp = TLegend( 0.53 , 0.65 , 0.9 , 0.9 )
lComp.SetTextSize   (0.035)
lComp.AddEntry( hSM , 'Slimmed Muons' , 'l' )
lComp.AddEntry( hSDM, 'Slimmed Displaced Muons'  , 'l' )
lComp.AddEntry( hCan,  'Packed Candidates (Muons)'  , 'l' )
lComp.AddEntry( hInvGC , 'Slimmed Muons'  , 'l' )
lComp.AddEntry( hInvGC , '(P > 0.2, #Delta v_{z} < 0.3 cm, Global muons)'  , '' )


hSM.DrawCopy()
hSDM.DrawCopy('same')
hCan.DrawCopy('same')
hInvGC.DrawCopy('same')
lComp.Draw()
cComp.Draw()
cComp.Print("Plots/Comparison_PC_SM_SDM_26_05.png")

hInvG = gROOT.FindObject('hInvG')
DrawHisto(hInvG, 19 )
####### J/Psi Fit from global
cZoomJPsi = TCanvas("cZoomJPsi", "cZoomJPsi", 1200, 800)
hZoomJPsi = hInvG.Clone('hZoomJPsi')
hZoomJPsi.SetTitle( 'J/#Psi(1S)' )
hZoomJPsi.GetXaxis().SetRangeUser( 2.85 , 3.35 )


range1 = 2.9
range2 = 3.25

fBWG.SetRange( range1, range2)

#fBWG.SetParameters(3.1 , 1. , -4. , 150. , 30000., 0.5 )
hZoomJPsi.Fit( fGplusG , 'R' ,"", range1 , range2)
hZoomJPsi.Fit( fBWG , 'R' ,"", range1 , range2)
sigmaJPsiGLB = np.abs( fBWG.GetParameter(0) - pdgJPsi ) / np.sqrt( pdgJPsiErr**2 + fBWG.GetParError(0)**2 )
sigmaJPsiGLBG = np.abs( fGplusG.GetParameter(0) - pdgJPsi ) / np.sqrt( pdgJPsiErr**2 + fGplusG.GetParError(0)**2 )
relDiffJPsiGLB = 100* (pdgJPsi-fBWG.GetParameter(0))/pdgJPsi
relDiffJPsiGLBG= 100* (pdgJPsi-fGplusG.GetParameter(0))/pdgJPsi

fBWGbckg = TF1('fBWGbckg', '[0]*x + [1]', xMin-0.05, xMax+0.05)
Background(fBWG.GetParameter(2), fBWG.GetParameter(3), fBWGbckg , fBWG.GetLineColor() )

lZoomJPsi = TLegend( 0.57 , 0.75 , 0.9 , 0.9 )
lZoomJPsi.SetTextSize   (0.03)
lZoomJPsi.AddEntry( hZoomJPsi    , 'Data' , 'f' )
lZoomJPsi.AddEntry( fBWG    , 'Convolution'  , 'l' )
lZoomJPsi.AddEntry( fBWGbckg   , 'Linear background'  , 'l' )


aveText1 = TPaveText(0.57, 0.57, 0.9, 0.75, "NDC")
aveText1.SetFillColor(0)
aveText1.SetTextSize(0.032)
aveText1.SetBorderSize(1)

aveText1.AddText(" ")
aveText1.AddText("Fit Parameters:")
aveText1.AddText("#mu :     (%.6f #pm %.6f) GeV" % (fBWG.GetParameter(0), fBWG.GetParError(0)))
aveText1.AddText("#sigma :     (%.5f #pm %.5f) GeV" % (np.abs(fBWG.GetParameter(5)), fBWG.GetParError(5)))
aveText1.AddText("#Gamma :     (%.5f #pm %.5f) GeV" % (np.abs(fBWG.GetParameter(1)), fBWG.GetParError(1)))

aveText1.SetTextAlign(10)
aveText1.SetTextFont(42) 

cZoomJPsi.cd()
hZoomJPsi.Draw('hist')

fBWG.DrawCopy( 'same' )
fBWGbckg.DrawCopy( 'same' )
aveText1.Draw()
lZoomJPsi.Draw()
cZoomJPsi.Draw()
cZoomJPsi.Print("Plots/GLB_JPSI_26_05.png")


cZoomPsi = TCanvas("cZoomPsi", "cZoomPsi", 1200, 800)
hZoomPsi = hInvG.Clone('hZoomPsi')
hZoomPsi.SetTitle( '#Psi(2S)' )
hZoomPsi.GetXaxis().SetRangeUser( 3.45 , 3.9 )


range1 = 3.5
range2 = 3.85
#fG.SetRange( range1, range2)
fBWG.SetRange( range1, range2)
fBWG.SetParameters(3.68 , 1. , -100. , 300. , 8000., 0.5 )
#fGplusG.SetParameters(3.68 , 0.1 , -200. , 500. , 9000., 9000., 1.)
hZoomPsi.Fit( fBWG , 'R' ,"", range1 , range2)
sigmaPsiGLB = np.abs( fBWG.GetParameter(0) - pdgPsi ) / np.sqrt( pdgPsiErr**2 + fBWG.GetParError(0)**2 )
relDiffPsiSM = 100* (pdgPsi-fBWG.GetParameter(0))/pdgPsi

lZoomPsi = TLegend( 0.57 , 0.75 , 0.9 , 0.9 )
lZoomPsi.SetTextSize   (0.03)
lZoomPsi.AddEntry( hZoomPsi    , 'Data' , 'f' )
lZoomPsi.AddEntry( fBWG        , 'Convolution'  , 'l' )
lZoomPsi.AddEntry( fBWGbckg    , 'Linear background'  , 'l' )


aveText1 = TPaveText(0.57, 0.57, 0.9, 0.75, "NDC")
aveText1.SetFillColor(0)
aveText1.SetTextSize(0.032)
aveText1.SetBorderSize(1)

aveText1.AddText(" ")
aveText1.AddText("Fit Parameters:")
aveText1.AddText("#mu :     (%.5f #pm %.5f) GeV" % (fBWG.GetParameter(0), fBWG.GetParError(0)))
aveText1.AddText("#sigma :     (%.4f #pm %.4f) GeV" % (np.abs(fBWG.GetParameter(5)), fBWG.GetParError(5)))
aveText1.AddText("#Gamma :     (%.4f #pm %.4f) GeV" % (np.abs(fBWG.GetParameter(1)), fBWG.GetParError(1)))

aveText1.SetTextAlign(10)
aveText1.SetTextFont(42) 
fBWGbckg = TF1('fBWGbckg', '[0]*x + [1]', range1-0.05, range2+0.05)
Background(fBWG.GetParameter(2), fBWG.GetParameter(3), fBWGbckg , fBWG.GetLineColor() )


cZoomPsi.cd()
hZoomPsi.Draw('hist')
fBWG.DrawCopy( 'same' )
fBWGbckg.DrawCopy( 'same' )
aveText1.Draw()
lZoomPsi.Draw()
cZoomPsi.Print("Plots/GLB_PSI_2S_26_05.png")

#### UPSILON FAMILY
range1 = 8.95
range2 = 10.65

fUps = TF1('fUps' , '[6] * TMath::BreitWigner(x, [0], [1]) * TMath::Gaus(x, [0], [11])  + [7] * TMath::Gaus(x, [2], [3] ) + [8] * TMath::Gaus(x, [4], [5] )+ [9]*x + [10] ' , range1 , range2 +0.05 )
fUps.SetParameters( 9.46 , 1. ,  10.02 , 0.1, 10.58, 0.2, 2400., 1000. , 1200. , -200., 500.)
fUps.SetParameter(11, 0.5)
fUps.SetParNames ( 'Mass1' , 'Width1' ,'Mass2' , 'Sigma2' ,'Mass3' , 'Sigma3' , 'Normalization1',  'Normalization2', 'Normalization3' , 'Slope' , 'Offset' )
fUps.SetLineColor ( 4 )
fUps.SetLineWidth ( 3 )

cUpsilon1 = TCanvas("cUpsilon1", "cUpsilon1", 1200, 800)
hUpsilon1 = hInvG.Clone('hUpsilon1')
hUpsilon1.SetTitle( '#varUpsilon(1S), #varUpsilon(2S), #varUpsilon(3S)' )
hUpsilon1.GetXaxis().SetRangeUser( 8.8 , 10.8 )

hUpsilon1.Fit( fUps, 'R' ,"", range1 , range2 )
sigmaUps1GLB = np.abs( fUps.GetParameter(0) - pdgUps1 ) / np.sqrt( pdgUps1Err**2 + fUps.GetParError(0)**2 )
relDiffUps1GLB = 100* (pdgUps1-fUps.GetParameter(0))/pdgUps1
sigmaUps2GLB = np.abs( fUps.GetParameter(2) - pdgUps2 ) / np.sqrt( pdgUps2Err**2 + fUps.GetParError(2)**2 )
relDiffUps2GLB = 100* (pdgUps2-fUps.GetParameter(2))/pdgUps2
sigmaUps3GLB = np.abs( fUps.GetParameter(4) - pdgUps3 ) / np.sqrt( pdgUps3Err**2 + fUps.GetParError(4)**2 )
relDiffUps3GLB = 100* (pdgUps3-fUps.GetParameter(4))/pdgUps3

fUpsBckg = TF1('fUpsBckg', '[0]*x + [1]', range1 - 0.05, range2 + 0.05 )
fGplusGbckg = TF1('fGplusGbckg', '[0]*x + [1]', range1-0.05, 9.8 )

Background(fUps.GetParameter(9), fUps.GetParameter(10), fUpsBckg , fUps.GetLineColor() )

lUpsilon1 = TLegend( 0.57 , 0.8 , 0.9 , 0.9 )
lUpsilon1.SetTextSize   (0.03)
lUpsilon1.AddEntry( hUpsilon1, 'Data'                            , 'f' )
lUpsilon1.AddEntry( fUps     , 'Convolution, Gaussians'  , 'l' )
lUpsilon1.AddEntry( fUpsBckg  , 'Linear background'           , 'l' )

aveText2 = TPaveText(0.57, 0.4, 0.9, 0.8, "NDC")
aveText2.SetFillColor(0)
aveText2.SetTextSize(0.03)
aveText2.SetBorderSize(1)

aveText2.AddText(" ")
aveText2.AddText("Fit Parameters:")
aveText2.AddText("#mu_{(1S)} :     (%.5f #pm %.5f) GeV" % (fUps.GetParameter(0), fUps.GetParError(0)))
aveText2.AddText("#sigma_{(1S)} :     (%.4f #pm %.4f) GeV" % (np.abs(fUps.GetParameter(11)), fUps.GetParError(11)))
aveText2.AddText("#Gamma_{(1S)} :     (%.4f #pm %.4f) GeV" % (np.abs(fUps.GetParameter(1)), fUps.GetParError(1)))
aveText2.AddText("#mu_{(2S)} :     (%.4f #pm %.4f) GeV" % (fUps.GetParameter(2), fUps.GetParError(2)))
aveText2.AddText("#sigma_{(2S)} :     (%.4f #pm %.4f) GeV" % (np.abs(fUps.GetParameter(3)), fUps.GetParError(3)))
aveText2.AddText("#mu_{(3S)} :     (%.4f #pm %.4f) GeV" % (fUps.GetParameter(4), fUps.GetParError(4)))
aveText2.AddText("#sigma_{(3S)} :     (%.4f #pm %.4f) GeV" % (np.abs(fUps.GetParameter(5)), fUps.GetParError(5)))

aveText2.SetTextAlign(10)
aveText2.SetTextFont(42) 

cUpsilon1.cd()
hUpsilon1.Draw('hist')

fG.DrawCopy( 'same' )
fGbckg.DrawCopy( 'same' )

fGplusG.DrawCopy( 'same' )
fGplusGbckg.DrawCopy( 'same' )

fUps.DrawCopy( 'same' )
fUpsBckg.DrawCopy( 'same' )
aveText2.Draw("same")
lUpsilon1.Draw()
cUpsilon1.Draw()
cUpsilon1.Print("Plots/GLB_Upsilon_26_05.png")

cInvG = TCanvas( 'cInvG' , 'cInvG' , 1000 , 600 )
hInvG.SetTitle( 'Invariant mass of #mu^{-}#mu^{+} pairs (Slimmed muons, global)' )
cInvG.SetLeftMargin(0.14)


cInvG.cd()
hInvG.Draw('hist')
cInvG.Draw()
cInvG.Print("Plots/MuMuInv_GLB_26_05.png")
'''
hExp = gROOT.FindObject('hExp')
cExp = TCanvas( 'cExp' , 'cExp' , 1000 , 600 )
DrawHisto(hExp, 19)
cExp.SetLeftMargin(0.14)


cExp.cd()
hExp.Draw('hist')
cExp.Draw()

'''

cInvMuMuK = TCanvas( 'cInvMuMuK' , 'cInvMuMuK' , 1000 , 600 )
hInvMuMuK = gROOT.FindObject('hInvMuMuK')

cInvMuMuK.SetLeftMargin(0.14)


DrawHisto(hInvMuMuK,19)
cInvMuMuK.cd()

hInvMuMuK.Draw('hist')
#fG.DrawCopy('same')

#fG2.DrawCopy('same')
#fGplusG.DrawCopy('same')
cInvMuMuK.Draw()
cInvMuMuK.Print("Plots/MuMuKInv_GLB_Muons_Kaonp_26_05.png")

hInvMuMuK = hInvMuMuK.Rebin(2)

cZoomB = TCanvas("cZoomB", "cZoomB", 1200, 800)
hZoomB = hInvMuMuK.Clone('hZoomB')
hZoomB.SetTitle( 'B^{+}' )
hZoomB.GetXaxis().SetRangeUser( 4.8 , 5.8 )


range1 = 5.05
range2 = 5.5
fG.SetRange( range1, range2)
fBWG.SetRange( range1, range2)
fBWG.SetParameters(5.3 , 0.1 , -100. , 200. , 500., 0.1 )
fG.SetParameters(5.3 , 0.1 , -100. , 500. , 750.)
#hZoomB.Fit( fBWG , 'R' ,"", range1 , range2)
hZoomB.Fit( fG , 'R' ,"", range1 , range2)
sigmaBp = np.abs( fG.GetParameter(0) - pdgBp ) / np.sqrt( pdgBpErr**2 + fG.GetParError(0)**2 )
relDiffBp = 100* (pdgBp-fG.GetParameter(0))/pdgBp

aveText3 = TPaveText(0.57, 0.57, 0.9, 0.75, "NDC")
aveText3.SetFillColor(0)
aveText3.SetTextSize(0.032)
aveText3.SetBorderSize(1)

aveText3.AddText(" ")
aveText3.AddText("Fit Parameters (Gaussian):")
aveText3.AddText("#mu :     (%.4f #pm %.4f) GeV" % (fG.GetParameter(0), fG.GetParError(0)))
aveText3.AddText("#sigma :     (%.4f #pm %.4f) GeV" % (np.abs(fG.GetParameter(1)), fG.GetParError(1)))
#aveText3.AddText("#Gamma :     (%.4f #pm %.4f) GeV" % (np.abs(fG.GetParameter(1)), fG.GetParError(1)))

aveText3.SetTextAlign(10)
aveText3.SetTextFont(42) 

fGbckg = TF1('fGbckg', '[0]*x + [1]', range1-0.05, range2+0.05)
Background(fG.GetParameter(2), fG.GetParameter(3), fGbckg , fG.GetLineColor() )

lZoomB = TLegend( 0.57 , 0.75 , 0.9 , 0.9 )
lZoomB.SetTextSize   (0.03)
lZoomB.AddEntry( hZoomB  , 'Data' , 'f' )
lZoomB.AddEntry( fG      , 'Gaussian'  , 'l' )
lZoomB.AddEntry( fGbckg  , 'Linear background'  , 'l' )

cZoomB.cd()
hZoomB.Draw('hist')
fG.DrawCopy( 'same' )
fGbckg.DrawCopy( 'same' )
aveText3.Draw()
lZoomB.Draw()
cZoomB.Print("Plots/B_meson_26_05.png")

cZoomX = TCanvas("cZoomX", "cZoomX", 1200, 800)
hZoomX = hInvMuMuK.Clone('hZoomX')
#hZoomX.SetTitle( '#chi_{c0}' )
hZoomX.GetXaxis().SetRangeUser( 4. , 5. )


range1 = 4.25
range2 = 4.65
fG.SetRange( range1, range2)
fBWG.SetRange( range1, range2)
fBWG.SetParameters(4.4 , 0.1 , -500. , 200. , 1000., 0.1 )
fG.SetParameters(4.4 , 0.1 , -500. , 500. , 1000.)
hZoomX.Fit( fBWG , 'R' ,"", range1 , range2)
hZoomX.Fit( fG , 'R' ,"", range1 , range2)
sigmaX = np.abs( fG.GetParameter(0) - pdgX ) / np.sqrt( pdgXErr**2 + fG.GetParError(0)**2 )

aveText4 = TPaveText(0.57, 0.57, 0.9, 0.75, "NDC")
aveText4.SetFillColor(0)
aveText4.SetTextSize(0.032)
aveText4.SetBorderSize(1)

aveText4.AddText(" ")
aveText4.AddText("Fit Parameters (Gaussian):")
aveText4.AddText("#mu :     (%.4f #pm %.4f) GeV" % (fG.GetParameter(0), fG.GetParError(0)))
aveText4.AddText("#sigma :     (%.4f #pm %.4f) GeV" % (np.abs(fG.GetParameter(1)), fG.GetParError(1)))
#aveText4.AddText("#Gamma :     (%.4f #pm %.4f) GeV" % (np.abs(fG.GetParameter(1)), fG.GetParError(1)))

aveText4.SetTextAlign(10)
aveText4.SetTextFont(42) 

fBWGbckg = TF1('fBWGbckg', '[0]*x + [1]', range1-0.05, range2+0.05)
Background(fBWG.GetParameter(2), fBWG.GetParameter(3), fBWGbckg , fBWG.GetLineColor() )
fGbckg = TF1('fGbckg', '[0]*x + [1]', range1-0.05, range2+0.05)
Background(fG.GetParameter(2), fG.GetParameter(3), fGbckg , fG.GetLineColor() )

lZoomX = TLegend( 0.57 , 0.75 , 0.9 , 0.9 )
lZoomX.SetTextSize   (0.03)
lZoomX.AddEntry( hZoomX    , 'Data' , 'f' )
lZoomX.AddEntry( fG        , 'Gaussian'  , 'l' )
lZoomX.AddEntry( fGbckg        , 'Linear background'  , 'l' )

cZoomX.cd()
hZoomX.Draw('hist')
fG.DrawCopy( 'same' )
fGbckg.DrawCopy( 'same' )
aveText4.Draw()
lZoomX.Draw()
cZoomX.Print("Plots/X_26_05.png")


cTest = TCanvas( 'cTest' , 'cTest' , 1000 , 600 )
hTest = gROOT.FindObject('hTest')

cTest.SetLeftMargin(0.14)
hTest = hTest.Rebin(8)

DrawHisto(hTest,19)
cTest.cd()

hTest.Draw('hist')
#fG.DrawCopy('same')

#fG2.DrawCopy('same')
#fGplusG.DrawCopy('same')
cTest.Draw()
cTest.Print("Plots/Test_GLB_Muons_Pion_26_05.png")

cMuMuKK = TCanvas( 'cMuMuKK' , 'cMuMuKK' , 1000 , 600 )
hMuMuKK = gROOT.FindObject('hMuMuKK')

cMuMuKK.SetLeftMargin(0.14)
hMuMuKK = hMuMuKK.Rebin(5)
DrawHisto(hMuMuKK,19)
cMuMuKK.cd()

hMuMuKK.Draw('hist')
#fG.DrawCopy('same')

#fG2.DrawCopy('same')
#fGplusG.DrawCopy('same')
cMuMuKK.Draw()
cMuMuKK.Print("Plots/MuMuKK19_06.png")

cKK = TCanvas( 'cKK' , 'cKK' , 1000 , 600 )
hKK = gROOT.FindObject('hKK')

cKK.SetLeftMargin(0.14)
hKK = hKK.Rebin(10)

DrawHisto(hKK,19)
cKK.cd()

hKK.Draw('hist')
#fG.DrawCopy('same')

#fG2.DrawCopy('same')
#fGplusG.DrawCopy('same')
cKK.Draw()
cKK.Print("Plots/KK19_06.png")

##### SIGMA


print("########################## ")
print( "J/Psi, Slimmed Muons: ", sigmaJPsiSM )
print( "J/Psi, Slimmed Muons, GLB: ", sigmaJPsiSMG )
print( "J/Psi, Slimmed Muons: ", sigmaJPsiGLB )
print( "J/Psi, Slimmed Muons, GLB: ", sigmaJPsiGLBG )
print( "Psi, Slimmed Muons, GLB: ", sigmaPsiGLB )
print( "Upsilon 1, Slimmed Muons, GLB: ", sigmaUps1GLB )
print( "Upsilon 2, Slimmed Muons, GLB: ", sigmaUps2GLB )
print( "Upsilon 3, Slimmed Muons, GLB: ", sigmaUps3GLB )
print( "B meson +, Slimmed Muons: ", sigmaBp )
print( "X , Slimmed Muons: ", sigmaX )
print("Relative mass difference in %")
print( "J/Psi, Slimmed Muons: ", relDiffJPsiSM)
print( "J/Psi, Slimmed Muons, GLB: ", relDiffJPsiSMG )
print( "J/Psi, Slimmed Muons: ", relDiffJPsiGLB)
print( "J/Psi, Slimmed Muons, GLB: ", relDiffJPsiGLBG )
print( "Psi, Slimmed Muons, GLB: ", relDiffPsiSM )
print( "Upsilon 1, Slimmed Muons, GLB: ", relDiffUps1GLB )
print( "Upsilon 2, Slimmed Muons, GLB: ", relDiffUps2GLB )
print( "Upsilon 3, Slimmed Muons, GLB: ", relDiffUps3GLB )
print( "B meson +, Slimmed Muons: ", relDiffBp )

#print( '',sigmaToGamma*fGplusG.GetParameter(1) )
#print( '',sigmaToGamma*fGplusG.GetParameter(2) )


input('press enter to exit')
