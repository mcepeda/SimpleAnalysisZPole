#!/usr/bin/env python
import ROOT 
from ROOT import TH1F,TFile

ROOT.gStyle.SetOptStat(0)

MESON="#rho (2#gamma)"
cuts="GEN" #2_0.1_2"

NGEN=100e6 # 8338055  #1000*1000 #1000 files  #976000+9.99e5 
xsec=1476.58*1000 # en fb
lumi=NGEN/xsec/ 1000. # 17 ab-1
scale= 1 # lumi*xsec/NGEN

print (scale)

#file=TFile("../makePlots/histos_KKM91p2_test_changeAtau_dm1.root")
#file=TFile("../makePlots/histos_KKM90_test_changeAtau_dm1.root")
#file=TFile("../makePlots/histos_PY8GENWI23_test_changeAtau_dm1.root")

file=TFile("../makePlots/histos_PY8GENWI23_LONG_changeAtau_dm1_02315_v2.root")
#file=TFile("../makePlots/histos_PY8GENWI23_LONG_changeAtau_dm1.root")
#file=TFile("../makePlots/histos_PY8GENWI23_LONG_changeAtau_dm1_1498_v2.root")


#fileOut=TFile("BINED_templates_KKMC91p2_GEN_Sept_1.root","RECREATE")
#fileOut=TFile("BINED_templates_KKMC90_GEN_Sept_1.root","RECREATE")
fileOut=TFile("BINED_templates_PY8WI23_GEN_LONG.root","RECREATE")
#GENOmegaCosTheta_M1_TAUMINUS


#var="GENOmegaCosThetaMeson"
var="GENOmegaCosTheta"
title="#omega_{#rho}"

rebin=1

#sample=["TAUMINUS","PosHel_TAUMINUS","NegHel_TAUMINUS"]
sample=["TAUMINUS","P1_TAUMINUS","M1_TAUMINUS"]
#sample2=["TAUPLUS","P1_TAUPLUS","M1_TAUPLUS"]
color=[ROOT.kBlack,ROOT.kGreen+2,ROOT.kRed]
sampleName=["SM","A_{#tau}=+1","A_{#tau}=-1"]

nBins=50  
binLength=100/nBins

for bin in range(0,nBins):

  print ("Hey!")

  c=ROOT.TCanvas("canvas","",800,800)
  leg=ROOT.TLegend(0.5,0.89,0.9,0.6)
  leg.SetFillStyle(0)
  leg.SetLineColor(0)
  leg.SetLineWidth(0)

  
  binIni=int(bin*binLength+1)
  binEnd=int((bin+1)*binLength)
  print (bin,binIni,binEnd)

  maxY=0

  histo_bin={}

  for i in range(0,len(sample)):
    print (var+"_"+sample[i])
    histo2D=file.Get(var+"_"+sample[i])
#    histo2D_2=file.Get(var+"_"+sample2[i])
#    histo2D.Add(histo2D_2)
    print ("histo_"+sample[i]+"_"+str(bin)) 
    histo_bin[i]=histo2D.ProjectionX("histo_"+sample[i]+"_"+str(bin),binIni,binEnd)
    histo_bin[i].Rebin(rebin)
    histo_bin[i].Scale(scale)
    histo_bin[i].SetXTitle(title)  
    histo_bin[i].SetLineWidth(2)
    histo_bin[i].SetLineColor(color[i])

    if sample[i]=="BG":
      histo_bin[i].SetFillColor(color[i])
      #histo_bin[i].SetFillStyle(3004)
      leg.AddEntry(histo_bin[i],sampleName[i],"f")
    else: 
      leg.AddEntry(histo_bin[i],sampleName[i],"l")

    print("Integral: ",sampleName, histo_bin[i].Integral())

    if maxY<histo_bin[i].GetMaximum():
      maxY=histo_bin[i].GetMaximum()
      
    if i==0:
      histo_bin[i].Draw("hist")
    else:
      histo_bin[i].Draw("hist,same")

  histo_bin[0].SetMaximum(maxY*1.2)
  leg.Draw()
  c.Draw()
  c.SaveAs("BINS/BIN_"+var+"_"+cuts+"_"+str(binIni)+"_"+str(binEnd)+".png")

  fileOut.cd()
  for i in range(0,len(sample)):
    histo_bin[i].Write()
  
  print ("oh?")


# Also for the full range
c = ROOT.TCanvas("canvas_full", "", 800, 800)
leg = ROOT.TLegend(0.5, 0.89, 0.9, 0.6)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.SetLineWidth(0)

maxY = 0

histo_full = {}

for i in range(0, len(sample)):

  histo2D = file.Get(var + "_" + sample[i])
  print("histo_" + sample[i] + "_full")
  histo_full[i] = histo2D.ProjectionX("histo_" + sample[i] + "_full", 1, 100)
  histo_full[i].Rebin(rebin)
  histo_full[i].Scale(scale)
  histo_full[i].SetXTitle(title)
  histo_full[i].SetLineWidth(2)
  histo_full[i].SetLineColor(color[i])

  if sample[i] == "BG":
    histo_full[i].SetFillColor(color[i])
    # histo_full[i].SetFillStyle(3004)
    leg.AddEntry(histo_full[i], sampleName[i], "f")
  else:
    leg.AddEntry(histo_full[i], sampleName[i], "l")

  print("Integral: ", sampleName, histo_full[i].Integral())

  if maxY < histo_full[i].GetMaximum():
    maxY = histo_full[i].GetMaximum()

  if i == 0:
    histo_full[i].Draw("hist")
  else:
    histo_full[i].Draw("hist,same")

histo_full[0].SetMaximum(maxY * 1.2)
leg.Draw()
c.Draw()
c.SaveAs("BINS/BIN_" + var + "_" + cuts + "_full.png")

fileOut.cd()
for i in range(0, len(sample)):
  histo_full[i].Write()


file.Close()
