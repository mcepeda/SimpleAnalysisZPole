import sys, os, math 
from array import array
import ROOT
from ROOT import TFile, TTree, TH1F, TH2F
import numpy as np
from pathlib import Path
import ctypes

import argparse
parser = argparse.ArgumentParser(
    description="Configure the analysis",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)


#parser.add_argument("-f", "--sample", default="../../testPlotsKK/SAVETREES/out_test_miniKKMCee2_mondayLONG_TREE.root")
#parser.add_argument("-o", "--outfile", default="histos_KKM91p2_test_changeAtau_dm1.root")

#parser.add_argument("-f", "--sample", default="../makeTrees/TREE_KKMCee_ONLYSIM_10M_90GeV_V2_partial.root")
#parser.add_argument("-o", "--outfile", default="histos_KKM90_test_changeAtau_dm1.root")

parser.add_argument("-f", "--sample", default="../makeTrees/TREE_p8_ee_Ztautau_ecm91_V2_partial.root")
parser.add_argument("-o", "--outfile", default="histos_PY8GENWI23_test_changeAtau_dm1.root")
parser.add_argument("-dP", "--decayP", default=1, type=int)
parser.add_argument("-dM", "--decayM", default=1, type=int)
args = parser.parse_args()
config = vars(args)
print(config)

filename = args.sample
outfilename = args.outfile


f = TFile.Open(filename)
tree=f.Get("outtree")

outfile = TFile.Open(outfilename, "RECREATE")

selectDecayP = args.decayP
selectDecayM = args.decayM


h = {}

def add1(name, nb, xmin, xmax):
    h[name] = ROOT.TH1F(name, "", nb, xmin, xmax)

def add2(name, nbx, xmin, xmax, nby, ymin, ymax):
    h[name] = ROOT.TH2F(name, "", nbx, xmin, xmax, nby, ymin, ymax)

def add1Variab(name, nbins):
    edges = array('d')

    for i in range(nbins + 1):
        t = float(i) / nbins          # 0 -> 1
        x = 0.5 * (1.0 - (1.0 - t)**2)
        edges.append(x)

    h[name] = ROOT.TH1F(name, "", nbins, edges)
    h[name].Sumw2()

add1("hEvents", 2, 0, 2)
add2(f"SpinWTvsSpinWThelApprox", 100, 0, 2, 100, 0, 2)
add1(f"SpinWT", 100, 0, 2)
add1(f"SpinWThelApprox", 100, 0, 2)

add1("GENZMass",500,80,95)
add1("GENZVisMass",500,0,95)
add1("GENZP",100,0,100)
add1("GENZTheta",100,0,3.2)
add1("GENZCosTheta",100,-1,1)
add1("GENZPhi",100,-3.2,3.2)
add1("GENZVisP",100,0,100)

add1("GENXCorr",100,-0.5,0.5)
add1("GENSqrtCorr",500,80,95)
add2("GENZMassVsSqrtCorr",500,80,95,500,80,95)

add1("GENCosThetaSimple",100,-1,1)
add1("GENCosThetaHat",100,-1,1)
add2("GENCosThetaHatVsSimple",100,-1,1,100,-1,1)

add1("GENSumPhotonE",100,0,50)
add1("GENZRecoil",500,80,95)
add1("GENZRecoil2",500,80,95)

add1("GENPhotonP",100,0,50)
add1("GENPhotonMom",50,-25,25)
add1("GENPhotonPISR",100,0,50)
add1("GENPhotonPFSR",100,0,50)
add1("GENPhotonDRISR",500,0,6.28)
add1("GENPhotonDRFSR",500,0,6.28)
add1("GENPhotonDR",500,0,6.28)
add2("GENPhotonPDRISR",100,0,50,500,0,6.28)
add2("GENPhotonPDRFSR",100,0,50,500,0,6.28)
add1("GENNPhoton",5,0,5)
add1("GENPhotonDRTau1",500,0,6.28)
add1("GENPhotonDRTau2",500,0,6.28)
add1("GENTausDR",500,0,6.28)

#for tag in ["", "_PosHel", "_NegHel", "_P1", "_M1"]:
#    add1("GENQX" + tag, 100, -1, 1)

for charge in ["TAUMINUS", "TAUPLUS"]:
    for tag in ["", "_PosHel", "_NegHel", "_P1", "_M1"]:
        add1(f"GENMesonEOverBeamE{tag}_{charge}", 50, 0, 1)
        add1(f"GENOmega{tag}_{charge}", 100, -1, 1)

        add2(f"GENOmegaCosTheta{tag}_{charge}", 100, -1, 1,100, -1, 1)
        add2(f"GENOmegaCosThetaHat{tag}_{charge}", 100, -1, 1,100, -1, 1)
        add2(f"GENOmegaCosThetaMeson{tag}_{charge}", 100, -1, 1,100, -1, 1)


        add1(f"GENCosTheta{tag}_{charge}", 100, -1, 1)
        add1(f"GENCosThetaHat{tag}_{charge}", 100, -1, 1)
        add1(f"GENCosThetaStar{tag}_{charge}", 100, -1, 1)
        add2(f"GENMesonEOverBeamE_CosTheta{tag}_{charge}", 10, 0, 1, 4, -1, 1)
        add1(f"GENMesonMass{tag}_{charge}", 100, 0, 2)
        add1(f"GENTauP{tag}_{charge}", 100, 0, 50)
        add1(f"GENTauTheta{tag}_{charge}", 40, -3.2, 3.2)
        add1(f"GENTauPhi{tag}_{charge}", 40, -3.2, 3.2)
        add1(f"GENTauDecayMode{tag}_{charge}", 31, -15, 15)
        add1(f"GENAbsZ{tag}_{charge}", 20, 0, 1)
        add1Variab(f"GENFZ{tag}_{charge}", 20)# , 0, 0.5) 
        add1Variab(f"GENFZ_BW{tag}_{charge}", 20)#, 0, 0.5)
        add1Variab(f"GENFZ_FW{tag}_{charge}", 20)#, 0, 0.5)
        add1(f"GENAbsZ_BW{tag}_{charge}", 20, 0, 1)
        add1(f"GENAbsZ_FW{tag}_{charge}", 20, 0, 1) 

        add1(f"GENFZFine{tag}_{charge}", 1000, 0, 0.5)
        add1(f"GENFZFine_BW{tag}_{charge}", 1000, 0, 0.5)
        add1(f"GENFZFine_FW{tag}_{charge}", 1000, 0, 0.5)






def build_p4(p, theta, phi, mass):
    p4 = ROOT.TLorentzVector()

    px = p * math.sin(theta) * math.cos(phi)
    py = p * math.sin(theta) * math.sin(phi)
    pz = p * math.cos(theta)

    p4.SetPxPyPzE(px, py, pz, math.sqrt(p*p + mass*mass))
    return p4

def dRAngle(p1,p2):
   dphi=abs(p1.Phi()-p2.Phi())
   if (dphi>math.pi) : dphi=2*math.pi-dphi
   dtheta=p1.Theta()-p2.Theta()
   dR=math.sqrt(dtheta*dtheta+dphi*dphi)
   return dR


totalEvents=0
selectedEvents=0

sumWeightsP1=0
sumWeightsM1=0
sumWeights=0

for ev in tree:

    if totalEvents % 100000 == 0:
        print(totalEvents)

#    if totalEvents> 1000000:
#        break    

    totalEvents += 1
    dmplus = getattr(ev, f"genTauID_plus")
    dmminus=getattr(ev, f"genTauID_minus") 
    
    if selectDecayM!=-999 and dmminus !=selectDecayM:
        continue

    if selectDecayP!=-999 and dmplus !=selectDecayP:
        continue

    if  dmminus not in (-11,-13,0,1):
        continue

    if dmplus not in (-11,-13,0,1):
        continue


    wP1 = ev.weight_P1_plus * ev.weight_P1_minus
    wM1 = ev.weight_M1_plus * ev.weight_M1_minus

    AeSM=0.1472 # 0.1498  #0.1472
    AtauSM=AeSM
    New_Atau=1
    New_Ae= AeSM # New_Atau # AeSM
    New_AtauM1=-1
    New_AeM1= AeSM # New_AtauM1 # AeSM
    AFBSM = 3/4 * AeSM * AtauSM
    AFBP1 = 3/4 * New_Ae  * New_Atau    # = 3/4
    AFBM1 = 3/4 * New_AeM1 * New_AtauM1 # = 3/4 * (-1)*(-1) = 3/4 

    Z=math.cos(ev.genTauTheta_minus) # this is the theta of the Tau, not the meson 

    cosThetaHat =  math.sin ( (ev.genTauTheta_plus-ev.genTauTheta_minus)/2 )/math.sin ( (ev.genTauTheta_minus+ev.genTauTheta_plus)/2 ) 

    cosThetaMeson=math.cos(ev.genMesonTheta_minus)

    Z=cosThetaHat

    z =  Z # cosThetaHat # cosTheta
    absz=abs(z)
    fz=absz/(1+absz*absz)
    signZ = +1 if z > 0 else -1

#    New_AeM1=-1
#    New_Ae=1

    PtauSM= - (AtauSM * (1+  Z*Z) + 2*AeSM*Z) / (1+Z*Z + 2*AeSM*AtauSM*Z)
    PtauP1 = - ( New_Atau   * (1+  Z*Z) + 2*New_Ae*Z) / (1+Z*Z + 2*New_Ae* New_Atau *Z)
    PtauM1 = - ( New_AtauM1   * (1+  Z*Z) + 2*New_AeM1*Z) / (1+Z*Z + 2*New_AeM1* New_AtauM1 *Z)

    angularSM = (3/8*(1 + Z*Z) + AFBSM * Z)
    angularP1 = (3/8*(1 + Z*Z) + AFBP1 * Z)
    angularM1 = (3/8*(1 + Z*Z) + AFBM1 * Z)

    decaySM = (1 + PtauSM * (ev.gen_w_plus + ev.gen_w_minus) + ev.gen_w_plus * ev.gen_w_minus)
    decayP1 = (1 + PtauP1 * (ev.gen_w_plus + ev.gen_w_minus) + ev.gen_w_plus * ev.gen_w_minus)
    decayM1 = (1 + PtauM1 * (ev.gen_w_plus + ev.gen_w_minus) + ev.gen_w_plus * ev.gen_w_minus)

    denW         = angularSM * decaySM
    weightTestP1 = (angularP1 * decayP1) / denW
    weightTestM1 = (angularM1 * decayM1) / denW 

    wP1=weightTestP1      
    wM1=weightTestM1

    tau_plus = build_p4(
        ev.genTauP_plus,
        ev.genTauTheta_plus,
        ev.genTauPhi_plus,
        ev.genTauM_plus,
    )
    
    tau_minus = build_p4(
        ev.genTauP_minus,
        ev.genTauTheta_minus,
        ev.genTauPhi_minus,
        ev.genTauM_minus,
    )
    
    meson_plus = build_p4(
        ev.genMesonP_plus,
        ev.genMesonTheta_plus,
        ev.genMesonPhi_plus,
        ev.genMesonM_plus,
    )
    
    meson_minus = build_p4(
        ev.genMesonP_minus,
        ev.genMesonTheta_minus,
        ev.genMesonPhi_minus,
        ev.genMesonM_minus,
    )

    Z_total = tau_plus + tau_minus
    Z_visible = meson_plus + meson_minus
    dRTaus=dRAngle(tau_plus,tau_minus)

    XCorr=2*abs(math.sin(ev.genTauTheta_plus+ev.genTauTheta_minus))/(math.sin(ev.genTauTheta_plus)+math.sin(ev.genTauTheta_minus)+abs(math.sin(ev.genTauTheta_plus+ev.genTauTheta_minus)))

    sqrtCorr=math.sqrt(91.188*91.188*(1-XCorr))

    h[f"GENZTheta"].Fill(Z_total.Theta())
    h[f"GENZCosTheta"].Fill(math.cos(Z_total.Theta()))    
    h[f"GENZPhi"].Fill(Z_total.Phi())
    h[f"GENZP"].Fill(Z_total.P())
    h[f"GENZMass"].Fill(Z_total.M())
    h[f"GENZVisMass"].Fill(Z_visible.M())
    h[f"GENZVisP"].Fill(Z_visible.P())
    h[f"GENTausDR"].Fill(dRTaus)

    h[f"GENXCorr"].Fill(XCorr)
    h[f"GENSqrtCorr"].Fill(sqrtCorr)
    h[f"GENZMassVsSqrtCorr"].Fill(Z_total.M(),sqrtCorr)

    h[f"GENCosThetaHat"].Fill(cosThetaHat)
    h[f"GENCosThetaSimple"].Fill(math.cos(ev.genTauTheta_minus))
    h[f"GENCosThetaHatVsSimple"].Fill(cosThetaHat,math.cos(ev.genTauTheta_minus))


#    if Z_total.M()<91.1 or Z_total.M()>91.2:
#    if Z_total.M()<89. or Z_total.M()>90.5:
#         continue

#    if abs(Z_total.M()-91.188)<0.050:
#        continue

#    if sqrtCorr>91:
#        continue


#    h[f"GENNPhoton"].Fill(ev.nPhotons)
#    photonSumE=0
#
#    photonSUMP4 = build_p4(0,0,0,0)
#    leadSUMP4 = build_p4(0,0,0,0)
#
#
#    if ev.nPhotons>0:
#       photon = build_p4(ev.genPhotonP, ev.genPhotonTheta,  ev.genPhotonPhi,0)
#
#       #if photon.E()<0.1:
#       # continue
#   
#       photonSUMP4+=photon
#     
#       if leadSUMP4.E()<photon.E():
#          leadSUMP4=photon
#
#       dRZPhoton = dRAngle(photon,Z_total)
#       dRPhotonTauP=dRAngle(photon,tau_plus)
#       dRPhotonTauM=dRAngle(photon,tau_minus)
#       h[f"GENPhotonP"].Fill(ev.genPhotonP)
#       h[f"GENPhotonMom"].Fill(ev.genPhotonMom)
#       if (abs(ev.genPhotonMom)==11):
#           h[f"GENPhotonPISR"].Fill(ev.genPhotonP)   
#           h[f"GENPhotonDRISR"].Fill(dRZPhoton)
#           h[f"GENPhotonPDRISR"].Fill(ev.genPhotonP,dRZPhoton)
#       elif (abs(ev.genPhotonMom)==15):
#           h[f"GENPhotonDRFSR"].Fill(dRZPhoton)  
#           h[f"GENPhotonPFSR"].Fill(ev.genPhotonP)   
#           h[f"GENPhotonPDRFSR"].Fill(ev.genPhotonP,dRZPhoton)
#       h[f"GENPhotonDRTau1"].Fill(dRZPhoton)
#       h[f"GENPhotonDRTau2"].Fill(dRZPhoton)
#       h[f"GENPhotonDR"].Fill(dRZPhoton)
#       photonSumE+=photon.E()
#
##    if photonSumE> 5:
##       continue 
#
#    h[f"GENSumPhotonE"].Fill(photonSumE)
#
#    totalP4 = ROOT.TLorentzVector() 
#    totalP4.SetPxPyPzE(0.0, 0.0, 0.0, 91.188)    
#    recoil = totalP4-photonSUMP4
#
#    mtautauRecoil=recoil.M()
#
#    mtautauRecoil2=math.sqrt( max(0,91.188*91.188-2*91.188*leadSUMP4.E()) )
#
#    h[f"GENZRecoil"].Fill(mtautauRecoil)
#    h[f"GENZRecoil2"].Fill(mtautauRecoil2)

   
    #print (Z_total.M(),sqrtCorr,mtautauRecoil,mtautauRecoil2,leadSUMP4.E(),photonSUMP4.E())


#    if mtautauRecoil<91.1 or mtautauRecoil>91.2:
#       continue 


    spinWT=ev.SpinWT
    spinWThelApprox=ev.SpinWThelApprox
    h[f"SpinWT"].Fill(spinWT)
    h[f"SpinWThelApprox"].Fill(spinWThelApprox)
    h[f"SpinWTvsSpinWThelApprox"].Fill(spinWT, spinWThelApprox)

    for suffix, charge in [("minus", "TAUMINUS"), ("plus",  "TAUPLUS")]:
        x = getattr(ev, f"genMesonE_{suffix}") / getattr(ev, f"genTauE_{suffix}")
        hel = getattr(ev, f"genTauHel_{suffix}")
        w = getattr(ev, f"gen_w_{suffix}")
        cosTheta = getattr(ev, f"gen_cos_theta_tau_{suffix}")
        cosThetaStar = getattr(ev, f"gen_cos_thetastar_{suffix}")
        dm = getattr(ev, f"genTauID_{suffix}")
        taucharge= getattr(ev, f"genTauQ_{suffix}") # we dont really need this, but just to be sure we are filling the right histograms
        mesonmass = getattr(ev, f"genMesonM_{suffix}")
        tauP=getattr(ev, f"genTauP_{suffix}")
        tauTheta=getattr(ev, f"genTauTheta_{suffix}")
        tauPhi=getattr(ev, f"genTauPhi_{suffix}")

 


        h[f"GENMesonEOverBeamE_{charge}"].Fill(x)
        h[f"GENOmega_{charge}"].Fill(w)
        h[f"GENCosTheta_{charge}"].Fill(cosTheta)
        h[f"GENCosThetaHat_{charge}"].Fill(cosThetaHat)        
        h[f"GENCosThetaStar_{charge}"].Fill(cosThetaStar)
        h[f"GENMesonEOverBeamE_CosTheta_{charge}"].Fill(x, cosTheta)
        h[f"GENMesonMass_{charge}"].Fill(mesonmass)
        h[f"GENTauP_{charge}"].Fill(tauP)
        h[f"GENTauTheta_{charge}"].Fill(tauTheta)
        h[f"GENTauPhi_{charge}"].Fill(tauPhi)
        h[f"GENTauDecayMode_{charge}"].Fill(dm)

        h[f"GENAbsZ_{charge}"].Fill(absz)
        h[f"GENFZ_{charge}"].Fill(fz)
        h[f"GENFZFine_{charge}"].Fill(fz)

        h[f"GENOmegaCosTheta_{charge}"].Fill(w,cosTheta)
        h[f"GENOmegaCosThetaHat_{charge}"].Fill(w,cosThetaHat)
        h[f"GENOmegaCosThetaMeson_{charge}"].Fill(w,cosThetaMeson)

        if (z>0):
          h[f"GENFZ_FW_{charge}"].Fill(fz)
          h[f"GENAbsZ_FW_{charge}"].Fill(absz) 
          h[f"GENFZFine_FW_{charge}"].Fill(fz)

        else:
          h[f"GENFZ_BW_{charge}"].Fill(fz)
          h[f"GENAbsZ_BW_{charge}"].Fill(absz)
          h[f"GENFZFine_BW_{charge}"].Fill(fz)

        if taucharge < 0:
            if (hel < 0):
               helTag = "NegHel"
            else:
               helTag = "PosHel"
        if taucharge > 0:
            if (hel > 0):
               helTag = "PosHel"
            else:
               helTag = "NegHel"

        h[f"GENMesonEOverBeamE_{helTag}_{charge}"].Fill(x)
        h[f"GENOmega_{helTag}_{charge}"].Fill(w)
        h[f"GENCosTheta_{helTag}_{charge}"].Fill(cosTheta)
        h[f"GENCosThetaHat_{helTag}_{charge}"].Fill(cosThetaHat)

        h[f"GENOmegaCosTheta_{helTag}_{charge}"].Fill(w,cosTheta)
        h[f"GENOmegaCosThetaHat_{helTag}_{charge}"].Fill(w,cosThetaHat)
        h[f"GENOmegaCosThetaMeson_{helTag}_{charge}"].Fill(w,cosThetaMeson)

        h[f"GENCosThetaStar_{helTag}_{charge}"].Fill(cosThetaStar)
        h[f"GENMesonEOverBeamE_CosTheta_{helTag}_{charge}"].Fill(x, cosTheta)
        h[f"GENMesonMass_{helTag}_{charge}"].Fill(mesonmass)

        h[f"GENAbsZ_{helTag}_{charge}"].Fill(absz)
        h[f"GENFZ_{helTag}_{charge}"].Fill(fz)
        h[f"GENFZFine_{helTag}_{charge}"].Fill(fz)

        if (z>0):
          h[f"GENFZ_FW_{helTag}_{charge}"].Fill(fz)
          h[f"GENAbsZ_FW_{helTag}_{charge}"].Fill(absz)
          h[f"GENFZFine_FW_{helTag}_{charge}"].Fill(fz)

        else:
          h[f"GENFZ_BW_{helTag}_{charge}"].Fill(fz)
          h[f"GENAbsZ_BW_{helTag}_{charge}"].Fill(absz)
          h[f"GENFZFine_BW_{helTag}_{charge}"].Fill(fz)


        h[f"GENTauP_{helTag}_{charge}"].Fill(tauP)
        h[f"GENTauTheta_{helTag}_{charge}"].Fill(tauTheta)
        h[f"GENTauPhi_{helTag}_{charge}"].Fill(tauPhi)
        h[f"GENTauDecayMode_{helTag}_{charge}"].Fill(dm)

        h[f"GENMesonEOverBeamE_P1_{charge}"].Fill(x, wP1)
        h[f"GENMesonEOverBeamE_M1_{charge}"].Fill(x, wM1)

        h[f"GENOmega_P1_{charge}"].Fill(w, wP1)
        h[f"GENOmega_M1_{charge}"].Fill(w, wM1)

        h[f"GENCosTheta_P1_{charge}"].Fill(cosTheta, wP1)
        h[f"GENCosTheta_M1_{charge}"].Fill(cosTheta, wM1)

        h[f"GENCosThetaHat_P1_{charge}"].Fill(cosThetaHat, wP1)
        h[f"GENCosThetaHat_M1_{charge}"].Fill(cosThetaHat, wM1)

        h[f"GENOmegaCosTheta_P1_{charge}"].Fill(w,cosTheta,wP1)
        h[f"GENOmegaCosTheta_M1_{charge}"].Fill(w,cosTheta,wM1)

        h[f"GENOmegaCosThetaMeson_P1_{charge}"].Fill(w,cosThetaMeson,wP1)
        h[f"GENOmegaCosThetaMeson_M1_{charge}"].Fill(w,cosThetaMeson,wM1)

        h[f"GENOmegaCosThetaHat_P1_{charge}"].Fill(w,cosThetaHat,wP1)
        h[f"GENOmegaCosThetaHat_M1_{charge}"].Fill(w,cosThetaHat,wM1)

        h[f"GENCosThetaStar_P1_{charge}"].Fill(cosThetaStar, wP1)
        h[f"GENCosThetaStar_M1_{charge}"].Fill(cosThetaStar, wM1)

        h[f"GENMesonEOverBeamE_CosTheta_P1_{charge}"].Fill(x, cosTheta, wP1)
        h[f"GENMesonEOverBeamE_CosTheta_M1_{charge}"].Fill(x, cosTheta, wM1)

        h[f"GENMesonMass_P1_{charge}"].Fill(mesonmass, wP1)
        h[f"GENMesonMass_M1_{charge}"].Fill(mesonmass, wM1)

        h[f"GENTauP_P1_{charge}"].Fill(tauP, wP1)
        h[f"GENTauP_M1_{charge}"].Fill(tauP, wM1)
        h[f"GENTauTheta_P1_{charge}"].Fill(tauTheta, wP1)
        h[f"GENTauTheta_M1_{charge}"].Fill(tauTheta, wM1)
        h[f"GENTauPhi_P1_{charge}"].Fill(tauPhi, wP1)
        h[f"GENTauPhi_M1_{charge}"].Fill(tauPhi, wM1)
        h[f"GENTauDecayMode_P1_{charge}"].Fill(dm, wP1)
        h[f"GENTauDecayMode_M1_{charge}"].Fill(dm, wM1)

        h[f"GENAbsZ_P1_{charge}"].Fill(absz,wP1)
        h[f"GENFZ_P1_{charge}"].Fill(fz,wP1)
        h[f"GENAbsZ_M1_{charge}"].Fill(absz,wM1)
        h[f"GENFZ_M1_{charge}"].Fill(fz,wM1)

        h[f"GENFZFine_P1_{charge}"].Fill(fz,wP1)
        h[f"GENFZFine_M1_{charge}"].Fill(fz,wM1)


        if (z>0):
          h[f"GENFZ_FW_P1_{charge}"].Fill(fz,wP1)
          h[f"GENFZ_FW_M1_{charge}"].Fill(fz,wM1)
          h[f"GENAbsZ_FW_P1_{charge}"].Fill(absz,wP1)
          h[f"GENAbsZ_FW_M1_{charge}"].Fill(absz,wM1)
          h[f"GENFZFine_FW_P1_{charge}"].Fill(fz,wP1)
          h[f"GENFZFine_FW_M1_{charge}"].Fill(fz,wM1)
        else:
          h[f"GENFZ_BW_P1_{charge}"].Fill(fz,wP1)
          h[f"GENFZ_BW_M1_{charge}"].Fill(fz,wM1)
          h[f"GENAbsZ_BW_P1_{charge}"].Fill(absz,wP1)
          h[f"GENAbsZ_BW_M1_{charge}"].Fill(absz,wM1)
          h[f"GENFZFine_BW_P1_{charge}"].Fill(fz,wP1)
          h[f"GENFZFine_BW_M1_{charge}"].Fill(fz,wM1)


    selectedEvents += 1

outfile.cd()

print (totalEvents)

h["hEvents"].Fill(0, totalEvents)
h["hEvents"].Fill(1, selectedEvents)

for hist in h.values():
    hist.Write()

outfile.Close()


