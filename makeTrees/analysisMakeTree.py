import sys, os, math
from array import array
import ROOT
from ROOT import TFile, TTree, TH1F, TH2F
import numpy as np
from podio import root_io
from pathlib import Path
import ctypes

#import edm4hep


import reco.tauReco
import reco.weightsPol
import reco.optimalVariabRho
import reco.optimalobservable2

import argparse
parser = argparse.ArgumentParser(
    description="Configure the analysis",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("-f", "--inputfile", default="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/KKMCee_ONLYSIM_10M_90GeV/out_sim_edm4hep_1.root") #Py8_ONLYSIM_2M_NOISRFSR")

parser.add_argument("-s", "--sample", default="KKMCee_ONLYSIM_10M_90GeV") #Py8_ONLYSIM_2M_NOISRFSR")
parser.add_argument("-fn", "--filename", default="out_sim_edm4hep") # can also be mini_sim_edm4hep_withhel")

parser.add_argument("-o", "--outfile", default="test_GEN_tree.root")
parser.add_argument("-hl", "--hel", default="False")

parser.add_argument("-w", "--Winter23", default="False")

args = parser.parse_args()
config = vars(args)
print(config)

fileOutName = args.outfile

saveHel=False 
if args.hel=="True":
   saveHel = True

winter23=False 
if args.Winter23=="True":
   winter23=True 
# bool(args.Winter23) #This does not work?

print (args.hel,saveHel)
print (args.Winter23,winter23)


genparts = "a"

if winter23:
   import cppyy
   edm4hep = cppyy.gbl.edm4hep
   genparts = "Particle"
elif not winter23:
   import edm4hep
   genparts = "MCParticles"
 

# ----------------------------------------------------------------------
# Input files
# ----------------------------------------------------------------------

file=args.inputfile
sample = args.sample
filetemplate=args.filename

filenames = []

if file!="":

   filenames.append(file)
   print ("Reading file", file)

else:

   path = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
   dir_path = path + "/" + sample
   filetemplate=args.filename
   badfiles = [-1]
   
   print(dir_path)
   
   nfiles = 100
   
   for i in range(1, nfiles + 1):
   
       if i in badfiles:
           continue
   
       filename = dir_path + "/" + filetemplate + "_{}.root".format(i)
   #    print(filename)
   
       my_file = Path(filename)
   
       if my_file.is_file():
           root_file = reco.tauReco.open_root_file(filename)
   
           if not root_file or root_file.IsZombie():
               continue
   
           filenames.append(filename)
   
   print("Read %d files" % len(filenames))

reader = root_io.Reader(filenames)

# ----------------------------------------------------------------------
# Output tree
# ----------------------------------------------------------------------

treeName = "outtree"

tau_vars = [
    "genTauID",
    "genTauP",
    "genMesonP",
    "genPionP",
    "genTauE",
    "genTauM",
    "genMesonE",
    "genMesonM",
    "genPionE",
    "genPionM",
    "genTauTheta",
    "genTauPhi",
    "genMesonTheta",
    "genMesonPhi",
    "genPionTheta",
    "genPionPhi",
    "gen_cos_thetastar",
    "gen_cos_psi",
    "gen_cos_beta",
    "gen_w",
    "weight_P1",
    "weight_M1",
    "weight_P1_simple",
    "weight_M1_simple",
    "gen_cos_theta_tau",
    "genTauQ",
    "genTauHel",
]

variabs = []

for var in tau_vars:
    variabs.append(var + "_plus")
    variabs.append(var + "_minus")

variabs += [
    "GenZMass",
    "GenZVisMass",
    "beamE",
    "SpinWT",
    "SpinWThelApprox",
    "nPhotons",
    "genPhotonP","genPhotonTheta","genPhotonPhi","genPhotonMom","genPhotonStatus", "weight_P1_corrangle","weight_M1_corrangle",
    "gen_cos_theta_hat","gen_cos_theta_tau"
]

outfile = ROOT.TFile(fileOutName, "RECREATE")
new_tree = ROOT.TTree(treeName, "processed variables")

branches = {}

for var in variabs:
    branches[var] = array("f", [-999.])
    new_tree.Branch(var, branches[var], "%s/F" % var)


def reset_branches(branches, default=-999.):
    for b in branches.values():
        b[0] = default


def safe_float(x, default=-999.):
    try:
        if x is None:
            return default
        if math.isnan(float(x)):
            return default
        return float(x)
    except Exception:
        return default


def make_pion_p4(genPion):
    p4 = ROOT.TLorentzVector()
    p4.SetXYZM(
        genPion.getMomentum().x,
        genPion.getMomentum().y,
        genPion.getMomentum().z,
        genPion.getMass()
    )
    return p4


def compute_tau_dict(
    genTauID,
    genTauQ,
    genTauHelicity,
    genTauP4,
    genMesonP4,
    genPionP4,
    beamE
):
    
    weight_P1_simple = reco.weightsPol.newAtau(
        genTauP4,
        genMesonP4,
        genTauID,
        +1
    )

    weight_M1_simple = reco.weightsPol.newAtau(
        genTauP4,
        genMesonP4,
        genTauID,
        -1
    )
    
    # this only makes sense for the rho decay mode, running for all for tests 
    if genTauID==1:
        (
        gen_cos_thetastar,
        gen_cos_psi,
        gen_cos_beta,
        gen_w,
        weight_P1,
        weight_M1,
        ) = reco.optimalVariabRho.wVariab(
        genTauP4,
        genMesonP4,
        genPionP4,
        beamE
        )
    else: 
       (
        gen_cos_thetastar,
        gen_cos_psi,
        gen_cos_beta,
        gen_w,
        weight_P1,
        weight_M1,
        )= (0,0,0,genMesonP4.E()/genTauP4.E(),weight_P1_simple,weight_M1_simple) 

    # Alternative, if needed:
    # (
    #     gen_cos_thetastar,
    #     gen_cos_psi,
    #     gen_cos_beta,
    #     gen_w,
    #     weight_P1,
    #     weight_M1,
    # ) = optimalobservable2.wVariab(
    #     genTauP4,
    #     genMesonP4,
    #     genPionP4,
    #     beamE
    # )



    gen_cos_theta_tau = math.cos(genTauP4.Theta())

    tau_dict = {
        "genTauID": genTauID,
        "genTauP": genTauP4.P(),
        "genMesonP": genMesonP4.P(),
        "genPionP": genPionP4.P(),

        "genTauE": genTauP4.E(),
        "genTauM": genTauP4.M(),
        "genMesonE": genMesonP4.E(),
        "genMesonM": genMesonP4.M(),
        "genPionE": genPionP4.E(),
        "genPionM": genPionP4.M(),

        "genTauTheta": genTauP4.Theta(),
        "genTauPhi": genTauP4.Phi(),
        "genMesonTheta": genMesonP4.Theta(),
        "genMesonPhi": genMesonP4.Phi(),
        "genPionTheta": genPionP4.Theta(),
        "genPionPhi": genPionP4.Phi(),

        "gen_cos_thetastar": gen_cos_thetastar,
        "gen_cos_psi": gen_cos_psi,
        "gen_cos_beta": gen_cos_beta,
        "gen_w": gen_w,

        "weight_P1": weight_P1,
        "weight_M1": weight_M1,
        "weight_P1_simple": weight_P1_simple,
        "weight_M1_simple": weight_M1_simple,

        "gen_cos_theta_tau": gen_cos_theta_tau,

        "genTauQ": genTauQ,
        "genTauHel": genTauHelicity,

    }

    return tau_dict


def fill_tau_branches(branches, tau_dict):
    q = tau_dict["genTauQ"]

    if q > 0:
        suffix = "plus"
    elif q < 0:
        suffix = "minus"
    else:
        print ("this makes no sense")
        return

    for var, value in tau_dict.items():
        branch_name = "%s_%s" % (var, suffix)

        if branch_name in branches:
            branches[branch_name][0] = safe_float(value)



# ----------------------------------------------------------------------
# Event loop
# ----------------------------------------------------------------------

totalEvents = 0
selectedEvents = 0

for event in reader.get("events"):

    if totalEvents % 10000 == 0:
        print(totalEvents)

    #print (totalEvents)

    totalEvents += 1

    reset_branches(branches)

    if saveHel==True: 
       SpinWT = event.get_parameter("SpinWT")
       SpinWThelApprox = event.get_parameter("SpinWThelApprox")
    else:
       SpinWT=-99
       SpinWThelApprox=-99

   # print (SpinWT)

    mc_particles = event.get(genparts)

    if len(mc_particles) < 2:
        continue

    beamE = mc_particles[0].getEnergy()


    countPhotons=0
    highestPho=-1 
    phoSaveP4=ROOT.TLorentzVector()
    phoMom=0
    phoStatusGen=-1
    for part in mc_particles:
        if part.getPDG()==22 and part.getGeneratorStatus() != 0:
            phoP4=ROOT.TLorentzVector()
            phoP4.SetXYZM(part.getMomentum().x,part.getMomentum().y,part.getMomentum().z,part.getMass())
            mothers = part.getParents()
            motherPDG = -999
            if len(mothers) > 0:
                  motherPDG = mothers[0].getPDG()
            if phoP4.P()>1:      
               #print (part.getPDG(), part.getGeneratorStatus(),part.getSimulatorStatus(),phoP4.P(),motherPDG,part.getMass())
               if phoP4.P()>highestPho:
                   highestPho=phoP4.P()
                   phoSaveP4=phoP4
                   phoMom=motherPDG
                   phoStatusGen=part.getGeneratorStatus()
               countPhotons+=1


    # GEN level taus
    genTaus = reco.tauReco.findAllGenTaus(mc_particles,saveHel)
    nGenTaus = len(genTaus)

    if nGenTaus < 2:
        continue




    # Optional decay-mode selection.
    #if selectDecay >= -999:
    #    if genTaus[0][1] != selectDecay:
    #        continue
    #    if genTaus[1][1] != selectDecay:
    #        continue

    try:
        ZGen = genTaus[0][3] + genTaus[1][3]
        ZVisGen = genTaus[0][0] + genTaus[1][0]

        GenZMass = ZGen.M()
        GenZVisMass = ZVisGen.M()

        # Tau 0
        genMesonP4_0 = genTaus[0][0]
        genTauID_0 = genTaus[0][1]
        genTauQ_0 = genTaus[0][2]
        genTauP4_0 = genTaus[0][3]
        genTauConst_0 = genTaus[0][6]
        if saveHel==True:
           genTauHelicity_0 = genTaus[0][7]
        else:
           genTauHelicity_0 = 0


        # Tau 1
        genMesonP4_1 = genTaus[1][0]
        genTauID_1 = genTaus[1][1]
        genTauQ_1 = genTaus[1][2]
        genTauP4_1 = genTaus[1][3]
        genTauConst_1 = genTaus[1][6]
        if saveHel==True:
           genTauHelicity_1 = genTaus[1][7]
        else:
           genTauHelicity_1 = 0
 
        if len(genTauConst_0) < 1:
            continue

        if len(genTauConst_1) < 1:
            continue

        genPion_0 = genTauConst_0[0]
        genPion_1 = genTauConst_1[0]

        if( (genTauID_0==1 and abs(genPion_0.getPDG())!=211)  ): 
             print ("Caution, constituents not ordered!! 0 is not a pion: ", genPion_0.getPDG())

        if( (genTauID_1==1 and abs(genPion_1.getPDG())!=211) ):
             print ("Caution, constituents not ordered!! 0 is not a pion: ", genPion_1.getPDG())


        genPionP4_0 = make_pion_p4(genPion_0)
        genPionP4_1 = make_pion_p4(genPion_1)

        tau_dict_0 = compute_tau_dict(
            genTauID_0,
            genTauQ_0,
            genTauHelicity_0,
            genTauP4_0,
            genMesonP4_0,
            genPionP4_0,
            beamE
        )

        tau_dict_1 = compute_tau_dict(
            genTauID_1,
            genTauQ_1,
            genTauHelicity_1,
            genTauP4_1,
            genMesonP4_1,
            genPionP4_1,
            beamE
        )
  
        #print(GenZVisMass)

#        if genTauQ_0==0:
#            print (nGenTaus,genTauQ_0,genTauID_0)
#            for c in range(0,len(genTauConst_0)):
#                print (make_pion_p4(genTauConst_0[c]).P(),genTauConst_0[c].getPDG())
#        if genTauQ_1==0:
#            print (nGenTaus,genTauQ_1,genTauID_1)
#            for c in range(0,len(genTauConst_1)):
#                print (make_pion_p4(genTauConst_1[c]).P(),genTauConst_1[c].getPDG())

        fill_tau_branches(branches, tau_dict_0)
        fill_tau_branches(branches, tau_dict_1)


################
        # Corrected Weight, only for rho 

        AeSM=0.1472 # corresponds to sin2thetaeff=0.2315, for pythia 
        AtauSM=AeSM
        New_Atau=1
        New_Ae= AeSM # New_Atau # AeSM
        New_AtauM1=-1
        New_AeM1= AeSM # New_AtauM1 # AeSM
        AFBSM = 3/4 * AeSM * AtauSM
        AFBP1 = 3/4 * New_Ae  * New_Atau    # = 3/4
        AFBM1 = 3/4 * New_AeM1 * New_AtauM1 # = 3/4 * (-1)*(-1) = 3/4
   
        genTauTheta_minus=tau_dict_1["genTauTheta"]
        genTauTheta_plus=tau_dict_0["genTauTheta"]
        gen_w_minus=tau_dict_1["gen_w"]
        gen_w_plus=tau_dict_0["gen_w"]
        genQ_minus=tau_dict_1["genTauQ"]

        if genQ_minus!=-1:
         genTauTheta_minus=tau_dict_0["genTauTheta"]
         genTauTheta_plus=tau_dict_1["genTauTheta"]
         gen_w_minus=tau_dict_0["gen_w"]
         gen_w_plus=tau_dict_1["gen_w"]
      

        gen_cos_theta_hat =  math.sin ( (genTauTheta_plus-genTauTheta_minus)/2 )/math.sin ( (genTauTheta_minus+genTauTheta_plus)/2 )
        #Z=gen_cos_theta_hat
    
        Z=math.cos(genTauTheta_minus) 

    #    New_AeM1=-1
    #    New_Ae=1
    
        PtauSM= - (AtauSM * (1+  Z*Z) + 2*AeSM*Z) / (1+Z*Z + 2*AeSM*AtauSM*Z)
        PtauP1 = - ( New_Atau   * (1+  Z*Z) + 2*New_Ae*Z) / (1+Z*Z + 2*New_Ae* New_Atau *Z)
        PtauM1 = - ( New_AtauM1   * (1+  Z*Z) + 2*New_AeM1*Z) / (1+Z*Z + 2*New_AeM1* New_AtauM1 *Z)
    
        angularSM = (3/8*(1 + Z*Z) + AFBSM * Z)
        angularP1 = (3/8*(1 + Z*Z) + AFBP1 * Z)
        angularM1 = (3/8*(1 + Z*Z) + AFBM1 * Z)
    
        decaySM = (1 + PtauSM * (gen_w_plus + gen_w_minus) + gen_w_plus * gen_w_minus)
        decayP1 = (1 + PtauP1 * (gen_w_plus + gen_w_minus) + gen_w_plus * gen_w_minus)
        decayM1 = (1 + PtauM1 * (gen_w_plus + gen_w_minus) + gen_w_plus * gen_w_minus)
    
        denW         = angularSM * decaySM
        weightTestP1 = (angularP1 * decayP1) / denW
        weightTestM1 = (angularM1 * decayM1) / denW
    
        weight_P1_corrangle=weightTestP1
        weight_M1_corrangle=weightTestM1
        branches["weight_M1_corrangle"][0]=safe_float(weight_M1_corrangle)
        branches["weight_P1_corrangle"][0]=safe_float(weight_P1_corrangle)

###############################3

        branches["gen_cos_theta_hat"][0]=gen_cos_theta_hat
        branches["gen_cos_theta_tau"][0]=genTauTheta_minus

        branches["GenZMass"][0] = safe_float(GenZMass)
        branches["GenZVisMass"][0] = safe_float(GenZVisMass)
        branches["beamE"][0] = safe_float(beamE)

        branches["SpinWT"][0] = safe_float(SpinWT)
        branches["SpinWThelApprox"][0] = safe_float(SpinWThelApprox)

        branches["nPhotons"][0] = countPhotons
        if countPhotons!=0:
            branches["genPhotonP"][0]= phoSaveP4.P()
            branches["genPhotonTheta"][0]= phoSaveP4.Theta()
            branches["genPhotonPhi"][0]= phoSaveP4.Phi()
            branches["genPhotonMom"][0]=phoMom
            branches["genPhotonStatus"][0]=phoStatusGen
        else:
            branches["genPhotonP"][0]= -1
            branches["genPhotonTheta"][0]=-1
            branches["genPhotonPhi"][0]= -1
            branches["genPhotonMom"][0]= -1

        new_tree.Fill()
        selectedEvents += 1

    except Exception as e:
        print("Skipping event due to error:", e)
        continue


print("Total events:", totalEvents)
print("Selected events:", selectedEvents)

outfile.cd()
new_tree.Write()
outfile.Close()

print("Wrote:", fileOutName)
