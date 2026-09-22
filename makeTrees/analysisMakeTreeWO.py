import sys, os, math
from array import array
import ROOT
from podio import root_io
from pathlib import Path

# Existing reconstruction
import reco.tauReco

# Olmo versions: place these in makeTrees/reco/
#   reco/weightsPolOLMO.py
#   reco/optimalVariabRhoOLMO.py
import reco.weightsPolOLMO as weightsPolWO
import reco.optimalVariabRhoOLMO as optimalVariabRhoWO


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

# Pythia default effective weak mixing angle used for fermion Z couplings.
# Pass this explicitly everywhere because the Olmo modules currently default
# to 0.2312 if sin_eff is omitted.
SIN2THETA_EFF = 0.2315

# For the clean closure test use the true tau- direction.
# Set to True only when explicitly studying the reconstructed/proxy angle.
USE_COSTHETAHAT_FOR_PRODUCTION = False


import argparse
parser = argparse.ArgumentParser(
    description="Configure the analysis (Olmo weight modules)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-f", "--inputfile",
    default="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
            "KKMCee_ONLYSIM_10M_90GeV/out_sim_edm4hep_1.root"
)
parser.add_argument("-s", "--sample", default="KKMCee_ONLYSIM_10M_90GeV")
parser.add_argument("-fn", "--filename", default="out_sim_edm4hep")
parser.add_argument("-o", "--outfile", default="test_GEN_tree_WO.root")
parser.add_argument("-hl", "--hel", default="False")
parser.add_argument("-w", "--Winter23", default="False")

args = parser.parse_args()
print(vars(args))

fileOutName = args.outfile

saveHel = (args.hel == "True")
winter23 = (args.Winter23 == "True")

print(args.hel, saveHel)
print(args.Winter23, winter23)


# ----------------------------------------------------------------------
# EDM4hep collection setup
# ----------------------------------------------------------------------

if winter23:
    import cppyy
    edm4hep = cppyy.gbl.edm4hep
    genparts = "Particle"
else:
    import edm4hep
    genparts = "MCParticles"


# ----------------------------------------------------------------------
# Input files
# ----------------------------------------------------------------------

file = args.inputfile
sample = args.sample
filetemplate = args.filename

filenames = []

if file != "":
    filenames.append(file)
    print("Reading file", file)

else:
    path = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
    dir_path = path + "/" + sample
    badfiles = [-1]

    print(dir_path)

    nfiles = 100

    for i in range(1, nfiles + 1):
        if i in badfiles:
            continue

        filename = dir_path + "/" + filetemplate + "_{}.root".format(i)
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
    "genPhotonP",
    "genPhotonTheta",
    "genPhotonPhi",
    "genPhotonMom",
    "genPhotonStatus",

    # Final full-event weights, analogous to the old corrangle branches.
    "weight_P1_corrangle",
    "weight_M1_corrangle",

    # Additional WO diagnostics:
    # joint tau-decay correlation only
    "weight_P1_joint",
    "weight_M1_joint",

    # production-angle / AFB factor only
    "weight_P1_angular",
    "weight_M1_angular",
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


def tau_pdg_from_charge(q):
    # tau- : PDG +15, charge -1
    # tau+ : PDG -15, charge +1
    if q < 0:
        return 15
    if q > 0:
        return -15
    raise ValueError("Tau with zero charge")


def compute_tau_dict(
    genTauID,
    genTauQ,
    genTauHelicity,
    genTauP4,
    genMesonP4,
    genPionP4,
    beamE
):
    """
    Compute per-tau observables/weights with Olmo's modules.

    weight_P1_simple / weight_M1_simple:
        channel-level single-tau weights from weightsPolOLMO.

    weight_P1 / weight_M1:
        for rho, the full omega-based weight returned by
        optimalVariabRhoOLMO.wVariab.
    """

    tau_pdg = tau_pdg_from_charge(genTauQ)

    # --------------------------------------------------------------
    # "Simple" single-tau weight from Olmo's weightsPol module.
    # --------------------------------------------------------------

    if genTauID in (0, 1, 10):
        weight_P1_simple = weightsPolWO.newAtau(
            genTauP4,
            genMesonP4,
            genTauID,
            +1,
            tau_pdg=tau_pdg,
            sin_eff=SIN2THETA_EFF
        )

        weight_M1_simple = weightsPolWO.newAtau(
            genTauP4,
            genMesonP4,
            genTauID,
            -1,
            tau_pdg=tau_pdg,
            sin_eff=SIN2THETA_EFF
        )

    elif genTauID in (-11, -13):
        # For leptonic tau decays genMesonP4 is the visible e/mu system.
        weight_P1_simple = weightsPolWO.newAtauLep(
            genMesonP4,
            genTauP4,
            beamE,
            +1,
            tau_pdg=tau_pdg,
            sin_eff=SIN2THETA_EFF
        )

        weight_M1_simple = weightsPolWO.newAtauLep(
            genMesonP4,
            genTauP4,
            beamE,
            -1,
            tau_pdg=tau_pdg,
            sin_eff=SIN2THETA_EFF
        )

    else:
        # New diagnostic IDs (eta, K, mixed strange modes, etc.) are kept in
        # the tree but are not assigned a physics weight here.
        weight_P1_simple = 1.0
        weight_M1_simple = 1.0

    # --------------------------------------------------------------
    # Full rho optimal variable.
    # --------------------------------------------------------------

    if genTauID == 1:
        (
            gen_cos_thetastar,
            gen_cos_psi,
            gen_cos_beta,
            gen_w,
            weight_P1,
            weight_M1,
        ) = optimalVariabRhoWO.wVariab(
            genTauP4,
            genMesonP4,
            genPionP4,
            beamE,
            tau_pdg=tau_pdg,
            sin_eff=SIN2THETA_EFF
        )

    else:
        # Preserve the old branch convention for non-rho modes.
        if genTauP4.E() > 0:
            gen_w = genMesonP4.E() / genTauP4.E()
        else:
            gen_w = -999.0

        gen_cos_thetastar = 0.0
        gen_cos_psi = 0.0
        gen_cos_beta = 0.0
        weight_P1 = weight_P1_simple
        weight_M1 = weight_M1_simple

    gen_cos_theta_tau = math.cos(genTauP4.Theta())

    return {
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


def fill_tau_branches(branches, tau_dict):
    q = tau_dict["genTauQ"]

    if q > 0:
        suffix = "plus"
    elif q < 0:
        suffix = "minus"
    else:
        print("Tau with zero charge; cannot assign plus/minus branch")
        return

    for var, value in tau_dict.items():
        branch_name = "%s_%s" % (var, suffix)

        if branch_name in branches:
            branches[branch_name][0] = safe_float(value)


def get_minus_plus(
    genTauP4_0, genTauP4_1,
    tau_dict_0, tau_dict_1
):
    """Return p4/dict ordered as tau-, tau+."""
    if tau_dict_0["genTauQ"] < 0:
        return genTauP4_0, tau_dict_0, genTauP4_1, tau_dict_1

    if tau_dict_1["genTauQ"] < 0:
        return genTauP4_1, tau_dict_1, genTauP4_0, tau_dict_0

    raise ValueError("Could not identify tau- in the event")


def compute_full_event_weights(
    tau_minus_p4,
    tau_minus_dict,
    tau_plus_p4,
    tau_plus_dict
):
    """
    Full rho-rho event reweighting using Olmo's two-tau joint weight.

    The Olmo joint function contains the correlated tau-decay factor:
      [1 + P_new (omega- + omega+) + omega- omega+]
      ------------------------------------------------
      [1 + P_SM  (omega- + omega+) + omega- omega+]

    The production-angle / AFB factor is kept separately, as in the original
    analysisMakeTree.py, and multiplied afterwards.

    Returns:
      full_P1, full_M1,
      joint_P1, joint_M1,
      angular_P1, angular_M1
    """

    # The joint omega expression is meaningful here for rho-rho.
    if tau_minus_dict["genTauID"] != 1 or tau_plus_dict["genTauID"] != 1:
        return (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    omega_minus = tau_minus_dict["gen_w"]
    omega_plus = tau_plus_dict["gen_w"]

    theta_minus = tau_minus_dict["genTauTheta"]
    theta_plus = tau_plus_dict["genTauTheta"]

    # Olmo's joint two-tau decay weights.
    joint_P1 = weightsPolWO.newAtauJoint(
        tau_minus_p4,
        omega_minus,
        omega_plus,
        +1,
        tau_pdg=15,
        sin_eff=SIN2THETA_EFF
    )

    joint_M1 = weightsPolWO.newAtauJoint(
        tau_minus_p4,
        omega_minus,
        omega_plus,
        -1,
        tau_pdg=15,
        sin_eff=SIN2THETA_EFF
    )

    # Production angle.
    cosThetaTauMinus = math.cos(theta_minus)

    denom_hat = math.sin((theta_minus + theta_plus) / 2.0)
    if abs(denom_hat) > 1e-12:
        cosThetaHat = math.sin((theta_plus - theta_minus) / 2.0) / denom_hat
    else:
        cosThetaHat = cosThetaTauMinus

    if USE_COSTHETAHAT_FOR_PRODUCTION:
        Z = cosThetaHat
    else:
        Z = cosThetaTauMinus

    # SM point defined by the same sin_eff passed to Olmo's modules.
    AeSM = weightsPolWO._compute_ae_sm(SIN2THETA_EFF)
    AtauSM = AeSM

    New_Atau_P1 = +1.0
    New_Atau_M1 = -1.0

    # Ae is held fixed while Atau is reweighted.
    AFB_SM = 3.0 / 4.0 * AeSM * AtauSM
    AFB_P1 = 3.0 / 4.0 * AeSM * New_Atau_P1
    AFB_M1 = 3.0 / 4.0 * AeSM * New_Atau_M1

    angular_SM = 3.0 / 8.0 * (1.0 + Z * Z) + AFB_SM * Z
    angular_P1_num = 3.0 / 8.0 * (1.0 + Z * Z) + AFB_P1 * Z
    angular_M1_num = 3.0 / 8.0 * (1.0 + Z * Z) + AFB_M1 * Z

    if abs(angular_SM) < 1e-12:
        angular_P1 = 1.0
        angular_M1 = 1.0
    else:
        angular_P1 = angular_P1_num / angular_SM
        angular_M1 = angular_M1_num / angular_SM

    full_P1 = angular_P1 * joint_P1
    full_M1 = angular_M1 * joint_M1

    return (
        full_P1,
        full_M1,
        joint_P1,
        joint_M1,
        angular_P1,
        angular_M1,
    )


# ----------------------------------------------------------------------
# Event loop
# ----------------------------------------------------------------------

totalEvents = 0
selectedEvents = 0

for event in reader.get("events"):

    if totalEvents % 10000 == 0:
        print(totalEvents)

    totalEvents += 1

    reset_branches(branches)

    if saveHel:
        SpinWT = event.get_parameter("SpinWT")
        SpinWThelApprox = event.get_parameter("SpinWThelApprox")
    else:
        SpinWT = -99
        SpinWThelApprox = -99

    mc_particles = event.get(genparts)

    if len(mc_particles) < 2:
        continue

    beamE = mc_particles[0].getEnergy()

    # --------------------------------------------------------------
    # Photon information
    # --------------------------------------------------------------

    countPhotons = 0
    highestPho = -1
    phoSaveP4 = ROOT.TLorentzVector()
    phoMom = 0
    phoStatusGen = -1

    for part in mc_particles:
        if part.getPDG() == 22 and part.getGeneratorStatus() != 0:
            phoP4 = ROOT.TLorentzVector()
            phoP4.SetXYZM(
                part.getMomentum().x,
                part.getMomentum().y,
                part.getMomentum().z,
                part.getMass()
            )

            mothers = part.getParents()
            motherPDG = -999
            if len(mothers) > 0:
                motherPDG = mothers[0].getPDG()

            if phoP4.P() > 1:
                if phoP4.P() > highestPho:
                    highestPho = phoP4.P()
                    phoSaveP4 = phoP4
                    phoMom = motherPDG
                    phoStatusGen = part.getGeneratorStatus()

                countPhotons += 1

    # --------------------------------------------------------------
    # GEN taus
    # --------------------------------------------------------------

    genTaus = reco.tauReco.findAllGenTaus(mc_particles, saveHel)
    nGenTaus = len(genTaus)

    if nGenTaus < 2:
        continue

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

        if saveHel:
            genTauHelicity_0 = genTaus[0][7]
        else:
            genTauHelicity_0 = 0

        # Tau 1
        genMesonP4_1 = genTaus[1][0]
        genTauID_1 = genTaus[1][1]
        genTauQ_1 = genTaus[1][2]
        genTauP4_1 = genTaus[1][3]
        genTauConst_1 = genTaus[1][6]

        if saveHel:
            genTauHelicity_1 = genTaus[1][7]
        else:
            genTauHelicity_1 = 0

        if len(genTauConst_0) < 1:
            continue
        if len(genTauConst_1) < 1:
            continue

        genPion_0 = genTauConst_0[0]
        genPion_1 = genTauConst_1[0]

        # Diagnostic only for events actually classified as rho.
        if genTauID_0 == 1 and abs(genPion_0.getPDG()) != 211:
            print(
                "Caution: tau 0 classified as rho but constituent 0 is",
                genPion_0.getPDG()
            )

        if genTauID_1 == 1 and abs(genPion_1.getPDG()) != 211:
            print(
                "Caution: tau 1 classified as rho but constituent 0 is",
                genPion_1.getPDG()
            )

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

        fill_tau_branches(branches, tau_dict_0)
        fill_tau_branches(branches, tau_dict_1)

        # ----------------------------------------------------------
        # WO full-event rho-rho weights
        # ----------------------------------------------------------

        (
            tau_minus_p4,
            tau_minus_dict,
            tau_plus_p4,
            tau_plus_dict
        ) = get_minus_plus(
            genTauP4_0,
            genTauP4_1,
            tau_dict_0,
            tau_dict_1
        )

        (
            weight_P1_corrangle,
            weight_M1_corrangle,
            weight_P1_joint,
            weight_M1_joint,
            weight_P1_angular,
            weight_M1_angular,
        ) = compute_full_event_weights(
            tau_minus_p4,
            tau_minus_dict,
            tau_plus_p4,
            tau_plus_dict
        )

        branches["weight_P1_corrangle"][0] = safe_float(weight_P1_corrangle)
        branches["weight_M1_corrangle"][0] = safe_float(weight_M1_corrangle)

        branches["weight_P1_joint"][0] = safe_float(weight_P1_joint)
        branches["weight_M1_joint"][0] = safe_float(weight_M1_joint)

        branches["weight_P1_angular"][0] = safe_float(weight_P1_angular)
        branches["weight_M1_angular"][0] = safe_float(weight_M1_angular)

        # ----------------------------------------------------------
        # Event-level branches
        # ----------------------------------------------------------

        branches["GenZMass"][0] = safe_float(GenZMass)
        branches["GenZVisMass"][0] = safe_float(GenZVisMass)
        branches["beamE"][0] = safe_float(beamE)

        branches["SpinWT"][0] = safe_float(SpinWT)
        branches["SpinWThelApprox"][0] = safe_float(SpinWThelApprox)

        branches["nPhotons"][0] = countPhotons

        if countPhotons != 0:
            branches["genPhotonP"][0] = phoSaveP4.P()
            branches["genPhotonTheta"][0] = phoSaveP4.Theta()
            branches["genPhotonPhi"][0] = phoSaveP4.Phi()
            branches["genPhotonMom"][0] = phoMom
            branches["genPhotonStatus"][0] = phoStatusGen
        else:
            branches["genPhotonP"][0] = -1
            branches["genPhotonTheta"][0] = -1
            branches["genPhotonPhi"][0] = -1
            branches["genPhotonMom"][0] = -1
            branches["genPhotonStatus"][0] = -1

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

