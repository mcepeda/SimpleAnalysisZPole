import ROOT
import math

def _clip(x, lo=-1.0, hi=1.0):
    return max(lo, min(hi, x))

def wVariab(genTauP4, genRhoP4, genPionP4, beamE, testAtau=0):

    mtau = 1.7769

    # -----------------------------
    # GEN: rho in tau rest frame
    # -----------------------------
    boostGEN = genTauP4.BoostVector()

    genRho_TauRes = ROOT.TLorentzVector()
    genRho_TauRes.SetPxPyPzE(
        genRhoP4.Px(), genRhoP4.Py(), genRhoP4.Pz(), genRhoP4.E()
    )
    genRho_TauRes.Boost(-boostGEN)

    # theta: rho direction in tau rest frame wrt tau lab/CM direction
    v_rho_tau = genRho_TauRes.Vect()
    v_tau_lab = genTauP4.Vect()

    gen_theta_Rho = v_rho_tau.Angle(v_tau_lab)
    z = math.cos(gen_theta_Rho)
    cos_theta = z

    # -----------------------------
    # GEN: charged pion in rho rest frame
    # -----------------------------
    genBoostRho = genRhoP4.BoostVector()

    genPion_RhoRes = ROOT.TLorentzVector()
    genPion_RhoRes.SetPxPyPzE(
        genPionP4.Px(), genPionP4.Py(), genPionP4.Pz(), genPionP4.E()
    )
    genPion_RhoRes.Boost(-genBoostRho)

    v_pi_rho = genPion_RhoRes.Vect()
    v_rho_lab = genRhoP4.Vect()

    gen_beta = v_pi_rho.Angle(v_rho_lab)
    gen_cosBeta = math.cos(gen_beta)

    # -----------------------------
    # Kinematic variables
    # -----------------------------
    x = genRhoP4.E() / genTauP4.E()
    sqrts = 2.0 * beamE
    mRho = genRhoP4.M()

    if mRho <= 0:
        raise ValueError("Bad rho mass")

    r = mtau * mtau / (mRho * mRho)

    arg = x*x - 4.0*mRho*mRho/(sqrts*sqrts)
    if arg <= 0:
        raise ValueError("Bad cosPsi denominator")

    cosPsi = (
        x * (mtau*mtau + mRho*mRho) - 2.0*mRho*mRho
    ) / (
        (mtau*mtau - mRho*mRho) * math.sqrt(arg)
    )

    cosPsi = _clip(cosPsi)
    anglePsi = math.acos(cosPsi)

    Cpsi  = 0.5 * (3.0*cosPsi*cosPsi - 1.0)
    Cbeta = 0.5 * (3.0*gen_cosBeta*gen_cosBeta - 1.0)

    sin_theta = math.sin(gen_theta_Rho)
    sin2psi = math.sin(2.0 * anglePsi)

    # -----------------------------
    # Decay density W(P) = f + P g
    # for tau -> rho nu -> pi pi0 nu
    # -----------------------------
    def W_decay(P):
        term_a = (
            2.0/3.0 * ((1.0 - P*z) - r*(1.0 + P*z))
            + r*(1.0 + P*z)
        )

        term_b = (
            -2.0/3.0
            * (
                ((1.0 - P*z) - r*(1.0 + P*z)) * Cpsi
                - 1.5 * math.sqrt(r) * P * sin2psi * sin_theta
            )
            * Cbeta
        )

        return term_a + term_b

    W_p = W_decay(+1.0)
    W_m = W_decay(-1.0)

    denom = W_p + W_m
    if abs(denom) < 1e-12:
        raise ValueError("Bad omega denominator")

    # Davier optimal observable omega = g/f
    w = (W_p - W_m) / denom

    # -----------------------------
    # SM polarization and Atau reweights
    # -----------------------------
    sin2theta_effective = 0.2315
    gv_ga = 1.0 - 4.0*sin2theta_effective

    Ae_sm = 2.0 * gv_ga / (1.0 + gv_ga*gv_ga)
    Atau_sm = Ae_sm
    Ae = Ae_sm

    costheta_tau = math.cos(genTauP4.Theta())

    def Ptau_from_Atau(Atau):
        return -(
            Atau * (1.0 + costheta_tau*costheta_tau)
            + 2.0 * Ae * costheta_tau
        ) / (
            1.0 + costheta_tau*costheta_tau
            + 2.0 * Ae * Atau * costheta_tau
        )

    Ptau_sm = Ptau_from_Atau(Atau_sm)

    den_sm = W_decay(Ptau_sm)
    if abs(den_sm) < 1e-12:
        raise ValueError("Bad SM reweighting denominator")

    Pnew_P1 = Ptau_from_Atau(+1.0)
    Pnew_M1 = Ptau_from_Atau(-1.0)

    weight_P1 = W_decay(Pnew_P1) / den_sm
    weight_M1 = W_decay(Pnew_M1) / den_sm

    if testAtau == 0:
        return (cos_theta, cosPsi, gen_cosBeta, w, weight_P1, weight_M1)
    else:
        Pnew_test = Ptau_from_Atau(testAtau)
        weight_Ptest = W_decay(Pnew_test) / den_sm
        return (cos_theta, cosPsi, gen_cosBeta, w,
                weight_P1, weight_M1, weight_Ptest)
    


def wVariabDanvier(genTauP4, genRhoP4, genPionP4, beamE, testAtau=0):

    mtau = 1.7769

    def clip(x):
        return max(-1.0, min(1.0, x))

    def makeTLV(p4):
        v = ROOT.TLorentzVector()
        v.SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.E())
        return v

    # ------------------------------------------------------------
    # Reconstruct pi0 and neutrino at generator level
    # tau -> rho nu, rho -> pi_ch pi0
    # ------------------------------------------------------------
    genPi0P4 = genRhoP4 - genPionP4
    genNuP4  = genTauP4 - genRhoP4

    # ------------------------------------------------------------
    # Boost decay products to tau rest frame
    # ------------------------------------------------------------
    boostTau = genTauP4.BoostVector()

    pi_ch_tau = makeTLV(genPionP4)
    pi0_tau   = makeTLV(genPi0P4)
    nu_tau    = makeTLV(genNuP4)
    rho_tau   = makeTLV(genRhoP4)

    pi_ch_tau.Boost(-boostTau)
    pi0_tau.Boost(-boostTau)
    nu_tau.Boost(-boostTau)
    rho_tau.Boost(-boostTau)

    # ------------------------------------------------------------
    # Diagnostic variables, kept compatible with your old function
    # ------------------------------------------------------------
    cos_theta = math.cos(rho_tau.Vect().Angle(genTauP4.Vect()))

    genBoostRho = genRhoP4.BoostVector()
    pion_rho = makeTLV(genPionP4)
    pion_rho.Boost(-genBoostRho)

    gen_cosBeta = math.cos(pion_rho.Vect().Angle(genRhoP4.Vect()))

    x = genRhoP4.E() / genTauP4.E()
    sqrts = 2.0 * beamE
    mRho = genRhoP4.M()

    arg = x*x - 4.0*mRho*mRho/(sqrts*sqrts)

    if arg > 0 and mRho > 0:
        cosPsi = (
            x * (mtau*mtau + mRho*mRho) - 2.0*mRho*mRho
        ) / (
            (mtau*mtau - mRho*mRho) * math.sqrt(arg)
        )
        cosPsi = clip(cosPsi)
    else:
        cosPsi = -999.0

    # ------------------------------------------------------------
    # Davier optimal observable omega = g/f
    #
    # For rho:
    # q = p(pi_ch) - p(pi0)
    # N = p(nu)
    #
    # H^mu = 2 (q.N) q^mu - q^2 N^mu
    #
    # In the tau rest frame, h = H_vec / H_0
    # omega = h · tau_direction
    # ------------------------------------------------------------
    q = pi_ch_tau - pi0_tau
    N = nu_tau

    q_dot_N = q.Dot(N)
    q2 = q.Dot(q)

    H = ROOT.TLorentzVector()
    H.SetPxPyPzE(
        2.0*q_dot_N*q.Px() - q2*N.Px(),
        2.0*q_dot_N*q.Py() - q2*N.Py(),
        2.0*q_dot_N*q.Pz() - q2*N.Pz(),
        2.0*q_dot_N*q.E()  - q2*N.E()
    )

    if abs(H.E()) < 1e-12:
        w = 0.0
    else:
        h_vec = H.Vect() * (1.0 / H.E())

        tau_axis = genTauP4.Vect()
        tau_axis = tau_axis.Unit()

        w = h_vec.Dot(tau_axis)

    # protect against tiny numerical excursions
    w = clip(w)

    # ------------------------------------------------------------
    # Reweights using Davier form:
    #
    # W(P) ∝ 1 + P omega
    #
    # weight(Pnew/Psm) = (1 + Pnew*w) / (1 + Psm*w)
    # ------------------------------------------------------------
    sin2theta_effective = 0.2315
    gv_ga = 1.0 - 4.0*sin2theta_effective

    Ae_sm = 2.0 * gv_ga / (1.0 + gv_ga*gv_ga)
    Atau_sm = Ae_sm
    Ae = Ae_sm

    costheta_tau = math.cos(genTauP4.Theta())

    def Ptau_from_Atau(Atau):
        return -(
            Atau * (1.0 + costheta_tau*costheta_tau)
            + 2.0 * Ae * costheta_tau
        ) / (
            1.0 + costheta_tau*costheta_tau
            + 2.0 * Ae * Atau * costheta_tau
        )

    Ptau_sm = Ptau_from_Atau(Atau_sm)
    Pnew_P1 = Ptau_from_Atau(+1.0)
    Pnew_M1 = Ptau_from_Atau(-1.0)

    den = 1.0 + Ptau_sm*w

    if abs(den) < 1e-12:
        weight_P1 = 1.0
        weight_M1 = 1.0
    else:
        weight_P1 = (1.0 + Pnew_P1*w) / den
        weight_M1 = (1.0 + Pnew_M1*w) / den

    return (cos_theta, cosPsi, gen_cosBeta, w, weight_P1, weight_M1)
