#!/usr/bin/env python3

import argparse
import ctypes
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# ============================================================
# Polarization fit
#
# Expected input histograms:
#
#   histo_TAUMINUS_<bin>      nominal SM distribution
#   histo_P1_TAUMINUS_<bin>   A_tau = +1 -> P_tau = -1
#   histo_M1_TAUMINUS_<bin>   A_tau = -1 -> P_tau = +1
#
# and the corresponding "_full" histograms.
#
# The P1/M1 histograms should be made with the decay-polarization
# weights, e.g.
#
#   wP1 = ev.weight_P1_joint
#   wM1 = ev.weight_M1_joint
# 
# Ie, ignoring the effect on cosTheta that comes from production
#
# IMPORTANT:
#   - P1/M1 templates are NOT normalized independently.
#   - One common normalization ("scale") floats in each omega fit.
#
# Two measurements are performed:
#
#   1) Inclusive omega fit:
#        omega(full range) -> P_tau -> A_tau = -P_tau
#
#   2) Angular fit:
#        omega in each cos(theta) bin -> P_tau(cos theta)
#        -> simultaneous fit for A_tau and A_e
# ============================================================


def ptau(z, Atau, Ae):
    return -(Atau*(1.0 + z*z) + 2.0*Ae*z) / (1.0 + z*z + 2.0*Ae*Atau*z)


def ptau_bin_average(zmin, zmax, Atau, Ae):
    """Exact production-weighted average P_tau in a finite cos(theta) bin."""
    I0 = (zmax - zmin) + (zmax**3 - zmin**3)/3.0
    I1 = zmax**2 - zmin**2
    return -(Atau*I0 + Ae*I1) / (I0 + Ae*Atau*I1)


def clone_hist(root_file, name, suffix):
    hist = root_file.Get(name)
    if not hist:
        raise RuntimeError("Missing histogram: " + name)

    out = hist.Clone(name + "_" + suffix)
    out.SetDirectory(0)
    return out


def get_parameter(minuit, index):
    value = ctypes.c_double(0.0)
    error = ctypes.c_double(0.0)
    minuit.GetParameter(index, value, error)
    return value.value, error.value


def fit_omega(hist_data, hist_minus, hist_plus, use_likelihood=True, print_level=1):
    """
    Fit the omega distribution with

      expected = scale * [
          (1-Ptau)/2 * T(Ptau=-1)
        + (1+Ptau)/2 * T(Ptau=+1)
      ]

    hist_minus = A_tau=+1 template = P_tau=-1
    hist_plus  = A_tau=-1 template = P_tau=+1

    No independent template normalization is applied.
    """

    nbins = hist_data.GetNbinsX()

    def fcn(npar, gin, fval, par, iflag):
        scale = par[0]
        pol = par[1]
        value = 0.0

        for i in range(1, nbins + 1):
            observed = hist_data.GetBinContent(i)
            n_minus = ((1.0 - pol)/2.0) * hist_minus.GetBinContent(i)
            n_plus = ((1.0 + pol)/2.0) * hist_plus.GetBinContent(i)
            expected = scale * (n_minus + n_plus)

            if expected <= 0.0:
                continue

            if use_likelihood:
                value += 2.0*(expected - observed*math.log(expected)) if observed > 0.0 else 2.0*expected
            else:
                value += (observed - expected)**2 / expected

        fval.value = value

    minuit = ROOT.TMinuit(2)
    minuit.SetFCN(fcn)
    minuit.SetPrintLevel(print_level)

    minuit.DefineParameter(0, "scale", 1.0, 1e-4, 0.0, 10.0)
    minuit.DefineParameter(1, "Ptau", -0.14, 1e-4, -1.0, 1.0)

    status = minuit.Command("MIGRAD")
    if status != 0:
        print("WARNING: omega MIGRAD did not converge")

    scale, scale_err = get_parameter(minuit, 0)
    pol, pol_err = get_parameter(minuit, 1)

    return scale, scale_err, pol, pol_err


def draw_omega_fit(hist_data, hist_minus, hist_plus, scale, pol, output_name, title=""):
    hist_minus_fit = hist_minus.Clone("hist_minus_fit_" + os.path.basename(output_name))
    hist_plus_fit = hist_plus.Clone("hist_plus_fit_" + os.path.basename(output_name))

    hist_minus_fit.Scale(scale*(1.0 - pol)/2.0)
    hist_plus_fit.Scale(scale*(1.0 + pol)/2.0)

    hist_fit = hist_minus_fit.Clone("hist_fit_" + os.path.basename(output_name))
    hist_fit.Add(hist_plus_fit)

    canvas = ROOT.TCanvas("c_" + os.path.basename(output_name), "", 800, 600)

    hist_data.SetTitle(title)
    hist_data.SetMarkerStyle(20)
    hist_data.SetMarkerSize(0.7)
    hist_data.SetLineColor(ROOT.kBlack)
    hist_data.GetXaxis().SetTitle("#omega_{#rho}")
    hist_data.GetYaxis().SetTitle("Events")

    ymax = max(hist_data.GetMaximum(), hist_fit.GetMaximum())
    hist_data.SetMaximum(1.20*ymax)

    hist_data.Draw("E")

    hist_fit.SetLineColor(ROOT.kRed)
    hist_fit.SetLineWidth(2)
    hist_fit.Draw("HIST SAME")

    hist_minus_fit.SetLineColor(ROOT.kBlue + 1)
    hist_minus_fit.SetLineStyle(2)
    hist_minus_fit.Draw("HIST SAME")

    hist_plus_fit.SetLineColor(ROOT.kGreen + 2)
    hist_plus_fit.SetLineStyle(2)
    hist_plus_fit.Draw("HIST SAME")

    legend = ROOT.TLegend(0.58, 0.68, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_data, "SM pseudo-data", "lep")
    legend.AddEntry(hist_fit, "Fit", "l")
    legend.AddEntry(hist_minus_fit, "P_{#tau}=-1", "l")
    legend.AddEntry(hist_plus_fit, "P_{#tau}=+1", "l")

    legend.Draw()

    canvas.SaveAs(output_name)
    canvas.Close()


def main(filename, filenameData,nBins=50, rebin=1, outdir="plots_pol", lumi=67.7,  use_likelihood=True, print_level=1):

    os.makedirs(outdir, exist_ok=True)

    root_file = ROOT.TFile.Open(filename)
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot open " + filename)

    root_file_data = ROOT.TFile.Open(filenameData)
    if not root_file_data or root_file_data.IsZombie():
        raise RuntimeError("Cannot open " + filenameData)


    bin_width = 2.0/nBins

    vectorPol = {}
    vectorPolError = {}
    vectorCosThetaMin = {}
    vectorCosThetaMax = {}


    # ========================================================
    # 1) INCLUSIVE OMEGA FIT
    # ========================================================

    hist_minus_full = clone_hist(root_file, "histo_P1_TAUMINUS_full", "minus_full")
    hist_plus_full = clone_hist(root_file, "histo_M1_TAUMINUS_full", "plus_full")
    hist_data_full = clone_hist(root_file_data, "histo_TAUMINUS_full", "data_full")

    if rebin > 1:
        hist_minus_full.Rebin(rebin)
        hist_plus_full.Rebin(rebin)
        hist_data_full.Rebin(rebin)

    print("\n============================================================")
    print("INCLUSIVE OMEGA FIT")
    print("============================================================")
    print("Raw integrals:")
    print("  data      = %.6f" % hist_data_full.Integral())
    print("  Ptau = -1 = %.6f" % hist_minus_full.Integral())
    print("  Ptau = +1 = %.6f" % hist_plus_full.Integral())

    scale_full, scale_full_err, Ptau_full, Ptau_full_err = fit_omega(
        hist_data_full,
        hist_minus_full,
        hist_plus_full,
        use_likelihood=use_likelihood,
        print_level=print_level
    )

    Atau_inclusive = -Ptau_full
    Atau_inclusive_err = Ptau_full_err

    print("\n============================================================")
    print("INCLUSIVE RESULT")
    print("Ptau = %.8f +- %.8f" % (Ptau_full, Ptau_full_err))
    print("Atau = %.8f +- %.8f" % (Atau_inclusive, Atau_inclusive_err))
    print("scale = %.8f +- %.8f" % (scale_full, scale_full_err))
    print("============================================================")

    draw_omega_fit(
        hist_data_full,
        hist_minus_full,
        hist_plus_full,
        scale_full,
        Ptau_full,
        os.path.join(outdir, "fit_omega_full.png"),
        "Inclusive #omega fit"
    )


    # ========================================================
    # 2) P_tau IN EACH cos(theta) BIN
    # ========================================================

    for ibin in range(nBins):

        zmin = -1.0 + ibin*bin_width
        zmax = zmin + bin_width
        bin_name = "_" + str(ibin)

        vectorCosThetaMin[ibin] = zmin
        vectorCosThetaMax[ibin] = zmax

        hist_minus = clone_hist(root_file, "histo_P1_TAUMINUS" + bin_name, "minus_" + str(ibin))
        hist_plus = clone_hist(root_file, "histo_M1_TAUMINUS" + bin_name, "plus_" + str(ibin))
        hist_data = clone_hist(root_file_data, "histo_TAUMINUS" + bin_name, "data_" + str(ibin))

        if rebin > 1:
            hist_minus.Rebin(rebin)
            hist_plus.Rebin(rebin)
            hist_data.Rebin(rebin)

        print("\n============================================================")
        print("BIN %d, cosTheta [%.4f, %.4f]" % (ibin, zmin, zmax))

        scale, scale_err, pol, pol_err = fit_omega(
            hist_data,
            hist_minus,
            hist_plus,
            use_likelihood=use_likelihood,
            print_level=print_level
        )

        vectorPol[ibin] = pol
        vectorPolError[ibin] = pol_err

        print("Ptau = %.8f +- %.8f" % (pol, pol_err))
        print("scale = %.8f +- %.8f" % (scale, scale_err))

        draw_omega_fit(
            hist_data,
            hist_minus,
            hist_plus,
            scale,
            pol,
            os.path.join(outdir, "fit_omega_bin_%02d.png" % ibin),
            "cos#theta bin %d" % ibin
        )


    # ========================================================
    # 3) COMBINED P_tau(cos theta) FIT FOR A_tau AND A_e
    # ========================================================

    def fcn_aeatau(npar, gin, fval, par, iflag):
        Atau = par[0]
        Ae = par[1]
        chi2 = 0.0
        val = 0.0

#        for bin in range(nBins):
#
#            ctheta = 0.5*(vectorCosThetaMin[bin] + vectorCosThetaMax[bin])
#
#            observed = vectorPol[bin]
#            observed_err = vectorPolError[bin]
#
#            Atau = par[0]
#            Ae = par[1]
#
#            expected = -(Atau*(1 + ctheta*ctheta) + 2*Ae*ctheta) / \
#                       (1 + ctheta*ctheta + 2*Ae*Atau*ctheta)
#
#            if observed_err > 0:
#               val += ((observed - expected)/observed_err)**2
#
#        fval.value = val

#       Alternative: do not assume narrow bins 
        for ibin in range(nBins):
            observed = vectorPol[ibin]
            error = vectorPolError[ibin]

            if error <= 0.0:
                continue
            expected = ptau_bin_average(vectorCosThetaMin[ibin], vectorCosThetaMax[ibin], Atau, Ae)
            chi2 += ((observed - expected)/error)**2
            fval.value = chi2




    minuit2 = ROOT.TMinuit(2)
    minuit2.SetFCN(fcn_aeatau)
    minuit2.SetPrintLevel(print_level)

    minuit2.DefineParameter(0, "Atau", 0.147, 1e-5, 0.0, 0.3)
    minuit2.DefineParameter(1, "Ae", 0.147, 1e-5, 0.0, 0.3)

    status2 = minuit2.Command("MIGRAD")
    if status2 != 0:
        print("WARNING: Atau/Ae MIGRAD did not converge")

    Atau, Atau_err = get_parameter(minuit2, 0)
    Ae, Ae_err = get_parameter(minuit2, 1)

    chi2 = 0.0
    nused = 0

    for ibin in range(nBins):
        observed = vectorPol[ibin]
        error = vectorPolError[ibin]

        if error <= 0.0:
            continue

        expected = ptau_bin_average(vectorCosThetaMin[ibin], vectorCosThetaMax[ibin], Atau, Ae)
        chi2 += ((observed - expected)/error)**2
        nused += 1

    ndf = nused - 2

    gv_ga = (1.0 - math.sqrt(1.0 - Ae*Ae))/Ae
    term1 = Ae*Ae/math.sqrt(1.0 - Ae*Ae)
    term2 = 1.0 - math.sqrt(1.0 - Ae*Ae)
    d_gv_ga_d_Ae = (term1 - term2)/(Ae*Ae)
    gv_ga_err = abs(d_gv_ga_d_Ae)*Ae_err

    sin2theta_eff = (1.0 - gv_ga)/4.0
    sin2theta_eff_err = gv_ga_err/4.0

    print("\n============================================================")
    print("ANGULAR Ptau(cosTheta) FIT")
    print("Atau = %.8f +- %.8f" % (Atau, Atau_err))
    print("Ae   = %.8f +- %.8f" % (Ae, Ae_err))

    if ndf > 0:
        print("chi2/ndf = %.4f / %d = %.4f" % (chi2, ndf, chi2/ndf))

    print("gv/ga = %.8f +- %.8f" % (gv_ga, gv_ga_err))
    print("sin2theta_eff = %.8f +- %.8f" % (sin2theta_eff, sin2theta_eff_err))
    print("============================================================")


    # ========================================================
    # Final P_tau(cos theta) plot
    # ========================================================

    graph = ROOT.TGraphErrors()
    graph.SetName("P_tau_vs_costheta")

    for ibin in range(nBins):
        zmin = vectorCosThetaMin[ibin]
        zmax = vectorCosThetaMax[ibin]
        x = 0.5*(zmin + zmax)
        xerr = 0.5*(zmax - zmin)

        graph.SetPoint(ibin, x, vectorPol[ibin])
        graph.SetPointError(ibin, xerr, vectorPolError[ibin])

    canvas = ROOT.TCanvas("canvas_pol", "", 800, 600)

    graph.SetMarkerStyle(20)
    graph.GetXaxis().SetTitle("cos #theta_{#tau}")
    graph.GetYaxis().SetTitle("P_{#tau}")
    graph.Draw("AP")

    fit_func = ROOT.TF1(
        "fit_func",
        "-([0]*(1+x*x)+2*[1]*x)/(1+x*x+2*[1]*[0]*x)",
        -1.0,
        1.0
    )

    fit_func.SetParameters(Atau, Ae)
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(2)
    fit_func.Draw("SAME")

    legend = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry("NULL","#sqrt{s}=91 GeV, %4.2f  fb^{-1}" %lumi ,"")
    legend.AddEntry(graph, "Extracted P_{#tau}", "pl")
    legend.AddEntry(fit_func, "Fit for A_{#tau}, A_{e}", "l")
    legend.Draw()

    canvas.SaveAs(os.path.join(outdir, "polarization_vs_costheta_fit.png"))
    canvas.SaveAs(os.path.join(outdir, "polarization_vs_costheta_fit.pdf"))


    # ========================================================
    # Summary
    # ========================================================

    summary = os.path.join(outdir, "fit_summary.txt")

    with open(summary, "w") as fout:
        fout.write("# Inclusive omega fit\n")
        fout.write("Ptau_inclusive  %.10f  %.10f\n" % (Ptau_full, Ptau_full_err))
        fout.write("Atau_inclusive  %.10f  %.10f\n" % (Atau_inclusive, Atau_inclusive_err))
        fout.write("scale_inclusive %.10f  %.10f\n" % (scale_full, scale_full_err))

        fout.write("\n# Angular Ptau(cosTheta) fit\n")
        fout.write("Atau_angular  %.10f  %.10f\n" % (Atau, Atau_err))
        fout.write("Ae_angular    %.10f  %.10f\n" % (Ae, Ae_err))
        fout.write("gv_ga         %.10f  %.10f\n" % (gv_ga, gv_ga_err))
        fout.write("sin2theta_eff %.10f  %.10f\n" % (sin2theta_eff, sin2theta_eff_err))
        fout.write("chi2          %.10f\n" % chi2)
        fout.write("ndf           %d\n" % ndf)

        fout.write("\n# Per-bin Ptau\n")
        fout.write("# bin  zmin  zmax  Ptau  error\n")

        for ibin in range(nBins):
            fout.write("%d  %.8f  %.8f  %.10f  %.10f\n" % (
                ibin,
                vectorCosThetaMin[ibin],
                vectorCosThetaMax[ibin],
                vectorPol[ibin],
                vectorPolError[ibin]
            ))

    root_file.Close()

    print("\nSummary written to:", summary)


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Tau polarization fits: inclusive omega -> Atau, and Ptau(cosTheta) -> Atau,Ae"
    )

    parser.add_argument("-i", "--input", default="BINED_templates_PY8WI23_GEN_LONG.root")
    parser.add_argument("-d", "--data", default="BINED_templates_PY8WI23_GEN_LONG.root")

    parser.add_argument("--nBins", type=int, default=50)
    parser.add_argument("--rebin", type=int, default=1)
    parser.add_argument("--chi2", action="store_true", help="Use chi2 instead of binned Poisson likelihood for omega fits")
    parser.add_argument("--print-level", type=int, default=1, help="TMinuit print level")
    parser.add_argument("-o", "--outdir", default="plots_pol")
    parser.add_argument("-l", "--lumi", type=float, default=67.7)


    args = parser.parse_args()

    main(
        filename=args.input,
        filenameData=args.data,
        nBins=args.nBins,
        rebin=args.rebin,
        outdir=args.outdir,
        use_likelihood=not args.chi2,
        print_level=args.print_level,
        lumi=args.lumi
    )

