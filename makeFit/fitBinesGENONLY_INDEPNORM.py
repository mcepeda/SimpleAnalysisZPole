#!/usr/bin/env python3

import argparse
import ctypes
import math
import os
from array import array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def ptau(z, Atau, Ae):
    return -(Atau*(1.0 + z*z) + 2.0*Ae*z) / (1.0 + z*z + 2.0*Ae*Atau*z)


def ptau_bin_average(zmin, zmax, Atau, Ae):
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


def normalize_shape(hist, suffix):
    out = hist.Clone(hist.GetName() + "_" + suffix)
    out.SetDirectory(0)
    integral = out.Integral(1, out.GetNbinsX())
    if integral <= 0.0:
        raise RuntimeError("Cannot normalize histogram %s: integral = %g" % (hist.GetName(), integral))
    out.Scale(1.0/integral)
    return out, integral


def get_parameter(minuit, index):
    value = ctypes.c_double(0.0)
    error = ctypes.c_double(0.0)
    minuit.GetParameter(index, value, error)
    return value.value, error.value


def get_covariance(minuit):
    matrix = array("d", [0.0]*4)
    minuit.mnemat(matrix, 2)
    return matrix[0], matrix[1], matrix[3]


def poisson_deviance(hist_data, hist_minus, hist_plus, NM1, NP1):
    value = 0.0
    nused = 0
    for i in range(1, hist_data.GetNbinsX() + 1):
        observed = hist_data.GetBinContent(i)
        expected = NM1*hist_minus.GetBinContent(i) + NP1*hist_plus.GetBinContent(i)
        if expected <= 0.0:
            continue
        value += 2.0*(expected - observed + observed*math.log(observed/expected)) if observed > 0.0 else 2.0*expected
        nused += 1
    return value, nused - 2


def fit_omega(hist_data, hist_minus, hist_plus, use_likelihood=True, print_level=1):
    """Fit two yields multiplying unit-area Ptau=-1/+1 shapes."""

    nbins = hist_data.GetNbinsX()
    Ndata = hist_data.Integral(1, nbins)
    if Ndata <= 0.0:
        raise RuntimeError("Data histogram has zero/negative integral")

    P0 = -0.14
    NM10 = 0.5*(1.0 - P0)*Ndata
    NP10 = 0.5*(1.0 + P0)*Ndata

    def fcn(npar, gin, fval, par, iflag):
        NM1 = par[0]
        NP1 = par[1]
        value = 0.0

        for i in range(1, nbins + 1):
            observed = hist_data.GetBinContent(i)
            expected = NM1*hist_minus.GetBinContent(i) + NP1*hist_plus.GetBinContent(i)
            if expected <= 0.0:
                continue
            if use_likelihood:
                value += 2.0*(expected - observed*math.log(expected)) if observed > 0.0 else 2.0*expected
            else:
                value += (observed - expected)**2/expected

        fval.value = value

    minuit = ROOT.TMinuit(2)
    minuit.SetFCN(fcn)
    minuit.SetPrintLevel(print_level)
#    minuit.DefineParameter(0, "NM1", NM10, max(1.0, math.sqrt(NM10)), 0.0, 2.0*Ndata)
#    minuit.DefineParameter(1, "NP1", NP10, max(1.0, math.sqrt(NP10)), 0.0, 2.0*Ndata)
    minuit.DefineParameter(0, "NM1", NM10, 1, 0.0, 2.0*Ndata)
    minuit.DefineParameter(1, "NP1", NP10, 1, 0.0, 2.0*Ndata)


    status = minuit.Command("MIGRAD")
    if status != 0:
        print("WARNING: omega MIGRAD did not converge")

    NM1, NM1_err = get_parameter(minuit, 0)
    NP1, NP1_err = get_parameter(minuit, 1)
    var_NM1, cov_NM1_NP1, var_NP1 = get_covariance(minuit)

    Nfit = NM1 + NP1
    pol = (NP1 - NM1)/Nfit if Nfit > 0.0 else 0.0

    dP_dNM1 = -2.0*NP1/(Nfit*Nfit)
    dP_dNP1 =  2.0*NM1/(Nfit*Nfit)
    pol_var = dP_dNM1*dP_dNM1*var_NM1 + dP_dNP1*dP_dNP1*var_NP1 + 2.0*dP_dNM1*dP_dNP1*cov_NM1_NP1
    pol_err = math.sqrt(max(0.0, pol_var))

    pol_from_NM1 = 1.0 - 2.0*NM1/Ndata
    pol_from_NP1 = 2.0*NP1/Ndata - 1.0
    delta_pol = pol_from_NP1 - pol_from_NM1
    yield_closure = Nfit/Ndata

    deviance, deviance_ndf = poisson_deviance(hist_data, hist_minus, hist_plus, NM1, NP1)

    return {
        "NM1": NM1, "NM1_err": NM1_err,
        "NP1": NP1, "NP1_err": NP1_err,
        "Ndata": Ndata, "Nfit": Nfit,
        "Ptau": pol, "Ptau_err": pol_err,
        "Ptau_from_NM1": pol_from_NM1,
        "Ptau_from_NP1": pol_from_NP1,
        "delta_Ptau": delta_pol,
        "yield_closure": yield_closure,
        "deviance": deviance, "deviance_ndf": deviance_ndf,
    }


def print_omega_result(result):
    print("NM1 = %.6f +- %.6f" % (result["NM1"], result["NM1_err"]))
    print("NP1 = %.6f +- %.6f" % (result["NP1"], result["NP1_err"]))
    print("Nfit/Ndata = %.8f" % result["yield_closure"])
    print("Ptau from ratio = %.8f +- %.8f" % (result["Ptau"], result["Ptau_err"]))
    print("Ptau from NM1/data = %.8f" % result["Ptau_from_NM1"])
    print("Ptau from NP1/data = %.8f" % result["Ptau_from_NP1"])
    print("Delta Ptau(NP1-NM1 diagnostic) = %.3e" % result["delta_Ptau"])
    if result["deviance_ndf"] > 0:
        print("omega deviance/ndf = %.3f / %d = %.3f" % (result["deviance"], result["deviance_ndf"], result["deviance"]/result["deviance_ndf"]))


def draw_omega_fit(hist_data, hist_minus, hist_plus, result, output_name, title=""):
    hist_minus_fit = hist_minus.Clone("hist_minus_fit_" + os.path.basename(output_name))
    hist_plus_fit = hist_plus.Clone("hist_plus_fit_" + os.path.basename(output_name))
    hist_minus_fit.Scale(result["NM1"])
    hist_plus_fit.Scale(result["NP1"])

    hist_fit = hist_minus_fit.Clone("hist_fit_" + os.path.basename(output_name))
    hist_fit.Add(hist_plus_fit)

    canvas = ROOT.TCanvas("c_" + os.path.basename(output_name), "", 800, 600)
    hist_data.SetTitle(title)
    hist_data.SetMarkerStyle(20)
    hist_data.SetMarkerSize(0.7)
    hist_data.SetLineColor(ROOT.kBlack)
    hist_data.GetXaxis().SetTitle("#omega_{#rho}")
    hist_data.GetYaxis().SetTitle("Events")
    hist_data.SetMaximum(1.20*max(hist_data.GetMaximum(), hist_fit.GetMaximum()))
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

    legend = ROOT.TLegend(0.55, 0.66, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_data, "Data", "lep")
    legend.AddEntry(hist_fit, "Fit", "l")
    legend.AddEntry(hist_minus_fit, "N_{-} T(P_{#tau}=-1)", "l")
    legend.AddEntry(hist_plus_fit, "N_{+} T(P_{#tau}=+1)", "l")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.032)
    latex.DrawLatex(0.16, 0.88, "P_{#tau} = %.4f #pm %.4f" % (result["Ptau"], result["Ptau_err"]))
    latex.DrawLatex(0.16, 0.83, "(N_{-}+N_{+})/N_{data} = %.4f" % result["yield_closure"])

    canvas.SaveAs(output_name)
    canvas.Close()


def main(filename, filenameData, nBins=50, rebin=1, outdir="plots_pol_indepnorm", lumi=67.7, lumiTemplate=67.7, use_likelihood=True, print_level=1):

    os.makedirs(outdir, exist_ok=True)

    root_file = ROOT.TFile.Open(filename)
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot open template file " + filename)

    root_file_data = ROOT.TFile.Open(filenameData)
    if not root_file_data or root_file_data.IsZombie():
        raise RuntimeError("Cannot open data file " + filenameData)

    print("============================================================")
    print("Template file :", filename)
    print("Data file     :", filenameData)
    print("Data lumi     :", lumi)
    print("Template lumi :", lumiTemplate)
    print("NOTE: template luminosity is not used in the omega fit.")
    print("      P1/M1 are normalized independently and used as shapes only.")
    print("============================================================")

    bin_width = 2.0/nBins

    vectorPol = {}
    vectorPolError = {}
    vectorNM1 = {}
    vectorNM1Error = {}
    vectorNP1 = {}
    vectorNP1Error = {}
    vectorPfromNM1 = {}
    vectorPfromNP1 = {}
    vectorDeltaPol = {}
    vectorYieldClosure = {}
    vectorCosThetaMin = {}
    vectorCosThetaMax = {}

    # ========================================================
    # 1) INCLUSIVE OMEGA FIT
    # ========================================================

    hist_minus_raw = clone_hist(root_file, "histo_P1_TAUMINUS_full", "minus_full_raw")
    hist_plus_raw = clone_hist(root_file, "histo_M1_TAUMINUS_full", "plus_full_raw")
    hist_data_full = clone_hist(root_file_data, "histo_TAUMINUS_full", "data_full")

    if rebin > 1:
        hist_minus_raw.Rebin(rebin)
        hist_plus_raw.Rebin(rebin)
        hist_data_full.Rebin(rebin)

    hist_minus_full, Iminus_full = normalize_shape(hist_minus_raw, "shape")
    hist_plus_full, Iplus_full = normalize_shape(hist_plus_raw, "shape")

    print("\n============================================================")
    print("INCLUSIVE OMEGA FIT: INDEPENDENT SHAPE NORMALIZATION")
    print("============================================================")
    print("Raw template integrals:")
    print("  data      = %.6f" % hist_data_full.Integral())
    print("  Ptau = -1 = %.6f" % Iminus_full)
    print("  Ptau = +1 = %.6f" % Iplus_full)

    result_full = fit_omega(hist_data_full, hist_minus_full, hist_plus_full, use_likelihood, print_level)
    Ptau_full = result_full["Ptau"]
    Ptau_full_err = result_full["Ptau_err"]
    Atau_inclusive = -Ptau_full
    Atau_inclusive_err = Ptau_full_err

    print_omega_result(result_full)
    print("Atau inclusive = %.8f +- %.8f" % (Atau_inclusive, Atau_inclusive_err))

    draw_omega_fit(hist_data_full, hist_minus_full, hist_plus_full, result_full,
                   os.path.join(outdir, "fit_omega_full.png"), "Inclusive #omega fit")

    # ========================================================
    # 2) P_tau IN EACH cos(theta) BIN
    # ========================================================

    for ibin in range(nBins):
        zmin = -1.0 + ibin*bin_width
        zmax = zmin + bin_width
        bin_name = "_" + str(ibin)

        vectorCosThetaMin[ibin] = zmin
        vectorCosThetaMax[ibin] = zmax

        hist_minus_raw = clone_hist(root_file, "histo_P1_TAUMINUS" + bin_name, "minus_%d_raw" % ibin)
        hist_plus_raw = clone_hist(root_file, "histo_M1_TAUMINUS" + bin_name, "plus_%d_raw" % ibin)
        hist_data = clone_hist(root_file_data, "histo_TAUMINUS" + bin_name, "data_%d" % ibin)

        if rebin > 1:
            hist_minus_raw.Rebin(rebin)
            hist_plus_raw.Rebin(rebin)
            hist_data.Rebin(rebin)

        hist_minus, Iminus = normalize_shape(hist_minus_raw, "shape")
        hist_plus, Iplus = normalize_shape(hist_plus_raw, "shape")

        print("\n============================================================")
        print("BIN %d, cosTheta [%.4f, %.4f]" % (ibin, zmin, zmax))
        print("Raw template integrals: Ptau=-1 %.6f, Ptau=+1 %.6f" % (Iminus, Iplus))

        result = fit_omega(hist_data, hist_minus, hist_plus, use_likelihood, print_level)

        vectorPol[ibin] = result["Ptau"]
        vectorPolError[ibin] = result["Ptau_err"]
        vectorNM1[ibin] = result["NM1"]
        vectorNM1Error[ibin] = result["NM1_err"]
        vectorNP1[ibin] = result["NP1"]
        vectorNP1Error[ibin] = result["NP1_err"]
        vectorPfromNM1[ibin] = result["Ptau_from_NM1"]
        vectorPfromNP1[ibin] = result["Ptau_from_NP1"]
        vectorDeltaPol[ibin] = result["delta_Ptau"]
        vectorYieldClosure[ibin] = result["yield_closure"]

        print_omega_result(result)

        draw_omega_fit(hist_data, hist_minus, hist_plus, result,
                       os.path.join(outdir, "fit_omega_bin_%02d.png" % ibin), "cos#theta bin %d" % ibin)

    # ========================================================
    # 3) COMBINED P_tau(cos theta) FIT FOR A_tau AND A_e
    # ========================================================

    def fcn_aeatau(npar, gin, fval, par, iflag):
        Atau = par[0]
        Ae = par[1]
        chi2 = 0.0

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
    minuit2.DefineParameter(0, "Atau", 0.147, 1e-6, 0.0, 0.3)
    minuit2.DefineParameter(1, "Ae", 0.147, 1e-6, 0.0, 0.3)

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

    fit_func = ROOT.TF1("fit_func", "-([0]*(1+x*x)+2*[1]*x)/(1+x*x+2*[1]*[0]*x)", -1.0, 1.0)
    fit_func.SetParameters(Atau, Ae)
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(2)
    fit_func.Draw("SAME")

    legend = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry("NULL", "#sqrt{s}=91 GeV, %4.2f fb^{-1}" % lumi, "")
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
        fout.write("# Input\n")
        fout.write("data_lumi       %.10f\n" % lumi)
        fout.write("template_lumi   %.10f\n" % lumiTemplate)
        fout.write("# Template luminosity is informational only in INDEPNORM mode\n")

        fout.write("\n# Inclusive omega fit: independent normalized shapes\n")
        fout.write("Ptau_inclusive       %.10f  %.10f\n" % (Ptau_full, Ptau_full_err))
        fout.write("Atau_inclusive       %.10f  %.10f\n" % (Atau_inclusive, Atau_inclusive_err))
        fout.write("NM1_inclusive        %.10f  %.10f\n" % (result_full["NM1"], result_full["NM1_err"]))
        fout.write("NP1_inclusive        %.10f  %.10f\n" % (result_full["NP1"], result_full["NP1_err"]))
        fout.write("Ndata_inclusive      %.10f\n" % result_full["Ndata"])
        fout.write("Nfit_inclusive       %.10f\n" % result_full["Nfit"])
        fout.write("yield_closure        %.10f\n" % result_full["yield_closure"])
        fout.write("Ptau_from_NM1        %.10f\n" % result_full["Ptau_from_NM1"])
        fout.write("Ptau_from_NP1        %.10f\n" % result_full["Ptau_from_NP1"])
        fout.write("delta_Ptau_yields    %.10e\n" % result_full["delta_Ptau"])
        fout.write("omega_deviance       %.10f\n" % result_full["deviance"])
        fout.write("omega_ndf            %d\n" % result_full["deviance_ndf"])

        fout.write("\n# Angular Ptau(cosTheta) fit\n")
        fout.write("Atau_angular         %.10f  %.10f\n" % (Atau, Atau_err))
        fout.write("Ae_angular           %.10f  %.10f\n" % (Ae, Ae_err))
        fout.write("gv_ga                %.10f  %.10f\n" % (gv_ga, gv_ga_err))
        fout.write("sin2theta_eff        %.10f  %.10f\n" % (sin2theta_eff, sin2theta_eff_err))
        fout.write("chi2                 %.10f\n" % chi2)
        fout.write("ndf                  %d\n" % ndf)

        fout.write("\n# Per-bin Ptau and fitted yields\n")
        fout.write("# bin zmin zmax Ptau err NM1 errNM1 NP1 errNP1 PfromNM1 PfromNP1 deltaP NfitOverNdata\n")

        for ibin in range(nBins):
            fout.write("%d %.8f %.8f %.10f %.10f %.6f %.6f %.6f %.6f %.10f %.10f %.3e %.10f\n" % (
                ibin, vectorCosThetaMin[ibin], vectorCosThetaMax[ibin], vectorPol[ibin], vectorPolError[ibin],
                vectorNM1[ibin], vectorNM1Error[ibin], vectorNP1[ibin], vectorNP1Error[ibin],
                vectorPfromNM1[ibin], vectorPfromNP1[ibin], vectorDeltaPol[ibin], vectorYieldClosure[ibin]
            ))

    root_file.Close()
    root_file_data.Close()
    print("\nSummary written to:", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tau polarization fit with independently normalized P1/M1 shapes")
    parser.add_argument("-i", "--input", default="BINED_templates_PY8WI23_GEN_LONG.root")
    parser.add_argument("-d", "--data", default="BINED_templates_PY8WI23_GEN_LONG.root")
    parser.add_argument("--nBins", type=int, default=50)
    parser.add_argument("--rebin", type=int, default=1)
    parser.add_argument("--chi2", action="store_true", help="Use chi2 instead of extended Poisson likelihood")
    parser.add_argument("--print-level", type=int, default=1)
    parser.add_argument("-o", "--outdir", default="plots_pol_indepnorm")
    parser.add_argument("-l", "--lumi", type=float, default=67.7)
    parser.add_argument("-lt", "--lumiTemplate", type=float, default=67.7)

    args = parser.parse_args()
    main(args.input, args.data, args.nBins, args.rebin, args.outdir, args.lumi, args.lumiTemplate,
         not args.chi2, args.print_level)

