#!/usr/bin/env python3

import argparse
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

SAMPLES = [
    ("TAUMINUS", "SM", ROOT.kBlack),
    ("P1_TAUMINUS", "A_{#tau}=+1", ROOT.kGreen + 2),
    ("M1_TAUMINUS", "A_{#tau}=-1", ROOT.kRed),
]


def get_hist2d(root_file, name):
    hist = root_file.Get(name)
    if not hist:
        raise RuntimeError("Missing histogram: " + name)
    return hist


def main(input_file, output_file, variable, nBins=20, rebin=1, plot_dir="BINS", make_plots=True):

    fin = ROOT.TFile.Open(input_file)
    if not fin or fin.IsZombie():
        raise RuntimeError("Cannot open input file: " + input_file)

    first = get_hist2d(fin, variable + "_" + SAMPLES[0][0])
    nY = first.GetNbinsY()

    if nY % nBins != 0:
        raise RuntimeError("Y axis has %d bins, not divisible by nBins=%d" % (nY, nBins))

    ybins_per_output_bin = nY // nBins

    fout = ROOT.TFile(output_file, "RECREATE")
    if not fout or fout.IsZombie():
        raise RuntimeError("Cannot create output file: " + output_file)

    if make_plots:
        os.makedirs(plot_dir, exist_ok=True)

    print("============================================================")
    print("Input      :", input_file)
    print("Output     :", output_file)
    print("Variable   :", variable)
    print("Y bins     :", nY)
    print("Output bins:", nBins)
    print("Y/bin      :", ybins_per_output_bin)
    print("Omega rebin:", rebin)
    print("============================================================")

    # --------------------------------------------------------
    # Angular bins
    # --------------------------------------------------------

    for ibin in range(nBins):

        bin_ini = ibin*ybins_per_output_bin + 1
        bin_end = (ibin + 1)*ybins_per_output_bin

        print("Angular bin %d: Y bins %d -> %d" % (ibin, bin_ini, bin_end))

        projected = []
        max_y = 0.0

        if make_plots:
            canvas = ROOT.TCanvas("canvas_%d" % ibin, "", 800, 800)
            legend = ROOT.TLegend(0.52, 0.66, 0.90, 0.89)
            legend.SetFillStyle(0)
            legend.SetLineColor(0)
            legend.SetLineWidth(0)

        for sample_tag, sample_label, color in SAMPLES:

            hist2d = get_hist2d(fin, variable + "_" + sample_tag)
            hist = hist2d.ProjectionX("histo_%s_%d" % (sample_tag, ibin), bin_ini, bin_end)
            hist.SetDirectory(0)

            if rebin > 1:
                hist.Rebin(rebin)

            hist.SetXTitle("#omega_{#rho}")
            hist.SetLineWidth(2)
            hist.SetLineColor(color)

            projected.append(hist)

            print("  %-16s integral = %.6f" % (sample_tag, hist.Integral()))

            if make_plots:
                max_y = max(max_y, hist.GetMaximum())
                legend.AddEntry(hist, sample_label, "l")

        if make_plots:
            projected[0].SetMaximum(1.2*max_y)
            projected[0].Draw("hist")

            for hist in projected[1:]:
                hist.Draw("hist same")

            legend.Draw()
            canvas.SaveAs(os.path.join(plot_dir, "%s_bin_%02d.png" % (variable, ibin)))
            canvas.Close()

        fout.cd()
        for hist in projected:
            hist.Write()

    # --------------------------------------------------------
    # Full angular range
    # --------------------------------------------------------

    print("Full range: Y bins 1 -> %d" % nY)

    projected = []
    max_y = 0.0

    if make_plots:
        canvas = ROOT.TCanvas("canvas_full", "", 800, 800)
        legend = ROOT.TLegend(0.52, 0.66, 0.90, 0.89)
        legend.SetFillStyle(0)
        legend.SetLineColor(0)
        legend.SetLineWidth(0)

    for sample_tag, sample_label, color in SAMPLES:

        hist2d = get_hist2d(fin, variable + "_" + sample_tag)
        hist = hist2d.ProjectionX("histo_%s_full" % sample_tag, 1, nY)
        hist.SetDirectory(0)

        if rebin > 1:
            hist.Rebin(rebin)

        hist.SetXTitle("#omega_{#rho}")
        hist.SetLineWidth(2)
        hist.SetLineColor(color)

        projected.append(hist)

        print("  %-16s integral = %.6f" % (sample_tag, hist.Integral()))

        if make_plots:
            max_y = max(max_y, hist.GetMaximum())
            legend.AddEntry(hist, sample_label, "l")

    if make_plots:
        projected[0].SetMaximum(1.2*max_y)
        projected[0].Draw("hist")

        for hist in projected[1:]:
            hist.Draw("hist same")

        legend.Draw()
        canvas.SaveAs(os.path.join(plot_dir, "%s_full.png" % variable))
        canvas.Close()

    fout.cd()
    for hist in projected:
        hist.Write()

    fout.Close()
    fin.Close()

    print("Wrote:", output_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Project omega histograms in cos(theta) bins")

    parser.add_argument("-i", "--input", required=True, help="Input ROOT file")
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file")
    parser.add_argument("-v", "--variable", default="GENOmegaCosThetaHat", help="2D variable")
    parser.add_argument("--nBins", type=int, default=20, help="Number of cos(theta) bins")
    parser.add_argument("--rebin", type=int, default=1, help="Rebin factor for omega axis")
    parser.add_argument("--plot-dir", default="BINS", help="Directory for diagnostic plots")
    parser.add_argument("--no-plots", action="store_true", help="Do not make projection plots")

    args = parser.parse_args()

    main(args.input, args.output, args.variable, args.nBins, args.rebin, args.plot_dir, not args.no_plots)

