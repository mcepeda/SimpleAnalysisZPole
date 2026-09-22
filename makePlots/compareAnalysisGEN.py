#!/usr/bin/env python3

import argparse
import math
import os
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def sanitize(name):
    return (
        name.replace("/", "__")
            .replace(" ", "_")
            .replace("#", "")
    )


def get_hist1(file_, name):
    obj = file_.Get(name)
    if not obj:
        return None
    if not obj.InheritsFrom("TH1"):
        return None
    if obj.InheritsFrom("TH2"):
        return None

    h = obj.Clone(name + "_clone_" + str(id(obj)))
    h.SetDirectory(0)
    return h


def get_common_th1_names(file1, file2):
    names = []

    for key in file1.GetListOfKeys():
        name = key.GetName()
        obj1 = file1.Get(name)

        if not obj1:
            continue
        if not obj1.InheritsFrom("TH1"):
            continue
        if obj1.InheritsFrom("TH2"):
            continue

        obj2 = file2.Get(name)
        if not obj2:
            continue
        if not obj2.InheritsFrom("TH1"):
            continue
        if obj2.InheritsFrom("TH2"):
            continue

        names.append(name)

    return sorted(names)


def get_meson_costheta_projections(file1, file2):
    """
    Build 1D cos(theta_meson) histograms from the Y projection of
    GENOmegaCosThetaMeson... TH2 histograms.

    The TH2 itself is never plotted.
    """
    projections = []

    prefix = "GENOmegaCosThetaMeson"

    for key in file1.GetListOfKeys():
        name = key.GetName()

        if not name.startswith(prefix):
            continue

        h2_1 = file1.Get(name)
        h2_2 = file2.Get(name)

        if not h2_1 or not h2_2:
            continue
        if not h2_1.InheritsFrom("TH2"):
            continue
        if not h2_2.InheritsFrom("TH2"):
            continue

        suffix = name[len(prefix):]
        outname = "GENCosThetaMeson" + suffix

        h1 = h2_1.ProjectionY(
            outname + "_file1",
            1,
            h2_1.GetNbinsX(),
            "e"
        )
        h2 = h2_2.ProjectionY(
            outname + "_file2",
            1,
            h2_2.GetNbinsX(),
            "e"
        )

        h1.SetDirectory(0)
        h2.SetDirectory(0)

        projections.append((outname, h1, h2))

    return sorted(projections, key=lambda x: x[0])


def normalize_hist(h):
    integral = h.Integral()
    if integral != 0:
        h.Scale(1.0 / integral)


def style_hist(h, color):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetLineWidth(2)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.55)


def compatible_binning(h1, h2):
    if h1.GetNbinsX() != h2.GetNbinsX():
        return False

    ax1 = h1.GetXaxis()
    ax2 = h2.GetXaxis()

    for i in range(1, h1.GetNbinsX() + 2):
        if abs(ax1.GetBinLowEdge(i) - ax2.GetBinLowEdge(i)) > 1e-9:
            return False

    return True


def is_ratio_observable(name):
    """
    Ratio plots requested for:
      - optimal variable omega
      - true tau cos(theta)
      - meson cos(theta)
      - cos(theta_hat)

    Do not select TH2 names here.
    """
    if name.startswith("GENOmega") and "CosTheta" not in name:
        return True

    if name == "GENCosThetaSimple":
        return True

    if name == "GENCosThetaHat":
        return True

    if name.startswith("GENCosTheta_"):
        return True

    if name.startswith("GENCosThetaHat_"):
        return True

    if name.startswith("GENCosThetaMeson"):
        return True

    return False


def ratio_hist(h_num, h_den, name):
    r = h_num.Clone(name)
    r.SetDirectory(0)
    r.Reset("ICES")

    # Explicit division: bins with zero reference content are left empty.
    for i in range(1, h_num.GetNbinsX() + 1):
        num = h_num.GetBinContent(i)
        den = h_den.GetBinContent(i)
        en = h_num.GetBinError(i)
        ed = h_den.GetBinError(i)

        if den == 0:
            r.SetBinContent(i, 0.0)
            r.SetBinError(i, 0.0)
            continue

        value = num / den

        # Independent-error propagation; useful for plotting only.
        # These two files may be statistically correlated.
        if num != 0:
            err = abs(value) * math.sqrt((en / num) ** 2 + (ed / den) ** 2)
        else:
            err = abs(en / den)

        r.SetBinContent(i, value)
        r.SetBinError(i, err)

    return r


def basic_metrics(h1, h2):
    i1 = h1.Integral()
    i2 = h2.Integral()

    integral_ratio = i2 / i1 if i1 != 0 else float("nan")

    return {
        "integral1": i1,
        "integral2": i2,
        "integral_ratio": integral_ratio,
        "mean1": h1.GetMean(),
        "mean2": h2.GetMean(),
        "rms1": h1.GetRMS(),
        "rms2": h2.GetRMS(),
    }


def draw_overlay(
    h1,
    h2,
    name,
    label1,
    label2,
    outdir,
    normalize=False,
    make_ratio=False,
    ratio_min=0.8,
    ratio_max=1.2,
    save_pdf=False
):
    if not compatible_binning(h1, h2):
        print("WARNING: incompatible binning, skipping:", name)
        return None

    # Preserve unscaled values for the summary.
    metrics = basic_metrics(h1, h2)

    p1 = h1.Clone(name + "_plot1")
    p2 = h2.Clone(name + "_plot2")
    p1.SetDirectory(0)
    p2.SetDirectory(0)

    if normalize:
        normalize_hist(p1)
        normalize_hist(p2)

    style_hist(p1, ROOT.kBlue + 1)
    style_hist(p2, ROOT.kRed + 1)

    if make_ratio:
        canvas = ROOT.TCanvas("c_" + sanitize(name), "", 800, 800)

        pad1 = ROOT.TPad("pad1_" + sanitize(name), "", 0.0, 0.30, 1.0, 1.0)
        pad2 = ROOT.TPad("pad2_" + sanitize(name), "", 0.0, 0.00, 1.0, 0.30)

        pad1.SetBottomMargin(0.025)
        pad2.SetTopMargin(0.035)
        pad2.SetBottomMargin(0.32)

        pad1.Draw()
        pad2.Draw()

        pad1.cd()

    else:
        canvas = ROOT.TCanvas("c_" + sanitize(name), "", 800, 700)
        canvas.cd()

    ymax = max(p1.GetMaximum(), p2.GetMaximum())
    if ymax <= 0:
        ymax = 1.0

    p1.SetMaximum(1.25 * ymax)
    p1.SetMinimum(0.0)

    ytitle = "Normalized entries" if normalize else "Entries"
    p1.GetYaxis().SetTitle(ytitle)

    title = name
    if normalize:
        title += " [shape normalized]"
    p1.SetTitle(title)

    p1.Draw("HIST E")
    p2.Draw("HIST E SAME")

    legend = ROOT.TLegend(0.62, 0.76, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(p1, label1, "l")
    legend.AddEntry(p2, label2, "l")
    legend.Draw()

    # Useful small validation annotation.
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.030)
    latex.DrawLatex(
        0.14, 0.84,
        "I_{2}/I_{1} = %.5g" % metrics["integral_ratio"]
        if math.isfinite(metrics["integral_ratio"])
        else "I_{2}/I_{1} = undefined"
    )

    if make_ratio:
        # Hide upper x labels.
        p1.GetXaxis().SetLabelSize(0.0)
        p1.GetXaxis().SetTitleSize(0.0)

        pad2.cd()

        r = ratio_hist(p2, p1, name + "_ratio")
        r.SetLineColor(ROOT.kBlack)
        r.SetMarkerColor(ROOT.kBlack)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0.55)

        r.SetTitle("")
        r.GetYaxis().SetTitle(label2 + "/" + label1)
        r.GetYaxis().SetRangeUser(ratio_min, ratio_max)
        r.GetYaxis().SetNdivisions(505)

        r.GetYaxis().SetTitleSize(0.085)
        r.GetYaxis().SetTitleOffset(0.48)
        r.GetYaxis().SetLabelSize(0.075)

        r.GetXaxis().SetTitle(p1.GetXaxis().GetTitle())
        r.GetXaxis().SetTitleSize(0.10)
        r.GetXaxis().SetTitleOffset(0.95)
        r.GetXaxis().SetLabelSize(0.085)

        r.Draw("E1")

        line = ROOT.TLine(
            r.GetXaxis().GetXmin(),
            1.0,
            r.GetXaxis().GetXmax(),
            1.0
        )
        line.SetLineStyle(2)
        line.Draw()

    canvas.Update()

    base = Path(outdir) / sanitize(name)
    canvas.SaveAs(str(base) + ".png")

    if save_pdf:
        canvas.SaveAs(str(base) + ".pdf")

    canvas.Close()

    return metrics


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare two analysisGENFromTREE.py output ROOT files. "
            "All common 1D histograms are overlaid; omega and production-angle "
            "observables additionally get ratio pads."
        )
    )

    parser.add_argument("file1", help="Reference ROOT file")
    parser.add_argument("file2", help="Comparison ROOT file")

    parser.add_argument("--label1", default="Reference")
    parser.add_argument("--label2", default="Alternative")
    parser.add_argument("--outdir", default="validation_plots")

    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize each histogram to unit area before plotting"
    )

    parser.add_argument("--ratio-min", type=float, default=0.8)
    parser.add_argument("--ratio-max", type=float, default=1.2)

    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Also save PDF versions"
    )

    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    f1 = ROOT.TFile.Open(args.file1)
    f2 = ROOT.TFile.Open(args.file2)

    if not f1 or f1.IsZombie():
        raise RuntimeError("Cannot open " + args.file1)

    if not f2 or f2.IsZombie():
        raise RuntimeError("Cannot open " + args.file2)

    common = get_common_th1_names(f1, f2)

    print("Found %d common 1D histograms" % len(common))

    summary_rows = []

    # --------------------------------------------------------------
    # All ordinary common TH1 histograms
    # --------------------------------------------------------------

    for name in common:
        h1 = get_hist1(f1, name)
        h2 = get_hist1(f2, name)

        if h1 is None or h2 is None:
            continue

        do_ratio = is_ratio_observable(name)

        metrics = draw_overlay(
            h1,
            h2,
            name,
            args.label1,
            args.label2,
            args.outdir,
            normalize=args.normalize,
            make_ratio=do_ratio,
            ratio_min=args.ratio_min,
            ratio_max=args.ratio_max,
            save_pdf=args.pdf
        )

        if metrics is not None:
            summary_rows.append((name, metrics, do_ratio))

    # --------------------------------------------------------------
    # Meson cos(theta):
    # derive 1D Y projections from the existing TH2 histograms.
    # No 2D histogram is drawn.
    # --------------------------------------------------------------

    meson_projections = get_meson_costheta_projections(f1, f2)

    print(
        "Built %d meson-cos(theta) 1D projections from TH2 inputs"
        % len(meson_projections)
    )

    for name, h1, h2 in meson_projections:
        h1.GetXaxis().SetTitle("cos #theta_{meson}")
        h2.GetXaxis().SetTitle("cos #theta_{meson}")

        metrics = draw_overlay(
            h1,
            h2,
            name,
            args.label1,
            args.label2,
            args.outdir,
            normalize=args.normalize,
            make_ratio=True,
            ratio_min=args.ratio_min,
            ratio_max=args.ratio_max,
            save_pdf=args.pdf
        )

        if metrics is not None:
            summary_rows.append((name, metrics, True))

    # --------------------------------------------------------------
    # Text summary
    # --------------------------------------------------------------

    summary_path = Path(args.outdir) / "validation_summary.tsv"

    with open(summary_path, "w") as out:
        out.write(
            "histogram\tintegral_file1\tintegral_file2\t"
            "integral_ratio_2_over_1\tmean_file1\tmean_file2\t"
            "rms_file1\trms_file2\tratio_plot\n"
        )

        for name, m, do_ratio in summary_rows:
            out.write(
                "%s\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%s\n"
                % (
                    name,
                    m["integral1"],
                    m["integral2"],
                    m["integral_ratio"],
                    m["mean1"],
                    m["mean2"],
                    m["rms1"],
                    m["rms2"],
                    str(do_ratio),
                )
            )

    print("Wrote plots to:", args.outdir)
    print("Wrote summary:", summary_path)

    f1.Close()
    f2.Close()


if __name__ == "__main__":
    main()
