#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import griddata
import ROOT as r

r.gROOT.SetBatch(True)


# -----------------------------------------------------
# Palette definition (CMS-style smooth gradients)
# -----------------------------------------------------
@dataclass
class Palette:
    stops: list
    reds: list
    greens: list
    blues: list
    contours: int = 80

PALETTES = {
    "blue-green": Palette(
        stops=[0.00, 0.30, 0.55, 0.75, 1.00],
        reds=[0.043, 0.157, 0.298, 0.510, 0.956],
        greens=[0.075, 0.275, 0.490, 0.725, 0.894],
        blues=[0.180, 0.478, 0.592, 0.478, 0.356],
        contours=120,
    )
}


def set_palette(pal):
    stops = r.std.vector("double")()
    reds = r.std.vector("double")()
    greens = r.std.vector("double")()
    blues = r.std.vector("double")()
    for s, rr, gg, bb in zip(pal.stops, pal.reds, pal.greens, pal.blues):
        stops.push_back(s)
        reds.push_back(rr)
        greens.push_back(gg)
        blues.push_back(bb)

    r.TColor.CreateGradientColorTable(
        len(pal.stops),
        stops.data(), reds.data(), greens.data(), blues.data(),
        pal.contours
    )
    r.gStyle.SetNumberContours(pal.contours)


# -----------------------------------------------------
# Load Combine scan points
# -----------------------------------------------------
def load_points(tree, args):
    xs, ys, zs = [], [], []

    safe_globals = {"abs": abs, "min": min, "max": max}

    for entry in tree:
        qexp = getattr(entry, "quantileExpected", -1)

        # Remove warmup / prefit points
        # if qexp < -0.4:
        #     continue

        x = getattr(entry, args.xvar)
        y = getattr(entry, args.yvar)
        d = getattr(entry, "deltaNLL")
        if d < 0: 
            d = 0

        ctx = {args.xvar: x, args.yvar: y, "deltaNLL": d}
        if args.cut and not eval(args.cut, safe_globals, ctx):
            continue

        xs.append(x)
        ys.append(y)
        zs.append(2 * d)

    return np.array(xs), np.array(ys), np.array(zs)


# -----------------------------------------------------
# Extract contour lines
# -----------------------------------------------------
def extract_contours(hist, level):
    tmp = hist.Clone()
    tmp.SetContour(1)
    tmp.SetContourLevel(0, level)

    canv = r.TCanvas("tmp_cont", "", 100, 100)
    tmp.Draw("CONT LIST")
    canv.Update()

    out = []
    conts = r.gROOT.GetListOfSpecials().FindObject("contours")
    if conts and conts.GetSize() > 0:
        lst = conts.At(0)
        for i in range(lst.GetSize()):
            out.append(lst.At(i).Clone())

    if conts:
        r.gROOT.GetListOfSpecials().Remove(conts)
    canv.Close()

    return out


# -----------------------------------------------------
# Draw interpolated scan
# -----------------------------------------------------
def draw_scan(grid_x, grid_y, grid_z, best, args):

    hist = r.TH2F("h2", "",
                  args.npx, grid_x.min(), grid_x.max(),
                  args.npy, grid_y.min(), grid_y.max())

    flat_x = grid_x.flatten()
    flat_y = grid_y.flatten()
    flat_z = np.minimum(grid_z.flatten(), args.zmax)

    for xx, yy, zz in zip(flat_x, flat_y, flat_z):
        hist.Fill(xx, yy, zz)

    hist.SetDirectory(0)

    # -----------------------------------------------------
    # Canvas
    # -----------------------------------------------------
    c = r.TCanvas("c", "", 900, 820)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.20)
    c.SetBottomMargin(0.14)
    c.SetTopMargin(0.08)
    c.SetTickx()
    c.SetTicky()

    # Remove ROOT frame border artifacts
    c.SetFrameLineWidth(2)
    c.SetFrameBorderMode(0)
    c.SetBorderSize(0)
    c.SetBorderMode(0)

    # -----------------------------------------------------
    # Draw heatmap
    # -----------------------------------------------------
    hist.Draw("COLZ")
    hist.SetLineWidth(0)

    hist.GetXaxis().SetTitle(args.xlabel or args.xvar)
    hist.GetYaxis().SetTitle(args.ylabel or args.yvar)
    hist.GetZaxis().SetTitle("-2 #Delta log L")

    hist.GetXaxis().SetTitleSize(0.055)
    hist.GetYaxis().SetTitleSize(0.055)
    hist.GetZaxis().SetTitleSize(0.055)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetZaxis().SetLabelSize(0.045)

    hist.GetZaxis().SetTitleOffset(1.0)

    # -----------------------------------------------------
    # Draw 1σ / 2σ contours
    # -----------------------------------------------------
    for level, style in [(2.30, 1), (5.99, 7)]:
        graphs = extract_contours(hist, level)
        for g in graphs:
            g.SetLineColor(r.kBlack)
            g.SetLineWidth(3)
            g.SetLineStyle(style)
            g.Draw("L SAME")

    # ==========================================================
    # 1. SM point FIRST (so it stays UNDER Best fit)
    # ==========================================================
    sm_x, sm_y = 0.0, 0.0

    # Outline (black)
    sm_outline = r.TGraph(1)
    sm_outline.SetPoint(0, sm_x, sm_y)
    sm_outline.SetMarkerStyle(33)
    sm_outline.SetMarkerSize(3.8)
    sm_outline.SetMarkerColor(r.kBlack)
    sm_outline.Draw("P SAME")

    # Red diamond
    sm = r.TGraph(1)
    sm.SetPoint(0, sm_x, sm_y)
    sm.SetMarkerStyle(33)
    sm.SetMarkerSize(3.4)
    sm.SetMarkerColor(r.kRed + 1)
    sm.Draw("P SAME")


    # ==========================================================
    # 2. Best-fit point LAST (so it stays ON TOP of SM)
    # ==========================================================
    # Black outline (below)
    bf_outline = r.TGraph(1)
    bf_outline.SetPoint(0, best[0], best[1])
    bf_outline.SetMarkerStyle(43)
    bf_outline.SetMarkerSize(2.0)
    bf_outline.SetMarkerColor(r.kGreen)
    bf_outline.Draw("P SAME")

    # White cross (above)
    bf = r.TGraph(1)
    bf.SetPoint(0, best[0], best[1])
    bf.SetMarkerStyle(43)
    bf.SetMarkerSize(1.5)
    bf.SetMarkerColor(r.kGreen)
    bf.Draw("P SAME")

    # -----------------------------------------------------
    # Crosshair reference lines (x=0, y=0)
    # -----------------------------------------------------
    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    ymin = hist.GetYaxis().GetXmin()
    ymax = hist.GetYaxis().GetXmax()

    # vertical line at x = 0
    vline = r.TLine(0.0, ymin, 0.0, ymax)
    vline.SetLineColorAlpha(r.kGray+2, 0.35)
    vline.SetLineStyle(3)
    vline.SetLineWidth(2)
    vline.Draw("SAME")

    # horizontal line at y = 0
    hline = r.TLine(xmin, 0.0, xmax, 0.0)
    hline.SetLineColorAlpha(r.kGray+2, 0.35)
    hline.SetLineStyle(3)
    hline.SetLineWidth(2)
    hline.Draw("SAME")
    
    # -----------------------------------------------------
    # Legend
    # -----------------------------------------------------
    leg = r.TLegend(0.60, 0.68, 0.78, 0.87)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.038)

    g1 = r.TGraph()
    g1.SetLineColor(r.kBlack)
    g1.SetLineWidth(3)
    g1.SetLineStyle(1)
    leg.AddEntry(g1, "1 #sigma CL", "l")

    g2 = r.TGraph()
    g2.SetLineColor(r.kBlack)
    g2.SetLineWidth(3)
    g2.SetLineStyle(7)
    leg.AddEntry(g2, "2 #sigma CL", "l")

    g3 = r.TGraph()
    g3.SetMarkerStyle(43)
    g3.SetMarkerColor(r.kGreen)
    g3.SetMarkerSize(2.8)
    leg.AddEntry(g3, "Best fit", "p")

    g4 = r.TGraph()
    g4.SetMarkerStyle(33)
    g4.SetMarkerColor(r.kRed + 1)
    g4.SetMarkerSize(3.4)
    leg.AddEntry(g4, "SM", "p")

    leg.Draw()

    # -----------------------------------------------------
    # CMS + Preliminary + Luminosity
    # -----------------------------------------------------
    texCMS = r.TLatex(0.14, 0.93, "CMS")
    texCMS.SetNDC()
    texCMS.SetTextFont(61)
    texCMS.SetTextSize(0.055)
    texCMS.Draw()

    texPre = r.TLatex(0.27, 0.93, "Preliminary")
    texPre.SetNDC()
    texPre.SetTextFont(52)
    texPre.SetTextSize(0.045)
    texPre.Draw()

    texLumi = r.TLatex(0.78, 0.93, f"{args.lumi:.1f} fb^{{-1}} (13 TeV)")
    texLumi.SetNDC()
    texLumi.SetTextAlign(31)
    texLumi.SetTextFont(42)
    texLumi.SetTextSize(0.045)
    texLumi.Draw()

    # -----------------------------------------------------
    # Save
    # -----------------------------------------------------
    c.SaveAs(f"{args.output}.png")
    c.SaveAs(f"{args.output}.pdf")



# -----------------------------------------------------
# Args
# -----------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--xvar", default="FT8")
    ap.add_argument("--yvar", default="FT9")
    ap.add_argument("--xlabel", default="")
    ap.add_argument("--ylabel", default="")
    ap.add_argument("--output", default="scan2D_griddata")
    ap.add_argument("--zmax", type=float, default=10)
    ap.add_argument("--cut", default="")
    ap.add_argument("--npx", type=int, default=200)
    ap.add_argument("--npy", type=int, default=200)
    ap.add_argument("--lumi", type=float, default=138)
    ap.add_argument("--palette", default="blue-green")
    return ap.parse_args()


# -----------------------------------------------------
# Main
# -----------------------------------------------------
def main():
    args = parse_args()

    set_palette(PALETTES[args.palette])
    r.gStyle.SetOptStat(0)

    f = r.TFile.Open(args.input)
    t = f.Get("limit")

    xs, ys, zs = load_points(t, args)

    best_idx = np.argmin(zs)
    best = (xs[best_idx], ys[best_idx])
    print("===================================")
    print("2D scan summary")
    print(f"  Best-fit point:")
    print(f"    {args.xvar} = {best[0]:.6g}")
    print(f"    {args.yvar} = {best[1]:.6g}")
    print(f"    -2ΔlogL     = {zs[best_idx]:.6g}")
    print("")
    print("  SM point:")
    print(f"    {args.xvar} = 0.0")
    print(f"    {args.yvar} = 0.0")
    print("===================================")
    xi = np.linspace(xs.min(), xs.max(), args.npx)
    yi = np.linspace(ys.min(), ys.max(), args.npy)
    grid_x, grid_y = np.meshgrid(xi, yi)

    grid_z = griddata((xs, ys), zs, (grid_x, grid_y), method="cubic")

    mask = np.isnan(grid_z)
    grid_z[mask] = args.zmax

    draw_scan(grid_x, grid_y, grid_z, best, args)


if __name__ == "__main__":
    main()
