#!/usr/bin/env python3
import ROOT as R
import os, re, math

R.gROOT.SetBatch(True)

# ========= 配置 =========
tag = "2018"
#tag = "2017"
#tag = "2016"
#tag = "2016APV"
#op_label  = "cT8"                         # AAC 进程后缀
hist_name = "eft_ZZ2l2nu"                 # 直方图基名
sm_file   = f"shapes-vbs-SR{tag}sm.root"  # SM 文件
aac_file  = f"aac_FT8_templates_{tag}.root"

# 要验证的 c 点
# c_points = [0.0, 0.1, -0.1, 0.2, -0.2, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0]
c_points = [
    0.0,
    0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
   -0.1,-0.2,-0.3,-0.5,-1.0,-2.0,-4.0,-8.0,-16.0,-32.0,-64.0,-128.0
]
# 要拼到大画布上的 c 点（只做 nominal）
grid_c_points = [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]

# 只做 nominal 验证？（True=只做 nominal）
ONLY_NOMINAL = True

# 输出目录
outdir = f"validation_plots_{tag}"
os.makedirs(outdir, exist_ok=True)

# ========= 工具 =========
def c_to_token(c: float) -> str:
    """把数值 c 转成文件名里的标记:
       +0.1 -> 0p1,  -0.1 -> m0p1,  +2.0 -> 2,  -2.0 -> m2
    """
    if abs(c) < 1e-12:
        return "0"  # c=0 单独处理
    sign = "m" if c < 0 else ""
    s = f"{abs(c):g}"
    s = s.replace(".", "p")
    return sign + s

def token_from_c(c: float) -> str:
    """把 c 值转换成文件名 token，含 0 处理。"""
    return "0" if abs(c) < 1e-12 else c_to_token(c)

def list_nuis_from_aac(tf):
    nuis = set()
    for key in tf.GetListOfKeys():
        n = key.GetName()
        m = (re.match(r"^(sm|Sm)_(.+)(Up|Down)$", n) or
             re.match(r"^(sm_lin_quad_.+|SmLinQuad)_(.+)(Up|Down)$", n) or
             re.match(r"^(quad_.+|Quad)_(.+)(Up|Down)$", n))
        if m:
            nuis.add(m.group(2))
    return sorted(nuis)

def get_hist(tfile, name):
    h = tfile.Get(name)
    if not h:
        return None
    h = h.Clone(f"cl_{name}")
    h.SetDirectory(0)
    return h

def fetch_components(tfaac, nuis=None, ud=None):
    suf = f"_{nuis}{ud}" if nuis and ud else ""
    names_try = [
        ("sm" + suf, f"sm_lin_quad" + suf, f"quad" + suf),
        ("Sm" + suf, "SmLinQuad" + suf, "Quad" + suf),
        ("sm" + suf, "sm_lin_quad_FT8" + suf, "quad_FT8" + suf),
    ]
    for sm_name, slq_name, q_name in names_try:
        h_sm  = get_hist(tfaac, sm_name)
        h_slq = get_hist(tfaac, slq_name)
        h_q   = get_hist(tfaac, q_name)
        if h_sm and h_slq and h_q:
            return h_sm, h_slq, h_q
    raise RuntimeError(f"Cannot find AAC components (nuis={nuis}, ud={ud}).")

def build_pred(tfaac, c, nuis=None, ud=None):
    h_sm, h_slq, h_q = fetch_components(tfaac, nuis, ud)
    pred = h_sm.Clone("pred")
    pred.Scale(1.0 - c)
    pred.Add(h_slq, c)
    pred.Add(h_q, c*c - c)
    pred.SetDirectory(0)
    return pred

def compute_metrics(h_data, h_pred):
    chi2, ndf = 0.0, 0
    for i in range(1, h_data.GetNbinsX()+1):
        y  = h_data.GetBinContent(i)
        ey = h_data.GetBinError(i)
        m  = h_pred.GetBinContent(i)
        if ey <= 0: continue
        chi2 += (y - m)**2 / (ey**2)
        ndf  += 1
    chi2_ndf = chi2 / ndf if ndf > 0 else float("nan")

    int_data = h_data.Integral()
    int_pred = h_pred.Integral()
    rel_int  = (int_data - int_pred) / int_pred if int_pred != 0 else float("inf")

    max_rel = 0.0
    for i in range(1, h_data.GetNbinsX()+1):
        y = h_data.GetBinContent(i)
        m = h_pred.GetBinContent(i)
        if m != 0:
            r = abs(y - m) / abs(m)
            if r > max_rel: max_rel = r
    return chi2_ndf, rel_int, max_rel

def styled(h, color=None, mstyle=20):
    if color is not None:
        h.SetLineColor(color)
        h.SetMarkerColor(color)
    h.SetMarkerStyle(mstyle)
    return h

def draw_compare(h_data, h_pred, title, outpng,
                 ratio_ylim=(0.6, 1.4),  # ratio 的 y 轴范围
                 legend_pos=(0.62, 0.72, 0.88, 0.88),  # (x1,y1,x2,y2)
                 grid_ratio=True):
    import ROOT as R

    # 全局样式
    R.gStyle.SetOptStat(0)

    # 计算统一的 y 轴上限（主图）
    ymax = max(h_data.GetMaximum(), h_pred.GetMaximum())
    h_pred.SetMaximum(ymax * 1.25)

    # 画布 & 两个 pad
    c = R.TCanvas("c","c",900,700)
    pad1 = R.TPad("pad1","main", 0,0.30, 1,1.00)
    pad2 = R.TPad("pad2","ratio",0,0.00, 1,0.30)
    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.06)
    pad2.SetBottomMargin(0.35)
    pad1.Draw(); pad2.Draw()

    # ------------ 主图 ------------
    pad1.cd()

    # 线/点样式
    h_pred.SetLineColor(R.kRed+1)
    h_pred.SetLineWidth(2)
    h_pred.SetTitle(title)
    h_pred.GetYaxis().SetTitle("Events")

    # 隐藏主图的 X 轴刻度与标题（只保留 ratio 的）
    h_pred.GetXaxis().SetLabelSize(0)
    h_pred.GetXaxis().SetTitleSize(0)

    # 字体尺寸（主图 Y 轴）
    h_pred.GetYaxis().SetTitleSize(0.05)
    h_pred.GetYaxis().SetLabelSize(0.045)

    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(1.0)
    h_data.SetLineColor(R.kBlack)
    h_data.SetMarkerColor(R.kBlack)

    # 先画预测，再叠数据
    h_pred.Draw("HIST")
    h_data.Draw("E1 SAME")

    # 图例
    x1,y1,x2,y2 = legend_pos
    leg = R.TLegend(x1,y1,x2,y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_data, "Original MC", "lep")
    leg.AddEntry(h_pred, "AAC pred",   "l")
    leg.Draw()

    # ------------ ratio ------------
    pad2.cd()
    ratio = h_data.Clone("ratio"); ratio.Divide(h_pred)
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("MC / pred")
    # ratio.GetXaxis().SetTitle("MT[GeV]")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetRangeUser(*ratio_ylim)

    # 轴样式（ratio）
    ratio.GetYaxis().SetTitleSize(0.11)
    ratio.GetYaxis().SetLabelSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.45)
    # ratio.GetXaxis().SetTitle(h_data.GetXaxis().GetTitle())
    ratio.GetXaxis().SetTitle("MT[GeV]")
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetLabelSize(0.10)

    ratio.Draw("E1")

    # y=1 参考线
    x_min = ratio.GetXaxis().GetXmin()
    x_max = ratio.GetXaxis().GetXmax()
    line = R.TLine(x_min,1.0,x_max,1.0)
    line.SetLineStyle(2)
    line.Draw()

    if grid_ratio:
        pad2.SetGridy(True)

    c.SaveAs(outpng)
    c.Close()

def arrange_validation_grid(png_infos, outpng, cols=3, outpdf=None):
    """把多个 png 摆放到一个大画布里。png_infos=[(label, path), ...]
       outpng 必填，outpdf 可选。
    """
    import ROOT as R

    if not png_infos:
        print("[INFO] no png found for grid composing, skip.")
        return

    first_img = R.TImage.Open(png_infos[0][1])
    if not first_img:
        print(f"[WARN] cannot open first image {png_infos[0][1]}, skip grid.")
        return
    w, h = first_img.GetWidth(), first_img.GetHeight()
    if w <= 0 or h <= 0:
        print(f"[WARN] invalid image size for {png_infos[0][1]}, skip grid.")
        return

    rows = math.ceil(len(png_infos) / cols)
    canvas = R.TCanvas("c_grid", "ValidationGrid", int(w * cols), int(h * rows))
    canvas.Divide(cols, rows)

    for idx, (label, path) in enumerate(png_infos):
        pad = canvas.cd(idx + 1)
        pad.SetFillColor(0)
        pad.SetMargin(0, 0, 0, 0)

        img = R.TImage.Open(path)
        if not img:
            print(f"[WARN] cannot open {path}, leave blank in grid.")
            continue
        img.Draw()

        latex = R.TLatex(0.05, 0.92, label)
        latex.SetNDC(True)
        latex.SetTextSize(0.06)
        latex.Draw()

    canvas.Update()
    canvas.SaveAs(outpng)
    if outpdf:
        try:
            from PIL import Image
            Image.open(outpng).convert("RGB").save(outpdf, "PDF")
        except Exception as e:
            print(f"[WARN] cannot convert {outpng} -> {outpdf}: {e}")
            outpdf = None
    canvas.Close()
    msg = f"[DONE] grid canvas saved => {outpng}"
    if outpdf:
        msg += f" & {outpdf}"
    print(msg)


# ========= 主程 =========
tfaac = R.TFile.Open(aac_file)
if not tfaac or tfaac.IsZombie():
    raise RuntimeError(f"Cannot open {aac_file}")

nuisances = []
if not ONLY_NOMINAL:
    nuisances = list_nuis_from_aac(tfaac)

for cval in c_points:
    if abs(cval) < 1e-12:
        data_file = sm_file
        tagname = "0"
    else:
        token = c_to_token(cval)
        data_file = f"shapes-vbs-SR{tag}FT8_{token}.root"
        tagname = token

    tfdata = R.TFile.Open(data_file)
    if not tfdata or tfdata.IsZombie():
        print(f"[WARN] cannot open {data_file}, skip this c")
        continue

    # -------- nominal ----------
    h_data = tfdata.Get(hist_name)
    if h_data:
        h_pred = build_pred(tfaac, cval, nuis=None, ud=None)
        ttl = f"FT8={cval:.3f} ; nominal"
        png = os.path.join(outdir, f"VBS_EFT_FT8_{tagname}.png")
        h_data = h_data.Clone(); h_data.SetDirectory(0)
        chi2_ndf, rel_int, max_rel = compute_metrics(h_data, h_pred)
        print(f"[METRIC] c={cval:+.3f}  nominal  chi2/ndf={chi2_ndf:.3f}  IntRel={rel_int:+.3%}  MaxBinRel={max_rel:.1%}")
        draw_compare(h_data, h_pred, ttl, png)

    # -------- systematics ----------
    if not ONLY_NOMINAL:
        for nuis in nuisances:
            for ud in ("Up","Down"):
                h_data_var = tfdata.Get(f"{hist_name}_{nuis}{ud}")
                if not h_data_var:
                    print(f"[INFO] missing {hist_name}_{nuis}{ud} in {data_file}, skip")
                    continue
                h_pred_var = build_pred(tfaac, cval, nuis=nuis, ud=ud)
                ttl = f"c={cval:.3f} ; {nuis}{ud}"
                png = os.path.join(outdir, f"cmp_{nuis}{ud}_c{tagname}.png")
                h_d = h_data_var.Clone(); h_d.SetDirectory(0)
                chi2_ndf, rel_int, max_rel = compute_metrics(h_d, h_pred_var)
                print(f"[METRIC] c={cval:+.3f}  {nuis}{ud}  chi2/ndf={chi2_ndf:.3f}  IntRel={rel_int:+.3%}  MaxBinRel={max_rel:.1%}")
                draw_compare(h_d, h_pred_var, ttl, png)

    tfdata.Close()

tfaac.Close()

# 组合 0~8 的 nominal 图到一个大画布里
grid_png_infos = []
for cv in grid_c_points:
    token = token_from_c(cv)
    png_path = os.path.join(outdir, f"VBS_EFT_FT8_{token}.png")
    if os.path.exists(png_path):
        grid_png_infos.append((f"FT8={cv:g}", png_path))
    else:
        print(f"[INFO] grid compose skip: {png_path} not found.")
if grid_png_infos:
    grid_out = os.path.join(outdir, "VBS_EFT_grid_0_to_8.png")
    grid_out_pdf = os.path.join(outdir, "VBS_EFT_grid_0_to_8.pdf")
    arrange_validation_grid(grid_png_infos, grid_out, cols=3, outpdf=grid_out_pdf)

print(f"[DONE] Plots => {outdir}")
