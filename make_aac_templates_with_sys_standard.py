#!/usr/bin/env python3
import ROOT as R
import re, glob
import numpy as np

R.gROOT.SetBatch(True)

# ============= 配置 =============
TAGS = ("2016APV", "2016", "2017", "2018")
OPERATORS = (
    "FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7",
    "FS0", "FS1", "FS2",
    "FT0", "FT1", "FT2", "FT3", "FT4", "FT5", "FT6", "FT7", "FT8", "FT9"
)

# OPERATORS = (
#     "FT8",
# )

# 拟合时使用的小 |c| 区间（按算符调整，默认兜底）
FIT_ABS_MAX = 128.0
FIT_ABS_MAX_BY_OP = {
    # FM / FS operators tend to have wider limits
    "FM0": 10.0,
    "FM1": 30.0,
    "FM2": 10.0,
    "FM3": 40.0,
    "FM4": 30.0,
    "FM5": 40.0,
    "FM7": 60.0,
    "FS0": 50.0,
    "FS1": 40.0,
    "FS2": 50.0,
    # FT operators are typically narrower; keep FT6 tighter
    "FT0": 2.0,
    "FT1": 2.0,
    "FT2": 5.0,
    "FT3": 5.0,
    "FT4": 10.0,
    "FT5": 5.0,
    "FT6": 8.0,
    "FT7": 15.0,
    "FT8": 5.0,
    "FT9": 8.0,
}
# 如需强制额外保留的点（例如 ±4），可写：FIT_KEEP = {-4.0, 4.0}
FIT_KEEP = set()
# 输入/输出
hist_name = "eft_ZZ2l2nu"                       # 输入直方图基名

# ============= 工具函数 =============
def parse_c_from_filename(operator, fn):
    m = re.search(rf"{re.escape(operator)}_([^./]+)\.root$", fn)
    if not m:
        raise RuntimeError(f"Cannot parse {operator} value from {fn}")
    s = m.group(1)
    neg = s.startswith("m")
    s2  = s[1:] if neg else s
    s2  = s2.replace("p",".")
    val = float(s2)
    return -val if neg else val

def hist_to_arrays(h):
    nb = h.GetNbinsX()
    y  = np.array([h.GetBinContent(i) for i in range(1, nb+1)], dtype=float)
    e2 = np.array([h.GetBinError(i)**2 for i in range(1, nb+1)], dtype=float)
    e2[e2 <= 0] = 1e-12
    return y, e2

def list_nuisances_from_sm(sm_tfile, hist_base):
    nuis = set()
    for key in sm_tfile.GetListOfKeys():
        name = key.GetName()
        m = re.match(rf"^{re.escape(hist_base)}_(.+)(Up|Down)$", name)
        if m:
            nuis.add(m.group(1))
    return sorted(nuis)

def fit_quadratic_per_bin(c_list, y_list, e2_list):
    cs  = np.array(c_list, dtype=float)
    nb  = y_list[0].shape[0]
    A   = np.vstack([np.ones_like(cs), cs, cs**2]).T  # (npt, 3)

    a = np.zeros(nb); b = np.zeros(nb); d = np.zeros(nb)
    Va = np.zeros(nb); Vb = np.zeros(nb); Vd = np.zeros(nb)
    Cab = np.zeros(nb); Cad = np.zeros(nb); Cbd = np.zeros(nb)

    for i in range(nb):
        y  = np.array([arr[i] for arr in y_list], dtype=float)
        e2 = np.array([arr[i] for arr in e2_list], dtype=float)
        w  = 1.0 / e2

        AtWA = A.T @ (w[:, None] * A)
        AtWy = A.T @ (w * y)

        coeff = np.linalg.solve(AtWA, AtWy)
        a[i], b[i], d[i] = coeff

        # 协方差矩阵（加权最小二乘的标准近似）
        cov = np.linalg.inv(AtWA)
        Va[i], Vb[i], Vd[i] = cov[0,0], cov[1,1], cov[2,2]
        Cab[i], Cad[i], Cbd[i] = cov[0,1], cov[0,2], cov[1,2]

    tiny_mask = (d < 0) & (d > -1e-10)
    d[tiny_mask] = 0.0
    return a, b, d, Va, Vb, Vd, Cab, Cad, Cbd

def write_three_hists(proto_hist, sm_name, quad_name, slq_name,
                      a, b, d, Va, Vb, Vd, Cab, Cad, Cbd, tfile_out):
    nb = proto_hist.GetNbinsX()
    h_sm   = proto_hist.Clone(sm_name);   h_sm.SetDirectory(0);   h_sm.Reset();   h_sm.Sumw2(True)
    h_quad = proto_hist.Clone(quad_name); h_quad.SetDirectory(0); h_quad.Reset(); h_quad.Sumw2(True)
    h_slq  = proto_hist.Clone(slq_name);  h_slq.SetDirectory(0);  h_slq.Reset();  h_slq.Sumw2(True)

    eps = 1e-12  # 防止0误差

    for i in range(1, nb+1):
        ai, bi, di = a[i-1], b[i-1], d[i-1]
        # Vai, Vbi, Vdi = Va[i-1], Vb[i-1], Vd[i-1]
        # Cabi, Cadi, Cbd_i = Cab[i-1], Cad[i-1], Cbd[i-1]

        #var_slq = Vai + Vbi + Vdi + 2*Cabi + 2*Cadi + 2*Cbd_i

        h_sm.SetBinContent(i,   ai)
        h_quad.SetBinContent(i, di)
        h_slq.SetBinContent(i,  ai + bi + di)

        # h_sm.SetBinError(i,   max(np.sqrt(max(Vai,    0.0)), eps))
        # h_quad.SetBinError(i, max(np.sqrt(max(Vdi,    0.0)), eps))
        # h_slq.SetBinError(i,  max(np.sqrt(max(var_slq,0.0)), eps))

        h_sm.SetBinError(i,   0.0)
        h_quad.SetBinError(i, 0.0)
        h_slq.SetBinError(i,  0.0)

    tfile_out.cd()
    h_sm.Write(); h_quad.Write(); h_slq.Write()

def select_pairs_for_fit(pairs, abs_max, keep=None):
    """pairs: [(c, y, e2), ...]  -> 仅保留 |c|<=abs_max 或在 keep 集合内的点"""
    keep = keep or set()
    sel = [(c,y,e2) for (c,y,e2) in pairs if (abs(c) <= abs_max) or (c in keep)]
    # 按 |c| 排序，便于打印
    sel = sorted(sel, key=lambda t: (abs(t[0]), t[0]))
    return sel

def build_templates_for(tag, operator):
    sm_proc  = "sm"
    slq_proc = f"sm_lin_quad_{operator}"
    quad_proc = f"quad_{operator}"

    sm_file  = f"shapes-vbs-SR{tag}sm.root"            # SM 文件（含形变）
    sig_glob = f"shapes-vbs-SR{tag}{operator}_*.root"  # 非零算符文件
    out_file = f"aac_{operator}_templates_{tag}.root"  # 输出

    print(f"[INFO] ===== tag={tag}, operator={operator} =====")

    # 读 SM 名义
    f_sm = R.TFile.Open(sm_file)
    if not f_sm or f_sm.IsZombie():
        print(f"[WARN] Cannot open {sm_file}, skip {tag}/{operator}")
        return

    h_sm_nom = f_sm.Get(hist_name)
    if not h_sm_nom:
        print(f"[WARN] Histogram {hist_name} not found in {sm_file}, skip {tag}/{operator}")
        f_sm.Close()
        return

    y_sm_nom, e2_sm_nom = hist_to_arrays(h_sm_nom)
    nbins = h_sm_nom.GetNbinsX()

    # 枚举全部算符样本（名义）
    sig_files = sorted(glob.glob(sig_glob))
    if not sig_files:
        print(f"[WARN] No files matched {sig_glob}, skip {tag}/{operator}")
        f_sm.Close()
        return

    pairs_all = []
    for fn in sig_files:
        try:
            c = parse_c_from_filename(operator, fn)
        except RuntimeError as exc:
            print(f"[WARN] {exc}, skip {fn}")
            continue

        if abs(c) < 1e-12:  # 跳过 0
            continue
        f = R.TFile.Open(fn)
        if not f or f.IsZombie():
            print(f"[WARN] cannot open {fn}, skip")
            continue
        h = f.Get(hist_name)
        if not h:
            print(f"[WARN] {hist_name} not in {fn}, skip")
            f.Close()
            continue
        if h.GetNbinsX() != nbins:
            print(f"[WARN] bin mismatch in {fn}, skip")
            f.Close()
            continue
        y, e2 = hist_to_arrays(h)
        pairs_all.append((c, y, e2))
        f.Close()

    if len(pairs_all) < 2:
        print(f"[WARN] {tag}/{operator}: 至少需要两个非零的点，当前 {len(pairs_all)} 个，跳过。")
        f_sm.Close()
        return

    # 选择用于名义拟合的子集（按算符范围）
    fit_abs_max = FIT_ABS_MAX_BY_OP.get(operator, FIT_ABS_MAX)
    pairs_fit = select_pairs_for_fit(pairs_all, fit_abs_max, FIT_KEEP)
    if len(pairs_fit) < 2:
        print(f"[WARN] {tag}/{operator}: 用于拟合的非零点不足（{len(pairs_fit)}），请调整 FIT_ABS_MAX/FIT_KEEP。")
        f_sm.Close()
        return
    print("[INFO] Fit (nominal) abs(c) <= ", fit_abs_max)
    print("[INFO] Fit (nominal) will use c points:", [round(c,6) for c,_,_ in pairs_fit])

    # 输出文件
    fout = R.TFile(out_file, "RECREATE")

    # ===== 名义：拟合并写出 Sm / Quad / SmLinQuad =====
    c_list  = [0.0] + [c for c,_,_ in pairs_fit]
    y_list  = [y_sm_nom] + [y for _,y,_ in pairs_fit]
    e2_list = [e2_sm_nom] + [e2 for _,_,e2 in pairs_fit]
    a_nom, b_nom, d_nom, Va, Vb, Vd, Cab, Cad, Cbd = fit_quadratic_per_bin(c_list, y_list, e2_list)
    write_three_hists(
        h_sm_nom,
        sm_name=sm_proc,
        quad_name=quad_proc,
        slq_name=slq_proc,
        a=a_nom, b=b_nom, d=d_nom,
        Va=Va, Vb=Vb, Vd=Vd, Cab=Cab, Cad=Cad, Cbd=Cbd,
        tfile_out=fout
    )
    print(f"[OK] Wrote nominal: {sm_proc}, {quad_proc}, {slq_proc}")

    # categories 顺序必须与 datacard 写的 process 顺序一致
    cats = [sm_proc, slq_proc, quad_proc]   # ["sm", "sm_lin_quad", "quad"]

    # 创建一个 3D histogram：X = bin index, Y = template 1, Z = template 2
    h_corr = R.TH3D(
        "histo_correlation",
        "correlated mc stat uncertainty; xbin; procY; procZ",
        nbins, 0.5, nbins + 0.5,        # X axis: bins
        len(cats), 0.5, len(cats) + 0.5, # Y axis: template categories
        len(cats), 0.5, len(cats) + 0.5  # Z axis
    )

    # 给 Y/Z 轴贴 bin label（可选，但调试很有用）
    for i, cat in enumerate(cats, start=1):
        h_corr.GetYaxis().SetBinLabel(i, cat)
        h_corr.GetZaxis().SetBinLabel(i, cat)

    # 逐 bin 填协方差矩阵
    for ib in range(1, nbins + 1):
        Va_i, Vb_i, Vd_i = Va[ib-1], Vb[ib-1], Vd[ib-1]
        Cab_i, Cad_i, Cbd_i = Cab[ib-1], Cad[ib-1], Cbd[ib-1]

        # 计算 Var/Cov
        VarS = Va_i
        VarQ = Vd_i
        VarL = Va_i + Vb_i + Vd_i + 2*(Cab_i + Cad_i + Cbd_i)
        CovSL = Va_i + Cab_i + Cad_i    # sm & sm_lin_quad
        CovSQ = Cad_i                   # sm & quad
        CovLQ = Vd_i + Cbd_i + Cad_i    # sm_lin_quad & quad

        # 写入对称矩阵
        def setcov(y,z,val):
            h_corr.SetBinContent(ib, y, z, val)
            h_corr.SetBinContent(ib, z, y, val)

        # 对角
        setcov(1,1,VarS)
        setcov(2,2,VarL)
        setcov(3,3,VarQ)

        # 非对角（协方差）
        setcov(1,2,CovSL)
        setcov(1,3,CovSQ)
        setcov(2,3,CovLQ)

    # 写到输出 ROOT
    fout.cd()
    h_corr.Write()

    # ===== 系统误差：逐形变（Up/Down）重建 =====
    nuis_list = list_nuisances_from_sm(f_sm, hist_name)
    print(f"[INFO] Found {len(nuis_list)} nuisances in SM file.")

    for nuis in nuis_list:
        for ud in ("Up","Down"):
            h_sm_var = f_sm.Get(f"{hist_name}_{nuis}{ud}")
            if not h_sm_var:
                print(f"[WARN] SM variation missing: {hist_name}_{nuis}{ud}, skip this side.")
                continue
            y_sm_var, e2_sm_var = hist_to_arrays(h_sm_var)

            # 收集该形变在各算符文件上的直方图
            pairs_var_all = []
            for fn in sig_files:
                try:
                    c = parse_c_from_filename(operator, fn)
                except RuntimeError:
                    continue
                if abs(c) < 1e-12:
                    continue
                f = R.TFile.Open(fn)
                if not f or f.IsZombie():
                    print(f"[WARN] {nuis}{ud} open failed @ c={c}")
                    continue
                h_var = f.Get(f"{hist_name}_{nuis}{ud}")
                if not h_var:
                    f.Close()
                    continue
                if h_var.GetNbinsX()!=nbins:
                    print(f"[WARN] bin mismatch in {fn} for {nuis}{ud}")
                    f.Close()
                    continue
                yv, e2v = hist_to_arrays(h_var)
                pairs_var_all.append((c, yv, e2v))
                f.Close()

            # 选择用于拟合的子集（同样只用小 |c|）
            pairs_var_fit = select_pairs_for_fit(pairs_var_all, fit_abs_max, FIT_KEEP)
            used = len(pairs_var_fit)
            if used < 2:
                print(f"[WARN] {tag}/{operator} NUIS={nuis}{ud}: 有效 c 点少于 2（used={used}），跳过该形变重建。")
                continue
            print(f"[INFO] Fit ({nuis}{ud}) uses c:", [round(c,6) for c,_,_ in pairs_var_fit])

            c_list_var  = [0.0] + [c for c,_,_ in pairs_var_fit]
            y_list_var  = [y_sm_var] + [y for _,y,_ in pairs_var_fit]
            e2_list_var = [e2_sm_var] + [e2 for _,_,e2 in pairs_var_fit]

            a_var, b_var, d_var, Va_var, Vb_var, Vd_var, Cab_var, Cad_var, Cbd_var = fit_quadratic_per_bin(c_list_var, y_list_var, e2_list_var)

            # 输出名必须符合 $PROCESS_$SYSTEMATIC 约定
            write_three_hists(
                h_sm_nom,
                sm_name=f"{sm_proc}_{nuis}{ud}",
                quad_name=f"{quad_proc}_{nuis}{ud}",
                slq_name=f"{slq_proc}_{nuis}{ud}",
                a=a_var, b=b_var, d=d_var,
                Va=Va_var, Vb=Vb_var, Vd=Vd_var, Cab=Cab_var, Cad=Cad_var, Cbd=Cbd_var,
                tfile_out=fout
            )
            print(f"[OK] Wrote {nuis}{ud}: {sm_proc}_{nuis}{ud}, {quad_proc}_{nuis}{ud}, {slq_proc}_{nuis}{ud}")

    fout.Close()
    f_sm.Close()
    print(f"[DONE] Output => {out_file}")


if __name__ == "__main__":
    for tag in TAGS:
        for operator in OPERATORS:
            build_templates_for(tag, operator)
