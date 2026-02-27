#!/usr/bin/env python3
"""Plot all nominal signal/background shapes on one canvas (CMS-style labels, tuned font sizes)."""
import argparse
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep


def parse_args():
    ap = argparse.ArgumentParser(description="Plot stacked backgrounds + FT8 signal overlay.")
    ap.add_argument(
        "-i",
        "--input",
        nargs="+",
        default=[
            "shapes-vbs-SR2016APVsm.root",
            "shapes-vbs-SR2016sm.root",
            "shapes-vbs-SR2017sm.root",
            "shapes-vbs-SR2018sm.root",
        ],
        help="Input ROOT files (backgrounds); same process names will be summed",
    )
    ap.add_argument("-o", "--output", default="nominals_all.png", help="Output image")
    ap.add_argument("--logy", action="store_true", help="Use log scale")
    ap.add_argument("--exclude", nargs="*", default=None, help="Exclude background processes")
    ap.add_argument("--signal", default="eft_ZZ2l2nu", help="Signal process to overlay")
    ap.add_argument("--sm-signal", default="eft_ZZ2l2nu", help="SM signal process from input file")
    ap.add_argument(
        "--ft8-file",
        nargs="+",
        default=[
            "shapes-vbs-SR2016APVFT8_4.root",
            "shapes-vbs-SR2016FT8_4.root",
            "shapes-vbs-SR2017FT8_4.root",
            "shapes-vbs-SR2018FT8_4.root",
        ],
        help="FT8=4 root files; same process names will be summed",
    )
    ap.add_argument("--ft8-tag", default="4", help="Label for FT8 signal")
    ap.add_argument(
        "--palette-yaml",
        default="/afs/cern.ch/user/l/lotan/smqawa/DCTools/config/input_UL_2018-fes.yaml",
        help="YAML file with group color definitions (kept for compatibility)",
    )

    # --- CMS label and font-size tuning (new) ---
    ap.add_argument("--lumi", default="138", help='Integrated luminosity in fb^-1, e.g. "138"')
    ap.add_argument("--com", type=float, default=13.0, help="Center-of-mass energy in TeV")
    ap.add_argument(
        "--cms-extra",
        default="Simulation Preliminary",
        help='CMS extra text (right plot style), e.g. "Simulation Preliminary" or "Preliminary"',
    )
    ap.add_argument("--cms-loc", type=int, default=0, help="0=UL, 1=UR, 2=LL, 3=LR (mplhep)")

    # match the right figure typography
    ap.add_argument("--cms-font", type=float, default=26, help="Font size for 'CMS'")
    ap.add_argument("--extra-font", type=float, default=20, help="Font size for extra text (Simulation Preliminary)")
    ap.add_argument("--lumi-font", type=float, default=20, help="Font size for lumi/com text")
    ap.add_argument("--tick-font", type=float, default=18, help="Tick label font size")
    ap.add_argument("--axis-font", type=float, default=26, help="Axis label font size (MT/Events)")
    ap.add_argument("--fig-w", type=float, default=10, help="Figure width")
    ap.add_argument("--fig-h", type=float, default=7, help="Figure height (7 recommended to match right plot)")
    ap.add_argument("--cms-pad", type=float, default=0.01, help="Padding above axes for CMS header")
    return ap.parse_args()


def list_keys(root_file):
    return [k.split(";")[0] for k in root_file.keys()]


def find_nominals(keys):
    procs = []
    for k in keys:
        if k == "data_obs":
            continue
        if k.endswith("Up") or k.endswith("Down"):
            continue
        procs.append(k)
    return sorted(set(procs))


def hist_values(root_file, name):
    h = root_file[name]
    values, edges = h.to_numpy()
    return values.astype(float), edges.astype(float)


def sum_processes(files, process):
    total_vals = None
    edges_ref = None
    for f in files:
        if process not in f:
            continue
        vals, edges = hist_values(f, process)
        if total_vals is None:
            total_vals = vals.copy()
            edges_ref = edges
        else:
            if not np.array_equal(edges_ref, edges):
                raise ValueError(f"Bin edges mismatch for {process}")
            total_vals += vals
    return total_vals, edges_ref


def load_colors_from_yaml(_path):
    # Hard-coded to match input_UL_2018-fes.yaml groups color scheme
    return {
        "WW": "#99B6F7",
        "WZ": "#46B3A5",
        "DY": "#F6D68D",
        "Top": "#2E6D92",
        "VVV": "#580C82",
        "ggVV": "#580CCC",
        "ZZ2l2nu": "#8E31A1",
        "ekw_NLO": "#FE4773",
    }


def apply_cms_header_sizes(ax, cms_font, extra_font, lumi_font):
    """
    mplhep.cms.label creates multiple Text objects in ax.texts.
    We post-adjust their font sizes to match the "right plot" style.
    """
    for t in ax.texts:
        s = t.get_text()

        # The "CMS" token is usually a standalone text
        if s.strip() == "CMS":
            t.set_fontsize(cms_font)
            t.set_fontweight("bold")
            continue

        # extra text like "Simulation Preliminary" is usually placed next to CMS
        if "Simulation" in s or "Preliminary" in s:
            # guard: avoid also catching the lumi line
            if "fb" not in s and "TeV" not in s:
                t.set_fontsize(extra_font)
                t.set_fontstyle("italic")  # like typical CMS extraText style
                continue

        # lumi/com text: contains fb^{-1} or TeV
        if "fb" in s or "TeV" in s:
            t.set_fontsize(lumi_font)
            continue


def main():
    args = parse_args()

    # CMS-like matplotlib style (fonts, ticks, etc.)
    plt.style.use(hep.style.CMS)

    input_paths = [Path(p) for p in args.input]
    for p in input_paths:
        if not p.exists():
            raise SystemExit(f"Input file not found: {p}")

    colors = load_colors_from_yaml(args.palette_yaml)

    with ExitStack() as stack:
        files = [stack.enter_context(uproot.open(str(p))) for p in input_paths]
        keys = sorted(set(k for f in files for k in list_keys(f)))
        procs = find_nominals(keys)

        if args.exclude:
            procs = [p for p in procs if p not in set(args.exclude)]

        # Remove signal(s) from background stack if present in input
        if args.signal in procs:
            procs = [p for p in procs if p != args.signal]
        if args.sm_signal in procs:
            procs = [p for p in procs if p != args.sm_signal]

        # Enforce background stacking order to match reference plot
        stack_order = ["WW", "WZ", "DY", "Top", "VVV", "ggVV", "ZZ2l2nu", "ekw_NLO"]
        ordered = [p for p in stack_order if p in procs]
        remaining = [p for p in procs if p not in ordered]
        procs = ordered + remaining

        fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h))

        # Typography tuned to match your "right" plot
        ax.tick_params(labelsize=args.tick_font)
        ax.xaxis.label.set_size(args.axis_font)
        ax.yaxis.label.set_size(args.axis_font)

        # Fallback palette if a process is missing
        palette = [
            "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
            "#9C755F", "#EDC948", "#76B7B2", "#FF9DA7", "#BAB0AC",
        ]
        color_idx = 0
        cumulative = None
        max_y = 0.0
        edges = None

        # ---- stack backgrounds (filled) ----
        for proc in procs:
            vals, edges = sum_processes(files, proc)
            if vals is None:
                continue
            if cumulative is None:
                cumulative = np.zeros_like(vals)
            bottom = cumulative.copy()
            cumulative = cumulative + vals
            max_y = max(max_y, float(np.max(cumulative)))

            x = edges
            y0 = np.r_[bottom, bottom[-1]]
            y1 = np.r_[cumulative, cumulative[-1]]

            facecolor = colors.get(proc)
            if not facecolor:
                facecolor = palette[color_idx % len(palette)]
                color_idx += 1

            ax.fill_between(
                x, y0, y1,
                step="post",
                alpha=0.6,
                linewidth=1.0,
                edgecolor="black",
                linestyle="-",
                facecolor=facecolor,
                label=proc,
            )

        # ---- axes ----
        ax.set_xlabel("MT [GeV]")
        ax.set_xlim(100, 1000)
        ax.set_ylabel("Events")
        ax.grid(False)
        if args.logy:
            ax.set_yscale("log")

        # ---- CMS header using mplhep, then post-adjust sizes ----
        hep.cms.label(
            ax=ax,
            label="Preliminary",  # e.g. "Simulation Preliminary"
            data=False,            # your plot has no data points here
            lumi=args.lumi,
            com=args.com,
            loc=args.cms_loc,
            pad=args.cms_pad,
        )
        apply_cms_header_sizes(ax, args.cms_font, args.extra_font, args.lumi_font)

        # ---- SM (FT8=0) signal line from background file (not stacked) ----
        if args.sm_signal in keys:
            sm_vals, sm_edges = sum_processes(files, args.sm_signal)
            if sm_vals is not None:
                x = sm_edges
                y = np.r_[sm_vals, sm_vals[-1]]
                max_y = max(max_y, float(np.max(y)))
                ax.step(
                    x, y,
                    where="post",
                    color="#E53935",
                    linestyle="--",
                    linewidth=2.0,
                    label="ekw_NLO (FT8=0)",
                    zorder=9,
                )

        # ---- FT8=4 signal from separate files ----
        ft8_paths = [Path(p) for p in args.ft8_file]
        for p in ft8_paths:
            if not p.exists():
                raise SystemExit(f"FT8 file not found: {p}")

        with ExitStack() as stack2:
            ft8_files = [stack2.enter_context(uproot.open(str(p))) for p in ft8_paths]
            keys2 = sorted(set(k for f2 in ft8_files for k in list_keys(f2)))
            if args.signal in keys2:
                vals, edges = sum_processes(ft8_files, args.signal)
                if vals is not None:
                    x = edges
                    y = np.r_[vals, vals[-1]]
                    max_y = max(max_y, float(np.max(y)))
                    ax.step(
                        x, y,
                        where="post",
                        color="#E53935",
                        linestyle="-",
                        linewidth=2.5,
                        label=f"{args.signal} (FT8={args.ft8_tag})",
                        zorder=10,
                    )

        # add headroom so curves don't touch the top frame
        if max_y > 0:
            ax.set_ylim(top=max_y * 1.25)

        # legend tuned similar to your reference
        ax.legend(fontsize=11, ncol=2, loc="upper right", frameon=False)

        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
