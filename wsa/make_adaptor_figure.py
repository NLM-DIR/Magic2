#!/usr/bin/env python3
"""
make_adaptor_figure.py  --  Magic2 adaptor-profile visualiser
=============================================================
Reads a Magic2 overhang TSF file and writes a two-panel PDF
(and PNG) showing the per-position base-frequency line plots
for Read 1 and Read 2.

Usage (standalone):
    python3 make_adaptor_figure.py <overhang.tsf>

Called from Magic2 C code:
    snprintf(cmd, sizeof(cmd),
             "python3 ~/ace/wsa/make_adaptor_figure.py %s", tsfPath);
    system(cmd);

Output:  <stem>.adaptor_histogram.pdf
         <stem>.adaptor_histogram.png
Dependencies: matplotlib, numpy
"""

import sys
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator

# ── colour palette (ACGT standard, colour-blind friendly) ────────────────────
COL   = {"A": "#4CAF50", "T": "#F44336", "G": "#FF9800", "C": "#2196F3",
         "N": "#BDBDBD"}
BASES = ["A", "T", "G", "C", "N"]

# ── TSF parser (unchanged) ────────────────────────────────────────────────────
def parse_tsf(path):
    blocks = []
    cur = None

    def new_block(label):
        return {"label": label, "adaptor_r1": "", "adaptor_r2": "",
                "r1": [], "r2": []}

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                m_a1 = re.search(r"Adaptor[_\s]*R1\s+(\S+)", line, re.I)
                m_a2 = re.search(r"Adaptor[_\s]*R2\s+(\S+)", line, re.I)
                if m_a1 or m_a2:
                    if cur is None:
                        cur = new_block("overhang"); blocks.append(cur)
                    if m_a1: cur["adaptor_r1"] = m_a1.group(1)
                    if m_a2: cur["adaptor_r2"] = m_a2.group(1)
                else:
                    m_label = re.search(
                        r"(3\s*prime|5\s*prime|overhang[^\t]*)", line, re.I)
                    label = m_label.group(1).strip() if m_label else "overhang"
                    cur = new_block(label); blocks.append(cur)
                continue
            if line.startswith("#"):
                continue
            if cur is None:
                cur = new_block("overhang"); blocks.append(cur)
            fields = line.split("\t")
            if len(fields) < 13:
                continue
            try:
                pos   = int(fields[1])
                pct_a = float(fields[9])
                pct_t = float(fields[10])
                pct_g = float(fields[11])
                pct_c = float(fields[12])
                pct_n = float(fields[13]) if len(fields) > 13 else 0.0
            except (ValueError, IndexError):
                continue
            row = (pos, pct_a, pct_t, pct_g, pct_c, pct_n)
            if fields[0].endswith("_R2"):
                cur["r2"].append(row)
            else:
                cur["r1"].append(row)

    for b in blocks:
        b["r1"].sort(key=lambda x: x[0])
        b["r2"].sort(key=lambda x: x[0])
    return [b for b in blocks if b["r1"] or b["r2"]]


def rows_to_arrays(rows):
    """Return (positions, dict base->np.array of pct)."""
    positions = np.array([r[0] for r in rows])
    pct = {b: np.array([rows[j][i+1] for j in range(len(rows))])
           for i, b in enumerate(BASES)}
    return positions, pct


# ── panel plotter ─────────────────────────────────────────────────────────────
def plot_panel(ax, rows, adaptor_seq, title, ylabel=True):
    if not rows:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=10, color="grey")
        ax.set_title(title, fontsize=9, fontweight="bold")
        return

    positions, pct = rows_to_arrays(rows)
    n = len(positions)

    for base in ["A", "T", "G", "C"]:          # skip N (near-zero, clutters)
        y = pct[base]
        color = COL[base]
        ax.plot(positions, y,
                color=color, linewidth=1.6,
                marker="o", markersize=4, markerfacecolor=color,
                markeredgewidth=0, zorder=3, label=base)

    # Inferred adaptor sequence just above the top axis
    if adaptor_seq:
        ax.text(positions[0] - 0.8, 103, "\u2192",
                ha="right", va="bottom", fontsize=6,
                color="#888888", style="italic",
                transform=ax.transData)
        for j, ch in enumerate(adaptor_seq[:n]):
            col = COL.get(ch.upper(), "#888888")
            ax.text(positions[j], 103, ch.upper(),
                    ha="center", va="bottom", fontsize=6.5,
                    fontweight="bold", color=col, zorder=5)

    ax.set_xlim(positions[0] - 0.8, positions[-1] + 0.8)
    ax.set_ylim(0, 110)
    ax.set_xticks(positions)
    ax.set_xticklabels(positions, fontsize=7)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(axis="y", labelsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold", pad=4)
    if ylabel:
        ax.set_ylabel("Base frequency (%)", fontsize=8)
    ax.set_xlabel("Overhang position (nt)", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    ax.grid(axis="x", linestyle=":", linewidth=0.3, alpha=0.4, zorder=0)


# ── figure builder ────────────────────────────────────────────────────────────
def make_figure(block, out_stem, run_label=""):
    has_r1 = bool(block["r1"])
    has_r2 = bool(block["r2"])
    ncols  = int(has_r1) + int(has_r2)
    if ncols == 0:
        print("  [make_adaptor_figure] block '%s': no data, skipping"
              % block["label"])
        return

    fig, axes = plt.subplots(1, ncols,
                             figsize=(5.5 * ncols + 0.5, 4.0),
                             gridspec_kw={"wspace": 0.30})
    if ncols == 1:
        axes = [axes]

    col = 0
    if has_r1:
        plot_panel(axes[col], block["r1"], block["adaptor_r1"],
                   title="(a)  Read 1", ylabel=True)
        col += 1
    if has_r2:
        plot_panel(axes[col], block["r2"], block["adaptor_r2"],
                   title="(b)  Read 2", ylabel=(ncols == 1))

    # Shared legend
    handles = [mpatches.Patch(color=COL[b], label=b) for b in "ATGC"]
    fig.legend(handles=handles, title="Base", title_fontsize=8,
               fontsize=8, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.04), frameon=False)

    run_str = ("run %s  \u2014  " % run_label) if run_label else ""
    fig.suptitle(
        "Per-position base-frequency of 3\u2032 unaligned overhangs\n"
        "(%s%s; composited after splitting by first overhang base)"
        % (run_str, block["label"]),
        fontsize=9, y=1.06)

    for ext in (".pdf", ".png"):
        fig.savefig(out_stem + ext, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("  [make_adaptor_figure] wrote %s.pdf  %s.png" % (out_stem, out_stem))


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 2:
        sys.exit(
            "Usage: python3 make_adaptor_figure.py <overhang.tsf>\n"
            "       PDF and PNG are written next to the TSF file."
        )
    tsf_path = sys.argv[1]
    if not os.path.isfile(tsf_path):
        sys.exit("[make_adaptor_figure] ERROR: file not found: %s" % tsf_path)

    base_stem = tsf_path[:-4] if tsf_path.endswith(".tsf") else tsf_path
    fname     = os.path.basename(base_stem)
    run_label = fname.split(".")[0] if "." in fname else fname

    print("[make_adaptor_figure] parsing %s" % tsf_path)
    blocks = parse_tsf(tsf_path)
    if not blocks:
        sys.exit("[make_adaptor_figure] ERROR: no data found in %s" % tsf_path)

    for i, block in enumerate(blocks):
        suffix = ".adaptor_histogram" if len(blocks) == 1 \
                 else (".adaptor_histogram_%d" % (i + 1))
        make_figure(block, base_stem + suffix, run_label=run_label)


if __name__ == "__main__":
    main()
