#!/usr/bin/env python3
"""
Make per‑MUTANT mutation maps for bacterial genomes from an Excel table.

Input Excel is expected to have columns (case-sensitive):
  - "#CHROM"   (chromosome/contig accession)
  - "MUTANT"   (mutant/sample identifier)
  - "POS"      (1-based genomic position of mutation)
  - "%REF"     (percentage of reference allele; may use comma decimal separator)

Other columns are ignored. The Excel workbook may have a single sheet or multiple; by default
we read the first sheet unless --sheet is provided.

For each MUTANT, this script produces a horizontal line for every chromosome (present in that MUTANT),
with short vertical ticks at every mutation position on that chromosome.
Tick color rule: red if %REF == 0 (after parsing), otherwise black.

Outputs (per mutant):
  figures/{MUTANT}_mutations.png
  figures/{MUTANT}_mutations.svg

Usage:
  python scripts/make_mutation_maps.py \
      --excel data/all_mutants_vcf_with_genes_14082025.xlsx \
      --outdir figures \
      --sheet all

"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from matplotlib.patches import Patch
import matplotlib.pyplot as plt


# -------------------------
# Data loading & cleaning
# -------------------------

REQUIRED_COLUMNS = ["#CHROM", "MUTANT", "POS", "%REF"]

@dataclass(frozen=True)
class PlotConfig:
    line_width: float = 2.0
    tick_height: float = 0.6
    tick_width: float = 0.8
    dpi: int = 200
    figsize_per_row: float = 0.45
    left_margin_in: float = 0.9
    right_margin_in: float = 0.3
    top_margin_in: float = 0.6
    bottom_margin_in: float = 0.8

    baseline_color: str = "#bfbfbf"
    bar_width_ratio: float = 0.22
    bar_edgecolor: Optional[str] = None

    # vertical spacing between chromosomes (1.0 = old spacing)
    row_gap: float = 2.0

    # gene label/legend tuning
    label_rotation: float = 45.0
    label_fontsize: float = 8.0
    leader_lines: bool = True
    leader_color: str = "0.3"
    leader_width: float = 0.6

    # label stacking (units are “row gaps”)
    label_base_offset: float = 0.70
    label_level_step: float = 0.35

    # horizontal anti-collision
    min_label_gap_frac: float = 0.03
    min_label_gap_abs: float = 1500.0

    # increase physical row spacing if many stacked levels
    row_spacing_scale_per_level: float = 0.18

    # legend sizing (inches)
    legend_inches_per_line: float = 0.18
    legend_inches_min: float = 0.7
    legend_inches_max: float = 3.0

    # NEW: put extra blank space between main plots and legend, and below legend
    legend_pad_in: float = 0.8            # gap between main plots and legend
    legend_bottom_margin_in: float = 0.6  # extra blank space below legend

    # NEW: allow a wider figure if needed
    fig_width_in: float = 12.0


def _parse_percent_ref(val: object) -> Optional[float]:
    """Parse %REF values that may use comma decimal separators.

    Returns a float in [0, 100] or None if parsing fails.
    """
    if val is None:
        return None
    if isinstance(val, (int, float)) and not pd.isna(val):
        return float(val)
    s = str(val).strip()
    if s == "" or s.lower() in {"nan", "na"}:
        return None
    # replace comma decimal separator
    s = s.replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return None


def load_mutation_table(excel_path: Path, sheet: Optional[str] = None) -> pd.DataFrame:
    """Load and normalize the mutation table from Excel.

    Ensures required columns exist, coerces types, and standardizes column names.
    """
    if sheet is None:
        df = pd.read_excel(excel_path)
    else:
        df = pd.read_excel(excel_path, sheet_name=sheet)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {excel_path}: {missing}.\n"
            f"Available columns: {list(df.columns)}"
        )

    # Coerce types
    df = df.copy()
    df["#CHROM"] = df["#CHROM"].astype(str)
    df["MUTANT"] = df["MUTANT"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")

    # Parse %REF with potential comma decimals
    df["PCT_REF"] = df["%REF"].map(_parse_percent_ref)

    # Parse %MUT if present; else derive as 100 - PCT_REF
    if "%MUT" in df.columns:
        df["PCT_MUT"] = df["%MUT"].map(_parse_percent_ref)
    else:
        df["PCT_MUT"] = np.nan

    # If PCT_MUT is missing but PCT_REF exists, fill with 100 - PCT_REF
    mask_fill = df["PCT_MUT"].isna() & (~df["PCT_REF"].isna())
    df.loc[mask_fill, "PCT_MUT"] = 100.0 - df.loc[mask_fill, "PCT_REF"]

    # Clamp to [0, 100]
    df["PCT_MUT"] = pd.to_numeric(df["PCT_MUT"], errors="coerce").clip(lower=0, upper=100)

    # drop rows with missing essential info
    df = df.dropna(subset=["POS"]).copy()
    df["POS"] = df["POS"].astype(int)

    return df



# -------------------------
# Plotting
# -------------------------

def _compute_axes_size(n_rows: int, cfg: PlotConfig) -> tuple[float, float]:
    height = max(cfg.bottom_margin_in + cfg.top_margin_in + n_rows * cfg.figsize_per_row, 2.0)
    # Width allocates room for map + bar; keep total ~10 in
    width = 10.0
    return width, height

def _is_hypothetical(product: Optional[str]) -> bool:
    """Return True if PRODUCT is (case-insensitive) 'hypothetical protein'."""
    if product is None or (isinstance(product, float) and pd.isna(product)):
        return False
    return str(product).strip().lower() == "hypothetical protein"


def compute_gene_label_positions(df_mutant: pd.DataFrame, cfg: PlotConfig) -> dict[str, list[dict]]:
    """
    For each chromosome, compute label positions for unique genes:
      text: 'PRODUCT (n)'  where n is the number of mutations in that gene on that chromosome,
      x: median POS across that gene's mutations on that chromosome,
      level: vertical stack level to reduce overlaps along x (greedy placement).

    Includes all genes (including 'hypothetical protein'). Requires non-null PRODUCT and GENE_ID.
    """
    mask = (~df_mutant["PRODUCT"].isna()) & (~df_mutant["GENE_ID"].isna())
    dfg = df_mutant.loc[mask, ["#CHROM", "GENE_ID", "PRODUCT", "POS"]].copy()
    out: dict[str, list[dict]] = {}

    if dfg.empty:
        return out

    # For PRODUCT per (chrom, gene), prefer the most frequent (mode), else first
    agg = dfg.groupby(["#CHROM", "GENE_ID"]).agg(
        n=("POS", "count"),
        x=("POS", "median"),
        product=("PRODUCT", lambda s: s.mode().iloc[0] if not s.mode().empty else s.iloc[0]),
    ).reset_index()

    xmax_global = float(df_mutant["POS"].max()) if not df_mutant.empty else 1.0
    min_gap = max(cfg.min_label_gap_frac * xmax_global, cfg.min_label_gap_abs)

    for chrom, sub in agg.groupby("#CHROM"):
        # left-to-right by position; for ties, show higher-count labels lower (earlier)
        sub = sub.sort_values(["x", "n"], ascending=[True, False]).reset_index(drop=True)

        levels_last_x: list[float] = []  # last x placed on each level
        labels: list[dict] = []
        for _, row in sub.iterrows():
            x = float(row["x"])
            prod = str(row["product"]).strip()
            text = f"{prod} ({int(row['n'])})"

            # place on the lowest level far enough from the last placed x
            placed_level = None
            for lvl, last_x in enumerate(levels_last_x):
                if abs(x - last_x) >= min_gap:
                    placed_level = lvl
                    levels_last_x[lvl] = x
                    break
            if placed_level is None:
                placed_level = len(levels_last_x)
                levels_last_x.append(x)

            labels.append({"x": x, "text": text, "level": placed_level})

        out[chrom] = labels

    return out


def _compute_dynamic_figsize(
    n_rows: int,
    cfg: PlotConfig,
    max_label_level: int,
    n_legend_lines: int
) -> tuple[float, float, float]:
    """
    Returns (width_in, total_height_in, legend_height_in).
    Height increases with row_gap, stacked labels, and legend length.
    """
    # physical row spacing scales with row_gap and additional stacking scale
    row_scale = cfg.row_gap * (1.0 + cfg.row_spacing_scale_per_level * max(0, max_label_level + 1))

    base_main_in = cfg.bottom_margin_in + cfg.top_margin_in + n_rows * cfg.figsize_per_row * row_scale

    # extra headroom above top for stacked labels (offsets are in “row gaps” → inches)
    label_headroom_in = 0.0
    if max_label_level >= 0:
        label_headroom_in = (
            (cfg.label_base_offset + max_label_level * cfg.label_level_step + 0.6)
            * cfg.row_gap * cfg.figsize_per_row
        )

    # legend height in inches
    legend_in = 0.0
    if n_legend_lines > 0:
        legend_in = min(cfg.legend_inches_max,
                        max(cfg.legend_inches_min, cfg.legend_inches_per_line * (n_legend_lines + 1)))

    total_h = max(2.0, base_main_in + label_headroom_in + legend_in)
    width = 10.0
    return width, total_h, legend_in


def compute_gene_legend_lines(df_mutant: pd.DataFrame) -> list[str]:
    """
    Build legend lines 'GENE_ID: PRODUCT' for all genes present
    (including 'hypothetical protein'), requiring non-null PRODUCT and GENE_ID.
    """
    mask = (~df_mutant["PRODUCT"].isna()) & (~df_mutant["GENE_ID"].isna())
    dfg = df_mutant.loc[mask, ["GENE_ID", "PRODUCT"]].copy()
    if dfg.empty:
        return []
    dfg = dfg.drop_duplicates(subset=["GENE_ID"]).sort_values("GENE_ID")
    return [f"{gid}: {prod}" for gid, prod in zip(dfg["GENE_ID"], dfg["PRODUCT"])]


def _mut_bin_color(v: Optional[float]) -> str:
    """
    Map %MUT to new color bins:
      100%         -> #d73027
      [80,100)     -> #fc8d59
      [60,80)      -> #fee090
      [40,60)      -> #e0f3f8
      [20,40)      -> #91bfdb
      [0,20)       -> #4575b4
    """
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "#808080"  # fallback if missing
    v = float(v)
    v = max(0.0, min(100.0, v))  # clamp
    if v == 100.0:
        return "#d73027"
    if 80.0 <= v < 100.0:
        return "#fc8d59"
    if 60.0 <= v < 80.0:
        return "#fee090"
    if 40.0 <= v < 60.0:
        return "#e0f3f8"
    if 20.0 <= v < 40.0:
        return "#91bfdb"
    # 0–20
    return "#4575b4"


def _mut_color_legend_handles() -> tuple[list[Patch], list[str]]:
    """
    Legend handles/labels for %MUT bins in ascending order:
      [0, 20), [20, 40), [40, 60), [60, 80), [80, 100), 100
    """
    colors = ["#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"]
    labels = ["[0, 20)", "[20, 40)", "[40, 60)", "[60, 80)", "[80, 100)", "100"]
    handles = [Patch(facecolor=c, edgecolor="none") for c in colors]
    return handles, labels



def plot_mutations_for_mutant(df_mutant: pd.DataFrame, mutant: str, outdir: Path, cfg: PlotConfig) -> None:
    """Create and save a mutation map figure for one MUTANT."""
    chrom_order = (
        df_mutant.groupby("#CHROM")["POS"].max().sort_values(ascending=False).index.tolist()
    )
    n = len(chrom_order)

    legend_lines = compute_gene_legend_lines(df_mutant)

    # Labels & stacking levels
    gene_labels_by_chrom = compute_gene_label_positions(df_mutant, cfg)
    max_level = -1
    for chrom in gene_labels_by_chrom:
        for item in gene_labels_by_chrom[chrom]:
            max_level = max(max_level, int(item["level"]))

    # --- Figure sizing pieces (inches)
    row_scale = cfg.row_gap * (1.0 + cfg.row_spacing_scale_per_level * max(0, max_level + 1))
    main_in = cfg.bottom_margin_in + cfg.top_margin_in + n * cfg.figsize_per_row * row_scale
    label_headroom_in = 0.0
    if max_level >= 0:
        label_headroom_in = (
            (cfg.label_base_offset + max_level * cfg.label_level_step + 0.6)
            * cfg.row_gap * cfg.figsize_per_row
        )
    main_in += label_headroom_in

    pad_in = max(0.4, cfg.legend_pad_in)
    legend_text_in = 0.0
    if legend_lines:
        legend_text_in = min(cfg.legend_inches_max,
                             max(cfg.legend_inches_min, cfg.legend_inches_per_line * (len(legend_lines) + 1)))
    legend_in = legend_text_in + cfg.legend_bottom_margin_in

    total_h = max(2.0, main_in + pad_in + legend_in)

    fig = plt.figure(figsize=(cfg.fig_width_in, total_h), dpi=cfg.dpi)

    # Grid: 3 rows -> [main | spacer | legend]
    bar_ratio = max(min(cfg.bar_width_ratio, 0.45), 0.05)
    gs = fig.add_gridspec(
        3, 2,
        width_ratios=[1.0 - bar_ratio, bar_ratio],
        height_ratios=[main_in, pad_in, legend_in],
        wspace=0.25, hspace=0.00
    )
    ax     = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_pad = fig.add_subplot(gs[1, :]); ax_pad.axis("off")
    ax_leg = fig.add_subplot(gs[2, :]); ax_leg.axis("off")

    # Y mapping with expanded row gap
    y_map = {chrom: idx * cfg.row_gap for idx, chrom in enumerate(chrom_order[::-1], start=1)}

    xmax = int(df_mutant["POS"].max()) if not df_mutant.empty else 1
    xmin = 0

    # Baselines
    for chrom in chrom_order:
        y = y_map[chrom]
        chrom_max = int(df_mutant.loc[df_mutant["#CHROM"] == chrom, "POS"].max())
        ax.hlines(y=y, xmin=xmin, xmax=chrom_max, linewidth=cfg.line_width, color=cfg.baseline_color)

    # Ticks colored by %MUT bins
    tick_half = (cfg.tick_height * cfg.row_gap) / 2.0
    for _, row in df_mutant.iterrows():
        y = y_map[row["#CHROM"]]
        x = int(row["POS"])
        color = _mut_bin_color(row.get("PCT_MUT", np.nan))
        ax.vlines(x=x, ymin=y - tick_half, ymax=y + tick_half, linewidth=cfg.tick_width, color=color)

    # ----- Right stacked bar — %MUT bins (ASCENDING order to match legend) -----
    # Bins: [0,20), [20,40), [40,60), [60,80), [80,100), {100}
    bin_edges        = [-1, 20, 40, 60, 80, 100, 101]    # right=False -> [a,b)
    bin_labels_asc   = [0, 1, 2, 3, 4, 5]                # ascending by edges
    order_for_plot   = [0, 1, 2, 3, 4, 5]                # low -> high
    colors_for_order = ["#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"]

    # use rows with non-NaN %MUT for bar counts
    df_bins = df_mutant.dropna(subset=["PCT_MUT"]).copy()
    df_bins["mut_bin"] = pd.cut(df_bins["PCT_MUT"], bins=bin_edges, labels=bin_labels_asc, right=False)

    # counts per chromosome per bin
    counts = (
        df_bins.groupby(["#CHROM", "mut_bin"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=chrom_order, columns=bin_labels_asc, fill_value=0)
    )

    y_vals = [y_map[c] for c in chrom_order]
    left = np.zeros(len(chrom_order), dtype=float)

    for i, color in zip(order_for_plot, colors_for_order):
        widths = counts[i].to_numpy() if i in counts.columns else np.zeros(len(chrom_order))
        ax_bar.barh(
            y=y_vals, width=widths, left=left,
            color=color, edgecolor=cfg.bar_edgecolor
        )
        left += widths

    # ----- Gene labels above the map (45°) + leader lines -----
    base_offset = cfg.label_base_offset * cfg.row_gap
    level_step  = cfg.label_level_step * cfg.row_gap
    top_extra   = 0.0 if max_level < 0 else (cfg.label_base_offset + max_level * cfg.label_level_step + 0.6) * cfg.row_gap

    # Axes formatting (left)
    y_top = n * cfg.row_gap + top_extra
    ax.set_xlim(xmin, xmax * 1.02)
    ax.set_ylim(0, y_top)
    ax.set_yticks([y_map[c] for c in chrom_order])
    ax.set_yticklabels(chrom_order)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("#CHROM")
    ax.set_title(f"Mutations in {mutant}")
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    for chrom in chrom_order:
        if chrom not in gene_labels_by_chrom:
            continue
        y_base = y_map[chrom]
        for item in gene_labels_by_chrom[chrom]:
            x = float(item["x"])
            lvl = int(item["level"])
            y_text = y_base + base_offset + lvl * level_step
            if cfg.leader_lines:
                ax.plot([x, x], [y_base + tick_half, y_text - 0.06 * cfg.row_gap],
                        color=cfg.leader_color, linewidth=cfg.leader_width, alpha=0.8)
            ax.text(
                x, y_text, item["text"],
                ha="center", va="bottom",
                rotation=cfg.label_rotation, rotation_mode="anchor",
                fontsize=cfg.label_fontsize
            )

    # Right axis formatting (no legend here; color legend is in the bottom row)
    ax_bar.set_ylim(0, y_top)
    ax_bar.set_yticks(y_vals)
    ax_bar.set_yticklabels([])
    ax_bar.set_xlabel("Mutations")
    xmax_bar = int(max(1, counts.sum(axis=1).max()))
    ax_bar.set_xlim(0, xmax_bar * 1.2)
    ax_bar.grid(axis="x", linestyle=":", alpha=0.3)

    # Bottom row: %MUT color legend + gene legend text
    handles, labels = _mut_color_legend_handles()
    ax_leg.legend(
        handles, labels,
        title="Mutant allele % (%MUT)",
        loc="upper left", frameon=False, ncol=6,
        handlelength=1.2, handletextpad=0.5, columnspacing=1.0, borderaxespad=0.0
    )

    if legend_lines:
        legend_text = "Genes shown (GENE_ID: PRODUCT):\n" + "\n".join(legend_lines)
        ax_leg.text(0.01, 0.60, legend_text, ha="left", va="top",
                    fontsize=cfg.label_fontsize, transform=ax_leg.transAxes)

    outdir.mkdir(parents=True, exist_ok=True)
    png_path = outdir / f"{mutant}_mutations.png"
    svg_path = outdir / f"{mutant}_mutations.svg"
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


# -------------------------
# Orchestration
# -------------------------

def make_all_plots(df: pd.DataFrame, outdir: Path, cfg: Optional[PlotConfig] = None) -> None:
    cfg = cfg or PlotConfig()
    for mutant, df_mut in df.groupby("MUTANT", sort=True):
        plot_mutations_for_mutant(df_mut, mutant, outdir, cfg)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Make per-MUTANT mutation maps from Excel.")
    p.add_argument("--excel", type=Path, required=True, help="Path to the Excel file with mutations.")
    p.add_argument("--sheet", type=str, default=None, help="Worksheet name; default: first sheet")
    p.add_argument("--outdir", type=Path, default=Path("figures"), help="Output directory for figures")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df = load_mutation_table(args.excel, sheet=args.sheet)
    make_all_plots(df, args.outdir)


if __name__ == "__main__":
    main()
