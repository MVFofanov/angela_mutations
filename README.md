# Angela Mutations — mutation map plots

This repository generates per-sample (per **MUTANT**) genome mutation maps from an Excel table. For each sample it draws:

- a **genomic map** with one horizontal line per `#CHROM` and vertical ticks for each mutation, **colored by `%MUT` bin** (see scale below);
- a **stacked bar chart** at right showing, per `#CHROM`, **counts of mutations binned by `%MUT`** (same colors, **left→right from low to high**);
- **gene labels** above the genomic map as **`PRODUCT (n)`** where `n` is the number of mutations in that gene on that chromosome;
- a bottom **legend** with the `%MUT` color scale **and** a list of **`GENE_ID: PRODUCT (n)`** for all genes (including *hypothetical protein*).

All figures are saved to `figures/` as both **PNG** and **SVG**.

---

## Quick start

From the repo root:

```bash
python3 scripts/make_mutation_maps.py \
  --excel data/all_mutants_vcf_with_genes_14082025.xlsx \
  --sheet all \
  --outdir figures
```

Output files will be created as figures/{MUTANT}_mutations.png and figures/{MUTANT}_mutations.svg.

## Input

Place your Excel file under data/. The script expects at least these columns (case-sensitive):

| Column    | Type   | Description                                                                                         |
| --------- | ------ | --------------------------------------------------------------------------------------------------- |
| `#CHROM`  | string | Chromosome/contig accession                                                                         |
| `MUTANT`  | string | Sample / mutant ID (one figure per unique value)                                                    |
| `POS`     | int    | 1-based genomic position of the mutation                                                            |
| `%REF`    | float  | Percent reads supporting the reference allele (comma decimals like `47,14` are accepted)            |
| `GENE_ID` | string | Gene identifier                                                                                     |
| `PRODUCT` | string | Product name                                                                                        |
| `%MUT`    | float  | *(optional)* Percent reads supporting the mutant allele. If missing, it’s computed as `100 − %REF`. |

Only the columns above are used; others are ignored.
If your workbook has multiple sheets, use --sheet (e.g. --sheet all).

## Output & plot semantics

For each MUTANT value, two files are written to the chosen --outdir (default figures/):

- {MUTANT}_mutations.png

- {MUTANT}_mutations.svg

Baselines: one gray horizontal line per #CHROM, length scaled to the maximum observed POS for that chromosome in the sample.
Ticks: vertical marks at each POS, colored by %MUT bin:

| %MUT label  | Numeric interval  | Color     |
| ----------- | ----------------- | --------- |
| `[0, 20)`   | `0 ≤ %MUT < 20`   | `#4575b4` |
| `[20, 40)`  | `20 ≤ %MUT < 40`  | `#91bfdb` |
| `[40, 60)`  | `40 ≤ %MUT < 60`  | `#e0f3f8` |
| `[60, 80)`  | `60 ≤ %MUT < 80`  | `#fee090` |
| `[80, 100)` | `80 ≤ %MUT < 100` | `#fc8d59` |
| `100`       | `%MUT = 100`      | `#d73027` |


Right bar: stacked counts per #CHROM, one stack per %MUT bin, ordered from low to high (same colors).
Gene labels: above each chromosome line at the gene’s median mutated position: PRODUCT (n) with simple stacking and optional leader lines.
Legend: bottom row shows the %MUT color legend (low→high) and a list of GENE_ID: PRODUCT (n) (includes hypothetical protein).

## CLI options
```bash
--excel  PATH   (required) Excel file path
--sheet  NAME   (optional) Worksheet name; default is the first sheet
--outdir PATH   (optional) Output directory; default: figures
```

Example:

```bash
python3 scripts/make_mutation_maps.py --excel data/all_mutants_vcf_with_genes_14082025.xlsx --sheet all --outdir figures
```

## Configuration knobs

Adjust these in PlotConfig inside scripts/make_mutation_maps.py:

- row_gap: vertical spacing between chromosomes (default 2.0, i.e. double spacing)

- label_rotation: gene label angle (default 45°)

- leader_lines: draw connectors from ticks to labels (default True)

- min_label_gap_frac / min_label_gap_abs: horizontal anti-collision between labels

- fig_width_in: figure width in inches (default 12.0)

- legend_pad_in, legend_bottom_margin_in: extra spacing to push the legend further down

- label_fontsize, label_base_offset, label_level_step: label readability/stacking

Re-run the command after tweaking.

## Dependencies

Python 3.8+

pandas, numpy, matplotlib, and an Excel engine (e.g. openpyxl for .xlsx)

Install (one way):

```bash
python3 -m pip install pandas numpy matplotlib openpyxl
```

Repo layout

```bash
.
├── data/                 # input Excel (ignored by git)
├── figures/              # output figures (ignored by git)
└── scripts/
    └── make_mutation_maps.py
```
