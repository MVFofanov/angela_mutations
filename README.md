# Angela Mutations — mutation map plots

This repository generates per‑sample (per **MUTANT**) genome mutation maps from an Excel table. For each sample it draws:

* a **genomic map** with one horizontal line per `#CHROM` and vertical ticks for each mutation (red if `%REF=0`, otherwise black), and
* a **stacked bar chart** summarizing counts per `#CHROM` (gray = present in reference+mutant, red = mutant‑only),
* **gene labels** above the genomic map (`GENE_ID (n)`), excluding rows where `PRODUCT` is *hypothetical protein*, plus a bottom legend `GENE_ID: PRODUCT`.

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

> Output files will be created as `figures/{MUTANT}_mutations.png` and `figures/{MUTANT}_mutations.svg`.

---

## Input

Place your Excel file at `data/`. The script expects at least these columns (case‑sensitive):

| Column    | Type   | Description                                                                              |
| --------- | ------ | ---------------------------------------------------------------------------------------- |
| `#CHROM`  | string | Chromosome/contig accession                                                              |
| `MUTANT`  | string | Sample / mutant ID (one figure per unique value)                                         |
| `POS`     | int    | 1‑based genomic position of the mutation                                                 |
| `%REF`    | float  | Percent reads supporting the reference allele (comma decimals like `47,14` are accepted) |
| `GENE_ID` | string | Gene identifier (used for labels/legend)                                                 |
| `PRODUCT` | string | Product name; labels exclude rows where value is *hypothetical protein*                  |

> Only the columns above are used; others are ignored.

If your workbook has multiple sheets, use `--sheet` to select one (e.g. `--sheet all`).

---

## Output

For each `MUTANT` value in the table, two files are written to the chosen `--outdir` (default `figures/`):

* `{MUTANT}_mutations.png`
* `{MUTANT}_mutations.svg`

### Plot semantics

* **Baselines**: one gray horizontal line per `#CHROM`, length scaled to the maximum observed `POS` for that chromosome in the sample.
* **Ticks**: vertical marks at each `POS`.

  * **red** = `%REF == 0` (mutant‑only)
  * **black** = otherwise (both ref+mut reads present)
* **Right bar**: stacked counts per `#CHROM` (gray = ref+mut, red = mutant‑only).
* **Gene labels**: `GENE_ID (n)` above the genome line at the gene’s median mutated position; simple stacking, with optional leader lines.
* **Legend**: bottom block listing `GENE_ID: PRODUCT` for all shown genes.

---

## CLI options

```text
--excel  PATH   (required) Excel file path
--sheet  NAME   (optional) Worksheet name; default is the first sheet
--outdir PATH   (optional) Output directory; default: figures
```

Example:

```bash
python3 scripts/make_mutation_maps.py --excel data/all_mutants_vcf_with_genes_14082025.xlsx --sheet all --outdir figures
```

---

## Configuration knobs

Most plot behavior can be tuned by editing `PlotConfig` in `scripts/make_mutation_maps.py`:

* `row_gap`: vertical spacing between chromosomes (default **2.0**, i.e. double spacing)
* `label_rotation`: rotates gene labels (default **45°**)
* `leader_lines`: draw connector lines from ticks to labels (default **True**)
* `min_label_gap_frac` / `min_label_gap_abs`: horizontal anti‑collision threshold
* `fig_width_in`: figure width in inches (default **12.0**)
* `legend_pad_in`, `legend_bottom_margin_in`: push the legend further down

After tweaking, re‑run the same command.

---

## Repo layout

```
.
├── data/                 # input Excel (ignored by git)
├── figures/              # output figures (ignored by git)
└── scripts/
    └── make_mutation_maps.py
```

---
