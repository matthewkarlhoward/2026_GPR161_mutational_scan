#!/bin/bash
# Rebuild every GPR161 panel in this repository from data/counts/ and data/scores/.
#
# Covers manuscript Figures 2, 3, 4 and 5 and Supplemental Figure 3. Figure 1 and
# Supplemental Figures 1-2 are the GPCRome screen and live in the companion
# repository, 2026_GPCRome_basal_cAMP.
#
# Figure 4 panels a and c are ChimeraX renders made by hand from
# data/structure/8smv.cif, coloured by the per-position values in
# data/processed/position_summary.tsv. Only panel b is drawn here.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

# --- scores -> the merged tables the figures read ---
for s in 01_merge_scores 02_position_summary 03_structure_annotation; do
  echo "=== pipeline/$s"
  "$PY" "pipeline/$s.py"
done

# --- the published panels ---
for s in fig1_landscape fig2_phenotype_classes fig3_structure fig4_human_variants; do
  echo "=== scripts/$s"
  "$PY" "scripts/$s.py"
done

# --- library and sequencing QC (Supplemental Figure 3) ---
for s in dms_qc_coverage dms_qc_scores dms_qc_panel; do
  echo "=== scripts/qc/$s"
  "$PY" "scripts/qc/$s.py"
done

echo
echo "Done."
echo "  Figure 2 panels b-c : figures/landscape_density_heat.pdf"
echo "  Figure 3 panels b-e : figures/phenotype_classes.pdf"
echo "  Figure 4 panel  b   : figures/segment_effects.pdf"
echo "                        (panels a, c: ChimeraX, by hand)"
echo "  Figure 5 panels a-b : figures/human_variants.pdf"
echo "  Supp. Figure 3      : figures/F_QC_panel.pdf"
