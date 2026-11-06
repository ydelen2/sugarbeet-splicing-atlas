#!/bin/bash
# Summary table of rMATS and SUPPA2 results for every comparison. Run on a login node after 07 and 08:
#   bash workflow/09_summarize.sh
# rMATS events count only when every replicate has at least MIN_EVENT_READS supporting reads.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env
MIN_EVENT_READS="${MIN_EVENT_READS:-20}"
# Output directory; give SUMMARY_DIR together with COMPARISONS_FILE to summarize a separate comparison set
# (for example the tissue comparisons) without touching the main tables
SUMMARY_DIR="${SUMMARY_DIR:-${RES_DIR}/summary}"; mkdir -p "$SUMMARY_DIR"
python "${SCRIPTS_DIR}/analysis/summarize_splicing.py" --comparisons "$COMPARISONS" \
  --rmats-dir "$RMATS_DIR" --suppa-dir "$SUPPA_DIR" --min-reads "$MIN_EVENT_READS" \
  --out "${SUMMARY_DIR}/splicing_summary.tsv"
echo "table: ${SUMMARY_DIR}/splicing_summary.tsv"
# Event-level agreement of the two tools (SE, A3, A5, RI, MX matched on splice-site coordinates), with a
# Benjamini-Hochberg correction of the SUPPA2 p-values over the testable events of each type and comparison
python "${SCRIPTS_DIR}/analysis/event_concordance.py" --comparisons "$COMPARISONS" \
  --suppa-dir "$SUPPA_DIR" --rmats-dir "$RMATS_DIR" --map "${FILT_DIR}/transcript_map.tsv" \
  --min-reads "$MIN_EVENT_READS" --out-prefix "${SUMMARY_DIR}/event_concordance"
# Reference gene, eggNOG description and cross-comparison recurrence for the high-confidence events
python "${SCRIPTS_DIR}/analysis/annotate_events.py" --events "${SUMMARY_DIR}/event_concordance_high_confidence.tsv" \
  --gff "${REF_DIR}/BeetSet-3.genes.2604.gff" --functional "${REF_DIR}/BeetSet-3.genes.2604.functional" \
  --out "${SUMMARY_DIR}/events_high_confidence_annotated.tsv"
# Coding consequences of the high-confidence events (CDS or UTR, reading frame, stop codons in retained introns)
python "${SCRIPTS_DIR}/analysis/event_consequences.py" --events "${SUMMARY_DIR}/events_high_confidence_annotated.tsv" \
  --gff "${REF_DIR}/BeetSet-3.genes.2604.gff" --genome "${REF_DIR}/RefBeet-3.0.fa" \
  --out "${SUMMARY_DIR}/events_high_confidence_consequences.tsv"
# Premature stop codons and NMD features of both forms of every event, built on the annotated gene models
python "${SCRIPTS_DIR}/analysis/nmd_features.py" --events "${SUMMARY_DIR}/events_high_confidence_annotated.tsv" \
  --consequences "${SUMMARY_DIR}/events_high_confidence_consequences.tsv" \
  --gff "${REF_DIR}/BeetSet-3.genes.2604.gff" --genome "${REF_DIR}/RefBeet-3.0.fa" \
  --out "${SUMMARY_DIR}/events_high_confidence_nmd.tsv"
# Splicing-regulator families and Pfam domains touched by the coding events (need the outputs of 10_domains.sbatch);
# the family tables of a separate comparison set go to SUMMARY_DIR, those of the main set stay in the domains directory
DOM="${RES_DIR}/domains"
FAM_DIR="$DOM"; [ "$SUMMARY_DIR" = "${RES_DIR}/summary" ] || FAM_DIR="$SUMMARY_DIR"
if [ -s "$DOM/pfam.domtblout" ]; then
  python "${SCRIPTS_DIR}/analysis/splicing_regulators.py" --domtbl "$DOM/pfam.domtblout" --proteins "$DOM/beet_proteins.fa" \
    --beet-ath "$DOM/beet_vs_ath.tsv" --ath-beet "$DOM/ath_vs_beet.tsv" \
    --families "${META_DIR}/curation/splicing_factor_families.tsv" \
    --events "${SUMMARY_DIR}/events_high_confidence_annotated.tsv" --suppa-dir "$SUPPA_DIR" --comparisons "$COMPARISONS" \
    --map "${FILT_DIR}/transcript_map.tsv" --out-prefix "${FAM_DIR}/splicing_regulators"
  python "${SCRIPTS_DIR}/analysis/domain_impact.py" --consequences "${SUMMARY_DIR}/events_high_confidence_consequences.tsv" \
    --gff "${REF_DIR}/BeetSet-3.genes.2604.gff" --domtbl "$DOM/pfam.domtblout" \
    --out "${SUMMARY_DIR}/events_high_confidence_domains.tsv"
else
  echo "splicing-regulator families and domain impact skipped: $DOM/pfam.domtblout not found (run 10_domains.sbatch)"
fi
# Abscisic acid pathway genes among the genes with events (terms in metadata/curation/aba_pathway_terms.tsv)
python "${SCRIPTS_DIR}/analysis/pathway_genes.py" --functional "${REF_DIR}/BeetSet-3.genes.2604.functional" \
  --terms "${META_DIR}/curation/aba_pathway_terms.tsv" --events "${SUMMARY_DIR}/events_high_confidence_annotated.tsv" \
  --suppa-dir "$SUPPA_DIR" --comparisons "$COMPARISONS" --map "${FILT_DIR}/transcript_map.tsv" --out-prefix "${SUMMARY_DIR}/aba_pathway"
# Abundance of the premature-stop form among the non-annotated forms of the coding events (SUPPA2 PSI over all runs)
PSI_ALL="${SUPPA_DIR}/psi/all_events.psi"
if [ -s "$PSI_ALL" ]; then
  python "${SCRIPTS_DIR}/analysis/ptc_form_abundance.py" --nmd "${SUMMARY_DIR}/events_high_confidence_nmd.tsv"     --consequences "${SUMMARY_DIR}/events_high_confidence_consequences.tsv" --psi "$PSI_ALL"     --out-prefix "${SUMMARY_DIR}/ptc_form_abundance"
else
  echo "PTC form abundance skipped: $PSI_ALL not found"
fi
# Expression of the NMD factors (UPF1, UPF2, UPF3, SMG7) in the main comparison set (needs 12_expression.sbatch)
IDMAP="${REF_DIR}/proteomes/UP000006548_3702.idmapping.gz"
if [ "$SUMMARY_DIR" = "${RES_DIR}/summary" ] && [ -d "${RES_DIR}/expression/de_results" ] && [ -s "$IDMAP" ]; then
  python "${SCRIPTS_DIR}/analysis/nmd_factor_expression.py" --beet-ath "$DOM/beet_vs_ath.tsv" --idmapping "$IDMAP"     --de-dir "${RES_DIR}/expression/de_results" --comparisons "$COMPARISONS" --out-prefix "${SUMMARY_DIR}/nmd_factors"
else
  echo "NMD factor expression skipped (main set only; needs ${RES_DIR}/expression/de_results and $IDMAP)"
fi
# GO enrichment of the genes with high-confidence events (term names from ref/go/go-basic.obo when present)
OBO="${REF_DIR}/go/go-basic.obo"
[ -s "$OBO" ] || echo "go-basic.obo not found in ${REF_DIR}/go, GO terms are reported without names"
python "${SCRIPTS_DIR}/analysis/go_enrichment.py" --functional "${REF_DIR}/BeetSet-3.genes.2604.functional" --obo "$OBO" \
  --events "${SUMMARY_DIR}/events_high_confidence_annotated.tsv" --suppa-dir "$SUPPA_DIR" --comparisons "$COMPARISONS" \
  --map "${FILT_DIR}/transcript_map.tsv" --out "${SUMMARY_DIR}/go_enrichment.tsv"
# Digest of all tables, the source of the numbers quoted in the manuscript
# The expression tables of 12_expression.sbatch belong to the main comparison set only
EXPR_OPT=""; [ "$SUMMARY_DIR" = "${RES_DIR}/summary" ] && EXPR_OPT="--expression ${RES_DIR}/expression"
python "${SCRIPTS_DIR}/analysis/results_digest.py" --summary "${SUMMARY_DIR}" --domains "$FAM_DIR" \
  --comparisons "$COMPARISONS" $EXPR_OPT > "${SUMMARY_DIR}/results_digest.txt"
echo "digest: ${SUMMARY_DIR}/results_digest.txt"
