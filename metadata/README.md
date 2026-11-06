# Metadata

## raw/

| File | Produced by | Content |
|---|---|---|
| `ena_runs_raw.tsv` | `01_fetch_ena_runs.py` | All public transcriptomic runs of Beta vulgaris in ENA with FASTQ links and MD5 sums |
| `ena_runs_raw.retrieved.txt` | `01_fetch_ena_runs.py` | Retrieval date and query, written when the table is rebuilt |
| `samples_candidate.tsv` | `02_filter_runs.py` | Runs passing the splice-analysis filter |
| `samples_candidate_attrs.tsv` | `03_fetch_biosample_attributes.py` | BioSample attributes per run |
| `run_aliases.tsv` | `04_fetch_run_aliases.py` | Submitter run, experiment, library and sample aliases |

The records in `raw/` were retrieved from ENA on 17 September 2026.

## Sample sheets

`samples_all.tsv` and `samples_pilot.tsv` are produced by `05_make_sample_sheets.py`.

| Column | Meaning |
|---|---|
| run | ENA run accession |
| study | BioProject accession |
| sample | BioSample accession |
| category | Broad experiment type (abiotic_cold, biotic_oomycete, development, ...) |
| condition | Treatment or control label used for comparisons |
| tissue | Controlled term derived from `tissue_reported` (leaf, root, taproot, seed, leaf+root, inflorescence, not reported) |
| genotype | Cultivar, line or genotype code as given by the submitter; codes such as GT1 are only meaningful within one BioProject |
| timepoint | Time point or developmental stage |
| rep | Biological replicate |
| use_das | `yes`: replicated design usable for differential splicing; `junction_only`: used for isoform and junction discovery only; `pending_curation`: labels still missing; `excluded`: not analyzed (reason in `curation_note`) |
| read_len | Mean read length per mate (bp) |
| gb | Sequenced bases per run (Gb) |
| instrument | Sequencing instrument |
| fq1, fq2 | HTTPS links to the paired FASTQ files |
| md5_1, md5_2 | MD5 checksums from ENA |
| sample_alias, library_name | Submitter aliases used for labeling |
| tissue_reported | Tissue exactly as written in the BioSample record (empty if not given) |
| curation_note | Reason for any manual override from `curation/run_overrides.tsv` (empty for untouched runs) |

## Manual curation

`curation/run_overrides.tsv` lists every manual change to the automatic labels (run, column, new value, reason). `05_make_sample_sheets.py` applies it and copies the reason into `curation_note`. Only runs with `use_das = yes` take part in the comparisons listed in `comparisons_pilot.tsv`. `comparisons_tissue.tsv` holds the leaf versus taproot comparisons of PRJNA602804 (same columns), run separately with `COMPARISONS_FILE` and summarized into their own directory with `SUMMARY_DIR`; they are not part of the 40 condition comparisons.

`curation/splicing_factor_families.tsv` lists the Arabidopsis thaliana UniProt entries (reference proteome UP000006548, retrieved 26 September 2026) that define the splicing-regulator families used by `scripts/analysis/splicing_regulators.py`: family, subfamily, UniProt entry name, gene name, the entry description copied from UniProt, and the source. SR subfamily labels follow the gene name prefixes (SR, RSZ, SC, SCL, RS, RS2Z); SR45 and SR45a are kept in a separate SR-like group. A sugar beet protein belongs to a family when its best DIAMOND hit in the Arabidopsis proteome is one of the listed entries; SR proteins are in addition tested by sequence with the plant definition of Barta, Kalyna and Reddy (2010, Plant Cell 22:2926): one or two RRMs followed by an RS domain of at least 50 residues with at least 20 percent of the residues in RS or SR dipeptides (five dipeptides in a 50-residue window; all 18 Arabidopsis SR proteins satisfy this reading).

## External gene lists

`external/calixto2018_cold_das_genes.txt` (2,442 AGI codes) and `external/calixto2018_expressed_genes.txt` (19,902 AGI codes) are the columns All_DAS_genes and Expressed_genes of sheet 1B of Supplemental Data Set 1 of Calixto et al. (2018), The Plant Cell 30:1424-1444, doi 10.1105/tpc.18.00177, taken unchanged (upper case, one code per line, sorted). The counts agree with sheet 1A of the same file and with the article. They are used by `scripts/analysis/cross_species_cold.py`.
