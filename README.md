# COVID-19 LAMP Primer Analysis

R scripts for analysing SARS-CoV-2 sequence diversity across multiple LAMP (Loop-mediated Isothermal Amplification) assays. The pipeline takes multiple sequence alignments (MSA) in FASTA format, counts duplicate sequences, maps strain codes, and extracts primer-binding regions for quality assessment.

---

## Scripts

### `covid_portfolio.R`
The main analysis pipeline. It:
- Reads three FASTA alignment files (one per assay)
- Builds a strain code ↔ sequence lookup table (Table 1)
- Counts the frequency of each unique sequence (Table 2)
- Maps each unique sequence back to all strain codes that share it (Table 3)
- Defines primer sets for all three assays (Harvard As1e, Color N, NEB Gene E1)
- Extracts the subsequence at each primer-binding site across all sequences

### `duplicate_counter_portfolio.R`
A standalone, reusable function (`duplicate_counter`) that wraps the sequence-frequency logic from the main script. Given a list of sequences and an output path, it writes a tab-delimited file containing each unique sequence and its count.

---

## Dependencies

Install the following R packages before running:

```r
install.packages(c("seqinr", "dplyr", "tidyr"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ampir")
```

| Package | Purpose |
|---|---|
| `seqinr` | Reading FASTA files |
| `dplyr` / `tidyr` | Data wrangling |
| `Biostrings` | DNA string manipulation and subsequence extraction |
| `ampir` | Antimicrobial peptide / sequence utilities |

---

## Input Files

| File | Description |
|---|---|
| `1_msa.fasta` | Multiple sequence alignment — Assay 1 (Harvard As1e) |
| `2_msa.fa` | Multiple sequence alignment — Assay 2 (Color N) |
| `3_msa.fa` | Multiple sequence alignment — Assay 3 (NEB Gene E1) |

---

## Outputs

| File | Contents |
|---|---|
| `table.strain_code_and_sequence_assay2.txt` | Strain code → sequence mapping (tab-delimited) |
| `table.NB_of_duplicates_assay2.txt` | Unique sequences and their frequency counts |
| `table.sequences_assay3.txt` | Raw sequence list for Assay 3 |
| `table.sequences_and_strain_code_assay2.txt` | Each unique sequence with all matching strain codes |
| `subset_sequences_HARVARD_As1e_FIPpart_2.txt` | Subsequences at the FIP part 2 primer binding site |

---

## Assays & Primers

Three LAMP assays are analysed, each targeting a different region of the SARS-CoV-2 genome:

| Assay | Target | Primers defined |
|---|---|---|
| Harvard As1e | — | F3, B3, FIP (×2), BIP (×2), LF, LB |
| Color N | N gene | F3, B3, FIP (×2), BIP (×2), LF, LB |
| NEB Gene E1 | E gene | F3, B3, FIP (×2), BIP (×2), LF, LB |
