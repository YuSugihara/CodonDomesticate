# NLRexpress Colab Workflow

[![Open NLRexpress In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/NLRexpress_Colab.ipynb)

This Colab notebook runs [NLRexpress](https://github.com/eliza-m/NLRexpress) and prepares optional handoff files for CodonDomesticate.

## What It Does

- accepts either protein or CDS input
- translates CDS input with Biopython before running NLRexpress
- runs NLRexpress on one sequence or many sequences
- downloads the original NLRexpress result folder as a zip archive
- extracts MHD motif hits without modifying the original NLRexpress files
- writes D-to-V mutation candidates such as `D489V` to a separate CSV file
- if CDS input is provided, writes a CodonDomesticate multiple-domestication CSV containing `name`, `sequence`, `aa_change`, and MHD candidate metadata
- can regenerate the CodonDomesticate batch table from data already present in the Colab session, asking for a matching CDS CSV only when no CDS table is available
- lets you set an MHD probability warning threshold before writing CodonDomesticate input, so low-confidence MHD calls are visible before mutation/domestication

## Recommended Workflows

### CDS Available

1. Open the NLRexpress Colab notebook.
2. Set the input type to `CDS`.
3. Run NLRexpress.
4. Download the original NLRexpress result archive.
5. Download the CodonDomesticate multiple-domestication CSV.
6. Open the CodonDomesticate notebook and upload that CSV for batch domestication/mutation.

### Protein Only

1. Open the NLRexpress Colab notebook.
2. Set the input type to `protein`.
3. Run NLRexpress.
4. Download the original NLRexpress result archive and MHD D-to-V candidate CSV.
5. Provide matching CDS sequences later if you want to perform CDS domestication.

[![Open CodonDomesticate In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/Domesticate_CDS_Colab.ipynb)

## Input Formats

### Single Sequence

Paste either one CDS sequence or one amino-acid sequence into the notebook and set the sequence type accordingly.

### Multiple Sequences

Upload either a FASTA file or a CSV file with these columns:

```csv
name,sequence
example_NLR_1,MA...
example_NLR_2,MA...
```

For CDS input, the `sequence` column should contain CDS sequences:

```csv
name,sequence
example_NLR_1,ATG...
example_NLR_2,ATG...
```

## MHD D-to-V Logic

NLRexpress reports MHD hits in the short output with a motif start position and motif sequence. The conserved D is treated as the last `D` in the motif sequence. The notebook converts that residue to valine and writes the result to a separate mutation-candidate CSV, producing an `aa_change` value like:

```text
D489V
```

This value is compatible with `domesticate_cds.py` and the CodonDomesticate batch notebook.

## Output Files

The notebook can download:

- the original NLRexpress result folder as a zip archive
- `single_mhd_dv_candidates.csv` or `multi_mhd_dv_candidates.csv`
- `single_codon_domesticate_input_with_mhd_dv.csv` or `multi_codon_domesticate_input_with_mhd_dv.csv` when CDS input is provided; these contain the `name`, `sequence`, and `aa_change` columns expected by CodonDomesticate batch mode, plus MHD candidate metadata such as motif position, probability, `mhd_probability_warning_threshold`, and `mhd_probability_below_threshold`

The original NLRexpress output files are not edited.

## Notes

- NLRexpress itself accepts protein FASTA input; this notebook translates CDS input with Biopython before running NLRexpress.
- CodonDomesticate accepts CDS input and optional `aa_change` values.
- If no MHD motif is detected, the notebook prints a warning and still writes a handoff CSV with an empty `aa_change` column.
- If some CDS records have no MHD candidate, or if a sequence has multiple MHD candidates, the notebook prints warnings and writes a name-check report CSV. For multiple candidates, the highest-probability MHD candidate is used.
- If an MHD candidate probability is below the configured warning threshold, the notebook prints a warning and records the low-confidence call in the name-check report. The candidate is not automatically removed; review the CSV before using it for CodonDomesticate.
- If only protein sequences are available, CDS domestication cannot be performed until matching CDS sequences are provided.
- The NLRexpress repository is GPL-3.0 licensed. This notebook downloads and runs NLRexpress but does not vendor its source code or modify its original result files.
- For publication or redistribution, cite NLRexpress as requested in the upstream repository.
