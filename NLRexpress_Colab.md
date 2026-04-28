# NLRexpress Colab Workflow

[![Open NLRexpress In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/NLRexpress_Colab.ipynb)

This Colab notebook runs [NLRexpress](https://github.com/eliza-m/NLRexpress) and prepares optional handoff files for CodonDomesticate.

## What It Does

- accepts either protein or CDS input
- translates CDS input with Biopython before running NLRexpress
- runs NLRexpress on one sequence or many sequences
- downloads the original NLRexpress result folder as a zip archive
- extracts MHD motif hits without modifying the original NLRexpress files
- writes D-to-V mutation candidates such as `D489V` to a temporary CSV in the Colab session
- can write a CodonDomesticate multiple-domestication CSV from data already present in the Colab session, asking for a matching CDS CSV only when no CDS table is available
- only writes the CodonDomesticate batch table when you run the optional handoff cell
- lets you set an MHD probability warning threshold in the optional batch-table regeneration cell on the same 0-100 scale as the NLRexpress `probability` column, so low-confidence MHD calls are visible before mutation/domestication

## Recommended Workflows

### CDS Available

1. Open the NLRexpress Colab notebook.
2. Set the input type to `CDS`.
3. Run NLRexpress.
4. Download the original NLRexpress result archive.
5. Review the MHD D-to-V candidates shown in the notebook or in the NLRexpress result archive.
6. Run the optional CodonDomesticate handoff cell if you want a batch domestication CSV.
7. Open the CodonDomesticate notebook and upload that CSV for batch domestication/mutation.

### Protein Only

1. Open the NLRexpress Colab notebook.
2. Set the input type to `protein`.
3. Run NLRexpress.
4. Download the original NLRexpress result archive.
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
- `codon_domesticate_input_with_mhd_dv.csv` from the optional batch-table regeneration cell; this contains the `name`, `sequence`, and `aa_change` columns expected by CodonDomesticate batch mode. When a probability threshold is set, it also includes `mhd_probability_warning_threshold` and `mhd_probability_higher_threshold`

The notebook also saves `single_mhd_dv_candidates.csv` or `multi_mhd_dv_candidates.csv` inside the Colab session as an intermediate table, but it does not download those files automatically. The same detailed NLRexpress information is available in the raw result archive.

The original NLRexpress output files are not edited.

## Notes

- NLRexpress itself accepts protein FASTA input; this notebook translates CDS input with Biopython before running NLRexpress.
- CodonDomesticate accepts CDS input and optional `aa_change` values.
- If no MHD motif is detected, the notebook prints a warning and still writes a handoff CSV with an empty `aa_change` column.
- If some CDS records have no MHD candidate, or if a sequence has multiple MHD candidates, the notebook prints warnings. For multiple candidates, the highest-probability MHD candidate is used.
- In the optional batch-table regeneration cell, if an MHD candidate probability is below the configured warning threshold, the notebook prints a warning. The threshold uses the same 0-100 scale as the NLRexpress `probability` column and defaults to `90`. The candidate is not automatically removed; review the CSV before using it for CodonDomesticate.
- The optional batch-table regeneration cell downloads only one file: `codon_domesticate_input_with_mhd_dv.csv`. Detailed NLRexpress information remains available in the raw result archive and intermediate MHD candidate CSV.
- If only protein sequences are available, CDS domestication cannot be performed until matching CDS sequences are provided.
- The NLRexpress repository is GPL-3.0 licensed. This notebook downloads and runs NLRexpress but does not vendor its source code or modify its original result files.
- For publication or redistribution, cite NLRexpress as requested in the upstream repository.
