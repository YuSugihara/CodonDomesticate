# NLRexpress Colab Workflow

[![Open NLRexpress In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/develop/notebooks/NLRexpress_Colab.ipynb)

This notebook runs [NLRexpress](https://github.com/eliza-m/NLRexpress) on Google Colab and prepares optional handoff files for CodonDomesticate.

## What It Does

- creates an isolated Colab `micromamba` environment with HMMER and the older Python/scikit-learn stack expected by NLRexpress
- clones the NLRexpress GitHub repository
- downloads the NLRexpress predictor models
- runs NLRexpress on one sequence or many sequences, accepting either protein or CDS input
- downloads the original NLRexpress result folder as a zip archive
- extracts MHD motif hits from the short output without modifying the original NLRexpress files
- writes D-to-V mutation candidates such as `D489V` to a separate CSV file
- if CDS input is provided, automatically writes a CodonDomesticate-ready CSV containing the CDS and `aa_change` values
- optionally merges those candidates with a separate CDS CSV when only protein input was used for NLRexpress

## Why This Is a Separate Notebook

Passing rich state directly between separate Colab notebooks is awkward. Files are the cleanest handoff. This notebook therefore produces CSV files that can be downloaded and then uploaded to the CodonDomesticate notebook.

For the MHD autoactive workflow, this is usually enough:

1. Run NLRexpress here on CDS or protein sequences.
2. Download the original NLRexpress result archive.
3. Download `multi_mhd_dv_candidates.csv` or, when CDS input is available, the automatically generated CodonDomesticate handoff CSV.
4. Open the CodonDomesticate notebook.
5. Upload the prepared CSV and run batch domestication/mutation.

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

### CDS Input and Optional CDS Handoff

If you provide CDS input, the notebook translates each CDS for NLRexpress and keeps the original CDS for the CodonDomesticate handoff file. If you provide protein input first, you can later upload a CDS CSV with matching names:

```csv
name,sequence
example_NLR_1,ATG...
example_NLR_2,ATG...
```

The new handoff CSV will add columns such as `aa_change`, `motif_start`, `motif_seq`, `d_position`, and `probability`. The original NLRexpress output files are not edited.

## MHD D-to-V Logic

NLRexpress reports MHD hits in the short output with a motif start position and motif sequence. The conserved D is treated as the last `D` in the motif sequence. The notebook converts that residue to valine and writes the result to a separate mutation-candidate CSV, producing an `aa_change` value like:

```text
D489V
```

This value is compatible with `domesticate_cds.py` and the CodonDomesticate batch notebook.

## Notes

- NLRexpress itself accepts protein FASTA input; this notebook translates CDS input before running NLRexpress.
- CodonDomesticate accepts CDS input and optional `aa_change` values. If only protein sequences are available, CDS domestication cannot be performed until matching CDS sequences are provided.
- If CDS input is provided and no MHD motif is detected, the notebook prints a warning and still writes a handoff CSV with an empty `aa_change` column.
- The NLRexpress repository is GPL-3.0 licensed. This notebook downloads and runs NLRexpress but does not vendor its source code or modify its original result files.
- For publication or redistribution, cite NLRexpress as requested in the upstream repository.
