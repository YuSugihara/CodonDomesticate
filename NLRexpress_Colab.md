# NLRexpress Colab Workflow

[![Open NLRexpress In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/NLRexpress_Colab.ipynb)

Use this notebook to run NLRexpress, review MHD candidates, and prepare a CodonDomesticate batch CSV.

## CDS Workflow

1. Run NLRexpress with CDS input.
2. Review the MHD candidates.
3. Run the optional CodonDomesticate handoff cell.
4. Download `codon_domesticate_input_with_mhd_mutations.csv`.
5. Upload that CSV in the CodonDomesticate notebook.

For CDS input, the NLRexpress result zip includes the matching CDS table. You can upload that zip later to regenerate the CodonDomesticate CSV without uploading a separate CDS file.

## Existing Result Zip

1. Run setup in the NLRexpress notebook.
2. Set `handoff_mhd_source = "upload_nlrexpress_zip"`.
3. Upload the NLRexpress result zip.
4. Download the CodonDomesticate CSV.

Older zips and protein-only runs still need a matching CDS CSV.

## Input

Use a single sequence, FASTA, or CSV:

```csv
name,sequence
example_NLR_1,ATG...
example_NLR_2,ATG...
```

## Output

The handoff file for CodonDomesticate is:

```text
codon_domesticate_input_with_mhd_mutations.csv
```

It contains `name`, `sequence`, and `aa_change`. Review `aa_change` before running CodonDomesticate.

[![Open CodonDomesticate In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/Domesticate_CDS_Colab.ipynb)
