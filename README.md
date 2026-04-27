# CodonDomesticate

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/Domesticate_CDS_Colab.ipynb)

CodonDomesticate removes restriction enzyme sites and user-defined DNA motifs from coding sequences by silent mutation. It can also introduce an optional amino-acid substitution, such as an MHD D-to-V mutation, while keeping the rest of the protein sequence unchanged.

## Features

- remove restriction enzyme sites from a CDS by silent mutation
- remove custom DNA motifs using the same silent-mutation strategy
- process one CDS or many CDS records from a CSV file
- optionally introduce one amino-acid substitution such as `D512V`
- use Kazusa codon usage frequencies through `python-codon-tables`
- run from the command line or in Google Colab

## Google Colab

Use the Colab notebook for the easiest workflow:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/main/notebooks/Domesticate_CDS_Colab.ipynb)

The notebook supports both single-CDS and batch-CSV domestication.

## NLRexpress Colab Workflow

A separate Colab notebook can run NLRexpress and prepare optional MHD D-to-V mutation handoff files for CodonDomesticate.

[![Open NLRexpress In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YuSugihara/CodonDomesticate/blob/develop/notebooks/NLRexpress_Colab.ipynb)

See `NLRexpress_Colab.md` for details.

## Install

```bash
pip install biopython python-codon-tables
```

If you do not want to use codon usage frequencies and any synonymous codon is acceptable, `python-codon-tables` is not required. Use `--no-codon-usage` in that case.

## Single CDS

```bash
python domesticate_cds.py single \
  --cds ATGGGTCTCTAA \
  --enzyme-names "BsaI;BpiI;Esp3I" \
  --organism "Nicotiana benthamiana"
```

You can also read the CDS from a text or FASTA file.

```bash
python domesticate_cds.py single \
  --cds input_cds.fasta \
  --enzyme-names "BsaI;BpiI;Esp3I" \
  --avoiding-motifs "AAAAAA;GCGGCCGC" \
  --output single_domesticated.csv
```

To avoid organism-name lookup through NCBI Entrez, provide the taxonomy ID directly.

```bash
python domesticate_cds.py single \
  --cds input_cds.txt \
  --taxonomy-id 4097 \
  --enzyme-names "BsaI;BpiI;Esp3I"
```

Minimal example without codon usage frequencies:

```bash
python domesticate_cds.py single \
  --cds ATGGGTCTCTAA \
  --enzyme-names BsaI \
  --no-codon-usage
```

## Multiple CDS

The input CSV must contain a `sequence` column. The `name` and `aa_change` columns are optional.

```csv
name,sequence,aa_change
example_1,ATGGGTCTCTAA,
example_2,ATGGGTCTCTAA,
```

Run batch domestication:

```bash
python domesticate_cds.py batch \
  --input-csv input_cds.csv \
  --output-csv domesticated_cds_results.csv \
  --enzyme-names "BsaI;BpiI;Esp3I" \
  --organism "Nicotiana benthamiana"
```

The output CSV preserves the original columns and adds:

- `domesticated_cds`
- `final_cds`
- `protein_before`
- `protein_after`
- `num_sites_before`
- `site_hits_before`
- `num_sites_after`
- `site_hits_after`
- `num_mutations`
- `mutations`
- `status`
- `error`

## Optional Amino-Acid Mutation

Use `--aa-change D512V` to first remove forbidden motifs by silent mutation and then introduce the requested amino-acid substitution.

```bash
python domesticate_cds.py single \
  --cds input_cds.fasta \
  --enzyme-names "BsaI;BpiI;Esp3I" \
  --aa-change D512V
```

For batch mode, add an `aa_change` column to the CSV. Leave the value empty for silent domestication only.

## Notes

- The CDS length must be a multiple of 3.
- The translated amino-acid sequence is preserved before and after silent domestication.
- If a forbidden motif cannot be removed by synonymous mutation, the batch output reports `status=error`.
- `python-codon-tables` can fetch Kazusa-derived codon usage tables by taxonomy ID. Internet access is required in Colab.
