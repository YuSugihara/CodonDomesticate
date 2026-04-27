# CodonDomesticate

This script extracts only the domestication-related workflow from
`CodonTransformer_YS_v2026_02_17.ipynb`.

It can:

- remove restriction enzyme sites from an existing CDS by silent mutation
- remove user-defined DNA motifs using the same silent-mutation strategy
- process multiple CDS records from a CSV file
- optionally introduce one amino-acid substitution such as `D512V`

It does not run codon optimization and does not use the CodonTransformer deep
learning model. Kazusa codon usage is fetched directly with
`python-codon-tables`, which is the same lightweight package used internally by
CodonTransformer for this task.

## Files

- `domesticate_cds.py`: command-line script and importable Python API for single-CDS and batch-CSV domestication
- `notebooks/Domesticate_CDS_Colab.ipynb`: Google Colab notebook containing both single-CDS and multiple-CDS workflows

## Install

```bash
pip install biopython python-codon-tables
```

If you do not want to use codon usage frequencies and any synonymous codon is
acceptable, `python-codon-tables` is not required. Use `--no-codon-usage` in
that case.

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

To avoid organism-name lookup through NCBI Entrez, provide the taxonomy ID
directly.

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

The input CSV must contain a `sequence` column. The `name` and `aa_change`
columns are optional.

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

Use `--aa-change D512V` to first remove forbidden motifs by silent mutation and
then introduce the requested amino-acid substitution.

```bash
python domesticate_cds.py single \
  --cds input_cds.fasta \
  --enzyme-names "BsaI;BpiI;Esp3I" \
  --aa-change D512V
```

For batch mode, add an `aa_change` column to the CSV. Leave the value empty for
silent domestication only.

## Google Colab

Open `Domesticate_CDS_Colab.ipynb` in Google Colab. The notebook includes:

- dependency installation
- shared domestication functions
- a single-CDS section
- a multiple-CDS CSV section

The Colab notebook installs `biopython` and `python-codon-tables`. It does not
install or import CodonTransformer.

## Notes

- The CDS length must be a multiple of 3.
- The translated amino-acid sequence is preserved before and after silent domestication.
- If a forbidden motif cannot be removed by synonymous mutation, the batch output reports `status=error`.
- `python-codon-tables` can fetch Kazusa-derived codon usage tables by taxonomy ID. Internet access is required in Colab.


## GitHub / Colab Source Download

The Colab notebook downloads `domesticate_cds.py` from GitHub at runtime. After
creating the GitHub repository, update the `github_raw_url` parameter in the
notebook if your repository URL is different from:

```text
https://raw.githubusercontent.com/YuSugihara/CodonDomesticate/main/domesticate_cds.py
```

A typical setup is:

```bash
git remote add origin https://github.com/YuSugihara/CodonDomesticate.git
git branch -M main
git push -u origin main
```
