#!/usr/bin/env python3
"""Domesticate CDS sequences with silent mutations.

This script removes forbidden DNA motifs from coding sequences while preserving
the translated amino-acid sequence. It can also apply one requested amino-acid
substitution after domestication.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import Entrez
from Bio.Restriction import Restriction
from Bio.Seq import Seq


CodonUsageByAA = Dict[str, Dict[str, float]]
RestrictionSites = Dict[str, str]


def normalize_dna(seq: str) -> str:
    """Return an uppercase A/C/G/T sequence from plain text or FASTA text."""
    lines = []
    for line in str(seq).splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        lines.append(line)
    dna = "".join(lines).upper().replace("U", "T").replace(" ", "")
    invalid = sorted(set(dna) - set("ACGT"))
    if invalid:
        raise ValueError(f"CDS contains invalid DNA bases: {''.join(invalid)}")
    if len(dna) % 3 != 0:
        raise ValueError("CDS length is not a multiple of 3.")
    return dna


def read_sequence_arg(value: str) -> str:
    """Read a CDS from a literal argument or from a file path."""
    path = Path(value)
    if path.exists():
        return normalize_dna(path.read_text())
    return normalize_dna(value)


def get_taxid_from_name(name: str, email: Optional[str] = None) -> Optional[int]:
    """Resolve an organism name to an NCBI taxonomy ID."""
    if email:
        Entrez.email = email
    handle = Entrez.esearch(db="taxonomy", term=name)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    if not ids:
        return None
    return int(ids[0])


def download_kazusa_codon_usage(
    organism: str,
    taxonomy_id: Optional[int] = None,
    email: Optional[str] = None,
) -> CodonUsageByAA:
    """Download Kazusa-format codon usage through python-codon-tables."""
    try:
        import python_codon_tables as pct
    except ImportError as exc:
        raise ImportError(
            "python-codon-tables is required for codon usage download. "
            "Install it with: pip install python-codon-tables"
        ) from exc

    taxid = taxonomy_id if taxonomy_id is not None else get_taxid_from_name(organism, email=email)
    if taxid is None:
        raise ValueError(f"Could not find taxonomy ID for organism: {organism}")

    return pct.get_codons_table(table_name=taxid)


def get_enzyme_restriction_site(enzyme_name: str) -> str:
    """Return the recognition sequence for a Bio.Restriction enzyme name."""
    try:
        enzyme_class = getattr(Restriction, enzyme_name)
    except AttributeError as exc:
        raise ValueError(f"Unknown restriction enzyme: {enzyme_name}") from exc

    site = enzyme_class.site
    if not site:
        raise ValueError(f"No restriction site found for enzyme: {enzyme_name}")
    return site.upper()


def get_restriction_sites_from_enzyme_names(enzyme_names: str) -> RestrictionSites:
    """Parse semicolon-separated enzyme names into enzyme -> site."""
    sites: RestrictionSites = {}
    if not enzyme_names:
        return sites

    for raw_name in enzyme_names.split(";"):
        name = raw_name.strip()
        if not name:
            continue
        sites[name] = get_enzyme_restriction_site(name)
    return sites


def get_avoiding_motifs_from_text(avoiding_motifs: str) -> RestrictionSites:
    """Parse semicolon-separated custom motifs into label -> motif."""
    motifs: RestrictionSites = {}
    if not avoiding_motifs:
        return motifs

    motif_index = 1
    for raw_motif in avoiding_motifs.split(";"):
        motif = raw_motif.strip().upper().replace("U", "T")
        if not motif:
            continue
        if any(base not in "ACGT" for base in motif):
            raise ValueError(f"Invalid avoiding motif '{raw_motif}'. Use only A/C/G/T.")
        motifs[f"motif_{motif_index}:{motif}"] = motif
        motif_index += 1
    return motifs


def build_forbidden_motifs(
    restriction_sites: Optional[RestrictionSites] = None,
    avoiding_motifs: Optional[RestrictionSites] = None,
) -> RestrictionSites:
    """Merge restriction sites and custom motifs into one dictionary."""
    forbidden: RestrictionSites = {}
    if restriction_sites:
        forbidden.update(restriction_sites)
    if avoiding_motifs:
        for label, motif in avoiding_motifs.items():
            safe_label = label
            suffix = 2
            while safe_label in forbidden:
                safe_label = f"{label}#{suffix}"
                suffix += 1
            forbidden[safe_label] = motif
    return forbidden


def codon_usage_from_kazusa(kazusa_codon_freqs: CodonUsageByAA) -> Dict[str, float]:
    """Convert AA -> {codon: freq} into codon -> freq."""
    codon_usage: Dict[str, float] = {}
    for codon_dict in kazusa_codon_freqs.values():
        for codon, freq in codon_dict.items():
            codon_usage[codon.upper()] = float(freq)
    return codon_usage


def _scan_forbidden_motifs(seq: str, forbidden_sites: RestrictionSites) -> Tuple[set, List[dict]]:
    """Scan sequence for motifs on forward and reverse-complement strands."""
    seq = seq.upper()
    seq_len = len(seq)
    site_positions = set()
    site_hits = []

    for name, site in forbidden_sites.items():
        if not site:
            continue
        site = site.upper()
        rc_site = str(Seq(site).reverse_complement()).upper()

        for strand, motif in (("forward", site), ("reverse", rc_site)):
            if strand == "reverse" and rc_site == site:
                continue
            start = 0
            while True:
                i = seq.find(motif, start)
                if i == -1:
                    break
                j = i + len(motif)
                site_positions.update(range(max(0, i), min(seq_len, j)))
                site_hits.append(
                    {
                        "enzyme": name,
                        "strand": strand,
                        "start": i + 1,
                        "end": j,
                        "site": motif,
                        "start0": i,
                        "end0": j,
                        "match_seq": seq[i:j],
                    }
                )
                start = i + 1

    return site_positions, site_hits


def scan_forbidden_motifs(seq: str, forbidden_sites: RestrictionSites) -> List[dict]:
    """Return all forbidden motif hits in a CDS."""
    _, hits = _scan_forbidden_motifs(normalize_dna(seq), forbidden_sites)
    return hits


def _attempt_fix_hit(
    seq: str,
    hit: dict,
    forbidden_sites: RestrictionSites,
    codon_usage_table: Optional[Dict[str, float]] = None,
    min_codon_freq: float = 0.2,
    seen_sequences: Optional[set] = None,
) -> Optional[dict]:
    """Try to destroy one motif with one silent mutation."""
    seq = seq.upper()
    start0 = hit["start0"]
    end0 = hit["end0"]
    old_subseq = hit["match_seq"]
    _, current_hits = _scan_forbidden_motifs(seq, forbidden_sites)
    current_hit_count = len(current_hits)
    best_candidate = None

    for pass_idx in (1, 2):
        for pos in range(start0, end0):
            codon_start = (pos // 3) * 3
            codon = seq[codon_start : codon_start + 3]
            aa = str(Seq(codon).translate())
            if aa == "*":
                continue

            for base in "ACGT":
                if base == seq[pos]:
                    continue

                codon_list = list(codon)
                codon_list[pos - codon_start] = base
                new_codon = "".join(codon_list)
                if str(Seq(new_codon).translate()) != aa:
                    continue

                new_freq = None
                if codon_usage_table is not None:
                    new_freq = codon_usage_table.get(new_codon.upper())
                if codon_usage_table is not None and pass_idx == 1:
                    if new_freq is None or new_freq < min_codon_freq:
                        continue

                new_seq = seq[:codon_start] + new_codon + seq[codon_start + 3 :]
                if new_seq[start0:end0] == old_subseq:
                    continue
                if seen_sequences is not None and new_seq in seen_sequences:
                    continue

                _, new_hits = _scan_forbidden_motifs(new_seq, forbidden_sites)
                new_hit_count = len(new_hits)
                if new_hit_count > current_hit_count:
                    continue

                score = new_freq if new_freq is not None else 0.0
                if best_candidate is None or (
                    new_hit_count,
                    -score,
                    pos,
                    base,
                ) < (
                    best_candidate["new_hit_count"],
                    -best_candidate["score"],
                    best_candidate["pos"],
                    best_candidate["base_new"],
                ):
                    best_candidate = {
                        "new_seq": new_seq,
                        "pos": pos,
                        "base_old": seq[pos],
                        "base_new": base,
                        "codon_start": codon_start,
                        "old_codon": codon,
                        "new_codon": new_codon,
                        "new_freq": new_freq,
                        "score": score,
                        "new_hit_count": new_hit_count,
                    }

        if best_candidate is not None:
            break

    return best_candidate


def domesticate_cds_with_silent_mutations(
    cds: str,
    forbidden_sites: RestrictionSites,
    kazusa_codon_freqs: Optional[CodonUsageByAA] = None,
    min_codon_freq: float = 0.2,
) -> Tuple[str, dict]:
    """Remove forbidden motifs from a CDS by silent mutations."""
    cds = normalize_dna(cds)
    codon_usage_table = None
    if kazusa_codon_freqs is not None:
        codon_usage_table = codon_usage_from_kazusa(kazusa_codon_freqs)

    protein_seq = str(Seq(cds).translate())
    _, site_hits_initial = _scan_forbidden_motifs(cds, forbidden_sites)

    seq = cds
    mutations = []
    max_steps = max(1, 10 * len(site_hits_initial))
    seen_sequences = {seq}

    for step in range(max_steps):
        _, site_hits = _scan_forbidden_motifs(seq, forbidden_sites)
        if not site_hits:
            break

        hit = site_hits[0]
        candidate = _attempt_fix_hit(
            seq,
            hit,
            forbidden_sites=forbidden_sites,
            codon_usage_table=codon_usage_table,
            min_codon_freq=min_codon_freq,
            seen_sequences=seen_sequences,
        )
        if candidate is None:
            raise ValueError(
                f"Could not domesticate forbidden motif '{hit['enzyme']}' "
                f"at {hit['start']}-{hit['end']} with silent mutations."
            )

        pos = candidate["pos"]
        old_codon = candidate["old_codon"]
        new_codon = candidate["new_codon"]
        original_freq = None
        new_freq = None
        if codon_usage_table is not None:
            original_freq = codon_usage_table.get(old_codon.upper())
            new_freq = codon_usage_table.get(new_codon.upper())

        mutations.append(
            {
                "position": pos + 1,
                "codon_index": (pos // 3) + 1,
                "codon_pos": (pos % 3) + 1,
                "original_base": candidate["base_old"],
                "new_base": candidate["base_new"],
                "original_codon": old_codon,
                "new_codon": new_codon,
                "original_freq": original_freq,
                "new_freq": new_freq,
            }
        )
        seq = candidate["new_seq"]
        seen_sequences.add(seq)
    else:
        _, remaining_hits = _scan_forbidden_motifs(seq, forbidden_sites)
        raise ValueError(
            "Reached the maximum number of domestication steps. "
            "The algorithm may be stuck in a loop. "
            f"Remaining forbidden motifs: {len(remaining_hits)}"
            f"{' (' + format_hits(remaining_hits[:5]) + ')' if remaining_hits else ''}."
        )

    domesticated = seq
    if str(Seq(domesticated).translate()) != protein_seq:
        raise ValueError("Amino acid sequence changed during domestication.")

    _, remaining_hits = _scan_forbidden_motifs(domesticated, forbidden_sites)
    if remaining_hits:
        raise ValueError("Some forbidden motifs remain after domestication.")

    report_data = {
        "protein": protein_seq,
        "input_cds": cds,
        "domesticated_cds": domesticated,
        "site_hits": site_hits_initial,
        "mutations": mutations,
        "remaining_hits": remaining_hits,
    }
    return domesticated, report_data


def _build_aa_to_codons() -> Dict[str, set]:
    aa_to_codons: Dict[str, set] = {}
    for b1 in "TCAG":
        for b2 in "TCAG":
            for b3 in "TCAG":
                codon = b1 + b2 + b3
                aa = str(Seq(codon).translate())
                aa_to_codons.setdefault(aa, set()).add(codon)
    return aa_to_codons


_AA_TO_CODONS = _build_aa_to_codons()


def mutate_amino_acid_with_constraints(
    cds: str,
    aa_change: str,
    forbidden_sites: RestrictionSites,
    codon_usage_by_aa: Optional[CodonUsageByAA] = None,
    min_codon_freq: float = 0.2,
) -> Tuple[str, dict]:
    """Apply a D456V-like amino-acid substitution while avoiding motifs."""
    cds = normalize_dna(cds)
    match = re.match(r"^([A-Z\*])(\d+)([A-Z\*])$", aa_change.strip())
    if not match:
        raise ValueError("aa_change must look like 'D456V'.")

    orig_aa = match.group(1)
    aa_index = int(match.group(2))
    new_aa = match.group(3)
    protein = str(Seq(cds).translate())
    if aa_index < 1 or aa_index > len(protein):
        raise ValueError(f"Amino acid index {aa_index} is out of range.")
    if protein[aa_index - 1] != orig_aa:
        raise ValueError(
            f"At position {aa_index}, protein has '{protein[aa_index - 1]}' "
            f"but aa_change expects '{orig_aa}'."
        )
    if new_aa not in _AA_TO_CODONS:
        raise ValueError(f"No codons found for amino acid '{new_aa}'.")

    codon_start = (aa_index - 1) * 3
    original_codon = cds[codon_start : codon_start + 3]

    def get_freq(aa: str, codon: str) -> Optional[float]:
        if codon_usage_by_aa is None:
            return None
        return codon_usage_by_aa.get(aa, {}).get(codon.upper())

    def search_candidates(freq_filtered: bool) -> List[dict]:
        candidates = []
        for codon in sorted(_AA_TO_CODONS[new_aa]):
            if codon == original_codon:
                continue
            freq = get_freq(new_aa, codon)
            if freq_filtered and codon_usage_by_aa is not None:
                if freq is None or freq < min_codon_freq:
                    continue
            new_cds = cds[:codon_start] + codon + cds[codon_start + 3 :]
            new_protein = str(Seq(new_cds).translate())
            if new_protein[: aa_index - 1] != protein[: aa_index - 1]:
                continue
            if new_protein[aa_index - 1] != new_aa:
                continue
            if new_protein[aa_index:] != protein[aa_index:]:
                continue
            _, sites_after = _scan_forbidden_motifs(new_cds, forbidden_sites)
            if sites_after:
                continue
            diffs = sum(c1 != c2 for c1, c2 in zip(original_codon, codon))
            candidates.append(
                {
                    "codon": codon,
                    "new_cds": new_cds,
                    "freq": freq,
                    "diffs": diffs,
                }
            )
        return candidates

    codon_usage_preferred_candidates = search_candidates(freq_filtered=True)
    used_codon_usage_preferred_pool = bool(codon_usage_preferred_candidates)
    candidates = codon_usage_preferred_candidates or search_candidates(freq_filtered=False)
    if not candidates:
        raise ValueError(
            f"Could not find a codon for {new_aa}{aa_index} that avoids forbidden motifs."
        )

    best = min(candidates, key=lambda c: (c["diffs"], -(c["freq"] or 0.0), c["codon"]))
    new_cds = best["new_cds"]
    new_codon = best["codon"]
    nucleotide_changes = []
    for i, (old_base, new_base) in enumerate(zip(original_codon, new_codon)):
        if old_base != new_base:
            nucleotide_changes.append(
                {
                    "position": codon_start + i + 1,
                    "from_base": old_base,
                    "to_base": new_base,
                }
            )

    _, sites_before = _scan_forbidden_motifs(cds, forbidden_sites)
    _, sites_after = _scan_forbidden_motifs(new_cds, forbidden_sites)
    codon_usage_checked = codon_usage_by_aa is not None
    new_freq_passes_min = None
    if codon_usage_checked:
        new_freq_passes_min = best["freq"] is not None and best["freq"] >= min_codon_freq
    info = {
        "aa_change": aa_change,
        "aa_index": aa_index,
        "original_aa": orig_aa,
        "new_aa": new_aa,
        "original_codon": original_codon,
        "new_codon": new_codon,
        "codon_index": aa_index,
        "nucleotide_changes": nucleotide_changes,
        "original_freq": get_freq(orig_aa, original_codon),
        "new_freq": best["freq"],
        "codon_usage_checked": codon_usage_checked,
        "min_codon_freq": min_codon_freq,
        "new_freq_passes_min": new_freq_passes_min,
        "used_codon_usage_preferred_pool": used_codon_usage_preferred_pool,
        "sites_before": sites_before,
        "sites_after": sites_after,
    }
    return new_cds, info


def domesticate_or_mutate_cds(
    cds: str,
    forbidden_sites: RestrictionSites,
    kazusa_codon_freqs: Optional[CodonUsageByAA] = None,
    min_codon_freq: float = 0.2,
    aa_change: Optional[str] = None,
) -> Tuple[str, dict]:
    """Domesticate a CDS and optionally apply one amino-acid substitution."""
    cds = normalize_dna(cds)
    domesticated_cds, dom_report = domesticate_cds_with_silent_mutations(
        cds=cds,
        forbidden_sites=forbidden_sites,
        kazusa_codon_freqs=kazusa_codon_freqs,
        min_codon_freq=min_codon_freq,
    )

    final_cds = domesticated_cds
    aa_change_info = None
    aa_change_str = str(aa_change).strip() if aa_change is not None else ""
    if aa_change_str:
        final_cds, aa_change_info = mutate_amino_acid_with_constraints(
            cds=domesticated_cds,
            aa_change=aa_change_str,
            forbidden_sites=forbidden_sites,
            codon_usage_by_aa=kazusa_codon_freqs,
            min_codon_freq=min_codon_freq,
        )

    return final_cds, {
        "input_cds": cds,
        "domesticated_cds": domesticated_cds,
        "final_cds": final_cds,
        "protein_before": str(Seq(cds).translate()),
        "protein_after": str(Seq(final_cds).translate()),
        "domestication_report": dom_report,
        "aa_change": aa_change_str or None,
        "aa_change_info": aa_change_info,
    }


def format_hits(hits: Iterable[dict]) -> str:
    return "; ".join(
        "{enzyme}:{strand}:{start}-{end}".format(
            enzyme=h["enzyme"],
            strand=h["strand"],
            start=h["start"],
            end=h["end"],
        )
        for h in hits
    )


def format_optional_float(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return f"{value:.2f}"


def format_codon_with_freq(codon: str, freq: Optional[float]) -> str:
    return f"{codon}({format_optional_float(freq)})"


def highlight_codon_base(codon: str, codon_pos: int) -> str:
    index = codon_pos - 1
    if index < 0 or index >= len(codon):
        return codon
    return f"{codon[:index]}[{codon[index]}]{codon[index + 1:]}"


def highlight_codon_bases(codon: str, codon_positions: Iterable[int]) -> str:
    indexes = {pos - 1 for pos in codon_positions}
    return "".join(f"[{base}]" if index in indexes else base for index, base in enumerate(codon))


def format_codon_change_with_freq(mutation: dict) -> str:
    return (
        f"{highlight_codon_base(mutation['original_codon'], mutation['codon_pos'])}"
        f"({format_optional_float(mutation.get('original_freq'))}) -> "
        f"{highlight_codon_base(mutation['new_codon'], mutation['codon_pos'])}"
        f"({format_optional_float(mutation.get('new_freq'))})"
    )


def format_mutations(mutations: Iterable[dict]) -> str:
    return "; ".join(
        (
            "pos{position}(codon{codon_index},pos{codon_pos}):"
            "{codon_change}"
        ).format(
            codon_change=format_codon_change_with_freq(m),
            **m,
        )
        for m in mutations
    )


def aa_change_nucleotide_mutations(aa_change_info: Optional[dict]) -> List[dict]:
    if not aa_change_info:
        return []
    mutations = []
    for change in aa_change_info.get("nucleotide_changes", []):
        position = change["position"]
        mutations.append(
            {
                "position": position,
                "codon_index": aa_change_info["codon_index"],
                "codon_pos": ((position - 1) % 3) + 1,
                "original_base": change["from_base"],
                "new_base": change["to_base"],
                "original_codon": aa_change_info["original_codon"],
                "new_codon": aa_change_info["new_codon"],
                "original_freq": aa_change_info.get("original_freq"),
                "new_freq": aa_change_info.get("new_freq"),
            }
        )
    return mutations


def format_aa_change_mutations(aa_change_info: Optional[dict]) -> str:
    if not aa_change_info:
        return ""
    nucleotide_changes = aa_change_info.get("nucleotide_changes", [])
    if not nucleotide_changes:
        return ""

    positions = [change["position"] for change in nucleotide_changes]
    codon_positions = [((position - 1) % 3) + 1 for position in positions]
    position_text = ",".join(f"pos{position}" for position in positions)
    codon_pos_text = ",".join(f"pos{codon_pos}" for codon_pos in codon_positions)
    original_codon = aa_change_info["original_codon"]
    new_codon = aa_change_info["new_codon"]
    codon_change = (
        f"{highlight_codon_bases(original_codon, codon_positions)}"
        f"({format_optional_float(aa_change_info.get('original_freq'))}) -> "
        f"{highlight_codon_bases(new_codon, codon_positions)}"
        f"({format_optional_float(aa_change_info.get('new_freq'))})"
    )
    return (
        f"{position_text}(codon{aa_change_info['codon_index']},{codon_pos_text}):"
        f"{codon_change}"
    )


def aa_change_codon_usage_status(aa_change_info: Optional[dict]) -> str:
    if not aa_change_info:
        return ""
    if not aa_change_info.get("codon_usage_checked"):
        return "not_checked"
    if aa_change_info.get("new_freq") is None:
        return "missing"
    if aa_change_info.get("new_freq_passes_min"):
        return "pass"
    return "below_min"


def print_report(info: dict) -> None:
    """Print a compact human-readable report."""
    dom_report = info["domestication_report"]
    site_hits = dom_report.get("site_hits", [])
    mutations = dom_report.get("mutations", [])

    print("Input CDS:")
    print(info["input_cds"])
    print("\nDomesticated CDS:")
    print(info["domesticated_cds"])
    if info["aa_change"]:
        print("\nFinal CDS after aa_change:")
        print(info["final_cds"])

    print("\nProtein before:")
    print(info["protein_before"])
    print("\nProtein after:")
    print(info["protein_after"])

    print("\nForbidden motifs before domestication:")
    if site_hits:
        for hit in site_hits:
            print(f"- {hit['enzyme']}: {hit['strand']} strand, {hit['start']}-{hit['end']}, site={hit['site']}")
    else:
        print("- none")

    print("\nSilent mutations:")
    if mutations:
        for mutation in mutations:
            print(
                "- pos {position} (codon {codon_index}, pos {codon_pos}): "
                "{codon_change}".format(
                    codon_change=format_codon_change_with_freq(mutation),
                    **mutation,
                )
            )
    else:
        print("- none")

    if info["aa_change_info"]:
        aa_info = info["aa_change_info"]
        print("\nAmino-acid mutation:")
        print(
            f"- {aa_info['aa_change']}: codon {aa_info['codon_index']} "
            f"{format_codon_with_freq(aa_info['original_codon'], aa_info.get('original_freq'))} -> "
            f"{format_codon_with_freq(aa_info['new_codon'], aa_info.get('new_freq'))}"
        )
        print(f"  Codon usage: {aa_change_codon_usage_status(aa_info)}")
        aa_change_mutations = format_aa_change_mutations(aa_info)
        if aa_change_mutations:
            print(f"  Nucleotide changes: {aa_change_mutations}")


def build_forbidden_from_args(args: argparse.Namespace) -> RestrictionSites:
    restriction_sites = get_restriction_sites_from_enzyme_names(args.enzyme_names)
    avoiding_sites = get_avoiding_motifs_from_text(args.avoiding_motifs)
    forbidden_sites = build_forbidden_motifs(restriction_sites, avoiding_sites)
    if not forbidden_sites:
        raise ValueError("No restriction enzymes or avoiding motifs were provided.")
    return forbidden_sites


def build_codon_usage_from_args(args: argparse.Namespace) -> Optional[CodonUsageByAA]:
    if args.no_codon_usage:
        return None
    if args.codon_usage_json:
        with open(args.codon_usage_json) as handle:
            return json.load(handle)
    return download_kazusa_codon_usage(
        organism=args.organism,
        taxonomy_id=args.taxonomy_id,
        email=args.entrez_email,
    )


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--organism",
        default="Nicotiana benthamiana",
        help="Organism name for Kazusa codon usage lookup.",
    )
    parser.add_argument(
        "--taxonomy-id",
        type=int,
        default=None,
        help="NCBI taxonomy ID. If omitted, the organism name is resolved through Entrez.",
    )
    parser.add_argument(
        "--entrez-email",
        default=None,
        help="Optional email address for NCBI Entrez requests.",
    )
    parser.add_argument(
        "--enzyme-names",
        default="BsaI;BpiI",
        help='Semicolon-separated restriction enzymes, e.g. "BsaI;BpiI".',
    )
    parser.add_argument(
        "--avoiding-motifs",
        default="",
        help='Semicolon-separated custom DNA motifs to remove, e.g. "AAAAAA;GCGGCCGC".',
    )
    parser.add_argument(
        "--min-codon-freq",
        type=float,
        default=0.2,
        help="Preferred minimum codon frequency for replacement codons.",
    )
    parser.add_argument(
        "--no-codon-usage",
        action="store_true",
        help="Do not download codon usage; use any synonymous codon.",
    )
    parser.add_argument(
        "--codon-usage-json",
        default=None,
        help="Optional Kazusa-format JSON file instead of downloading codon usage.",
    )


def run_single(args: argparse.Namespace) -> int:
    cds = read_sequence_arg(args.cds)
    forbidden_sites = build_forbidden_from_args(args)
    kazusa_codon_freqs = build_codon_usage_from_args(args)
    final_cds, info = domesticate_or_mutate_cds(
        cds=cds,
        forbidden_sites=forbidden_sites,
        kazusa_codon_freqs=kazusa_codon_freqs,
        min_codon_freq=args.min_codon_freq,
        aa_change=args.aa_change,
    )

    if args.output:
        silent_mutations = info["domestication_report"]["mutations"]
        aa_change_mutations = aa_change_nucleotide_mutations(info["aa_change_info"])
        row = {
            "input_cds": info["input_cds"],
            "domesticated_cds": info["domesticated_cds"],
            "final_cds": final_cds,
            "protein_before": info["protein_before"],
            "protein_after": info["protein_after"],
            "num_sites_before": len(info["domestication_report"]["site_hits"]),
            "site_hits_before": format_hits(info["domestication_report"]["site_hits"]),
            "num_sites_after": len(scan_forbidden_motifs(final_cds, forbidden_sites)),
            "num_silent_mutations": len(silent_mutations),
            "num_aa_change_mutations": len(aa_change_mutations),
            "num_mutations": len(silent_mutations) + len(aa_change_mutations),
            "mutations": format_mutations(silent_mutations),
            "aa_change": info["aa_change"] or "",
            "aa_change_mutations": format_aa_change_mutations(info["aa_change_info"]),
            "aa_change_codon_usage_status": aa_change_codon_usage_status(info["aa_change_info"]),
        }
        with open(args.output, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        print(f"Wrote: {args.output}")
    else:
        print_report(info)
    return 0


def run_batch(args: argparse.Namespace) -> int:
    forbidden_sites = build_forbidden_from_args(args)
    kazusa_codon_freqs = build_codon_usage_from_args(args)
    with open(args.input_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        dataset = list(reader)
        input_fieldnames = reader.fieldnames or []

    if "sequence" not in input_fieldnames:
        raise ValueError("Input CSV must contain a 'sequence' column.")

    output_rows = []
    for index, row in enumerate(dataset):
        name = row.get("name") or f"row_{index}"
        aa_change = str(row.get("aa_change", "") or "").strip()

        out = dict(row)
        sites_before = None
        protein_before = ""
        try:
            cds = normalize_dna(str(row["sequence"]))
            protein_before = str(Seq(cds).translate())
            sites_before = scan_forbidden_motifs(cds, forbidden_sites)
            final_cds, info = domesticate_or_mutate_cds(
                cds=cds,
                forbidden_sites=forbidden_sites,
                kazusa_codon_freqs=kazusa_codon_freqs,
                min_codon_freq=args.min_codon_freq,
                aa_change=aa_change,
            )
            hits_after = scan_forbidden_motifs(final_cds, forbidden_sites)
            dom_report = info["domestication_report"]
            silent_mutations = dom_report["mutations"]
            aa_change_mutations = aa_change_nucleotide_mutations(info["aa_change_info"])
            mutation_details = format_mutations(silent_mutations)
            aa_change_mutation_details = format_aa_change_mutations(info["aa_change_info"])
            out.update(
                {
                    "domesticated_cds": info["domesticated_cds"],
                    "final_cds": final_cds,
                    "protein_before": info["protein_before"],
                    "protein_after": info["protein_after"],
                    "num_sites_before": len(dom_report["site_hits"]),
                    "site_hits_before": format_hits(dom_report["site_hits"]),
                    "num_sites_after": len(hits_after),
                    "site_hits_after": format_hits(hits_after),
                    "num_silent_mutations": len(silent_mutations),
                    "num_aa_change_mutations": len(aa_change_mutations),
                    "num_mutations": len(silent_mutations) + len(aa_change_mutations),
                    "mutations": mutation_details,
                    "aa_change_mutations": aa_change_mutation_details,
                    "aa_change_codon_usage_status": aa_change_codon_usage_status(info["aa_change_info"]),
                    "status": "ok",
                    "error": "",
                }
            )
            print(f"[ok] {name}")
            print(f"  Silent mutations: {mutation_details or 'none'}")
            print(f"  Amino-acid mutation changes: {aa_change_mutation_details or 'none'}")
        except Exception as exc:
            out.update(
                {
                    "domesticated_cds": "",
                    "final_cds": "",
                    "protein_before": protein_before,
                    "protein_after": "",
                    "num_sites_before": len(sites_before) if sites_before is not None else "",
                    "site_hits_before": format_hits(sites_before or []),
                    "num_sites_after": "",
                    "site_hits_after": "",
                    "num_silent_mutations": "",
                    "num_aa_change_mutations": "",
                    "num_mutations": "",
                    "mutations": "",
                    "aa_change_mutations": "",
                    "aa_change_codon_usage_status": "",
                    "status": "error",
                    "error": str(exc),
                }
            )
            print(f"[error] {name}: {exc}", file=sys.stderr)
        output_rows.append(out)

    output_fieldnames = list(input_fieldnames)
    for fieldname in [
        "domesticated_cds",
        "final_cds",
        "protein_before",
        "protein_after",
        "num_sites_before",
        "site_hits_before",
        "num_sites_after",
        "site_hits_after",
        "num_silent_mutations",
        "num_aa_change_mutations",
        "num_mutations",
        "mutations",
        "aa_change_mutations",
        "aa_change_codon_usage_status",
        "status",
        "error",
    ]:
        if fieldname not in output_fieldnames:
            output_fieldnames.append(fieldname)

    with open(args.output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote: {args.output_csv}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Domesticate CDS sequences by removing forbidden motifs with silent mutations.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    single = subparsers.add_parser("single", help="Domesticate one CDS.")
    single.add_argument(
        "--cds",
        required=True,
        help="CDS sequence, FASTA text, or path to a file containing the CDS.",
    )
    single.add_argument(
        "--aa-change",
        default="",
        help="Optional amino-acid substitution like D512V after domestication.",
    )
    single.add_argument("--output", default=None, help="Optional output CSV path.")
    add_common_args(single)
    single.set_defaults(func=run_single)

    batch = subparsers.add_parser("batch", help="Domesticate multiple CDSs from CSV.")
    batch.add_argument("--input-csv", required=True, help="Input CSV with columns: name, sequence, optional aa_change.")
    batch.add_argument("--output-csv", required=True, help="Output CSV path.")
    add_common_args(batch)
    batch.set_defaults(func=run_batch)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
