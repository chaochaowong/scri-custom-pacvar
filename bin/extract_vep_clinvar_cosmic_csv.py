#!/usr/bin/env python3
"""Extract pathogenic ClinVar and COSMIC VEP annotations from a VCF."""

import argparse
import csv
import gzip
import re
import sys


CORE_CSQ_FIELDS = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "Existing_variation",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "CANONICAL",
    "MANE_SELECT",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "MAX_AF_POPS",
]

OUTPUT_FIELDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "ALT_INDEX",
    "QUAL",
    "FILTER",
    "HGVSg",
] + CORE_CSQ_FIELDS + [
    "SAMPLE",
    "GT",
    "GQ",
    "DP",
    "AD",
    "AD_REF",
    "AD_ALT",
    "REF_ALLELE_FRACTION",
    "ALT_ALLELE_FRACTION",
    "VAF",
    "PS",
]

COSMIC_RE = re.compile(r"(?:^|&)COS(?:V|M|N)\d+(?:$|&)", re.IGNORECASE)
PATHOGENIC_TERMS = {
    "pathogenic",
    "likely_pathogenic",
    "conflicting_classifications_of_pathogenicity",
    "conflicting_interpretations_of_pathogenicity",
}


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="VEP-annotated VCF or VCF.GZ")
    parser.add_argument("--clinvar-output", required=True, help="ClinVar CSV")
    parser.add_argument("--cosmic-output", required=True, help="COSMIC CSV")
    parser.add_argument(
        "--include-non-carriers",
        action="store_true",
        help="Also report annotations for ALT alleles absent from the sample GT",
    )
    return parser.parse_args()


def minimal_vep_allele(ref, alt):
    """Return VEP's CSQ Allele representation for a VCF REF/ALT pair."""
    if len(ref) != len(alt) and ref and alt and ref[0] == alt[0]:
        return alt[1:] or "-"
    return alt


def matching_alt_indices(ref, alts, vep_allele):
    exact = [i for i, alt in enumerate(alts, 1) if alt == vep_allele]
    if exact:
        return exact
    minimal = [
        i for i, alt in enumerate(alts, 1)
        if minimal_vep_allele(ref, alt) == vep_allele
    ]
    if minimal:
        return minimal
    return [1] if len(alts) == 1 else []


def make_hgvsg(chrom, pos, ref, alt):
    """Create a readable VCF-derived genomic HGVS expression."""
    start = int(pos)
    r, a = ref, alt
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]
        start += 1
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]

    prefix = f"{chrom}:g."
    if len(r) == 1 and len(a) == 1:
        return f"{prefix}{start}{r}>{a}"
    if r and not a:
        end = start + len(r) - 1
        locus = str(start) if start == end else f"{start}_{end}"
        return f"{prefix}{locus}del"
    if a and not r:
        return f"{prefix}{start - 1}_{start}ins{a}"
    if r and a:
        end = start + len(r) - 1
        locus = str(start) if start == end else f"{start}_{end}"
        return f"{prefix}{locus}delins{a}"
    return ""


def split_info(info):
    result = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = True
    return result


def value_at(values, index):
    return values[index] if index < len(values) else ""


def numeric_fraction(numerator, denominator):
    try:
        denominator = float(denominator)
        if denominator == 0:
            return ""
        return f"{float(numerator) / denominator:.6g}"
    except (TypeError, ValueError):
        return ""


def genotype_has_alt(gt, alt_index):
    alleles = re.split(r"[/|]", gt)
    return str(alt_index) in alleles


def build_row(fixed, annotation, csq_index, sample, format_map, alt, alt_index):
    chrom, pos, record_id, ref, _, qual, filt = fixed[:7]
    ad = format_map.get("AD", "").split(",")
    format_vaf = format_map.get("VAF", "").split(",")
    ad_ref = value_at(ad, 0)
    ad_alt = value_at(ad, alt_index)
    called_depth = ""
    try:
        called_depth = sum(float(x) for x in ad if x not in {"", "."})
    except ValueError:
        pass

    row = {
        "CHROM": chrom,
        "POS": pos,
        "ID": record_id,
        "REF": ref,
        "ALT": alt,
        "ALT_INDEX": alt_index,
        "QUAL": qual,
        "FILTER": filt,
        "HGVSg": make_hgvsg(chrom, pos, ref, alt),
        "SAMPLE": sample,
        "GT": format_map.get("GT", ""),
        "GQ": format_map.get("GQ", ""),
        "DP": format_map.get("DP", ""),
        "AD": format_map.get("AD", ""),
        "AD_REF": ad_ref,
        "AD_ALT": ad_alt,
        "REF_ALLELE_FRACTION": numeric_fraction(ad_ref, called_depth),
        "ALT_ALLELE_FRACTION": numeric_fraction(ad_alt, called_depth),
        "VAF": value_at(format_vaf, alt_index - 1),
        "PS": format_map.get("PS", ""),
    }
    for field in CORE_CSQ_FIELDS:
        index = csq_index.get(field)
        row[field] = value_at(annotation, index) if index is not None else ""
    return row


def main():
    args = parse_args()
    csq_fields = None
    sample_names = []
    clinvar_rows = cosmic_rows = unmatched_annotations = noncarrier_rows = 0
    unmatched_examples = []

    with open_text(args.vcf) as source, \
            open(args.clinvar_output, "w", newline="") as clinvar_handle, \
            open(args.cosmic_output, "w", newline="") as cosmic_handle:
        clinvar_writer = csv.DictWriter(
            clinvar_handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore"
        )
        cosmic_writer = csv.DictWriter(
            cosmic_handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore"
        )
        clinvar_writer.writeheader()
        cosmic_writer.writeheader()

        for line in source:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r"Format: ([^\">]+)", line)
                if not match:
                    sys.exit("Could not parse the VEP CSQ Format declaration")
                csq_fields = match.group(1).split("|")
                csq_index = {name: i for i, name in enumerate(csq_fields)}
                missing = {"Allele", "Existing_variation", "CLIN_SIG"} - set(csq_index)
                if missing:
                    sys.exit(f"Required CSQ fields missing: {sorted(missing)}")
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            if csq_fields is None:
                sys.exit("VCF has no parseable VEP CSQ header")

            columns = line.rstrip("\n").split("\t")
            fixed = columns[:8]
            ref = columns[3]
            alts = columns[4].split(",")
            info = split_info(columns[7])
            if "CSQ" not in info:
                continue

            format_keys = columns[8].split(":") if len(columns) > 8 else []
            sample_columns = columns[9:] or [""]
            names = sample_names or [""]

            for csq_text in info["CSQ"].split(","):
                annotation = csq_text.split("|")
                clin_sig = value_at(annotation, csq_index["CLIN_SIG"])
                existing = value_at(annotation, csq_index["Existing_variation"])
                clin_tokens = {
                    token.lower() for token in clin_sig.split("&") if token
                }
                is_clinvar = bool(clin_tokens & PATHOGENIC_TERMS)
                is_cosmic = bool(COSMIC_RE.search(existing))
                if not (is_clinvar or is_cosmic):
                    continue

                vep_allele = value_at(annotation, csq_index["Allele"])
                alt_indices = matching_alt_indices(ref, alts, vep_allele)
                if not alt_indices:
                    unmatched_annotations += 1
                    if len(unmatched_examples) < 10:
                        unmatched_examples.append(
                            f"{columns[0]}:{columns[1]} "
                            f"{ref}>{','.join(alts)} CSQ_Allele={vep_allele}"
                        )
                    continue

                for sample, sample_text in zip(names, sample_columns):
                    sample_values = sample_text.split(":")
                    format_map = dict(zip(format_keys, sample_values))
                    for alt_index in alt_indices:
                        if (
                            not args.include_non_carriers
                            and not genotype_has_alt(
                                format_map.get("GT", ""), alt_index
                            )
                        ):
                            noncarrier_rows += 1
                            continue
                        row = build_row(
                            fixed, annotation, csq_index, sample, format_map,
                            alts[alt_index - 1], alt_index
                        )
                        if is_clinvar:
                            clinvar_writer.writerow(row)
                            clinvar_rows += 1
                        if is_cosmic:
                            cosmic_writer.writerow(row)
                            cosmic_rows += 1

    print(f"ClinVar rows: {clinvar_rows}", file=sys.stderr)
    print(f"COSMIC rows: {cosmic_rows}", file=sys.stderr)
    print(f"Non-carrier annotation rows skipped: {noncarrier_rows}", file=sys.stderr)
    print(
        f"Matching annotations with unresolved ALT: {unmatched_annotations}",
        file=sys.stderr,
    )
    for example in unmatched_examples:
        print(f"  unresolved: {example}", file=sys.stderr)


if __name__ == "__main__":
    main()
