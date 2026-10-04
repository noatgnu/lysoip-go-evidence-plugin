#!/usr/bin/env python3
"""goatools GO enrichment test for the expected subcellular-compartment GO term(s)."""

import argparse
import contextlib
import csv
import gzip
import io
import math
import shutil
import sys
import urllib.request
from pathlib import Path

OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GAF_URL = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
USER_AGENT = "Mozilla/5.0 (compatible; cauldron-lysoip-go-evidence-plugin)"

DEFAULT_GO_TERMS = "GO:0005764,GO:0005765,GO:0043202,GO:0007041"


def fetch(url, path):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request) as response, open(path, "wb") as out_file:
        shutil.copyfileobj(response, out_file)


def ensure_go_data(go_data_dir):
    go_data_dir.mkdir(parents=True, exist_ok=True)
    obo_path = go_data_dir / "go-basic.obo"
    gaf_path = go_data_dir / "goa_human.gaf"
    gaf_gz_path = go_data_dir / "goa_human.gaf.gz"

    if not obo_path.exists():
        print(f"Downloading {OBO_URL} ...", file=sys.stderr)
        fetch(OBO_URL, obo_path)

    if not gaf_path.exists():
        print(f"Downloading {GAF_URL} ...", file=sys.stderr)
        fetch(GAF_URL, gaf_gz_path)
        with gzip.open(gaf_gz_path, "rb") as source, gaf_path.open("wb") as destination:
            shutil.copyfileobj(source, destination)
        gaf_gz_path.unlink()

    return obo_path, gaf_path


def expand_with_descendants(godag, term_ids):
    expanded = set(term_ids)
    for term_id in term_ids:
        term = godag.get(term_id)
        if term is not None:
            expanded |= term.get_all_children()
    return expanded


def main():
    parser = argparse.ArgumentParser(description="GO enrichment evidence for the expected subcellular compartment")
    parser.add_argument("--gene_lookup_file", required=True, help="Any TSV with protein and gene columns")
    parser.add_argument("--expected_go_terms", default=DEFAULT_GO_TERMS)
    parser.add_argument("--go_data_dir", default=".go_data")
    parser.add_argument("--output_folder", required=True)
    args = parser.parse_args()

    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # @step: Fetching GO reference data
    obo_path, gaf_path = ensure_go_data(Path(args.go_data_dir))

    from goatools.anno.gaf_reader import GafReader
    from goatools.go_enrichment import GOEnrichmentStudy
    from goatools.obo_parser import GODag

    # @step: Loading GO ontology and annotations
    with contextlib.redirect_stdout(io.StringIO()):
        godag = GODag(str(obo_path))
        id2gos = GafReader(str(gaf_path)).get_id2gos(namespace="CC")

    expected_terms = {t.strip() for t in args.expected_go_terms.split(",") if t.strip()}
    expanded_terms = expand_with_descendants(godag, expected_terms)

    # @step[id=load_genes]: Loading protein/gene list
    gene_by_protein = {}
    with open(args.gene_lookup_file) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_by_protein[row["protein"]] = row["gene"]

    study_accessions = [protein for protein in gene_by_protein if protein in id2gos]

    # @step-if[id=has_annotations,from=load_genes]: Any study proteins have GO annotations?
    if not study_accessions:
        # @step[from=has_annotations:no]: Writing empty result
        with open(output_folder / "go_evidence.tsv", "w", newline="") as f:
            csv.writer(f, delimiter="\t").writerow(["protein", "gene", "go_evidence", "p_fdr_bh", "matched_terms"])
        print("No proteins had GO annotations in the reference data; wrote an empty result.", file=sys.stderr)
        return

    # @step[from=has_annotations:yes]: Running GO enrichment study
    with contextlib.redirect_stdout(io.StringIO()):
        study = GOEnrichmentStudy(
            list(id2gos.keys()), id2gos, godag, propagate_counts=True, alpha=0.05, methods=["fdr_bh"]
        )
        enrichment_results = study.run_study(study_accessions)

    best_qvalue_by_accession = {}
    matched_terms_by_accession = {}
    for result in enrichment_results:
        if result.enrichment != "e" or result.GO not in expanded_terms:
            continue
        for accession in result.study_items:
            if accession not in best_qvalue_by_accession or result.p_fdr_bh < best_qvalue_by_accession[accession]:
                best_qvalue_by_accession[accession] = result.p_fdr_bh
            matched_terms_by_accession.setdefault(accession, []).append(result.GO)

    # @step: Writing GO evidence scores
    rows = []
    for protein, gene in gene_by_protein.items():
        qvalue = best_qvalue_by_accession.get(protein)
        if qvalue is None:
            continue
        value = -math.log10(max(qvalue, 1e-300))
        rows.append({
            "protein": protein,
            "gene": gene,
            "go_evidence": value,
            "p_fdr_bh": qvalue,
            "matched_terms": ";".join(sorted(matched_terms_by_accession[protein])),
        })

    with open(output_folder / "go_evidence.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["protein", "gene", "go_evidence", "p_fdr_bh", "matched_terms"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} GO evidence scores", file=sys.stderr)
    print("GO evidence scoring complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
