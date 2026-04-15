#!/usr/bin/env python3
"""
Find Pisum sativum orthologs for Medicago truncatula genes using OrthoFinder output.

This script maps genes from a genelist (with R108 ortholog annotations) to their
Pisum sativum orthologs by tracing through OrthoFinder orthogroups.

Strategy:
    Genelist A17 gene -> R108 ortholog -> OrthoFinder orthogroup -> Pisum ortholog

Usage:
    python find_orthologs.py --genelist genelist.tsv --orthofinder Results_Jan15/ --output orthologs.tsv

Requirements:
    - OrthoFinder 3.x output directory
    - Genelist TSV with columns: Mtv5r1.7, acronym, description, R108_orthologs_final
    - Pisum FASTA file (for gene ID extraction)
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def load_pisum_headers(fasta_path):
    """Load Pisum FASTA headers to extract gene IDs and symbols."""
    headers = {}
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split(' ', 1)
                protein_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                gene_match = re.search(r'gene:(PSAT_LOCUS\d+)', description)
                gene_symbol_match = re.search(r'gene_symbol:(\S+)', description)
                headers[protein_id] = {
                    'description': description,
                    'gene_id': gene_match.group(1) if gene_match else '',
                    'gene_symbol': gene_symbol_match.group(1) if gene_symbol_match else ''
                }
    return headers


def load_single_copy(orthofinder_dir):
    """Load list of single-copy orthogroups."""
    sc_file = orthofinder_dir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
    with open(sc_file) as f:
        return set(line.strip() for line in f if line.strip())


def load_orthogroups(orthofinder_dir):
    """Load orthogroups and create R108 -> orthogroup mapping."""
    og_file = orthofinder_dir / "Orthogroups" / "Orthogroups.tsv"
    orthogroups = {}
    r108_to_og = {}

    with open(og_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        species_cols = header[1:]

        for row in reader:
            og = row[0]
            orthogroups[og] = {}
            for i, species in enumerate(species_cols):
                proteins = [p.strip() for p in row[i+1].split(',') if p.strip()]
                orthogroups[og][species] = proteins
                if 'R108' in species:
                    for p in proteins:
                        r108_to_og[p] = og

    return orthogroups, species_cols, r108_to_og


def load_r108_pisum_pairwise(orthofinder_dir):
    """Load direct R108 -> Pisum ortholog mappings from pairwise file."""
    r108_to_pisum = defaultdict(list)

    # Find the pairwise file
    ortho_dir = orthofinder_dir / "Orthologues" / "Orthologues_Medicago_truncatula_R108"
    if not ortho_dir.exists():
        return r108_to_pisum

    for f in ortho_dir.glob("*Pisum*.tsv"):
        with open(f, 'r') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                for col in row:
                    if 'R108' in col:
                        r108_proteins = [p.strip() for p in row[col].split(',') if p.strip()]
                    elif 'Pisum' in col:
                        pisum_proteins = [p.strip() for p in row[col].split(',') if p.strip()]
                for r108 in r108_proteins:
                    r108_to_pisum[r108].extend(pisum_proteins)
        break

    return r108_to_pisum


def find_orthologs(genelist_path, orthofinder_dir, pisum_fasta, r108_prefix="medtr.R108.gnmHiC_1.ann1."):
    """Main function to find Pisum orthologs for all genes in genelist."""

    print("Loading data...")
    pisum_headers = load_pisum_headers(pisum_fasta)
    single_copy = load_single_copy(orthofinder_dir)
    orthogroups, species_cols, r108_to_og = load_orthogroups(orthofinder_dir)
    r108_to_pisum = load_r108_pisum_pairwise(orthofinder_dir)

    print(f"  Pisum proteins: {len(pisum_headers)}")
    print(f"  Orthogroups: {len(orthogroups)}")
    print(f"  Single-copy OGs: {len(single_copy)}")
    print(f"  R108 index: {len(r108_to_og)}")

    # Load genelist
    genelist = []
    with open(genelist_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genelist.append(row)

    print(f"  Genelist genes: {len(genelist)}")

    # Process each gene
    results = []
    stats = {'matched': 0, 'not_matched': 0, 'has_pisum': 0, 'no_pisum': 0}

    for gene in genelist:
        gene_id = gene.get('Mtv5r1.7', '')
        acronym = gene.get('acronym', '')
        description = gene.get('description', '')
        r108_ortholog = gene.get('R108_orthologs_final', '').strip()

        if r108_ortholog:
            # Convert genelist R108 ID to OrthoFinder format
            of_r108_id = f"{r108_prefix}{r108_ortholog}"
            og = r108_to_og.get(of_r108_id, '')

            if og:
                stats['matched'] += 1

                # Get Pisum orthologs from orthogroup
                pisum_from_og = orthogroups[og].get('Pisum_sativum_JI2822', [])
                pisum_from_pairwise = r108_to_pisum.get(of_r108_id, [])
                all_pisum = list(set(pisum_from_og + pisum_from_pairwise))

                is_single_copy = og in single_copy
                a17_proteins = orthogroups[og].get('Medicago_truncatula_A17', [])

                if all_pisum:
                    stats['has_pisum'] += 1
                else:
                    stats['no_pisum'] += 1

                # Get Pisum gene details
                pisum_genes = []
                pisum_symbols = []
                for p in all_pisum[:5]:
                    if p in pisum_headers:
                        info = pisum_headers[p]
                        if info['gene_id']:
                            pisum_genes.append(info['gene_id'])
                        if info['gene_symbol']:
                            pisum_symbols.append(info['gene_symbol'])

                results.append({
                    'gene_id': gene_id,
                    'acronym': acronym,
                    'description': description[:60],
                    'r108_ortholog': r108_ortholog,
                    'orthogroup': og,
                    'single_copy': is_single_copy,
                    'a17_proteins': a17_proteins[:3],
                    'pisum_proteins': all_pisum[:5],
                    'pisum_genes': pisum_genes,
                    'pisum_symbols': pisum_symbols,
                    'status': 'FOUND' if all_pisum else 'NO_PISUM'
                })
            else:
                stats['not_matched'] += 1
                results.append({
                    'gene_id': gene_id, 'acronym': acronym, 'description': description[:60],
                    'r108_ortholog': r108_ortholog, 'orthogroup': 'NOT_IN_OG',
                    'single_copy': False, 'a17_proteins': [], 'pisum_proteins': [],
                    'pisum_genes': [], 'pisum_symbols': [], 'status': 'R108_NOT_FOUND'
                })
        else:
            stats['not_matched'] += 1
            results.append({
                'gene_id': gene_id, 'acronym': acronym, 'description': description[:60],
                'r108_ortholog': '', 'orthogroup': 'NO_R108',
                'single_copy': False, 'a17_proteins': [], 'pisum_proteins': [],
                'pisum_genes': [], 'pisum_symbols': [], 'status': 'NO_R108_IN_GENELIST'
            })

    return results, stats


def write_output(results, output_path):
    """Write results to TSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'A17_Gene_ID', 'Acronym', 'Description', 'R108_Ortholog',
            'Orthogroup', 'SingleCopy_1to1to1',
            'Pisum_Protein_IDs', 'Pisum_Gene_IDs', 'Pisum_Gene_Symbols',
            'A17_Protein_IDs', 'Status'
        ])
        for r in results:
            writer.writerow([
                r['gene_id'], r['acronym'], r['description'], r['r108_ortholog'],
                r['orthogroup'], 'Yes' if r['single_copy'] else 'No',
                ';'.join(r['pisum_proteins']), ';'.join(r['pisum_genes']),
                ';'.join(r['pisum_symbols']), ';'.join(r['a17_proteins']), r['status']
            ])


def main():
    parser = argparse.ArgumentParser(description='Find Pisum orthologs for Medicago genes')
    parser.add_argument('--genelist', required=True, help='Path to genelist TSV')
    parser.add_argument('--orthofinder', required=True, help='Path to OrthoFinder Results directory')
    parser.add_argument('--pisum-fasta', required=True, help='Path to Pisum protein FASTA')
    parser.add_argument('--output', required=True, help='Output TSV path')
    parser.add_argument('--r108-prefix', default='medtr.R108.gnmHiC_1.ann1.',
                        help='Prefix for R108 IDs in OrthoFinder (default: medtr.R108.gnmHiC_1.ann1.)')

    args = parser.parse_args()

    orthofinder_dir = Path(args.orthofinder)

    results, stats = find_orthologs(
        args.genelist,
        orthofinder_dir,
        args.pisum_fasta,
        args.r108_prefix
    )

    write_output(results, args.output)

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total genes: {len(results)}")
    print(f"Matched via R108: {stats['matched']}")
    print(f"  - With Pisum ortholog: {stats['has_pisum']}")
    print(f"  - NO Pisum ortholog: {stats['no_pisum']}")
    print(f"Not matched: {stats['not_matched']}")
    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
