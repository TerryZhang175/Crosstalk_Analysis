#!/usr/bin/env python3
"""
Calculate minimum Euclidean distances between phosphorylation and glutathionylation sites
in the same proteins using AlphaFold PDB structures.

For each protein with both modification types, find all site pairs and calculate
the minimum distance considering all atoms in both residues.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from collections import defaultdict
import re

# Try to import BioPython for PDB parsing
try:
    from Bio.PDB import PDBParser, Selection
    from Bio.PDB.NeighborSearch import NeighborSearch
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("BioPython not available. Installing...")
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    from Bio.PDB import PDBParser, Selection
    from Bio.PDB.NeighborSearch import NeighborSearch
    BIOPYTHON_AVAILABLE = True

def load_modification_sites():
    """
    Load phosphorylation and glutathionylation sites.
    
    Returns:
        dict: {uniprot_id: {'phospho': [sites], 'gluta': [sites]}}
    """
    # Load phosphorylation sites
    phospho_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/phosphorylation_sites.csv')
    
    # Load glutathionylation sites
    gluta_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/cleaned_glutathionylation_sites.csv')
    
    # Group by UniProt ID
    protein_sites = defaultdict(lambda: {'phospho': [], 'gluta': []})
    
    # Add phosphorylation sites
    for _, row in phospho_df.iterrows():
        uniprot_id = row['UniProt_ID']
        site = int(row['Phosphorylation_Site'])
        protein_sites[uniprot_id]['phospho'].append(site)
    
    # Add glutathionylation sites
    for _, row in gluta_df.iterrows():
        uniprot_id = row['UniProt_ID']
        site = int(row['Glutathione_Site'])
        protein_sites[uniprot_id]['gluta'].append(site)
    
    # Keep only proteins with both modification types
    overlapping_proteins = {}
    for uniprot_id, sites in protein_sites.items():
        if sites['phospho'] and sites['gluta']:
            # Remove duplicates and sort
            sites['phospho'] = sorted(list(set(sites['phospho'])))
            sites['gluta'] = sorted(list(set(sites['gluta'])))
            overlapping_proteins[uniprot_id] = sites
    
    return overlapping_proteins

def parse_pdb_structure(pdb_file):
    """
    Parse PDB structure and extract residue information.
    
    Args:
        pdb_file: Path to PDB file
    
    Returns:
        dict: {residue_number: {'atoms': [(x, y, z), ...], 'residue_name': str}}
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('protein', pdb_file)
        
        residues_data = {}
        
        # Get the first model and chain (AlphaFold structures have single model/chain)
        model = structure[0]
        chain = list(model.get_chains())[0]  # Get first chain
        
        for residue in chain.get_residues():
            # Skip water and heteroatoms
            if residue.get_id()[0] != ' ':
                continue
                
            res_num = residue.get_id()[1]  # Residue number
            res_name = residue.get_resname()
            
            # Get all atom coordinates
            atoms = []
            for atom in residue.get_atoms():
                coord = atom.get_coord()
                atoms.append((float(coord[0]), float(coord[1]), float(coord[2])))
            
            if atoms:  # Only add if residue has atoms
                residues_data[res_num] = {
                    'atoms': atoms,
                    'residue_name': res_name
                }
        
        return residues_data
    
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return {}

def calculate_minimum_distance(residue1_atoms, residue2_atoms):
    """
    Calculate minimum Euclidean distance between all atom pairs in two residues.
    
    Args:
        residue1_atoms: List of (x, y, z) coordinates for residue 1
        residue2_atoms: List of (x, y, z) coordinates for residue 2
    
    Returns:
        float: Minimum distance between any two atoms
    """
    min_distance = float('inf')
    
    for atom1 in residue1_atoms:
        for atom2 in residue2_atoms:
            # Calculate Euclidean distance
            distance = np.sqrt(
                (atom1[0] - atom2[0])**2 + 
                (atom1[1] - atom2[1])**2 + 
                (atom1[2] - atom2[2])**2
            )
            if distance < min_distance:
                min_distance = distance
    
    return min_distance

def process_protein_distances(uniprot_id, sites_data, pdb_dir):
    """
    Process a single protein to calculate all pairwise distances.
    
    Args:
        uniprot_id: UniProt accession ID
        sites_data: Dict with 'phospho' and 'gluta' site lists
        pdb_dir: Directory containing PDB files
    
    Returns:
        list: List of distance results for this protein
    """
    # Find PDB file (try both v4 and v3)
    pdb_file_v4 = pdb_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
    pdb_file_v3 = pdb_dir / f"AF-{uniprot_id}-F1-model_v3.pdb"
    
    pdb_file = None
    if pdb_file_v4.exists():
        pdb_file = pdb_file_v4
    elif pdb_file_v3.exists():
        pdb_file = pdb_file_v3
    
    if not pdb_file:
        print(f"No PDB file found for {uniprot_id}")
        return []
    
    # Parse structure
    residues_data = parse_pdb_structure(pdb_file)
    if not residues_data:
        print(f"Failed to parse structure for {uniprot_id}")
        return []
    
    results = []
    
    # Calculate distances between all phospho-gluta pairs
    for phospho_site in sites_data['phospho']:
        for gluta_site in sites_data['gluta']:
            
            # Skip if same site (shouldn't happen, but just in case)
            if phospho_site == gluta_site:
                continue
            
            # Check if both sites exist in structure
            if phospho_site not in residues_data or gluta_site not in residues_data:
                # Record missing residues
                missing = []
                if phospho_site not in residues_data:
                    missing.append(f"phospho_{phospho_site}")
                if gluta_site not in residues_data:
                    missing.append(f"gluta_{gluta_site}")
                
                results.append({
                    'uniprot_id': uniprot_id,
                    'phospho_site': phospho_site,
                    'gluta_site': gluta_site,
                    'distance_angstrom': None,
                    'phospho_residue': None,
                    'gluta_residue': None,
                    'status': f"Missing residues: {', '.join(missing)}"
                })
                continue
            
            # Get residue data
            phospho_residue = residues_data[phospho_site]
            gluta_residue = residues_data[gluta_site]
            
            # Calculate minimum distance
            min_distance = calculate_minimum_distance(
                phospho_residue['atoms'],
                gluta_residue['atoms']
            )
            
            results.append({
                'uniprot_id': uniprot_id,
                'phospho_site': phospho_site,
                'gluta_site': gluta_site,
                'distance_angstrom': round(min_distance, 3),
                'phospho_residue': phospho_residue['residue_name'],
                'gluta_residue': gluta_residue['residue_name'],
                'status': 'calculated'
            })
    
    return results

def main():
    # Paths
    pdb_dir = Path('/SSD/Terry/Project/Gluta_phospho/alphafold_structures')
    output_dir = Path('/SSD/Terry/Project/Gluta_phospho/distance_analysis')
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("Phosphorylation-Glutathionylation Site Distance Analysis")
    print("=" * 70)
    
    # Load modification sites
    print("Loading modification sites...")
    protein_sites = load_modification_sites()
    
    print(f"Found {len(protein_sites)} proteins with both phosphorylation and glutathionylation sites")
    
    # Show summary statistics
    total_phospho = sum(len(sites['phospho']) for sites in protein_sites.values())
    total_gluta = sum(len(sites['gluta']) for sites in protein_sites.values())
    total_pairs = sum(len(sites['phospho']) * len(sites['gluta']) for sites in protein_sites.values())
    
    print(f"Total phosphorylation sites: {total_phospho}")
    print(f"Total glutathionylation sites: {total_gluta}")
    print(f"Total site pairs to analyze: {total_pairs}")
    
    # Process each protein
    print(f"\nProcessing proteins...")
    print("-" * 70)
    
    all_results = []
    proteins_processed = 0
    proteins_with_pdb = 0
    
    for i, (uniprot_id, sites_data) in enumerate(protein_sites.items(), 1):
        print(f"[{i}/{len(protein_sites)}] Processing {uniprot_id} "
              f"({len(sites_data['phospho'])} phospho, {len(sites_data['gluta'])} gluta sites)...")
        
        results = process_protein_distances(uniprot_id, sites_data, pdb_dir)
        
        if results:
            all_results.extend(results)
            proteins_with_pdb += 1
            
            # Show sample results for first few proteins
            if i <= 3:
                successful_results = [r for r in results if r['status'] == 'calculated']
                if successful_results:
                    print(f"  Sample distances: {[r['distance_angstrom'] for r in successful_results[:3]]}")
        
        proteins_processed += 1
        
        # Save intermediate results every 50 proteins
        if i % 50 == 0:
            df_temp = pd.DataFrame(all_results)
            df_temp.to_csv(output_dir / f"distances_intermediate_{i}.csv", index=False)
    
    # Create final results DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Save detailed results
    output_file = output_dir / "phospho_gluta_distances.csv"
    results_df.to_csv(output_file, index=False)
    
    # Generate summary statistics
    successful_results = results_df[results_df['status'] == 'calculated'].copy()
    
    summary_stats = {
        'total_proteins': len(protein_sites),
        'proteins_with_pdb': proteins_with_pdb,
        'total_site_pairs': len(results_df),
        'successful_calculations': len(successful_results),
        'failed_calculations': len(results_df) - len(successful_results),
        'distance_statistics': {}
    }
    
    if len(successful_results) > 0:
        distances = successful_results['distance_angstrom'].values
        summary_stats['distance_statistics'] = {
            'mean_distance': float(np.mean(distances)),
            'median_distance': float(np.median(distances)),
            'std_distance': float(np.std(distances)),
            'min_distance': float(np.min(distances)),
            'max_distance': float(np.max(distances)),
            'percentiles': {
                '25th': float(np.percentile(distances, 25)),
                '75th': float(np.percentile(distances, 75)),
                '90th': float(np.percentile(distances, 90)),
                '95th': float(np.percentile(distances, 95))
            }
        }
        
        # Distance distribution
        distance_ranges = [
            (0, 5, "Very close (0-5 Å)"),
            (5, 10, "Close (5-10 Å)"),
            (10, 20, "Medium (10-20 Å)"),
            (20, 50, "Far (20-50 Å)"),
            (50, float('inf'), "Very far (>50 Å)")
        ]
        
        distribution = {}
        for min_d, max_d, label in distance_ranges:
            count = len(distances[(distances >= min_d) & (distances < max_d)])
            distribution[label] = {
                'count': int(count),
                'percentage': float(count / len(distances) * 100)
            }
        
        summary_stats['distance_distribution'] = distribution
    
    # Save summary
    summary_file = output_dir / "distance_analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    # Print final summary
    print("\n" + "=" * 70)
    print("Distance Analysis Complete")
    print("=" * 70)
    print(f"Results saved to: {output_file}")
    print(f"Summary saved to: {summary_file}")
    print(f"\nProcessing Summary:")
    print(f"  Total proteins: {summary_stats['total_proteins']}")
    print(f"  Proteins with PDB: {summary_stats['proteins_with_pdb']}")
    print(f"  Total site pairs: {summary_stats['total_site_pairs']}")
    print(f"  Successful calculations: {summary_stats['successful_calculations']}")
    print(f"  Failed calculations: {summary_stats['failed_calculations']}")
    
    if summary_stats['distance_statistics']:
        stats = summary_stats['distance_statistics']
        print(f"\nDistance Statistics:")
        print(f"  Mean: {stats['mean_distance']:.2f} Å")
        print(f"  Median: {stats['median_distance']:.2f} Å")
        print(f"  Range: {stats['min_distance']:.2f} - {stats['max_distance']:.2f} Å")
        
        print(f"\nDistance Distribution:")
        for label, data in summary_stats['distance_distribution'].items():
            print(f"  {label}: {data['count']} pairs ({data['percentage']:.1f}%)")

if __name__ == "__main__":
    main()