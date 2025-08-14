#!/usr/bin/env python3
"""
Calculate distances between acetylation and glutathionylation sites using AlphaFold structures
Similar to phosphorylation and ubiquitination analysis
"""

import pandas as pd
import numpy as np
import os
from typing import List, Dict, Tuple, Optional
import json
from pathlib import Path

def parse_pdb_structure(pdb_file: str) -> Dict[int, np.ndarray]:
    """
    Parse PDB file and extract CA coordinates for each residue
    """
    coordinates = {}
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    # Extract residue number and coordinates
                    residue_num = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    coordinates[residue_num] = np.array([x, y, z])
    
    except Exception as e:
        print(f"Error parsing PDB file {pdb_file}: {e}")
        return {}
    
    return coordinates

def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    """
    Calculate Euclidean distance between two 3D coordinates
    """
    return np.linalg.norm(coord1 - coord2)

def analyze_protein_distances(uniprot_id: str, acetyl_sites: List[int], gluta_sites: List[int], 
                            alphafold_dir: str) -> List[dict]:
    """
    Calculate all pairwise distances between acetylation and glutathionylation sites for a protein
    """
    results = []
    
    # Try to find AlphaFold structure
    pdb_file = os.path.join(alphafold_dir, f'AF-{uniprot_id}-F1-model_v4.pdb')
    
    if not os.path.exists(pdb_file):
        print(f"Warning: No AlphaFold structure found for {uniprot_id}")
        return results
    
    # Parse structure
    coordinates = parse_pdb_structure(pdb_file)
    
    if not coordinates:
        print(f"Warning: Failed to parse structure for {uniprot_id}")
        return results
    
    # Calculate all pairwise distances
    successful_calculations = 0
    total_pairs = len(acetyl_sites) * len(gluta_sites)
    
    for acetyl_site in acetyl_sites:
        for gluta_site in gluta_sites:
            # Check if both sites have coordinates
            if acetyl_site in coordinates and gluta_site in coordinates:
                distance = calculate_distance(coordinates[acetyl_site], coordinates[gluta_site])
                
                results.append({
                    'uniprot_id': uniprot_id,
                    'acetyl_site': acetyl_site,
                    'gluta_site': gluta_site,
                    'distance_angstrom': distance
                })
                successful_calculations += 1
            else:
                missing_sites = []
                if acetyl_site not in coordinates:
                    missing_sites.append(f"acetyl_{acetyl_site}")
                if gluta_site not in coordinates:
                    missing_sites.append(f"gluta_{gluta_site}")
                
                print(f"Warning: Missing coordinates for {uniprot_id}: {', '.join(missing_sites)}")
    
    if successful_calculations > 0:
        print(f"Success: {uniprot_id} - {successful_calculations}/{total_pairs} distance calculations")
    
    return results

def get_residue_type(coordinates: Dict[int, np.ndarray], pdb_file: str, site: int) -> Optional[str]:
    """
    Get the residue type (amino acid) at a specific site from PDB file
    """
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    residue_num = int(line[22:26].strip())
                    if residue_num == site:
                        residue_type = line[17:20].strip()
                        return residue_type
    except:
        pass
    
    return None

def main():
    print("=== Calculate Acetylation-Glutathionylation Site Distances ===\n")
    
    base_dir = '/SSD/Terry/Project/Gluta_phospho'
    alphafold_dir = f'{base_dir}/alphafold_structures'
    
    # Load acetylation-glutathionylation overlap data
    acetyl_file = f'{base_dir}/acetylation_glutathionylation_overlap.csv'
    acetyl_df = pd.read_csv(acetyl_file)
    
    # Load glutathionylation data
    gluta_file = f'{base_dir}/cleaned_glutathionylation_sites.csv'
    gluta_df = pd.read_csv(gluta_file)
    
    print(f"Loaded acetylation sites: {len(acetyl_df):,}")
    print(f"Loaded glutathionylation sites: {len(gluta_df):,}")
    
    # Get proteins with both modifications
    overlap_proteins = set(acetyl_df['UniProt_ID']).intersection(set(gluta_df['UniProt_ID']))
    print(f"Proteins with both modifications: {len(overlap_proteins):,}")
    
    # Prepare data for distance calculations
    all_results = []
    processed_proteins = 0
    failed_proteins = 0
    
    print(f"\n🔄 Starting distance calculation...")
    
    for uniprot_id in overlap_proteins:
        # Get acetylation sites for this protein
        protein_acetyl = acetyl_df[acetyl_df['UniProt_ID'] == uniprot_id]
        acetyl_sites = protein_acetyl['Acetylation_Site'].tolist()
        
        # Get glutathionylation sites for this protein
        protein_gluta = gluta_df[gluta_df['UniProt_ID'] == uniprot_id]
        gluta_sites = protein_gluta['Glutathione_Site'].tolist()
        
        # Calculate distances
        protein_results = analyze_protein_distances(uniprot_id, acetyl_sites, gluta_sites, alphafold_dir)
        
        if protein_results:
            all_results.extend(protein_results)
            processed_proteins += 1
        else:
            failed_proteins += 1
        
        # Progress update
        if (processed_proteins + failed_proteins) % 50 == 0:
            print(f"Progress: {processed_proteins + failed_proteins}/{len(overlap_proteins)} proteins processed")
    
    print(f"\n📊 Calculation complete:")
    print(f"  Successfully processed proteins: {processed_proteins:,}")
    print(f"  Failed proteins: {failed_proteins:,}")
    print(f"  Total distance calculations: {len(all_results):,}")
    
    if all_results:
        # Convert to DataFrame
        results_df = pd.DataFrame(all_results)
        
        # Save raw results
        raw_output_file = f'{base_dir}/acetyl_gluta_distances_raw.csv'
        results_df.to_csv(raw_output_file, index=False)
        print(f"\n💾 Raw results saved to: {raw_output_file}")
        
        # Filter by residue types (Lys for acetylation, Cys for glutathionylation)
        print(f"\n🧬 Filtering by residue type...")
        
        filtered_results = []
        residue_validation_count = 0
        
        for _, row in results_df.iterrows():
            uniprot_id = row['uniprot_id']
            acetyl_site = int(row['acetyl_site'])
            gluta_site = int(row['gluta_site'])
            
            # Get PDB file path
            pdb_file = os.path.join(alphafold_dir, f'AF-{uniprot_id}-F1-model_v4.pdb')
            
            if os.path.exists(pdb_file):
                # Parse coordinates to get residue types
                coordinates = parse_pdb_structure(pdb_file)
                acetyl_residue = get_residue_type(coordinates, pdb_file, acetyl_site)
                gluta_residue = get_residue_type(coordinates, pdb_file, gluta_site)
                
                # Acetylation typically occurs on Lysine (LYS), Glutathionylation on Cysteine (CYS)
                if acetyl_residue == 'LYS' and gluta_residue == 'CYS':
                    filtered_results.append({
                        'uniprot_id': uniprot_id,
                        'acetyl_site': acetyl_site,
                        'acetyl_residue': acetyl_residue,
                        'gluta_site': gluta_site,
                        'gluta_residue': gluta_residue,
                        'distance_angstrom': row['distance_angstrom']
                    })
                    residue_validation_count += 1
        
        print(f"  Valid acetylation(Lys)-glutathionylation(Cys) pairs: {len(filtered_results):,}")
        
        if filtered_results:
            # Save filtered results
            filtered_df = pd.DataFrame(filtered_results)
            filtered_output_file = f'{base_dir}/acetyl_gluta_distances_filtered.csv'
            filtered_df.to_csv(filtered_output_file, index=False)
            print(f"  Filtered results saved to: {filtered_output_file}")
            
            # Statistics
            print(f"\n📈 Distance statistics:")
            distances = filtered_df['distance_angstrom']
            print(f"  Minimum distance: {distances.min():.3f} Å")
            print(f"  Maximum distance: {distances.max():.3f} Å")
            print(f"  Average distance: {distances.mean():.3f} Å")
            print(f"  Median distance: {distances.median():.3f} Å")
            
            # Count close interactions (≤10Å)
            close_interactions = filtered_df[filtered_df['distance_angstrom'] <= 10.0]
            print(f"\n🎯 Close interactions (≤10Å):")
            print(f"  Number of close interactions: {len(close_interactions):,}")
            print(f"  Percentage of close interactions: {len(close_interactions)/len(filtered_df)*100:.1f}%")
            print(f"  Number of proteins involved: {close_interactions['uniprot_id'].nunique():,}")
            
            if len(close_interactions) > 0:
                # Save close interactions
                close_output_file = f'{base_dir}/acetyl_gluta_close_interactions.csv'
                close_interactions.to_csv(close_output_file, index=False)
                print(f"  Close interactions saved to: {close_output_file}")
        
        else:
            print(f"⚠️ No valid acetylation(Lys)-glutathionylation(Cys) pairs found")
    
    else:
        print(f"❌ No distances were successfully calculated")

if __name__ == "__main__":
    main()