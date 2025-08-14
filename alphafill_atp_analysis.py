#!/usr/bin/env python3
"""
AlphaFill-based ATP/ADP binding site analysis and distance calculation to glutathionylation sites
"""

import os
import sys
import json
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import PDB
import urllib.request
import time
from typing import Dict, List, Tuple, Optional

def download_alphafill_structure(uniprot_id: str, output_dir: str) -> Optional[str]:
    """
    Download AlphaFill-enriched structure from the AlphaFill database
    """
    # AlphaFill API endpoint
    base_url = "https://alphafill.eu/v1"
    
    # Try to download the AlphaFill structure
    alphafill_url = f"{base_url}/aff/{uniprot_id}"
    output_file = os.path.join(output_dir, f"AlphaFill_{uniprot_id}.cif")
    
    try:
        response = requests.get(alphafill_url, timeout=30)
        if response.status_code == 200:
            # Check if the response contains structure data
            if len(response.content) > 100:  # Basic check for valid content
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                print(f"✓ Downloaded AlphaFill structure for {uniprot_id}")
                return output_file
        else:
            print(f"✗ No AlphaFill data available for {uniprot_id} (status: {response.status_code})")
    except Exception as e:
        print(f"✗ Error downloading AlphaFill for {uniprot_id}: {e}")
    
    return None

def parse_cif_for_ligands(cif_file: str) -> Dict:
    """
    Parse CIF file to extract ATP/ADP ligand information
    """
    ligands = {
        'ATP': [],
        'ADP': [],
        'other_nucleotides': []
    }
    
    nucleotide_ids = ['ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'CTP', 'CDP', 'UTP', 'UDP']
    
    try:
        with open(cif_file, 'r') as f:
            lines = f.readlines()
        
        in_atom_site = False
        for line in lines:
            if line.startswith('_atom_site.'):
                in_atom_site = True
            elif in_atom_site and line.startswith('HETATM'):
                parts = line.split()
                if len(parts) > 5:
                    comp_id = parts[5]  # compound ID
                    if comp_id in nucleotide_ids:
                        atom_info = {
                            'atom': parts[3],
                            'x': float(parts[10]),
                            'y': float(parts[11]),
                            'z': float(parts[12]),
                            'comp_id': comp_id
                        }
                        
                        if comp_id == 'ATP':
                            ligands['ATP'].append(atom_info)
                        elif comp_id == 'ADP':
                            ligands['ADP'].append(atom_info)
                        else:
                            ligands['other_nucleotides'].append(atom_info)
    except Exception as e:
        print(f"Error parsing CIF file: {e}")
    
    return ligands

def get_ligand_heavy_atoms(ligands: Dict, ligand_type: str) -> List[Tuple[float, float, float]]:
    """
    Extract heavy atom coordinates from ligand data
    """
    coords = []
    if ligand_type in ligands and ligands[ligand_type]:
        for atom in ligands[ligand_type]:
            # Skip hydrogen atoms
            if not atom['atom'].startswith('H'):
                coords.append((atom['x'], atom['y'], atom['z']))
    return coords

def parse_pdb_structure(pdb_file: str):
    """
    Parse PDB structure file
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    return structure

def get_cysteine_sg_coords(structure, cys_position: int) -> Optional[Tuple[float, float, float]]:
    """
    Get the coordinates of the sulfur atom (SG) of a cysteine at given position
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == cys_position:
                    if residue.resname == 'CYS':
                        if 'SG' in residue:
                            atom = residue['SG']
                            return tuple(atom.coord)
    return None

def calculate_minimum_distance(point: Tuple[float, float, float], 
                              coords_list: List[Tuple[float, float, float]]) -> float:
    """
    Calculate minimum distance from a point to a list of coordinates
    """
    if not coords_list:
        return float('inf')
    
    min_dist = float('inf')
    for coord in coords_list:
        dist = np.sqrt(sum((a - b)**2 for a, b in zip(point, coord)))
        if dist < min_dist:
            min_dist = dist
    
    return min_dist

def analyze_protein(uniprot_id: str, glutathione_sites: List[int], 
                   alphafold_dir: str, alphafill_dir: str) -> Dict:
    """
    Analyze a single protein for ATP/ADP binding and distances to glutathionylation sites
    """
    results = {
        'uniprot_id': uniprot_id,
        'has_atp': False,
        'has_adp': False,
        'has_other_nucleotides': False,
        'distances': []
    }
    
    # Download AlphaFill structure
    alphafill_file = download_alphafill_structure(uniprot_id, alphafill_dir)
    
    if not alphafill_file:
        return results
    
    # Parse ligands from AlphaFill CIF
    ligands = parse_cif_for_ligands(alphafill_file)
    
    # Check what ligands are present
    results['has_atp'] = len(ligands['ATP']) > 0
    results['has_adp'] = len(ligands['ADP']) > 0
    results['has_other_nucleotides'] = len(ligands['other_nucleotides']) > 0
    
    # If we have ATP or ADP, calculate distances
    if results['has_atp'] or results['has_adp']:
        # Load the original AlphaFold structure
        alphafold_pdb = os.path.join(alphafold_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
        
        if not os.path.exists(alphafold_pdb):
            print(f"Warning: AlphaFold PDB not found for {uniprot_id}")
            return results
        
        structure = parse_pdb_structure(alphafold_pdb)
        
        # Get ATP/ADP heavy atom coordinates
        atp_coords = get_ligand_heavy_atoms(ligands, 'ATP')
        adp_coords = get_ligand_heavy_atoms(ligands, 'ADP')
        
        # Calculate distances for each glutathionylation site
        for gluta_site in glutathione_sites:
            sg_coord = get_cysteine_sg_coords(structure, gluta_site)
            
            if sg_coord:
                distance_info = {
                    'gluta_site': gluta_site,
                    'distance_to_atp': None,
                    'distance_to_adp': None
                }
                
                if atp_coords:
                    distance_info['distance_to_atp'] = calculate_minimum_distance(sg_coord, atp_coords)
                
                if adp_coords:
                    distance_info['distance_to_adp'] = calculate_minimum_distance(sg_coord, adp_coords)
                
                results['distances'].append(distance_info)
            else:
                print(f"Warning: Could not find SG atom for Cys{gluta_site} in {uniprot_id}")
    
    return results

def main():
    # Create directories
    alphafill_dir = "/SSD/Terry/Project/Gluta_phospho/alphafill_structures"
    alphafold_dir = "/SSD/Terry/Project/Gluta_phospho/alphafold_structures"
    os.makedirs(alphafill_dir, exist_ok=True)
    
    # Load glutathionylation sites
    gluta_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/cleaned_glutathionylation_sites.csv')
    
    # Group by UniProt ID
    gluta_sites_by_protein = gluta_df.groupby('UniProt_ID')['Glutathione_Site'].apply(list).to_dict()
    
    # Test with proteins that have close phospho-gluta interactions
    # These are more likely to have nucleotide binding sites
    close_interaction_proteins = ['P50991', 'P68371', 'P27348', 'P00558', 'P04406', 
                                  'P07310', 'P04075', 'Q04447', 'P13929', 'P50461']
    
    # Filter for proteins we have glutathionylation data for
    test_proteins = [p for p in close_interaction_proteins if p in gluta_sites_by_protein][:10]
    
    print(f"Testing AlphaFill ATP/ADP analysis on {len(test_proteins)} proteins...")
    print("=" * 60)
    
    all_results = []
    
    for uniprot_id in test_proteins:
        print(f"\nAnalyzing {uniprot_id}...")
        gluta_sites = gluta_sites_by_protein[uniprot_id]
        print(f"  Glutathionylation sites: {gluta_sites}")
        
        result = analyze_protein(uniprot_id, gluta_sites, alphafold_dir, alphafill_dir)
        all_results.append(result)
        
        # Report findings
        if result['has_atp'] or result['has_adp']:
            print(f"  ✓ Found nucleotides: ATP={result['has_atp']}, ADP={result['has_adp']}")
            for dist_info in result['distances']:
                if dist_info['distance_to_atp'] is not None:
                    print(f"    Cys{dist_info['gluta_site']} -> ATP: {dist_info['distance_to_atp']:.2f} Å")
                if dist_info['distance_to_adp'] is not None:
                    print(f"    Cys{dist_info['gluta_site']} -> ADP: {dist_info['distance_to_adp']:.2f} Å")
        else:
            print(f"  ✗ No ATP/ADP found")
        
        time.sleep(1)  # Be nice to the server
    
    # Save results
    results_file = '/SSD/Terry/Project/Gluta_phospho/alphafill_atp_analysis_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print(f"Results saved to {results_file}")
    
    # Summary statistics
    proteins_with_atp = sum(1 for r in all_results if r['has_atp'])
    proteins_with_adp = sum(1 for r in all_results if r['has_adp'])
    proteins_with_nucleotides = sum(1 for r in all_results if r['has_atp'] or r['has_adp'])
    
    print(f"\nSummary:")
    print(f"  Proteins with ATP: {proteins_with_atp}/{len(test_proteins)}")
    print(f"  Proteins with ADP: {proteins_with_adp}/{len(test_proteins)}")
    print(f"  Proteins with any nucleotide: {proteins_with_nucleotides}/{len(test_proteins)}")

if __name__ == "__main__":
    main()