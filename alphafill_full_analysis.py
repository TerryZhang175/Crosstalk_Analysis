#!/usr/bin/env python3
"""
AlphaFill-based ATP/ADP binding site analysis for all proteins
Calculates distances to glutathionylation sites
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
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def download_alphafill_structure(uniprot_id: str, output_dir: str) -> Optional[str]:
    """
    Download AlphaFill-enriched structure from the AlphaFill database
    """
    # Check if already downloaded
    output_file = os.path.join(output_dir, f"AlphaFill_{uniprot_id}.cif")
    if os.path.exists(output_file):
        logger.debug(f"AlphaFill structure already exists for {uniprot_id}")
        return output_file
    
    # AlphaFill API endpoint
    base_url = "https://alphafill.eu/v1"
    alphafill_url = f"{base_url}/aff/{uniprot_id}"
    
    try:
        response = requests.get(alphafill_url, timeout=30)
        if response.status_code == 200:
            # Check if the response contains structure data
            if len(response.content) > 100:  # Basic check for valid content
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                logger.info(f"✓ Downloaded AlphaFill structure for {uniprot_id}")
                return output_file
        else:
            logger.debug(f"No AlphaFill data available for {uniprot_id} (status: {response.status_code})")
    except Exception as e:
        logger.error(f"Error downloading AlphaFill for {uniprot_id}: {e}")
    
    return None

def parse_cif_for_ligands(cif_file: str) -> Dict:
    """
    Parse CIF file to extract ATP/ADP ligand information
    """
    ligands = {
        'ATP': [],
        'ADP': [],
        'AMP': [],
        'GTP': [],
        'GDP': [],
        'other_nucleotides': []
    }
    
    nucleotide_ids = ['ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'CTP', 'CDP', 'UTP', 'UDP', 
                      'ANP', 'ACP', 'AGS']  # Include non-hydrolyzable analogs
    
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
                        
                        if comp_id in ligands:
                            ligands[comp_id].append(atom_info)
                        else:
                            ligands['other_nucleotides'].append(atom_info)
    except Exception as e:
        logger.error(f"Error parsing CIF file: {e}")
    
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
        'has_gtp': False,
        'has_gdp': False,
        'has_amp': False,
        'has_other_nucleotides': False,
        'nucleotide_types': [],
        'distances': [],
        'status': 'no_alphafill'
    }
    
    # Download AlphaFill structure
    alphafill_file = download_alphafill_structure(uniprot_id, alphafill_dir)
    
    if not alphafill_file:
        return results
    
    # Parse ligands from AlphaFill CIF
    ligands = parse_cif_for_ligands(alphafill_file)
    
    # Check what ligands are present
    results['has_atp'] = len(ligands.get('ATP', [])) > 0
    results['has_adp'] = len(ligands.get('ADP', [])) > 0
    results['has_gtp'] = len(ligands.get('GTP', [])) > 0
    results['has_gdp'] = len(ligands.get('GDP', [])) > 0
    results['has_amp'] = len(ligands.get('AMP', [])) > 0
    results['has_other_nucleotides'] = len(ligands.get('other_nucleotides', [])) > 0
    
    # List nucleotide types found
    nucleotide_types = []
    for nuc_type in ['ATP', 'ADP', 'GTP', 'GDP', 'AMP']:
        if ligands.get(nuc_type):
            nucleotide_types.append(nuc_type)
    results['nucleotide_types'] = nucleotide_types
    
    # If we have any nucleotides, calculate distances
    if any([results['has_atp'], results['has_adp'], results['has_gtp'], 
            results['has_gdp'], results['has_amp']]):
        
        results['status'] = 'has_nucleotides'
        
        # Load the original AlphaFold structure
        alphafold_pdb = os.path.join(alphafold_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
        
        if not os.path.exists(alphafold_pdb):
            logger.warning(f"AlphaFold PDB not found for {uniprot_id}")
            results['status'] = 'no_alphafold_pdb'
            return results
        
        structure = parse_pdb_structure(alphafold_pdb)
        
        # Get nucleotide heavy atom coordinates
        nucleotide_coords = {}
        for nuc_type in ['ATP', 'ADP', 'GTP', 'GDP', 'AMP']:
            coords = get_ligand_heavy_atoms(ligands, nuc_type)
            if coords:
                nucleotide_coords[nuc_type] = coords
        
        # Calculate distances for each glutathionylation site
        for gluta_site in glutathione_sites:
            sg_coord = get_cysteine_sg_coords(structure, gluta_site)
            
            if sg_coord:
                distance_info = {
                    'gluta_site': gluta_site,
                    'distances_to_nucleotides': {}
                }
                
                # Calculate distance to each nucleotide type
                for nuc_type, coords in nucleotide_coords.items():
                    distance = calculate_minimum_distance(sg_coord, coords)
                    distance_info['distances_to_nucleotides'][nuc_type] = distance
                
                # Find minimum distance across all nucleotides
                if distance_info['distances_to_nucleotides']:
                    min_distance = min(distance_info['distances_to_nucleotides'].values())
                    distance_info['min_distance'] = min_distance
                    distance_info['closest_nucleotide'] = min(
                        distance_info['distances_to_nucleotides'].items(), 
                        key=lambda x: x[1]
                    )[0]
                
                results['distances'].append(distance_info)
            else:
                logger.warning(f"Could not find SG atom for Cys{gluta_site} in {uniprot_id}")
    else:
        results['status'] = 'no_nucleotides'
    
    return results

def process_batch(protein_batch: List[Tuple[str, List[int]]], 
                 alphafold_dir: str, alphafill_dir: str) -> List[Dict]:
    """
    Process a batch of proteins
    """
    results = []
    for uniprot_id, gluta_sites in protein_batch:
        result = analyze_protein(uniprot_id, gluta_sites, alphafold_dir, alphafill_dir)
        results.append(result)
        time.sleep(0.5)  # Be nice to the server
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
    
    # Get list of proteins that have AlphaFold structures
    alphafold_files = os.listdir(alphafold_dir)
    alphafold_proteins = set()
    for file in alphafold_files:
        if file.startswith('AF-') and file.endswith('.pdb'):
            uniprot_id = file.replace('AF-', '').replace('-F1-model_v4.pdb', '')
            alphafold_proteins.add(uniprot_id)
    
    # Filter for proteins we have both glutathionylation data and AlphaFold structures
    proteins_to_analyze = [(uid, sites) for uid, sites in gluta_sites_by_protein.items() 
                           if uid in alphafold_proteins]
    
    logger.info(f"Found {len(proteins_to_analyze)} proteins with both glutathionylation sites and AlphaFold structures")
    
    # Process all proteins
    all_results = []
    batch_size = 10
    
    for i in range(0, len(proteins_to_analyze), batch_size):
        batch = proteins_to_analyze[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(proteins_to_analyze) + batch_size - 1)//batch_size}")
        
        batch_results = process_batch(batch, alphafold_dir, alphafill_dir)
        all_results.extend(batch_results)
        
        # Save intermediate results
        if (i + batch_size) % 50 == 0:
            intermediate_file = f'/SSD/Terry/Project/Gluta_phospho/alphafill_analysis_intermediate_{i+batch_size}.json'
            with open(intermediate_file, 'w') as f:
                json.dump(all_results, f, indent=2)
            logger.info(f"Saved intermediate results ({len(all_results)} proteins)")
    
    # Save final results
    results_file = '/SSD/Terry/Project/Gluta_phospho/alphafill_atp_analysis_complete.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    # Create summary CSV
    summary_data = []
    for result in all_results:
        if result['distances']:
            for dist_info in result['distances']:
                if 'min_distance' in dist_info:
                    summary_data.append({
                        'uniprot_id': result['uniprot_id'],
                        'gluta_site': dist_info['gluta_site'],
                        'closest_nucleotide': dist_info['closest_nucleotide'],
                        'distance_angstrom': dist_info['min_distance'],
                        'has_atp': result['has_atp'],
                        'has_adp': result['has_adp'],
                        'nucleotide_types': ','.join(result['nucleotide_types'])
                    })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv('/SSD/Terry/Project/Gluta_phospho/alphafill_nucleotide_distances.csv', index=False)
        
        # Create close interaction subset (<= 10 Angstrom)
        close_interactions = summary_df[summary_df['distance_angstrom'] <= 10]
        close_interactions.to_csv('/SSD/Terry/Project/Gluta_phospho/alphafill_close_nucleotide_interactions.csv', index=False)
    
    # Print summary statistics
    logger.info(f"\n{'=' * 60}")
    logger.info(f"Analysis Complete!")
    logger.info(f"Total proteins analyzed: {len(all_results)}")
    
    proteins_with_nucleotides = sum(1 for r in all_results if r['status'] == 'has_nucleotides')
    proteins_with_atp = sum(1 for r in all_results if r['has_atp'])
    proteins_with_adp = sum(1 for r in all_results if r['has_adp'])
    
    logger.info(f"Proteins with nucleotides: {proteins_with_nucleotides}")
    logger.info(f"Proteins with ATP: {proteins_with_atp}")
    logger.info(f"Proteins with ADP: {proteins_with_adp}")
    
    if summary_data:
        close_count = len(close_interactions)
        unique_proteins_close = close_interactions['uniprot_id'].nunique()
        logger.info(f"Site pairs with distance <= 10Å: {close_count}")
        logger.info(f"Unique proteins with close interactions: {unique_proteins_close}")

if __name__ == "__main__":
    main()