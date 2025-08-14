#!/usr/bin/env python3
"""
Calculate distances between ubiquitination and glutathionylation sites
Following the same approach as phospho-gluta analysis
"""

import pandas as pd
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser
import os
import sys
from typing import List, Tuple, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_pdb_structure(pdb_file: str):
    """Parse PDB structure using BioPython"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        return structure
    except Exception as e:
        logger.error(f"Error parsing {pdb_file}: {e}")
        return None

def get_residue_atoms(structure, position: int, expected_residue: str = None) -> List[Tuple[float, float, float]]:
    """
    Get all heavy atom coordinates for a residue at given position
    """
    coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == position:
                    # Check residue type if specified
                    if expected_residue and residue.resname != expected_residue:
                        logger.warning(f"Expected {expected_residue} at position {position}, found {residue.resname}")
                        continue
                    
                    # Get all heavy atoms (non-hydrogen)
                    for atom in residue:
                        if not atom.name.startswith('H'):
                            coords.append(tuple(atom.coord))
                    break
    
    return coords

def calculate_minimum_distance(atoms1: List[Tuple[float, float, float]], 
                              atoms2: List[Tuple[float, float, float]]) -> float:
    """
    Calculate minimum distance between two sets of atoms
    """
    if not atoms1 or not atoms2:
        return float('inf')
    
    min_distance = float('inf')
    
    for atom1 in atoms1:
        for atom2 in atoms2:
            distance = np.sqrt(sum((a - b)**2 for a, b in zip(atom1, atom2)))
            if distance < min_distance:
                min_distance = distance
    
    return min_distance

def analyze_protein_distances(uniprot_id: str, ubiq_sites: List[int], gluta_sites: List[int], 
                            alphafold_dir: str) -> List[dict]:
    """
    Calculate all pairwise distances between ubiquitination and glutathionylation sites
    """
    results = []
    
    # Load AlphaFold structure
    pdb_file = os.path.join(alphafold_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
    
    if not os.path.exists(pdb_file):
        logger.warning(f"PDB file not found for {uniprot_id}")
        return results
    
    structure = parse_pdb_structure(pdb_file)
    if structure is None:
        return results
    
    # Calculate all pairwise distances
    for ubiq_site in ubiq_sites:
        for gluta_site in gluta_sites:
            
            # Get atom coordinates (we'll validate residue types later)
            ubiq_atoms = get_residue_atoms(structure, ubiq_site)
            gluta_atoms = get_residue_atoms(structure, gluta_site)
            
            if ubiq_atoms and gluta_atoms:
                distance = calculate_minimum_distance(ubiq_atoms, gluta_atoms)
                
                results.append({
                    'uniprot_id': uniprot_id,
                    'ubiq_site': ubiq_site,
                    'gluta_site': gluta_site,
                    'distance_angstrom': distance,
                    'status': 'calculated'
                })
            else:
                logger.warning(f"Could not find atoms for {uniprot_id} ubiq:{ubiq_site} gluta:{gluta_site}")
                results.append({
                    'uniprot_id': uniprot_id,
                    'ubiq_site': ubiq_site,
                    'gluta_site': gluta_site,
                    'distance_angstrom': None,
                    'status': 'missing_atoms'
                })
    
    return results

def main():
    # Set up directories
    alphafold_dir = '/SSD/Terry/Project/Gluta_phospho/alphafold_structures'
    
    # Load datasets
    logger.info("Loading datasets...")
    
    # Load proteins ready for analysis
    ready_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/ubiq_gluta_proteins_with_af.csv')
    ready_proteins = set(ready_df['uniprot_id'].unique())
    
    # Load site data
    gluta_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/glutathionylation_sites_overlap.csv')
    ubiq_df = pd.read_csv('/SSD/Terry/Project/Gluta_phospho/ubiquitination_sites_overlap.csv')
    
    # Filter for proteins with AlphaFold structures
    gluta_df = gluta_df[gluta_df['uniprot_id'].isin(ready_proteins)]
    ubiq_df = ubiq_df[ubiq_df['uniprot_id'].isin(ready_proteins)]
    
    logger.info(f"Analyzing {len(ready_proteins)} proteins with AlphaFold structures")
    logger.info(f"Glutathionylation sites: {len(gluta_df)}")
    logger.info(f"Ubiquitination sites: {len(ubiq_df)}")
    
    # Group sites by protein
    gluta_sites_by_protein = gluta_df.groupby('uniprot_id')['gluta_site'].apply(list).to_dict()
    ubiq_sites_by_protein = ubiq_df.groupby('uniprot_id')['ubiq_site'].apply(list).to_dict()
    
    # Calculate distances
    all_results = []
    total_proteins = len(ready_proteins)
    
    for i, protein_id in enumerate(sorted(ready_proteins)):
        logger.info(f"Processing {i+1}/{total_proteins}: {protein_id}")
        
        ubiq_sites = ubiq_sites_by_protein.get(protein_id, [])
        gluta_sites = gluta_sites_by_protein.get(protein_id, [])
        
        if ubiq_sites and gluta_sites:
            logger.info(f"  Ubiq sites: {len(ubiq_sites)}, Gluta sites: {len(gluta_sites)} -> {len(ubiq_sites) * len(gluta_sites)} pairs")
            
            protein_results = analyze_protein_distances(protein_id, ubiq_sites, gluta_sites, alphafold_dir)
            all_results.extend(protein_results)
        else:
            logger.warning(f"  No sites found for {protein_id}")
    
    # Save results
    if all_results:
        results_df = pd.DataFrame(all_results)
        
        output_file = '/SSD/Terry/Project/Gluta_phospho/ubiq_gluta_distances_raw.csv'
        results_df.to_csv(output_file, index=False)
        
        # Generate summary
        successful_calcs = results_df[results_df['status'] == 'calculated']
        logger.info(f"\nDistance calculation completed!")
        logger.info(f"Total site pairs analyzed: {len(results_df)}")
        logger.info(f"Successful calculations: {len(successful_calcs)}")
        logger.info(f"Failed calculations: {len(results_df) - len(successful_calcs)}")
        logger.info(f"Results saved to: {output_file}")
        
        if len(successful_calcs) > 0:
            logger.info(f"Distance range: {successful_calcs['distance_angstrom'].min():.2f} - {successful_calcs['distance_angstrom'].max():.2f} Å")
            logger.info(f"Mean distance: {successful_calcs['distance_angstrom'].mean():.2f} Å")
            
            # Count close interactions
            close_interactions = successful_calcs[successful_calcs['distance_angstrom'] <= 10.0]
            logger.info(f"Close interactions (≤10Å): {len(close_interactions)}")
    
    else:
        logger.error("No results generated!")

if __name__ == "__main__":
    main()