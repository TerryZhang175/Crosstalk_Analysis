#!/usr/bin/env python3
"""
Analyze multi-PTM crosstalk mechanisms around the same glutathionylation sites
Find Cys sites that have multiple PTM partners nearby (acetylation, phosphorylation, ubiquitination)
"""

import pandas as pd
import numpy as np

def load_all_crosstalk_data():
    """
    Load all PTM crosstalk data
    """
    base_dir = '/SSD/Terry/Project/Gluta_phospho'
    
    # Load phosphorylation-glutathionylation data
    phospho_file = f'{base_dir}/Final_Results/02_glutathionylation_phosphorylation_interactions.csv'
    phospho_df = pd.read_csv(phospho_file)
    
    # Load ubiquitination-glutathionylation data  
    ubiq_file = f'{base_dir}/Final_Results/01_ubiquitination_glutathionylation_interactions.csv'
    ubiq_df = pd.read_csv(ubiq_file)
    
    # Load acetylation-glutathionylation data
    acetyl_file = f'{base_dir}/Final_Results/03_acetylation_glutathionylation_interactions.csv'
    acetyl_df = pd.read_csv(acetyl_file)
    
    return phospho_df, ubiq_df, acetyl_df

def standardize_columns(phospho_df, ubiq_df, acetyl_df):
    """
    Standardize column names across all datasets
    """
    # Phosphorylation data
    phospho_clean = phospho_df.rename(columns={
        'uniprot_id': 'UniProt_ID',
        'protein_name': 'Protein_Name',
        'phospho_site': 'PTM_Site',
        'gluta_site': 'Gluta_Site',
        'distance_angstrom': 'Distance_Angstrom'
    })
    phospho_clean['PTM_Type'] = 'Phosphorylation'
    
    # Ubiquitination data
    ubiq_clean = ubiq_df.rename(columns={
        'Ubiquitination_Site': 'PTM_Site',
        'Glutathionylation_Site': 'Gluta_Site'
    })
    ubiq_clean['PTM_Type'] = 'Ubiquitination'
    
    # Acetylation data
    acetyl_clean = acetyl_df.rename(columns={
        'Acetylation_Site': 'PTM_Site',
        'Glutathionylation_Site': 'Gluta_Site'
    })
    acetyl_clean['PTM_Type'] = 'Acetylation'
    
    # Select common columns
    common_cols = ['UniProt_ID', 'Protein_Name', 'PTM_Site', 'Gluta_Site', 'Distance_Angstrom', 'PTM_Type']
    
    phospho_clean = phospho_clean[common_cols]
    ubiq_clean = ubiq_clean[common_cols]
    acetyl_clean = acetyl_clean[common_cols]
    
    return phospho_clean, ubiq_clean, acetyl_clean

def find_multi_ptm_sites(phospho_df, ubiq_df, acetyl_df):
    """
    Find glutathionylation sites that have multiple PTM partners nearby
    """
    # Combine all PTM data
    all_ptm_df = pd.concat([phospho_df, ubiq_df, acetyl_df], ignore_index=True)
    
    # Create unique site identifier
    all_ptm_df['Protein_Gluta_Site'] = all_ptm_df['UniProt_ID'] + '_C' + all_ptm_df['Gluta_Site'].astype(str)
    
    # Group by protein and glutathionylation site
    multi_ptm_analysis = []
    
    for protein_gluta_site in all_ptm_df['Protein_Gluta_Site'].unique():
        site_data = all_ptm_df[all_ptm_df['Protein_Gluta_Site'] == protein_gluta_site]
        
        uniprot_id = site_data.iloc[0]['UniProt_ID']
        protein_name = site_data.iloc[0]['Protein_Name']
        gluta_site = site_data.iloc[0]['Gluta_Site']
        
        # Count different PTM types
        ptm_types = site_data['PTM_Type'].unique()
        ptm_count = len(ptm_types)
        
        if ptm_count >= 2:  # Multi-PTM sites
            # Get statistics for each PTM type
            ptm_details = {}
            total_interactions = 0
            min_distance = float('inf')
            
            for ptm_type in ptm_types:
                ptm_subset = site_data[site_data['PTM_Type'] == ptm_type]
                ptm_details[ptm_type] = {
                    'count': len(ptm_subset),
                    'min_distance': ptm_subset['Distance_Angstrom'].min(),
                    'sites': list(ptm_subset['PTM_Site'].unique())
                }
                total_interactions += len(ptm_subset)
                min_distance = min(min_distance, ptm_subset['Distance_Angstrom'].min())
            
            multi_ptm_analysis.append({
                'UniProt_ID': uniprot_id,
                'Protein_Name': protein_name,
                'Gluta_Site': gluta_site,
                'PTM_Types': ', '.join(sorted(ptm_types)),
                'PTM_Count': ptm_count,
                'Total_Interactions': total_interactions,
                'Min_Distance': min_distance,
                'Phospho_Count': ptm_details.get('Phosphorylation', {}).get('count', 0),
                'Phospho_Min_Dist': ptm_details.get('Phosphorylation', {}).get('min_distance', None),
                'Phospho_Sites': ', '.join(map(str, ptm_details.get('Phosphorylation', {}).get('sites', []))),
                'Ubiq_Count': ptm_details.get('Ubiquitination', {}).get('count', 0),
                'Ubiq_Min_Dist': ptm_details.get('Ubiquitination', {}).get('min_distance', None),
                'Ubiq_Sites': ', '.join(map(str, ptm_details.get('Ubiquitination', {}).get('sites', []))),
                'Acetyl_Count': ptm_details.get('Acetylation', {}).get('count', 0),
                'Acetyl_Min_Dist': ptm_details.get('Acetylation', {}).get('min_distance', None),
                'Acetyl_Sites': ', '.join(map(str, ptm_details.get('Acetylation', {}).get('sites', [])))
            })
    
    return pd.DataFrame(multi_ptm_analysis)

def analyze_triple_crosstalk(multi_ptm_df):
    """
    Analyze sites with triple PTM crosstalk
    """
    triple_sites = multi_ptm_df[multi_ptm_df['PTM_Count'] == 3]
    
    print(f"🌟 Analysis of Triple PTM Crosstalk Sites:")
    print(f"  Glutathionylation sites with triple crosstalk: {len(triple_sites)}")
    
    if len(triple_sites) > 0:
        print(f"  Average total interactions: {triple_sites['Total_Interactions'].mean():.1f}")
        print(f"  Average minimum distance: {triple_sites['Min_Distance'].mean():.2f} Å")
        
        # Show top triple crosstalk sites
        triple_sorted = triple_sites.sort_values('Total_Interactions', ascending=False)
        print(f"\n🏆 Top 10 Triple Crosstalk Sites:")
        
        for idx, row in triple_sorted.head(10).iterrows():
            print(f"  {row['UniProt_ID']} C{row['Gluta_Site']} ({row['Protein_Name'][:30]}...):")
            print(f"    Phosphorylation: {row['Phospho_Count']} (min {row['Phospho_Min_Dist']:.2f}Å)")
            print(f"    Ubiquitination: {row['Ubiq_Count']} (min {row['Ubiq_Min_Dist']:.2f}Å)")
            print(f"    Acetylation: {row['Acetyl_Count']} (min {row['Acetyl_Min_Dist']:.2f}Å)")
            print(f"    Total Interactions: {row['Total_Interactions']}")
            print()
    
    return triple_sites

def analyze_dual_crosstalk(multi_ptm_df):
    """
    Analyze sites with dual PTM crosstalk
    """
    dual_sites = multi_ptm_df[multi_ptm_df['PTM_Count'] == 2]
    
    print(f"🔄 Analysis of Dual PTM Crosstalk Sites:")
    print(f"  Glutathionylation sites with dual crosstalk: {len(dual_sites)}")
    
    if len(dual_sites) > 0:
        # Analyze different dual combinations
        dual_combinations = dual_sites['PTM_Types'].value_counts()
        
        print(f"\n📊 Dual Crosstalk Combination Distribution:")
        for combination, count in dual_combinations.items():
            percentage = count / len(dual_sites) * 100
            print(f"  {combination}: {count} sites ({percentage:.1f}%)")
        
        # Show top dual crosstalk sites for each combination
        for combination in dual_combinations.head(3).index:
            combo_sites = dual_sites[dual_sites['PTM_Types'] == combination]
            combo_sorted = combo_sites.sort_values('Total_Interactions', ascending=False)
            
            print(f"\n🎯 Top 5 {combination} Sites:")
            for idx, row in combo_sorted.head(5).iterrows():
                print(f"  {row['UniProt_ID']} C{row['Gluta_Site']} ({row['Protein_Name'][:25]}...): "
                      f"{row['Total_Interactions']} interactions, min {row['Min_Distance']:.2f}Å")
    
    return dual_sites

def find_spatial_clusters(all_ptm_df, distance_threshold=15):
    """
    Find spatial clusters where multiple PTMs are close to the same Cys site
    """
    print(f"\n🎪 Analysis of Spatial Clusters (≤{distance_threshold}Å):")
    
    # Filter for close interactions
    close_interactions = all_ptm_df[all_ptm_df['Distance_Angstrom'] <= distance_threshold]
    
    # Group by protein and glutathionylation site
    spatial_clusters = []
    
    for protein_gluta_site in close_interactions['Protein_Gluta_Site'].unique():
        site_data = close_interactions[close_interactions['Protein_Gluta_Site'] == protein_gluta_site]
        
        if len(site_data['PTM_Type'].unique()) >= 2:  # Multi-PTM clusters
            uniprot_id = site_data.iloc[0]['UniProt_ID']
            protein_name = site_data.iloc[0]['Protein_Name']
            gluta_site = site_data.iloc[0]['Gluta_Site']
            
            cluster_info = {
                'UniProt_ID': uniprot_id,
                'Protein_Name': protein_name,
                'Gluta_Site': gluta_site,
                'PTM_Types': ', '.join(sorted(site_data['PTM_Type'].unique())),
                'PTM_Count': len(site_data['PTM_Type'].unique()),
                'Close_Interactions': len(site_data),
                'Avg_Distance': site_data['Distance_Angstrom'].mean(),
                'Min_Distance': site_data['Distance_Angstrom'].min(),
                'Max_Distance': site_data['Distance_Angstrom'].max()
            }
            
            spatial_clusters.append(cluster_info)
    
    spatial_df = pd.DataFrame(spatial_clusters)
    
    if len(spatial_df) > 0:
        print(f"  Found {len(spatial_df)} spatial multi-PTM clusters")
        
        # Statistics by PTM count
        for ptm_count in sorted(spatial_df['PTM_Count'].unique(), reverse=True):
            count = len(spatial_df[spatial_df['PTM_Count'] == ptm_count])
            print(f"  {ptm_count}-fold PTM clusters: {count}")
        
        # Show top spatial clusters
        spatial_sorted = spatial_df.sort_values(['PTM_Count', 'Close_Interactions'], ascending=[False, False])
        print(f"\n🌟 Top 10 Spatial PTM Clusters:")
        
        for idx, row in spatial_sorted.head(10).iterrows():
            print(f"  {row['UniProt_ID']} C{row['Gluta_Site']} ({row['Protein_Name'][:25]}...):")
            print(f"    PTM Types: {row['PTM_Types']}")
            print(f"    Close Interactions: {row['Close_Interactions']} (avg {row['Avg_Distance']:.2f}Å)")
    
    return spatial_df

def main():
    print("=== Analysis of Multi-PTM Crosstalk Mechanisms ===\n")
    
    base_dir = '/SSD/Terry/Project/Gluta_phospho'
    
    # Load all crosstalk data
    print("📖 Loading all PTM crosstalk data...")
    phospho_df, ubiq_df, acetyl_df = load_all_crosstalk_data()
    
    print(f"  Phosphorylation-glutathionylation interactions: {len(phospho_df):,}")
    print(f"  Ubiquitination-glutathionylation interactions: {len(ubiq_df):,}")
    print(f"  Acetylation-glutathionylation interactions: {len(acetyl_df):,}")
    
    # Standardize data
    print(f"\n🔄 Standardizing data format...")
    phospho_clean, ubiq_clean, acetyl_clean = standardize_columns(phospho_df, ubiq_df, acetyl_df)
    
    # Find multi-PTM sites
    print(f"\n🔍 Finding multi-PTM crosstalk sites...")
    multi_ptm_df = find_multi_ptm_sites(phospho_clean, ubiq_clean, acetyl_clean)
    
    if len(multi_ptm_df) > 0:
        # Sort by interaction count
        multi_ptm_df = multi_ptm_df.sort_values(['PTM_Count', 'Total_Interactions'], ascending=[False, False])
        
        print(f"📊 Multi-PTM Crosstalk Statistics:")
        print(f"  Glutathionylation sites with multi-PTM crosstalk found: {len(multi_ptm_df):,}")
        
        # Statistics by PTM count
        ptm_counts = multi_ptm_df['PTM_Count'].value_counts().sort_index(ascending=False)
        for ptm_count, count in ptm_counts.items():
            percentage = count / len(multi_ptm_df) * 100
            print(f"  {ptm_count}-fold PTM crosstalk: {count} sites ({percentage:.1f}%)")
        
        # Analyze triple crosstalk
        triple_sites = analyze_triple_crosstalk(multi_ptm_df)
        
        # Analyze dual crosstalk  
        dual_sites = analyze_dual_crosstalk(multi_ptm_df)
        
        # Find spatial clusters
        all_ptm_combined = pd.concat([phospho_clean, ubiq_clean, acetyl_clean], ignore_index=True)
        all_ptm_combined['Protein_Gluta_Site'] = all_ptm_combined['UniProt_ID'] + '_C' + all_ptm_combined['Gluta_Site'].astype(str)
        
        spatial_clusters = find_spatial_clusters(all_ptm_combined)
        
        # Save results
        output_file = f'{base_dir}/multi_ptm_crosstalk_analysis.csv'
        multi_ptm_df.to_csv(output_file, index=False)
        print(f"\n💾 Multi-PTM analysis results saved to: {output_file}")
        
        if len(triple_sites) > 0:
            triple_output = f'{base_dir}/triple_ptm_crosstalk_sites.csv'
            triple_sites.to_csv(triple_output, index=False)
            print(f"💾 Triple PTM sites saved to: {triple_output}")
        
        if len(spatial_clusters) > 0:
            spatial_output = f'{base_dir}/spatial_ptm_clusters.csv'
            spatial_clusters.to_csv(spatial_output, index=False)
            print(f"💾 Spatial PTM clusters saved to: {spatial_output}")
        
        # Generate summary report
        summary_report = f"""
# Multi-PTM Crosstalk Mechanism Analysis Report

## Key Findings
- Glutathionylation sites with multi-PTM crosstalk: {len(multi_ptm_df):,}
- Triple PTM crosstalk sites: {len(triple_sites) if len(triple_sites) > 0 else 0}
- Dual PTM crosstalk sites: {len(dual_sites) if len(dual_sites) > 0 else 0}
- Spatial PTM clusters: {len(spatial_clusters) if len(spatial_clusters) > 0 else 0}

## Potential Regulatory Mechanisms
1. **Competitive Regulation**: Multiple PTMs may compete for the same regulatory region.
2. **Synergistic Regulation**: Different PTMs may work together to regulate protein function.
3. **Sequential Regulation**: PTMs may occur in a specific order.
4. **Spatial Confinement**: Spatial proximity influences PTM interactions.

## Biological Significance
- Multi-PTM crosstalk sites may be key regulatory nodes.
- These sites could play a hub role in cellular signal transduction.
- They provide important clues for understanding complex PTM regulatory networks.
"""
        
        summary_file = f'{base_dir}/multi_ptm_mechanisms_summary.md'
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(summary_report)
        
        print(f"📋 Mechanism analysis report saved to: {summary_file}")
        
    else:
        print("❌ No multi-PTM crosstalk sites found")
    
    print(f"\n✅ Multi-PTM crosstalk mechanism analysis complete!")

if __name__ == "__main__":
    main()
