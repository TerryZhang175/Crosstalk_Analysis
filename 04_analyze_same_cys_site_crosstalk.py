#!/usr/bin/env python3
"""
Analyze proteins where the SAME cysteine site has both phospho and ubiquitination crosstalk
"""

import pandas as pd
import os

def main():
    print("=== 同一半胱氨酸位点的双重交互分析 ===\n")
    
    base_dir = '/SSD/Terry/Project/Gluta_phospho'
    
    # 1. Load phospho-gluta interactions
    phospho_gluta_file = f'{base_dir}/Final_Results/02_glutathionylation_phosphorylation_interactions.csv'
    phospho_gluta_df = pd.read_csv(phospho_gluta_file)
    
    # Filter for close interactions only
    if 'distance_angstrom' in phospho_gluta_df.columns:
        phospho_gluta_df = phospho_gluta_df[phospho_gluta_df['distance_angstrom'] <= 10.0]
    elif 'Distance_Angstrom' in phospho_gluta_df.columns:
        phospho_gluta_df = phospho_gluta_df[phospho_gluta_df['Distance_Angstrom'] <= 10.0]
    
    print(f"Loading phospho-gluta interactions: {len(phospho_gluta_df)} close interactions")
    
    # 2. Load ubiq-gluta interactions
    ubiq_gluta_file = f'{base_dir}/Final_Results/01_ubiquitination_glutathionylation_interactions.csv'
    ubiq_gluta_df = pd.read_csv(ubiq_gluta_file)
    
    print(f"Loading ubiq-gluta interactions: {len(ubiq_gluta_df)} close interactions")
    
    # 3. Standardize column names
    # For phospho-gluta
    if 'uniprot_id' in phospho_gluta_df.columns:
        phospho_gluta_df = phospho_gluta_df.rename(columns={
            'uniprot_id': 'UniProt_ID',
            'phospho_site': 'Phospho_Site', 
            'gluta_site': 'Gluta_Site',
            'distance_angstrom': 'Distance'
        })
    
    # For ubiq-gluta
    if 'Ubiquitination_Site' in ubiq_gluta_df.columns:
        ubiq_gluta_df = ubiq_gluta_df.rename(columns={
            'Ubiquitination_Site': 'Ubiq_Site',
            'Glutathionylation_Site': 'Gluta_Site',
            'Distance_Angstrom': 'Distance'
        })
    
    # 4. Find proteins with the same Cys site in both interactions
    print("\n🔍 查找同一半胱氨酸位点的双重交互...\n")
    
    # Create unique keys for protein-cys combinations
    phospho_gluta_df['protein_cys'] = phospho_gluta_df['UniProt_ID'] + '_C' + phospho_gluta_df['Gluta_Site'].astype(str)
    ubiq_gluta_df['protein_cys'] = ubiq_gluta_df['UniProt_ID'] + '_C' + ubiq_gluta_df['Gluta_Site'].astype(str)
    
    # Find common protein-cys combinations
    phospho_cys_set = set(phospho_gluta_df['protein_cys'].unique())
    ubiq_cys_set = set(ubiq_gluta_df['protein_cys'].unique())
    
    same_cys_sites = phospho_cys_set.intersection(ubiq_cys_set)
    
    print(f"📊 统计结果:")
    print(f"  磷酸化-谷胱甘肽化的独特Cys位点: {len(phospho_cys_set):,}")
    print(f"  泛素化-谷胱甘肽化的独特Cys位点: {len(ubiq_cys_set):,}")
    print(f"  同一Cys位点参与两种交互: {len(same_cys_sites):,}")
    print(f"  重合率: {len(same_cys_sites)/min(len(phospho_cys_set), len(ubiq_cys_set))*100:.1f}%")
    
    if len(same_cys_sites) > 0:
        # 5. Analyze these special sites
        print(f"\n🎯 同一Cys位点双重交互的详细分析:")
        
        detailed_results = []
        
        for protein_cys in same_cys_sites:
            protein_id, cys_site = protein_cys.split('_C')
            cys_site = int(cys_site)
            
            # Get phospho interactions for this Cys
            phospho_matches = phospho_gluta_df[phospho_gluta_df['protein_cys'] == protein_cys]
            
            # Get ubiq interactions for this Cys
            ubiq_matches = ubiq_gluta_df[ubiq_gluta_df['protein_cys'] == protein_cys]
            
            # Get protein name
            protein_name = "Unknown"
            if 'Protein_Name' in phospho_matches.columns and len(phospho_matches) > 0:
                protein_name = phospho_matches.iloc[0]['Protein_Name']
            elif 'Protein_Name' in ubiq_matches.columns and len(ubiq_matches) > 0:
                protein_name = ubiq_matches.iloc[0]['Protein_Name']
            
            # Collect all phospho sites interacting with this Cys
            phospho_sites = phospho_matches['Phospho_Site'].unique() if 'Phospho_Site' in phospho_matches.columns else []
            phospho_distances = phospho_matches['Distance'].values if 'Distance' in phospho_matches.columns else []
            
            # Collect all ubiq sites interacting with this Cys
            ubiq_sites = ubiq_matches['Ubiq_Site'].unique() if 'Ubiq_Site' in ubiq_matches.columns else []
            ubiq_distances = ubiq_matches['Distance'].values if 'Distance' in ubiq_matches.columns else []
            
            detailed_results.append({
                'UniProt_ID': protein_id,
                'Protein_Name': protein_name,
                'Cys_Site': cys_site,
                'Phospho_Partners': len(phospho_sites),
                'Phospho_Sites': ','.join(map(str, phospho_sites)),
                'Min_Phospho_Distance': min(phospho_distances) if len(phospho_distances) > 0 else None,
                'Ubiq_Partners': len(ubiq_sites),
                'Ubiq_Sites': ','.join(map(str, ubiq_sites)),
                'Min_Ubiq_Distance': min(ubiq_distances) if len(ubiq_distances) > 0 else None,
                'Total_Interactions': len(phospho_sites) + len(ubiq_sites)
            })
        
        # Create DataFrame and sort by total interactions
        results_df = pd.DataFrame(detailed_results)
        results_df = results_df.sort_values('Total_Interactions', ascending=False)
        
        # Save results
        output_file = f'{base_dir}/Final_Results/same_cys_dual_crosstalk.csv'
        results_df.to_csv(output_file, index=False)
        print(f"\n📁 详细结果已保存到: {output_file}")
        
        # Show top results
        print(f"\n🏆 Top 10 同一Cys位点双重交互:")
        top_10 = results_df.head(10)
        
        for idx, row in top_10.iterrows():
            print(f"\n{row['UniProt_ID']} - {row['Protein_Name']}")
            print(f"  Cys位点: C{row['Cys_Site']}")
            print(f"  磷酸化伙伴: {row['Phospho_Partners']}个 (最短距离: {row['Min_Phospho_Distance']:.1f}Å)")
            print(f"  泛素化伙伴: {row['Ubiq_Partners']}个 (最短距离: {row['Min_Ubiq_Distance']:.1f}Å)")
        
        # 6. Statistical summary
        print(f"\n📈 统计总结:")
        
        # Count proteins (not sites)
        unique_proteins = results_df['UniProt_ID'].nunique()
        print(f"  涉及的独特蛋白质数: {unique_proteins}")
        print(f"  平均每个Cys位点的磷酸化伙伴: {results_df['Phospho_Partners'].mean():.1f}")
        print(f"  平均每个Cys位点的泛素化伙伴: {results_df['Ubiq_Partners'].mean():.1f}")
        
        # Check if same phospho/ubiq sites interact with multiple Cys
        print(f"\n🔬 交互模式分析:")
        
        # Sites with both close phospho and ubiq partners
        very_close = results_df[(results_df['Min_Phospho_Distance'] <= 5) & 
                                (results_df['Min_Ubiq_Distance'] <= 5)]
        print(f"  极近距离双重交互 (两种都≤5Å): {len(very_close)}个Cys位点")
        
        # Hub sites (multiple partners)
        hub_sites = results_df[results_df['Total_Interactions'] >= 5]
        print(f"  枢纽Cys位点 (≥5个交互伙伴): {len(hub_sites)}个")
        
        # 7. Biological interpretation
        print(f"\n💡 生物学意义:")
        print(f"  1. {len(same_cys_sites)}个Cys位点同时参与两种调控")
        print(f"  2. 这些位点可能是关键的调控开关")
        print(f"  3. 氧化还原状态可能同时影响磷酸化信号和泛素化降解")
        print(f"  4. 这些Cys位点可能是药物靶点的理想候选")
        
    else:
        print(f"\n❌ 没有发现同一Cys位点同时参与两种交互")
    
    # 8. Extended analysis - check if they're on same proteins but different Cys
    print(f"\n" + "="*60)
    print(f"扩展分析：同一蛋白质的不同Cys位点")
    print("="*60)
    
    # Get proteins with both types of interactions
    phospho_proteins = set(phospho_gluta_df['UniProt_ID'].unique())
    ubiq_proteins = set(ubiq_gluta_df['UniProt_ID'].unique())
    both_proteins = phospho_proteins.intersection(ubiq_proteins)
    
    print(f"有两种交互的蛋白质数: {len(both_proteins)}")
    
    if len(both_proteins) > 0:
        # Check how many have same vs different Cys sites
        proteins_same_cys = set([pc.split('_')[0] for pc in same_cys_sites])
        proteins_diff_cys = both_proteins - proteins_same_cys
        
        print(f"  同一Cys位点: {len(proteins_same_cys)}个蛋白质")
        print(f"  不同Cys位点: {len(proteins_diff_cys)}个蛋白质")
        
        if len(proteins_same_cys) > 0:
            print(f"\n✨ 关键发现：")
            print(f"  {len(proteins_same_cys)}个蛋白质的同一Cys位点是双重调控节点！")
            print(f"  这占所有双重交互蛋白质的 {len(proteins_same_cys)/len(both_proteins)*100:.1f}%")

if __name__ == "__main__":
    main()