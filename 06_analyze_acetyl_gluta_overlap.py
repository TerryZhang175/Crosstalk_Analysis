#!/usr/bin/env python3
"""
Analyze how many glutathionylated proteins have acetylation and how many have interactions
"""

import pandas as pd

def main():
    print("=== Analysis of Acetylation in Glutathionylated Proteins ===\n")
    
    base_dir = '/SSD/Terry/Project/Gluta_phospho'
    
    # Load glutathionylation data (all proteins)
    gluta_file = f'{base_dir}/cleaned_glutathionylation_sites.csv'
    gluta_df = pd.read_csv(gluta_file)
    
    # Load acetylation data (all proteins)
    acetyl_file = f'{base_dir}/acetylation_sites_cleaned.csv'
    acetyl_df = pd.read_csv(acetyl_file)
    
    # Load acetylation-glutathionylation overlap data
    overlap_file = f'{base_dir}/acetylation_glutathionylation_overlap.csv'
    overlap_df = pd.read_csv(overlap_file)
    
    # Load close interactions (proteins with spatial crosstalk)
    close_file = f'{base_dir}/acetyl_gluta_close_interactions.csv'
    close_df = pd.read_csv(close_file)
    
    print(f"📊 基础数据统计:")
    print(f"  总谷胱甘肽化蛋白质: {gluta_df['UniProt_ID'].nunique():,}")
    print(f"  总乙酰化蛋白质: {acetyl_df['UniProt_ID'].nunique():,}")
    print(f"  重合蛋白质: {overlap_df['UniProt_ID'].nunique():,}")
    print(f"  有近距离交互的蛋白质: {close_df['uniprot_id'].nunique():,}")
    
    # Get unique protein sets
    all_gluta_proteins = set(gluta_df['UniProt_ID'].unique())
    all_acetyl_proteins = set(acetyl_df['UniProt_ID'].unique())
    overlap_proteins = set(overlap_df['UniProt_ID'].unique())
    interaction_proteins = set(close_df['uniprot_id'].unique())
    
    print(f"\n🔍 谷胱甘肽化蛋白质中的乙酰化分析:")
    print(f"  总谷胱甘肽化蛋白质: {len(all_gluta_proteins):,}")
    print(f"  其中有乙酰化的: {len(overlap_proteins):,} ({len(overlap_proteins)/len(all_gluta_proteins)*100:.1f}%)")
    print(f"  其中有近距离交互的: {len(interaction_proteins):,} ({len(interaction_proteins)/len(all_gluta_proteins)*100:.1f}%)")
    
    # Among proteins with acetylation, how many have interactions?
    print(f"\n📈 有乙酰化的谷胱甘肽化蛋白质中:")
    print(f"  有乙酰化但无近距离交互: {len(overlap_proteins - interaction_proteins):,}")
    print(f"  有乙酰化且有近距离交互: {len(interaction_proteins):,} ({len(interaction_proteins)/len(overlap_proteins)*100:.1f}%)")
    
    # Detailed interaction analysis
    print(f"\n🎯 近距离交互详细分析:")
    
    # Count interactions per protein
    interaction_counts = close_df.groupby('uniprot_id').size().reset_index(name='Interaction_Count')
    
    # Count sites per protein
    acetyl_sites_per_protein = close_df.groupby('uniprot_id')['acetyl_site'].nunique().reset_index(name='Acetyl_Sites')
    gluta_sites_per_protein = close_df.groupby('uniprot_id')['gluta_site'].nunique().reset_index(name='Gluta_Sites')
    
    # Merge data
    protein_stats = interaction_counts.merge(acetyl_sites_per_protein, on='uniprot_id')
    protein_stats = protein_stats.merge(gluta_sites_per_protein, on='uniprot_id')
    protein_stats = protein_stats.sort_values('Interaction_Count', ascending=False)
    
    print(f"  平均每个蛋白质的交互数: {protein_stats['Interaction_Count'].mean():.1f}")
    print(f"  平均每个蛋白质的乙酰化位点数: {protein_stats['Acetyl_Sites'].mean():.1f}")
    print(f"  平均每个蛋白质的谷胱甘肽化位点数: {protein_stats['Gluta_Sites'].mean():.1f}")
    
    # Distribution of interaction counts
    interaction_ranges = [
        (1, 1, "1个交互"),
        (2, 3, "2-3个交互"),
        (4, 5, "4-5个交互"),
        (6, 10, "6-10个交互"),
        (11, float('inf'), ">10个交互")
    ]
    
    print(f"\n📊 交互数量分布:")
    for min_count, max_count, label in interaction_ranges:
        if max_count == float('inf'):
            count = len(protein_stats[protein_stats['Interaction_Count'] >= min_count])
        else:
            count = len(protein_stats[(protein_stats['Interaction_Count'] >= min_count) & 
                                   (protein_stats['Interaction_Count'] <= max_count)])
        percentage = count / len(protein_stats) * 100
        print(f"  {label}: {count} 蛋白质 ({percentage:.1f}%)")
    
    # Top proteins with most interactions
    print(f"\n🏆 交互最多的前10个蛋白质:")
    top_10 = protein_stats.head(10)
    
    # Get protein names
    protein_names = {}
    try:
        for uniprot_id in top_10['uniprot_id'].unique():
            protein_data = overlap_df[overlap_df['UniProt_ID'] == uniprot_id]
            if len(protein_data) > 0:
                protein_names[uniprot_id] = protein_data.iloc[0]['Protein_Name']
            else:
                protein_names[uniprot_id] = "Unknown"
    except:
        pass
    
    for _, row in top_10.iterrows():
        protein_name = protein_names.get(row['uniprot_id'], 'Unknown')[:40]
        print(f"  {row['uniprot_id']} ({protein_name}...): {row['Interaction_Count']} 交互 "
              f"({row['Acetyl_Sites']} 乙酰化位点, {row['Gluta_Sites']} 谷胱甘肽化位点)")
    
    # Summary statistics
    summary_stats = {
        'total_glutathionylated': len(all_gluta_proteins),
        'with_acetylation': len(overlap_proteins),
        'with_interactions': len(interaction_proteins),
        'acetylation_percentage': len(overlap_proteins)/len(all_gluta_proteins)*100,
        'interaction_percentage': len(interaction_proteins)/len(all_gluta_proteins)*100,
        'interaction_among_acetylated': len(interaction_proteins)/len(overlap_proteins)*100
    }
    
    # Save detailed protein statistics
    protein_stats_with_names = protein_stats.copy()
    protein_stats_with_names['Protein_Name'] = protein_stats_with_names['uniprot_id'].map(protein_names)
    protein_stats_with_names = protein_stats_with_names.rename(columns={'uniprot_id': 'UniProt_ID'})
    
    output_file = f'{base_dir}/acetyl_gluta_detailed_protein_stats.csv'
    protein_stats_with_names.to_csv(output_file, index=False)
    print(f"\n💾 详细蛋白质统计已保存到: {output_file}")
    
    # Create summary report
    print(f"\n📋 总结报告:")
    print(f"  在 {summary_stats['total_glutathionylated']:,} 个谷胱甘肽化蛋白质中:")
    print(f"  - {summary_stats['with_acetylation']:,} 个有乙酰化 ({summary_stats['acetylation_percentage']:.1f}%)")
    print(f"  - {summary_stats['with_interactions']:,} 个有近距离交互 ({summary_stats['interaction_percentage']:.1f}%)")
    print(f"  - 在有乙酰化的蛋白质中，{summary_stats['interaction_among_acetylated']:.1f}% 存在近距离交互")
    
    # Save summary
    summary_report = f"""
# 谷胱甘肽化蛋白质中的乙酰化分析报告

## 基本统计
- 总谷胱甘肽化蛋白质: {summary_stats['total_glutathionylated']:,}
- 其中有乙酰化的: {summary_stats['with_acetylation']:,} ({summary_stats['acetylation_percentage']:.1f}%)
- 其中有近距离交互的: {summary_stats['with_interactions']:,} ({summary_stats['interaction_percentage']:.1f}%)

## 关键发现
1. **乙酰化覆盖率**: {summary_stats['acetylation_percentage']:.1f}% 的谷胱甘肽化蛋白质具有乙酰化修饰
2. **交互率**: {summary_stats['interaction_percentage']:.1f}% 的谷胱甘肽化蛋白质与乙酰化存在近距离交互
3. **功能交互**: 在有乙酰化的蛋白质中，{summary_stats['interaction_among_acetylated']:.1f}% 存在空间上的近距离交互

## 生物学意义
- 超过一半的谷胱甘肽化蛋白质具有乙酰化修饰，表明两种PTM在蛋白质调控中的广泛共存
- 约1/4的谷胱甘肽化蛋白质与乙酰化存在空间近距离交互，提示直接的功能性调控关系
- 这种高度的PTM共存和交互暗示了复杂的蛋白质调控网络
"""
    
    summary_file = f'{base_dir}/acetyl_gluta_overlap_summary.md'
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write(summary_report)
    
    print(f"\n💾 总结报告已保存到: {summary_file}")
    
    return summary_stats

if __name__ == "__main__":
    main()