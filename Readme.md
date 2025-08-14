# **PTM Crosstalk Analysis Toolkit**

## **Overview**

This toolkit provides a set of Python scripts for analyzing **protein post-translational modification (PTM) crosstalk** with a focus on **glutathionylation** and its spatial relationships with other PTMs — including **phosphorylation**, **ubiquitination**, and **acetylation**.
 It also includes **AlphaFill-based analyses** to investigate **ATP/ADP binding proximity** to glutathionylated cysteine sites.

The toolkit leverages:

- **AlphaFold** structural models for distance measurements
- **BioPython** for PDB/CIF parsing
- **Pandas/Numpy** for data processing
- **Custom filtering** to focus on biologically relevant residue types

------

## **Directory Structure**

```
project_root/
│
├── 01_calculate_phospho_gluta_distances.py
├── 02_calculate_ubiq_gluta_distances.py
├── 03_calculate_acetyl_gluta_distances.py
├── 04_analyze_same_cys_site_crosstalk.py
├── 05_analyze_multi_ptm_crosstalk.py
├── 06_analyze_acetyl_gluta_overlap.py
├── alphafill_atp_analysis.py
├── alphafill_full_analysis.py
└── data/  # Expected CSV and structural data
```

------

## **Script Descriptions**

### **1. Distance Calculations**

- **`01_calculate_phospho_gluta_distances.py`**
   Calculates minimum Euclidean distances between phosphorylation and glutathionylation sites within the same protein, using **all atoms** of both residues.
- **`02_calculate_ubiq_gluta_distances.py`**
   Calculates distances between ubiquitination and glutathionylation sites, using heavy atom coordinates from AlphaFold models.
- **`03_calculate_acetyl_gluta_distances.py`**
   Calculates distances between acetylation and glutathionylation sites, filters for biologically relevant residue types (Lys for acetylation, Cys for glutathionylation), and identifies close interactions (≤ 10 Å).

------

### **2. Crosstalk Analyses**

- **`04_analyze_same_cys_site_crosstalk.py`**
   Identifies cases where the **same cysteine site** is involved in both **phospho-gluta** and **ubiq-gluta** interactions, calculates overlap statistics, and classifies interaction patterns.
- **`05_analyze_multi_ptm_crosstalk.py`**
   Integrates phospho-, ubiq-, and acetyl-gluta datasets to:
  - Identify **multi-PTM sites**
  - Classify into **dual** or **triple crosstalk**
  - Detect **spatial PTM clusters** within a distance threshold
  - Generate summary reports
- **`06_analyze_acetyl_gluta_overlap.py`**
   Quantifies overlap between acetylation and glutathionylation proteins, computes interaction statistics, and produces summary reports on coverage, proximity, and network implications.

------

### **3. AlphaFill ATP/ADP Binding Analyses**

- **`alphafill_atp_analysis.py`**
   Downloads AlphaFill-enriched structures for selected proteins, identifies ATP/ADP binding sites, and calculates distances from cysteine SG atoms at glutathionylation sites.
- **`alphafill_full_analysis.py`**
   Extends ATP/ADP analysis to all available proteins, supports multiple nucleotide types (ATP, ADP, GTP, GDP, AMP), and outputs:
  - Distance summaries
  - Closest nucleotide type per site
  - Interaction subsets (≤ 10 Å)

------

## **Data Requirements**

The scripts expect pre-processed CSV datasets, including:

- **PTM site lists** (`*_sites.csv`)
- **Overlap datasets** (e.g., `acetylation_glutathionylation_overlap.csv`)
- **AlphaFold PDB files** (`AF-<UniProt_ID>-F1-model_v4.pdb`)
- **AlphaFill CIF files** (downloaded automatically if missing)

**Default data directory:**
 `/SSD/Terry/Project/Gluta_phospho/`

------

## **Usage**

Each script is standalone and can be run as:

```bash
python script_name.py
```

**Examples:**

```bash
python 03_calculate_acetyl_gluta_distances.py
python 05_analyze_multi_ptm_crosstalk.py
python alphafill_full_analysis.py
```

------

## **Outputs**

Depending on the script, results include:

- Raw and filtered **distance tables** (`*_distances_raw.csv`, `*_distances_filtered.csv`)
- **Close interaction subsets** (≤ 10 Å)
- **Crosstalk analysis summaries** (`*_summary.md`)
- **AlphaFill nucleotide distance tables**
- Intermediate CSV/JSON files for large-scale runs

------

## **Key Features**

- Handles multiple PTM types with unified analysis pipelines
- Integrates **structural biology** with **PTM network analysis**
- Generates **both detailed per-protein results** and **high-level summaries**
- Flexible to adapt to different PTM datasets or distance thresholds

------

## **Dependencies**

- Python 3.8+
- `pandas`
- `numpy`
- `biopython`
- `requests`
- `pathlib`

Install via:

```bash
pip install pandas numpy biopython requests
```

------

## **Citation**

If using this toolkit for publication, please cite the relevant **AlphaFold**, **AlphaFill**, and **PTM dataset** sources, along with this analysis pipeline.

------


Do you want me to also create a **flowchart figure** showing how the different scripts connect in your analysis workflow? That would make the README even clearer.
