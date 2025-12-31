# Electrical-_Field-assisted_PCR

# DNA Sequence Analysis for E-PCR Results

This project provides a suite of MATLAB scripts designed to analyze and visualize DNA sequencing data (FASTQ format) obtained from E-PCR (Emulsion PCR). The tools focus on sequence alignment, quality control, and detailed error profiling (insertions, deletions, and substitutions).

## 📂 File Structure

| Filename | Description | Key Outputs |
| --- | --- | --- |
| `Fig_4C_i_and_4D_matlabcode.m` | Analyzes sequence length distribution and position-specific base ratios. | DNA length histogram, Stacked bar chart of base distribution. |
| `Fig_4C_ii_matlabcode.m` | Performs statistical analysis on Quality Scores (Q-scores) for specific DNA lengths. | Q-score Mean and Standard Deviation (saved as .csv). |
| `Fig_4E_matlabcode.m` | Calculates the frequency of sequence variations relative to a template. | Counts for Insertion (+1/+2), Deletion (-1/-2), and Substitution (1/2). |

---

## 🚀 Key Features

### 1. Data Preprocessing & Alignment

* **FASTQ Loading**: Loads sequence (`S`) and quality (`Q`) data using the `fastqread` function.
* **Directional Correction**: Automatically detects and converts sequences to their reverse complement based on specific primer sequences to ensure consistent alignment.
* **Quality Filtering**: Filters for high-quality reads where every base in a 100bp sequence exceeds a defined Q-score threshold.

### 2. Analytical Capabilities

* **Length Distribution**: Generates a histogram to visualize the variation in DNA fragment lengths.
* **Base Composition Analysis**: Computes the percentage of A, T, C, and G at every nucleotide position (1–100 bp) and displays the results in a stacked bar chart.
* **Error Profiling**:
* **Insertions (+1, +2)**: Identifies sequences longer than the target length that match the reverse primer.
* **Deletions (-1, -2)**: Identifies sequences shorter than the target length that match the reverse primer.
* **Substitutions (1, 2)**: Detects single or double nucleotide mismatches in sequences of the correct length.



### 3. Data Export

* Saves statistical results (Mean and Standard Deviation of Q-scores) directly to CSV files for further reporting.

---

## 🛠 Usage & Requirements

### Prerequisites

* **MATLAB** (with **Bioinformatics Toolbox** for `fastqread`).
* Input files `E-PCR_1.fastq` and `E-PCR_2.fastq` must be present in the working directory.

### Standard Workflow

1. **Initial Overview**: Run `Fig_4C_i_and_4D_matlabcode.m` to view the length histogram and base distribution.
2. **Quality Statistics**: Run `Fig_4C_ii_matlabcode.m`. Enter the target DNA length (e.g., `100`) when prompted to generate Q-score CSV files.
3. **Error Analysis**: Run `Fig_4E_matlabcode.m` to calculate specific mutation rates (InDels and Substitutions).

---

## 🧬 Reference Sequences

* **Target Length**: 100 bp
* **Forward Primer**: `ATGCACTCAGACGTACAACG`
* **Reverse Primer**: `CGAAGGCGTGAGGTGTAATG` or `CATTACACCTCACGCCTTCG`

## ⚠️ Notes

* The script `Fig_4C_ii_matlabcode.m` requires external helper functions (e.g., `evertComplementary`, `PrimerSorting`) to be available in the MATLAB path.
* Primer matching allows for a small margin of error (up to 2-5 mismatches depending on the script).
