# Electrical-_Field-assisted_PCR

# Figure 4
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
2. **Quality Statistics**: Run `Fig_4C_ii_matlabcode.m`. Enter the target DNA length (`100`) when prompted to generate Q-score CSV files.
4. **Error Analysis**: Run `Fig_4E_matlabcode.m` to calculate specific mutation rates (InDels and Substitutions).

---

## 🧬 Reference Sequences

* **Target Length**: 100 bp
* **Forward Primer**: `ATGCACTCAGACGTACAACG`
* **Reverse Primer**: `CATTACACCTCACGCCTTCG`

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Figure 5
📋 Overview
The pipeline takes raw FASTQ reads, filters them based on a high-quality Phred score threshold, ensures all sequences are oriented correctly using primer sequences, and finally visualizes the percentage of A, T, C, and G at each nucleotide position.

🛠 Features
Data Integration: Merges reads from multiple FASTQ files into a single dataset for analysis.

Quality Control (QC): Filters out low-quality reads by calculating the mean Phred quality score (threshold set to 35).

Sequence Standardization: Automatically detects and corrects sequence orientation (reverse complement) by matching against known primer sequences.

Positional Analysis: Calculates the frequency of each nucleotide (A, T, C, G) across a specific sequence length (e.g., 90 bp).

Visualization: Generates a histogram for read length distribution and a stacked bar chart for base composition.

🚀 Workflow Details
1. Data Loading and Quality Filtering
The script reads controlled DNA synthesis_1.fastq and controlled DNA synthesis_2.fastq.

It calculates the Phred quality score by subtracting 33 from the ASCII quality characters.

Reads with a mean quality score of less than 35 are discarded.

2. Orientation Correction
Defined primers:

F_primers: CGACACTGGCTACCAGTACG

R_primers: GTACACCTCACGCCTTCACC

The script checks if the first 20 base pairs of a read match the reverse complement of the primer (allowing up to 2 mismatches).

If a match is found, the entire sequence is converted to its reverse complement to ensure uniform directionality.

3. Base Composition Analysis
The script isolates reads with a specific length of 90 bp.

It iterates through each position (1 to 90) and counts the occurrences of 'A', 'T', 'C', and 'G'.

The counts are converted into percentages for relative distribution analysis.

📊 Outputs
Read Length Histogram: Displays the distribution of sequence lengths after quality filtering.

Stacked Bar Chart: A visual representation of the nucleotide distribution at each position from 1 to 90.

Legend: Blue (A), Orange (T), Yellow (C), Purple (G) — Note: Colors may vary based on MATLAB's default color map.

⚙️ Requirements
MATLAB

Bioinformatics Toolbox (required for fastqread)

Input FASTQ files: controlled DNA synthesis_1.fastq and controlled DNA synthesis_2.fastq

## ⚠️ Notes

* The script `Fig_4C_ii_matlabcode.m` requires external helper functions (e.g., `evertComplementary`, `PrimerSorting`) to be available in the MATLAB path.
* Primer matching allows for a small margin of error (up to 2-5 mismatches depending on the script).
