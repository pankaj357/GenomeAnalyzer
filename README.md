
# 🧬 Advanced Genome Analyzer

A mini bioinformatics project written in Python to analyze bacterial genome sequences. This tool calculates GC content, codon usage, GC skew, detects open reading frames (ORFs), and computes Codon Adaptation Index (CAI) for predicted coding sequences.

---

## 🚀 Features

- ✅ GC Content Calculation (whole genome)
- ✅ Codon Usage Analysis (frequency of each codon)
- ✅ GC-Skew Plot Generation (detect replication origin)
- ✅ ORF Detection (start and stop codons in 3 frames)
- ✅ Codon Adaptation Index (CAI) Calculation

---

## 📦 Requirements

Install dependencies using:

```bash
pip install biopython matplotlib
```

---

## 📂 File Structure

```
advanced-genome-analyzer/
├── advanced_genome_analyzer.py  # Main script
├── example_bacteria.fna         # Example genome file (FASTA)
├── gc_skew_plot.png             # Output GC skew plot
└── README.md                    # Project documentation
```

---

## ▶️ Usage

1. Place your genome `.fna` file in the folder and name it `example_bacteria.fna` (or change filename in the script).
2. Run the script:

```bash
python advanced_genome_analyzer.py
```

3. Output:
   - Prints GC content, codon usage, and ORFs in terminal
   - Saves GC skew plot as `gc_skew_plot.png`

---

## 📊 Output Example

```bash
Analyzing NZ_CP011113.1
Total length: 4641652 bp
GC Content: 50.79%

Top 10 most used codons:
GAA: 12412
GAG: 11230
GCT: 11058
...

Generating GC skew plot...
Plot saved as 'gc_skew_plot.png'

Finding ORFs (min 100 bp)...
Found 257 ORFs.
ORF 1: 202-685 (483 bp)
...

Calculating Codon Adaptation Index (CAI)...
ORF 1: CAI = 0.73
...
```

---

## 📌 Author

**Pankaj**  
Bioinformatics Researcher  
Indian Agricultural Research Institute (IARI), New Delhi  
Contact: [ft.pank@gmail.com](mailto:ft.pank@gmail.com)

---

## 📜 License

This project is open-source and available under the MIT License.
