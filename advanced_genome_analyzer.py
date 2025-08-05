from Bio import SeqIO
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
import matplotlib.pyplot as plt
from collections import defaultdict, Counter

def calculate_gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2) if len(seq) > 0 else 0

def codon_usage(seq):
    codon_count = defaultdict(int)
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            codon_count[codon] += 1
    return codon_count

def gc_skew(seq, window=1000):
    skew = []
    for i in range(0, len(seq)-window, window):
        window_seq = seq[i:i+window].upper()
        g = window_seq.count('G')
        c = window_seq.count('C')
        if g + c > 0:
            skew.append((i, (g - c)/(g + c)))
    return skew

def plot_gc_skew(skew_data, output_file="gc_skew_plot.png"):
    positions, skews = zip(*skew_data)
    plt.figure(figsize=(10, 4))
    plt.plot(positions, skews)
    plt.title("GC Skew Plot")
    plt.xlabel("Position")
    plt.ylabel("GC Skew")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def find_orfs(seq, min_len=100):
    orfs = []
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    seq = seq.upper()
    for frame in range(3):
        i = frame
        while i < len(seq) - 3:
            codon = seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(seq)-3, 3):
                    stop = seq[j:j+3]
                    if stop in stop_codons:
                        orf = seq[i:j+3]
                        if len(orf) >= min_len:
                            orfs.append((i, j+3, orf))
                        break
                i = j
            else:
                i += 3
    return orfs

def calculate_cai(seq_list, ref_seqs):
    cai = CodonAdaptationIndex()
    ref_seqs_str = [str(seq) for seq in ref_seqs]
    cai.generate_index(ref_seqs_str)
    return [cai.cai_for_gene(str(seq)) for seq in seq_list]

def analyze_genome(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq.upper()
        print(f"Analyzing {record.id}")
        print(f"Total length: {len(seq)} bp")
        print(f"GC Content: {calculate_gc_content(seq)}%")

        codons = codon_usage(seq)
        print("\nTop 10 most used codons:")
        for codon, count in Counter(codons).most_common(10):
            print(f"{codon}: {count}")

        print("\nGenerating GC skew plot...")
        skew = gc_skew(seq)
        plot_gc_skew(skew)
        print("Plot saved as 'gc_skew_plot.png'")

        print("\nFinding ORFs (min 100 bp)...")
        orfs = find_orfs(seq)
        print(f"Found {len(orfs)} ORFs.")
        if orfs:
            print("First 5 ORFs (start-end, length):")
            for i, (start, end, orf_seq) in enumerate(orfs[:5]):
                print(f"ORF {i+1}: {start}-{end} ({end - start} bp)")

        print("\nCalculating Codon Adaptation Index (CAI) for ORFs...")
        orf_seqs = [orf_seq for _, _, orf_seq in orfs]
        if len(orf_seqs) > 5:
            ref_seqs = orf_seqs[:5]
            cai_scores = calculate_cai(orf_seqs, ref_seqs)
            print("CAI Scores for first 5 ORFs:")
            for i, score in enumerate(cai_scores[:5]):
                print(f"ORF {i+1}: CAI = {round(score, 4)}")
        else:
            print("Not enough ORFs to calculate CAI.")

        print("\nAnalysis complete.")
        print("="*50)

# Change the filename here to your downloaded .fna genome file
if __name__ == "__main__":
    fasta_file = "example_bacteria.fna"
    analyze_genome(fasta_file)

