import pysam
import matplotlib.pyplot as plt
import argparse

def calculate_unique_molecules(bam_file):
    """Calculate number of unique molecules for an increasing number of processed reads."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    unique_molecules = set()
    unique_counts = []
    num_reads = []
    saturation_index = 0
    for i,read in enumerate(bam):
        # molecule is defined as (UMI, chromosome, start pos, end pos)
        molecule = (read.query_name.split(":")[-1], read.reference_name, read.reference_start, read.reference_end)
        # sets are unique, so it won't add a duplicate molecule
        unique_molecules.add(molecule)
        # print some stats every 100k reads
        if i % 100000 == 0:
            print(saturation_index)
            print(len(unique_molecules))
            print(i)
             # append number of unique molecules at each iteration
            unique_counts.append(len(unique_molecules))
            # append number of total reads at each iteration
            num_reads.append(i+1)
            saturation_index = 1 - (len(unique_molecules) / (i+1))
        
    bam.close()
    return unique_counts, num_reads

def plot_unique_vs_total_reads(read_counts, unique_counts):
    """Plot number of unique molecules vs. total processed reads."""
    plt.plot(read_counts, unique_counts)
    plt.plot(range(len(unique_counts)), range(len(unique_counts)), 'k--', label="x=y")
    plt.xlabel("Total Reads Processed")
    plt.ylabel("Unique Molecules Recovered")
    plt.title("Unique Molecules vs. Total Reads Processed")
    plt.grid(True)
    plt.show()
    plt.savefig('saturation_test.png')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot number of unique molecules vs. total processed reads from a BAM file.")
    parser.add_argument("--bam", type=str, help="Path to the input BAM file.")
    args = parser.parse_args()
    read_counts, unique_counts = calculate_unique_molecules(args.bam)
    plot_unique_vs_total_reads(read_counts, unique_counts)