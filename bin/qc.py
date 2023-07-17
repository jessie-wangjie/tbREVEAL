import argparse

def convert_csv_to_interval_list(csv_file, output_file):
    intervals = []  # List to store interval data

    with open(csv_file, 'r') as f:
        # Read and process each line
        for i,line in enumerate(f):
            if i == 0:
                continue
            row = line.strip().split(',')
            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            name = row[3]

            intervals.append((chromosome, start, end, name))  # Store interval data

    # Sort intervals by chromosome and start position
    sorted_intervals = sorted(intervals, key=lambda x: (chromosome_sort_key(x[0]), x[1]))

    with open(output_file, 'w') as outfile:
        # Write header
        outfile.write('@HD\tVN:1.0\tSO:coordinate\n')

        # Write sorted intervals to the output file
        for interval in sorted_intervals:
            chromosome, start, end, name = interval
            outfile.write(f"chr{chromosome}\t{start}\t{end}\t+\t{name}\n")

    print(f"Conversion completed. Interval list file saved as: {output_file}")


def chromosome_sort_key(chromosome):
    # Custom sorting key function for chromosome names
    if chromosome.isdigit():  # If chromosome is a number
        return int(chromosome)
    elif chromosome == 'X':
        return 23
    elif chromosome == 'Y':
        return 24
    else:
        return chromosome


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert CSV to Picard Interval List format.')
    parser.add_argument('--metadata_file', help='Path to the input CSV file')
    parser.add_argument('--output_file', help='Path to the output Interval List file')
    args = parser.parse_args()

    convert_csv_to_interval_list(args.metadata_file, args.output_file)
