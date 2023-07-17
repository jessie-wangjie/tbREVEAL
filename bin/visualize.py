import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def create_heatmap(input_file, output_file):
    # Read the input CSV file into a pandas DataFrame
    data = pd.read_csv(input_file, delimiter=',')
    data = data[~data['Target'].str.startswith('AA')]
    data = data.loc[~(data.iloc[:, 1:] == 0).all(axis=1)]


    # Extract the target column and remove it from the data
    targets = data['Target']
    data = data.drop('Target', axis=1)
    data = data.drop('median', axis=1)

    data = data.sort_values(by='averages', ascending=False)
    print(len(data['averages']))
    data = data.drop('averages', axis=1)
    targets = targets.loc[data.index]

    # Create a heatmap using seaborn
    plt.figure(figsize=(10, 12))  # Set the figure size
    sns.heatmap(data, cmap='YlGnBu', fmt='.2f', cbar=True)

    # Set the axis labels and title
    plt.xlabel('Replicates')
    plt.ylabel('Cryptic Site')  # Update the y-axis label
    plt.title('Heatmap of Data')

    # Rotate the y-axis labels for better visibility if needed
    # plt.yticks(range(len(data)), targets)

    plt.yticks(np.arange(len(data)) + 0.5, targets)

    # Save the plot to a PNG file
    plt.savefig(output_file)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Create a heatmap from a CSV file and save it to a PNG file.')
    parser.add_argument('--input_file', type=str, help='Path to the input CSV file')
    parser.add_argument('--output_file', type=str, help='Path to the output PNG file')
    args = parser.parse_args()

    # Call the create_heatmap function with the input and output files
    create_heatmap(args.input_file, args.output_file)
