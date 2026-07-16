import numpy as np
import pandas as pd
import gzip
import sys
import os

# Check for input argument
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <input_file.gz>")
    sys.exit(1)

input_file = sys.argv[1]

# Load only the first two columns from gzipped data
with gzip.open(input_file, "rt") as f:
    df = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1], names=["Index", "NUMT"])

# Parameters
num_bootstraps = 1000
bootstrap_means = []

# Bootstrap resampling
for _ in range(num_bootstraps):
    sample = df.sample(frac=1, replace=True)
    bootstrap_means.append(sample["NUMT"].mean())

# Generate output filename based on input file name
output_file = os.path.splitext(input_file)[0] + "_bootstrap_frequencies.txt"

# Save results
pd.DataFrame(bootstrap_means, columns=["Bootstrap_Frequency"]).to_csv(output_file, index=False)

print(f"Bootstrap results saved to {output_file}")
