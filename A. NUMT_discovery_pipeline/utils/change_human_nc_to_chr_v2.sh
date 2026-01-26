#!/bin/bash

set -ue

# Change human NC chromosome IDs into human-readable chr numbers.

input_file=$1

# Check if the input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Extract the file name without extension and the file extension separately
base_name="${input_file%.*}"
extension="${input_file##*.}"

# Create an inline mapping in awk
awk '
BEGIN {
    # Define the mapping
    map["NC_060925.1"] = "chr1";
    map["NC_060926.1"] = "chr2";
    map["NC_060927.1"] = "chr3";
    map["NC_060928.1"] = "chr4";
    map["NC_060929.1"] = "chr5";
    map["NC_060930.1"] = "chr6";
    map["NC_060931.1"] = "chr7";
    map["NC_060932.1"] = "chr8";
    map["NC_060933.1"] = "chr9";
    map["NC_060934.1"] = "chr10";
    map["NC_060935.1"] = "chr11";
    map["NC_060936.1"] = "chr12";
    map["NC_060937.1"] = "chr13";
    map["NC_060938.1"] = "chr14";
    map["NC_060939.1"] = "chr15";
    map["NC_060940.1"] = "chr16";
    map["NC_060941.1"] = "chr17";
    map["NC_060942.1"] = "chr18";
    map["NC_060943.1"] = "chr19";
    map["NC_060944.1"] = "chr20";
    map["NC_060945.1"] = "chr21";
    map["NC_060946.1"] = "chr22";
    map["NC_060947.1"] = "chrX";
    map["NC_060948.1"] = "chrY";
    map["NC_000001.11"] = "chr1";
    map["NC_000002.12"] = "chr2";
    map["NC_000003.12"] = "chr3";
    map["NC_000004.12"] = "chr4";
    map["NC_000005.10"] = "chr5";
    map["NC_000006.12"] = "chr6";
    map["NC_000007.14"] = "chr7";
    map["NC_000008.11"] = "chr8";
    map["NC_000009.12"] = "chr9";
    map["NC_000010.11"] = "chr10";
    map["NC_000011.10"] = "chr11";
    map["NC_000012.12"] = "chr12";
    map["NC_000013.11"] = "chr13";
    map["NC_000014.9"] = "chr14";
    map["NC_000015.10"] = "chr15";
    map["NC_000016.10"] = "chr16";
    map["NC_000017.11"] = "chr17";
    map["NC_000018.10"] = "chr18";
    map["NC_000019.10"] = "chr19";
    map["NC_000020.11"] = "chr20";
    map["NC_000021.9"] = "chr21";
    map["NC_000022.11"] = "chr22";
    map["NC_000023.11"] = "chrX";
    map["NC_000024.10"] = "chrY";
    map["NC_012920.1"] = "chrM";
    map["chrUn_NT"] = "chrUn_NT";
    map["chrUn_NW"] = "chrUn_NW";
}
{
    OFS="\t";
    for (i=1; i<=NF; i++) {
        for (key in map) {
            if ($i ~ key) {
                sub(key, map[key], $i);
            }
        }
    }
    print;
}' "$input_file" > "${base_name}.noNC.${extension}"

