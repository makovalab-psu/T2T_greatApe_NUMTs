#!/bin/bash

set -uex

# Variables for file paths
MEGACC_PATH="/Applications/MEGA11.app/Contents/MacOS/megacc"
NEWICK_FILE="tree_flanks500_0.8.nwk"
TEMPLATE_FILE="/Users/edmundo/Documents/GitHub/numts_primate_t2t/konkel/longer_flanks/tree_plot_template.mao"  # Use the full path to your MEGA .mao template file
OUTPUT_IMAGE="tree_plot.png"

# Check if the .mao file exists
if [[ -f "$TEMPLATE_FILE" ]]; then
    echo "$TEMPLATE_FILE exists."
else
    echo "$TEMPLATE_FILE does not exist! Exiting."
    exit 1
fi

# Run MEGA's Command Line Interface to plot the tree
$MEGACC_PATH -t "$TEMPLATE_FILE" -d "$NEWICK_FILE" -o "$OUTPUT_IMAGE"

# Check if the output image was successfully generated
if [[ -f "$OUTPUT_IMAGE" ]]; then
    echo "Tree plotted and saved as $OUTPUT_IMAGE"
else
    echo "Failed to generate the tree plot."
fi
