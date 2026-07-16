#!/bin/bash

# Check if the flank size is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <flank_size>"
    exit 1
fi

# Assign the flank size from the argument
FLANKS=$1

# Define filename suffixes
FILENAME_SUFFIXES=("0.2" "0.5" "0.8")

##############
## Tree-QMC ##
# Use Tree-QMC to compute the tree (from Erin Molloy's group).
for suffix in "${FILENAME_SUFFIXES[@]}"; do
 rm qmctree_flanks${FLANKS}_${suffix}.nwk
 tree-qmc -i results_${FLANKS}bp/species_matrix_flanks${FLANKS}_${suffix}.nexus --bp --root B._orangutan --support -o qmctree_flanks${FLANKS}_${suffix}.nwk
done

###############
## PAUP route ##
# Create the Nexus script
{
echo "#NEXUS"
echo ""
echo "begin paup;"
for suffix in "${FILENAME_SUFFIXES[@]}"; do
    echo "    execute results_${FLANKS}bp/species_matrix_flanks${FLANKS}_${suffix}.nexus;"
    echo "    outgroup 5-6;"
    echo "    hsearch;"
    echo "    showTrees / outRoot=monophyletic;"
    echo "    saveTrees format=Newick brLens=yes root=yes trees=all file=pauptree_flanks${FLANKS}_${suffix}.nwk;"
    echo ""
done
echo "end;"
} > temp_script.nexus

# Run the PAUP script with the generated Nexus file
paup_tool/paup4a168_osx -n < temp_script.nexus


##########
## Plot ##

# Then plot using web-based tool IcyTree or MEGA11:
# https://icytree.org

