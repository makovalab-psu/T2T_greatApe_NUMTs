#!/bin/bash

for SAMPLE in 'hprc5'; do

    ###############
    ## PAUP route ##
    # Create the Nexus script
    {
    echo "#NEXUS"
    echo ""
    echo "begin paup;"
    echo "    execute ${SAMPLE}.species_matrix_flanks500_0.5.nexus;"
    # echo "    outgroup 5-6;"
    echo "    hsearch;"
    echo "    showTrees / outRoot=monophyletic;"
    echo "    saveTrees format=Newick brLens=yes root=yes trees=all file=pauptree_${SAMPLE}.nwk;"
    echo ""
    echo "end;"
    } > temp_script.nexus

    # Run the PAUP script with the generated Nexus file
    paup_tool/paup4a168_osx -n < temp_script.nexus

done

##########
## Plot ##

# Then plot using web-based tool IcyTree or MEGA11:
# https://icytree.org

