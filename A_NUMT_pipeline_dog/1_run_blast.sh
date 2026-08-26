#!/bin/bash
set -eu

# Generic NUMT discovery runner (reads assemblies from assemblies.txt and uses config.sh)
script_dir=$(cd $(dirname "$0") && pwd)
cd "$script_dir"

if [[ -f config.sh ]]; then source config.sh; else echo "Please edit config.sh and set REFDIR and other variables."; exit 1; fi

# read assemblies and run BLASTn for each
while read -r ASS || [[ -n "$ASS" ]]; do
  ASS=$(echo "$ASS" | sed 's/#.*//g' | xargs)
  [[ -z "$ASS" ]] && continue
  echo "Processing assembly: $ASS"

  for HAPLOTYPE in "${HAPLOTYPES[@]}"; do
    echo " Haplotype: $HAPLOTYPE"

    # locate nuclear FASTA (allow .fasta/.fa/.fna)
    NUCL=$(ls ${REFDIR}/*${ASS}*${HAPLOTYPE}* 2>/dev/null | head -n1 || true)
    if [[ -z "$NUCL" ]]; then
      NUCL=$(ls ${REFDIR}/*${ASS}*.fasta ${REFDIR}/*${ASS}*.fa ${REFDIR}/*${ASS}*.fna 2>/dev/null | head -n1 || true)
    fi
    if [[ -z "$NUCL" ]]; then
      echo "Cannot find nuclear FASTA for ${ASS} in ${REFDIR}. Skipping."; continue
    fi

    # locate mt FASTA
    MT=$(ls ${REFDIR}/mt_genomes/${ASS}* 2>/dev/null | head -n1 || true)
    if [[ -z "$MT" ]]; then
      echo "Cannot find mt FASTA for ${ASS} in ${REFDIR}/mt_genomes. Provide mt or create doubled mt using utils/double_fasta.py"; continue
    fi

    mkdir -p "$OUTDIR"
    NUMT_TAB=${OUTDIR}/blast.${ASS}.${HAPLOTYPE}.tab
    NUMT_BED=${OUTDIR}/blast.${ASS}.${HAPLOTYPE}.bed

    echo "Running blastn: DB=$NUCL QUERY=$MT -> $NUMT_TAB"
    blastn -db $NUCL -query $MT -outfmt '7' -task blastn -evalue 0.0001 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 > "$NUMT_TAB"

    # convert blast to bed (script expects one arg)
    python3 2_process_blast.py "$NUMT_TAB" | bedtools sort | bedtools merge -c 4,5,6 -o distinct > "$NUMT_BED"

    # flank and extract sequences
    if [[ -f ${NUCL}.fai ]]; then CHR_SIZES=${NUCL}.fai; else CHR_SIZES=""; fi
    bedtools flank -b 500 -i "$NUMT_BED" -g ${NUCL}.fai > ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.bed
    bedtools getfasta -name -fi $NUCL -bed ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.bed -fo ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.fasta

    # concatenate flank sequences (writes concat.<filename>)
    python3 utils/concat_fasta.py ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.fasta

  done

done < "$ASSEMBLIES_FILE"

# Produce numt_set with unique IDs
mkdir -p numt_set
while read -r ASS || [[ -n "$ASS" ]]; do
  ASS=$(echo "$ASS" | sed 's/#.*//g' | xargs)
  [[ -z "$ASS" ]] && continue
  for H in "${HAPLOTYPES[@]}"; do
    if [[ -f ${OUTDIR}/blast.${ASS}.${H}.bed ]]; then
      awk -v prefix="${ASS}_${H}_N" 'BEGIN { OFS = "\t" } { print $0, prefix i++ }' ${OUTDIR}/blast.${ASS}.${H}.bed > numt_set/numts.${ASS}.${H}.bed
    fi
  done
done < "$ASSEMBLIES_FILE"

echo "Pipeline finished. Results in $OUTDIR and numt_set/"