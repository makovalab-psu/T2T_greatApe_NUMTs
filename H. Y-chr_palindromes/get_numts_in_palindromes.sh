#!/bin/bash

set -eu

# Get sex chromosome NUMTs.
for ASS in mPanPan1 mPanTro3 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1 CHM13 ; do 
 cat ../numts_from_laptop_blast/numt_set/numts.${ASS}.pri.bed | grep "^chr[XY]" | sed 's/_hsaY//g' | sed 's/_hsaX//g' | sed 's/_mat//g' | sed 's/_pat//g' | sed 's/_hap1//g' | sed 's/_hap2//g' > sex_numts/${ASS}.sex.numts.bed 
done

# Fix coordinate order in palindrome files.
awk '{if ($2 > $3) {tmp=$2; $2=$3; $3=tmp} print}' OFS="\t" Palindromes_apes_XY_v2_lifted/archive/mPanTro3_palindromes_v1.1_to_v2.0_lifted.bed > Palindromes_apes_XY_v2_lifted/mPanTro3_palindromes_v1.1_to_v2.0_lifted.fixed.bed 
awk '{if ($2 > $3) {tmp=$2; $2=$3; $3=tmp} print}' OFS="\t" Palindromes_apes_XY_v2_lifted/archive/mPanPan1_palindromes_v1.1_to_v2.0_lifted.bed > Palindromes_apes_XY_v2_lifted/mPanPan1_palindromes_v1.1_to_v2.0_lifted.fixed.bed  

# Get NUMT overlaps with palindromes by species.
mkdir -p overlaps
for ASS in mPanPan1 mPanTro3 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1 CHM13 ; do bedtools intersect -wo -a sex_numts/${ASS}*.bed -b Palindromes_apes_XY_v2_lifted/${ASS}*.bed > overlaps/${ASS}.overlaps.bed ; done

# Merge overlapping palindrome overlaps (same NUMT, multiple palindrome arms).
for ASS in mPanPan1 mPanTro3 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do bedtools merge -i overlaps/${ASS}.overlaps.bed -c 12 -o distinct > overlaps/${ASS}.merged.overlaps.bed ; done

# Keep what does not overlap with palindromes.
for ASS in mPanPan1 mPanTro3 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1 CHM13 ; do \
bedtools intersect -b Palindromes_apes_XY_v2_lifted/${ASS}*.bed -a sex_numts/${ASS}.sex.numts.bed -loj | grep 'chrY' | grep '\.\t\-1\t\-1\t\.' > overlaps/${ASS}.nonOverlaps.bed
done


# Merge overlapping palindrome overlaps (same NUMT, multiple palindrome arms).
for ASS in mPanPan1 mPanTro3 mGorGor1 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do bedtools merge -i overlaps/${ASS}.nonOverlaps.bed  > overlaps/${ASS}.merged.nonOverlaps.bed ; done

