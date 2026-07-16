#!/bin/bash

set -eu

# Add assembly name to chromsosome. Join dicts.
for ASS in CHM13 mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do  \
	cat dict.v1_to_v2.${ASS}.bed | \
	awk -v ass=$ASS '{OFS="\t"}{print ass":"$0}'\
		done | bedtools sort > joint.dict.v1_to_v2.bed

