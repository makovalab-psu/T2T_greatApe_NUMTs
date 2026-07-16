#!/bin/bash
set -eu

matrix="hprc5.species_matrix_flanks500_0.5.tsv"

# Print columns
echo ""
head -n1 "${matrix}" | cut -f2-12

# Grouping frequencies
# Only found in human T2T, show HPRC5 groupings.
echo ""
echo "# NUMT only found in human T2T, show HPRC5 groupings."
grep -E $'.*_numt[0-9]+\t0\t0\t1\t0\t0\t0' "${matrix}" | cut -f2-12 | sort | uniq -c | sort -hr

# Found in all T2T except human T2T, show HPRC5 groupings.
echo ""
echo "# NUMT found in all T2T except human T2T, show HPRC5 groupings."
grep -E $'.*_numt[0-9]+\t1\t1\t0\t1\t1\t1' "${matrix}" | cut -f2-12 | sort | uniq -c | sort -hr

# Found in all T2T except human T2T, show HPRC5 groupings.
echo ""
echo "# NUMT found in all T2T, show HPRC5 groupings."
grep -E $'.*_numt[0-9]+\t1\t1\t1\t1\t1\t1' "${matrix}" | cut -f2-12 | sort | uniq -c | sort -hr

# Found in CBHG but not in orangutans, show HPRC5 groupings.
echo ""
echo "# NUMT found in all T2T, show HPRC5 groupings."
grep -E $'.*_numt[0-9]+\t1\t1\t1\t1\t0\t0' "${matrix}" | cut -f2-12 | sort | uniq -c | sort -hr