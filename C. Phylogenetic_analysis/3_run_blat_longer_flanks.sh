#!/bin/bash
set -eu
echo "Requires BLAT v36"
mkdir -p log

# FLANK_SIZE=1000
# MIN_SCORE=1300

FLANK_SIZE=750
MIN_SCORE=975 #Flank_size * 1.3

# # Loop through each assembly as both query and database
# for Q in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do 
#   for DB in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do

# Loop through each assembly as both query and database
for Q in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do 
  for DB in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do

    # Define paths
    QUERY="concat.flanks${FLANK_SIZE}.${Q}.fasta"
    DATABASE="/nfs/brubeck.bx.psu.edu/scratch4/makova_lab/downloads/primate_T2T/assemblies/v2/${DB}.pri.cur.*[0-9].fasta"
    if [[ $DB == "CHM13" ]]; then DATABASE=/nfs/brubeck.bx.psu.edu/scratch6/makova_lab/downloads/Reference_Genomes/human/CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna; fi
    OUTPUT="out.${Q}.to.${DB}.flanks${FLANK_SIZE}.psl"

    # Generate unique job script name
    JOB_SCRIPT="../log/temp_${Q}_to_${DB}.flanks${FLANK_SIZE}.job"

    # Create the job script
    cat <<EOF > $JOB_SCRIPT
#!/bin/bash
#SBATCH --job-name=blat_${Q}_to_${DB}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --time=0-80:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ejt89@psu.edu
#SBATCH --chdir=/nfs/brubeck.bx.psu.edu/scratch4/edmundo/numts_t2t/bin/konkel_loftus/longer_flanks
#SBATCH --output=/nfs/brubeck.bx.psu.edu/scratch4/edmundo/numts_t2t/bin/konkel_loftus/log/ejt89-blat%j.out
#SBATCH --error=/nfs/brubeck.bx.psu.edu/scratch4/edmundo/numts_t2t/bin/konkel_loftus/log/ejt89-blat%j.err

# Conda environment
export CONDA_ENVS_PATH=/galaxy/home/ejt89/.conda/envs
source /galaxy/home/ejt89/anaconda3/etc/profile.d/conda.sh
conda activate /galaxy/home/ejt89/.conda/envs/blat
conda info --envs
conda list

# Run BLAT
time blat -minScore=${MIN_SCORE} $DATABASE $QUERY $OUTPUT
EOF

    # Submit the job
    sbatch $JOB_SCRIPT

  done
done