#!/bin/bash

set -eu


out_path="/scratch/ejt89/apes_output/alignments"


for fq_path in $(cat list_paths_fq.txt | tail -n1) ; do

bai_prefix=$( echo $fq_path | cut -d'/' -f5-8 | sed 's/\//_/g' | cut -d'_' -f1-4 )
bai_path="${out_path}/${bai_prefix}.bam.bai"

echo "#Started job for $bai_path"

cat <<EOF > temp.job
#!/bin/bash
#SBATCH --job-name=single_align_snakemake
#SBATCH --account=kdm16_sc_default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm/single_%j.out
#SBATCH --error=logs/slurm/single_%j.err

# 1. Initialize micromamba (using your proven settings)
export MAMBA_ROOT_PREFIX="/storage/home/ejt89"
eval "\$(\$MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"

# 2. Activate the bioinfo environment
micromamba activate bioinfo

# 3. Clear any stale locks from previous crashes.
snakemake --unlock

# 4. Run 1 sample at a time, each using 4 threads.
snakemake --cores 4 --latency-wait 60 ${bai_path}

EOF

# Submit job.
sbatch temp.job

# Wait for scheduler to breathe.
sleep 0.5

done

