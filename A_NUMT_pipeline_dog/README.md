NUMT discovery pipeline (dog-ready)

This folder contains a standalone, repurposable NUMT discovery pipeline adapted from the original project. It’s minimal and generic so it can be pointed at a dog reference.

Quick steps:
1. Edit config.sh to point REFDIR to your reference FASTA directory and set ASSEMBLY name(s).
2. Provide mt_sizes.tsv with mitochondrial genome sizes for the assemblies (or set MT_SIZE env var).
3. Run: bash 1_run_blast.sh

Files:
- config.sh  : configuration (REFDIR, ASSEMBLIES_FILE, HAPLOTYPES)
- assemblies.txt : one assembly name per line (example: canFam3)
- mt_sizes.tsv : tab-separated assembly<TAB>mt_size (example provided)
- 1_run_blast.sh : runs BLASTn, produces BED/flanks
- 2_process_blast.py : converts BLAST output to BED (uses mt_sizes.tsv)
- 3_process_bed.py : post-process BEDs
- utils/ : helper scripts (concat_fasta.py, double_fasta.py, bed_to_bigBed.sh)

Notes:
- Update REFDIR in config.sh to where your assembly FASTA and mt FASTA live.
- The pipeline expects a doubled mt FASTA for correct breakpoint handling; double_fasta.py is provided.
- This is a first-pass genericization; confirm paths and parameters before large runs.