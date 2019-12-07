# snakemake execution
# uncomment the last half of each command and modify as necessary to use with a slurm cluster

conda activate snakemake # conda environment with py>=3.6, snakemake


cd 1.initial-k-sweep
snakemake -s biobombe_initial_ksweep.snakefile --configfile ksweep_config.yaml --conda-prefix ../  --use-conda -k --rerun-incomplete --jobs 5 #--cluster "sbatch -t 1:00:00 -N 1 -n 1 -c 1 -p bmm --mem=15gb"
cd -

cd 2.sequential-compression
snakemake -s biobombe_sequential_compression.snakefile --configfile sequential_config.yaml --jobs 5 -k --rerun-incomplete --conda-prefix ../  --use-conda #--cluster "sbatch -t 1:00:00 -N 1 -n 1 -c 1 -p bmm --mem=15gb" 
cd -

