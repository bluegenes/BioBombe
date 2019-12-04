# Use snakemake to run BioBombe!

I've added a snakemake version of BioBombe's sequential compression. 

I wanted to try [BioBombe](https://github.com/greenelab/BioBombe) out on transcriptome data sequenced as part of the Marine Microbial Eukaryotic Sequencing Project (MMETSP). That work is in progress over in [2019-burgers-shrooms](https://github.com/bluegenes/2019-burgers-shrooms).

Underlying scripts were (minorly) modified from the versions found in the [BioBombe repo](https://github.com/greenelab/BioBombe) (G. Way and C. Greene).

BioBombe Preprint:[Sequential compression across latent space dimensions enhances gene expression signatures Way, G.P., Zietz, M., Himmelstein, D.S., Greene, C.S. biorXiv preprint (2019) doi:10.1101/573782](https://www.biorxiv.org/content/10.1101/573782v2)


## To run:

These scripts rely on snakemake workflow software and conda package management, and require conda (e.g. miniconda), snakemake and python>=3.6 to be installed.

Use the `--use-conda` flag to allow snakemake to install the `biobombe` environment (from the `environment.yml` file) for you and use it for each of the `biobombe` steps.
These steps take _a while_, so I've left the `--dryrun` on the end of the following commands. This will check that the workflow can be executed and show you which 
files will be generated. Remove this flag to actually run these steps. Use `--cores 4` to run on, say, 4 cores, and check out [this info](https://hackmd.io/K0FWjvlYQbCQ1gi-llgpPg) 
for running snakemake on a cluster job scheduling system (e.g. slurm).

Each snakefile can be run entirely independently.

Clone this repo:
```
git clone https://github.com/bluegenes/BioBombe
```


Run the initial ksweep:

```
cd 1.initial-k-sweep
snakemake -s biobombe_initial_ksweep.snakefile --use-conda --configfile ksweep_config.yaml --dryrun
```

Run sequential compression:

```
cd 2.sequential-compression
snakemake -s biobombe_sequential_compression.snakefile --use-conda --configfile sequential_config.yaml --dryrun
```

Modify the `ksweep_config.yaml` and `sequential_config.yaml` files to run this on your favorite dataset.

Modify the parameter sweep files found in the `config` folder to change the parameters swept over. If using tsv, watch out for tabs vs spaces, they can trip pandas up.

Contact: Tessa Pierce, @bluegenes on GitHub.
