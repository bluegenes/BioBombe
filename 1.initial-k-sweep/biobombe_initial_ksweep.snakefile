"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda --configfile ksweep_config.yaml
"""

import os
import pandas as pd
from scripts.biobombe_snakemake_utils import read_params
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

# read in adage params
adage_paramsfile=config.get('adage_paramsfile', "config/initial_z_parameter_sweep_adage_MMETSP.tsv")
adage_params = read_params(adage_paramsfile)
adage_params['sweep_values'] = adage_params['sweep_values'].str.split(',')
adage_paramsD = adage_params.to_dict()

# read in tybalt params
tybalt_paramsfile=config.get('tybalt_paramsfile', "config/initial_z_parameter_sweep_tybalt_MMETSP.tsv")
tybalt_params = read_params(tybalt_paramsfile)
tybalt_params['sweep_values'] = tybalt_params['sweep_values'].str.split(',')
tybalt_paramsD = tybalt_params.to_dict()

SAMPLE = config.get("sample_basename", "haptophyta_orthogroup")
data_link = config.get("quant_data_httplink", "https://osf.io/ek9nu/download")

rule all:
    input: 
        f"figures/viz_results/{SAMPLE}_z_parameter_final_loss_adage.png", f"figures/viz_results/{SAMPLE}_z_parameter_final_loss_tybalt.png",
        expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= SAMPLE, learning_rate = tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa = tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']), 
        expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= SAMPLE, learning_rate = adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity = adage_paramsD['sweep_values']['sparsity'], noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
        
        # use this instead for specific param combos version
        #"figures/viz_results/{sample}_z_parameter_final_loss_dae.png", "figures/viz_results/{sample}_z_parameter_final_loss_vae.png"

rule download_data:
    input: HTTP.remote(data_link)
    output: "data/{sample}.quant.tsv" 
    log: "logs/download_{sample}.log"
    shell: "mv {input} {output} 2> {log}"

rule preprocess_data:
    input:
        "data/{sample}.quant.tsv"
    output:
        processed = "data/{sample}.processed.tsv.gz", # only needed file
        #train = "data/{sample}.train.processed.tsv.gz",
        #test = "data/{sample}.test.processed.tsv.gz",
        #mad = "data/{sample}.mad.processed.tsv.gz",
        #mad_test = "data/{sample}.mad.test10.processed.tsv.gz",
        #mad_train = "data/{sample}.mad.train90.processed.tsv.gz",
    params:
        outdir = "data"
    conda:
        "../environment.yml"
    shell:
        """
        python scripts/process_expression_data.py {input} --mad --output_folder {params.outdir}
        """

# not needed - scaling happens in train_models script
rule preprocess_data_scale:
    input:
        "data/{sample}.quant.tsv"
    output:
        train = "data/{sample}.train.processed.zeroone.tsv.gz",
        test = "data/{sample}.test.processed.zeroone.tsv.gz",
        mad = "data/{sample}.mad.processed.zeroone.tsv.gz"
    params:
        outdir = "data"
    conda:
        "../environment.yml"
    shell:
        """
        python scripts/process_expression_data.py {input} --mad --output_folder {params.outdir} --scale --scale_method "min_max"
        """


rule run_adage:
    input: 
        expand("data/{sample}.processed.tsv.gz", sample = SAMPLE)
    output: 
        "results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv",
    conda: '../environment.yml'
    shell:
       """
        export KERAS_BACKEND=tensorflow
       python scripts/adage.py  --input_data {input} --learning_rate {wildcards.learning_rate} --batch_size {wildcards.batch_size} --epochs {wildcards.epochs} --sparsity {wildcards.sparsity} --noise {wildcards.noise} --output_filename {output} --num_components {wildcards.num_components} --subset_mad_genes 8000 --scale
        """

rule summarize_paramsweep_adage:
    input:
        all_adage_runs = expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= SAMPLE, learning_rate = adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity = adage_paramsD['sweep_values']['sparsity'], noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
    output: "results/adage/paramsweep_summary.txt"
    params: 
        results_dir = directory("results/adage")
    shell:
        """ 
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule run_tybalt:
    input: 
        expand("data/{sample}.processed.tsv.gz", sample = SAMPLE)
    output:
        "results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv"
    conda: '../environment.yml'
    shell:
       """
        export KERAS_BACKEND=tensorflow
       python scripts/vae.py    --input_data {input} --learning_rate {wildcards.learning_rate} --batch_size {wildcards.batch_size} --epochs {wildcards.epochs} --kappa {wildcards.kappa} --output_filename {output} --num_components {wildcards.num_components} --subset_mad_genes 8000 --scale
        """

rule summarize_paramsweep_tybalt:
    input: 
        all_tybalt_runs = expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= SAMPLE, learning_rate = tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa = tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']),
    params:
        results_dir = directory("results/tybalt")
    output: "results/tybalt/paramsweep_summary.txt"
    shell:
        """ 
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """


rule visualize_paramsweep:
    input: 
        adage = "results/adage/paramsweep_summary.txt", 
        tybalt = "results/tybalt/paramsweep_summary.txt",
    params:
        dataset_name = SAMPLE,
    output: 
        vae_loss_png="figures/viz_results/{sample}_z_parameter_final_loss_tybalt.png",
        vae_training_png="figures/viz_results/{sample}_z_parameter_training_tybalt.png",
        vae_best_params="results/{sample}_best_params_tybalt.tsv",
        vae_best_model_png="results/{sample}_best_model_tybalt.png",

        dae_loss_png="figures/viz_results/{sample}_z_parameter_final_loss_adage.png",
        dae_best_params="results/{sample}_best_params_adage.tsv",
        dae_best_model_png="figures/viz_results/{sample}_bestmodel_adage.png"
    log:
        "logs/{sample}_visualize_paramsweep_adage_tybalt.log"
    params:
        dataset_name = SAMPLE,
    conda:
        "../environment.yml"
    script:
        "scripts/visualize-parameter-sweep.R"


### code below NOT being used right now, but *can* be used to run individual dae/vae with specified parameter combos


def generate_dae_combos(w):
    dae_combos = []
    for component, zparams in zsweep_paramsD.items():
        learning_rate = zparams["dae_lr"]
        batch_size = zparams["dae_batch_size"]
        epochs = zparams["dae_epochs"]
        sparsity = zparams["dae_sparsity"]
        noise = zparams["dae_noise"]
        dae_combos.append(f"results/dae/{w.sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}.tsv")
    return dae_combos


def generate_vae_combos(w):
    vae_combos = []
    for component, zparams in zsweep_paramsD.items():
        learning_rate = zparams["vae_lr"]
        batch_size = zparams["vae_batch_size"]
        epochs = zparams["vae_epochs"]
        kappa = zparams["vae_kappa"]
        vae_combos.append(f"results/vae/{w.sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv")
    return vae_combos

rule run_dae:
    input: 
        "data/{sample}.scaled.tsv"
    output: 
        "results/dae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}.tsv"
    conda: '../environment.yml'
    shell:
       """
       python scripts/dae.py --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --sparsity {wildcards.sparsity}
                        --noise {wildcards.noise}
                        --output_filename {output}
                        --subset_mad_genes
        """

rule summarize_paramsweep_dae:
    input: generate_dae_combos
    output: "results/dae/{sample}_paramsweep_summary.txt"
    params: 
        results_dir = directory("results/dae")
    shell:
        """
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule run_vae:
    input:
        "data/{sample}.scaled.tsv"
    output:
        "results/vae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv"
    conda: '../environment.yml'
    shell:
       """
       python scripts/vae.py    --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --kappa {wildcards.kappa}
                        --output_filename {output}
                        --subset_mad_genes
        """

rule summarize_paramsweep_vae:
    input: generate_vae_combos
    params:
        results_dir = directory("results/vae")
    output:
        "results/vae/{sample}_paramsweep_summary.txt"
    shell:
        """
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule visualize_dae_vae_paramsweep:
    input:
        dae = "results/dae/{sample}_paramsweep_summary.txt",
        vae = "results/vae/{sample}_paramsweep_summary.txt",
    output:
        vae_loss_png="figures/viz_results/{sample}_z_parameter_final_loss_vae.png",
        vae_training_png="figures/viz_results/{sample}_z_parameter_training_vae.png",
        vae_best_params="results/{sample}_best_params_vae.tsv",

        dae_loss_png="figures/viz_results/{sample}_z_parameter_final_loss_dae.png",
        dae_best_params="results/{sample}_best_params_dae.tsv",
        dae_best_model_png="figures/viz_results/{sample}_bestmodel_dae.png"
    log:
        "logs/{sample}_visualize_dae_vae_paramsweep.log"
    params:
        dataset_name = SAMPLE, 
    script:
        "scripts/visualize-parameter-sweep.R"

