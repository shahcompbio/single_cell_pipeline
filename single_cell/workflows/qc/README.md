## pseudobulkQC

Snakemake pipeline to generate QC plots for DLP pseudobulk pipeline. Pipeline generates plots with [scgenome](https://github.com/shahcompbio/scgenome/tree/master/scgenome) and then compiles a report using Rmarkdown.

### Input
To run you will need a csv file that looks like `metadata/fitness_pseudobulk_test.csv`. In particular, it will need the following columns `datatag`, `sample_id`, `jira_id`, `library_id`. `jira_id` is the pseudobulk ticket. (If some of this is superfluous or you have a different csv file with different columns you might want to tweak the top of the `Snakefile`).

Then change the `config.yaml` file to point to the appropriate file.

### Running on juno

To run on juno you can do the following (from the directory containing the `Snakefile`). log files will get written to the `logs` directory. You will need to launch from within a python environment that has scgenome and snakemake installed. The html compilation runs through a singularity image that snakemake will download automatically the first time you run.
```
mkdir logs
module load singularity

CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

vepdata="/work/shah/reference/vep/"
genomeref="/work/shah/reference/genomes/GRCh37-lite/"
tantalusdir="/work/shah/tantalus/"

snakemake --jobs 50 \
  --use-singularity \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --singularity-args "--bind ${vepdata} --bind ${genomeref} --bind ${tantalusdir}" --keep-going
```

### Output
If everything works you should see a html file for each sample in `reports/`.
