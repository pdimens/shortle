# shortle
A Pretty Ok Short-Read Assembly Workflow


## How to use me
1. You'll need some reads
- Preferably paired-end
- Preferably one file each for foward and reverse reads

2. You'll need conda/mamba (recommended: mamba) and to install the main conda environment, provided by `shortle.yaml` in the main directory of this repository. The remaining software is provided in `workflow/envs/redundans.yaml`, which will be used by Snakemake, so you don't need to worry about it.
```bash
conda env create -f shortle.yaml
```
This command will create an environment called `shortle`. Activate it with:
```bash
conda activate shortle
```

3. You'll need to modify the first few lines that specify the forward/reverse read information in `Snakefile` to match your data

4. You'll need to run it once with the specific command to end it prematurely (at the purge_haplotigs stage):
```bash
snakemake -j 20 purge_haplotigs/mapped2.bam.gencov
```
where `j` is how many cores you are willing to reserve for the Snakemake workflow.

5. Inspect `purge_haplotigs/mapped2.bam.histogram.png` to figure out your cutoffs as described [here](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/): ![purge haplotigs histogram](https://bitbucket.org/repo/Ej8Mz7/images/1039246939-coverage_histogram.png)

6. Modify `rule purge_haplotigs_suspects` to match your chosen cutoffs and rerun the entire workflow:
```bash
snakemake -j 20
```

7. Your resulting assembly will be in `polish_2`
