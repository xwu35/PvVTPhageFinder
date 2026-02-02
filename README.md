# PvVTPhageFinder

PvVTPhageFinder is largely similar to [VTPhageFinder](https://github.com/xwu35/VTPhageFinder), with modifications to enable phage identification from viral tagging data where the host genome does not contain any prophage regions.

## Set up environment

### Clone the repository 

```bash
git clone https://github.com/xwu35/PvVTPhageFinder.git
```

### Install Snakemake

PvVTPhageFinder is built for Snakemake version 7. Version 8 and above introduce breaking changes and deprecations and have not been tested. It may not function correctly with newer versions. Please install Snakemake version 7 using the script below.

```bash
cd PvVTPhageFinder

# option 1 using conda
conda env create -n snakemake -f snakemake_env.yml

# option 2 using mamba if it's installed
mamba env create -n snakemake -f snakemake_env.yml
```

### Download snakemake profile

The profile is required to run the workflow on HPC. Skip this step if you already have a SLURM profile in `~/.config/snakemake`.

```bash
# download the profile
git clone https://github.com/xwu35/slurm

# move the profile to the right directory
mv slurm ~/.config/snakemake 
```

## Sample information table

The sample information table should look like this:

| sample   | R1                                    | R2                                    |
|----------|---------------------------------------|---------------------------------------|
| Fp22_10A | Baldridge_10A_Fp22_CD_R1_001.fastq.gz | Baldridge_10A_Fp22_CD_R2_001.fastq.gz | 
| Fp22_3C  | Baldridge_3C_Fp22_HHC_R1_001.fastq.gz | Baldridge_3C_Fp22_HHC_R2_001.fastq.gz | 
| Fp22_7B  | Baldridge_7B_Fp22_CD_R1_001.fastq.gz  | Baldridge_7B_Fp22_CD_R2_001.fastq.gz  | 
| Fp22_8G  | Baldridge_8G_Fp22_CD_R1_001.fastq.gz  | Baldridge_8G_Fp22_CD_R2_001.fastq.gz  | 

## Usage

PvVTPhageFinder supports two mapping software options (bowtie2 and minimap2) and three assembler options (megahit, metaspades and spades_sc), with minimap2 and spades_sc (single cell mode) used by default. Detailed usage information can be viewed using the -h or --help flags `python PvVTPhageFinder.py -h`.

Do not run the analysis on the login node. Submit it as a sbatch job. See `run_pvvtphagefinder.sh` for an example, or check the HTCF usage guide here (https://github.com/xwu35/baldridge_lab/blob/main/HTCF.md). 

A dry-run can be performed to check which rules will be executed and which files will be produced by specifying `--dryrun`.

```bash
conda activate snakemake

python PvVTPhageFinder.py \
    --reads_dir /path/to/dir/containing/sequences \
    --sample_info /path/to/sample_info/file \
    --output_dir output_dir \
    --reference_genome /path/to/reference/genome/fasta/file 
```

### Specific steps

Specific steps can be run using the `--step` flag. 

- **fastqc**: QC on raw reads
- **preprocess**: run fastqc and trim the reads
- **assemble**: all steps (fastqc, preprocess, assemble trimmed reads into contigs and remove contigs aligned to the host genome with >= 95% ANI over 85% AF)

PvVTPhageFinder runs all steps by default.

## Output description

- **Quality control results**: output_dir/reads_processing/fastqc
- **Trimmed reads used for assembly**: output_dir/reads_processing/filtered_reads
- **Read counts**: output_dir/reads_processing/reads_statistics/{number_of_reads_removed_at_each_step.txt,reads_composition_barplot.svg}
- **Retained contig sequences**: output_dir/check_contig_contamination/no_host_contig_sequences


