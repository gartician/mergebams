# mergebam

Use these scripts to merge multiple bam files and then generate bigwig tracks; this is appropriate to merge replicate samples into one file. Please use the `mergebams_srun.sh` to submit jobs to a SLURM system for optimal performance, otherwise use the `mergebams.sh`. 

# Input

Edit the tab-separated file `mergebam.config` to specify the merge.

```
bam merged_bam
/path/to/control1.bam   control_merged.bam
/path/to/control2.bam   control_merged.bam
/path/to/control3.bam   control_merged.bam
/path/to/treated1.bam   treated_merged.bam
/path/to/treated2.bam   treated_merged.bam
/path/to/treated3.bam   treated_merged.bam
```

And then specify which flags to use. The following can be brought up `bash mergebams_srun.sh`:

```
usage: mergebams.sh [-h] [-c CONFIG] [-p CORES] [-o OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        BAM file merging scheme
  -p CORES, --cores CORES
                        Number of cores to use (default = 4)
  -o OUTDIR, --outdir OUTDIR
                        Automatically creates and export merged bams + bigwigs
                        to this directory (default is current working
                        directory). Do not end this option with a '/'
```

Please don't specify the output directory in the config file itself. Instead, use the `-o` flag!

Example command:

```
bash mergebams_srun.sh -c mergebams.config -p 20 -o my_output
```

# Output

The output of this script is replicate-merged bam file and the bigwig file. The exact commands used for these transformations are logged into the `mergebam.log` file.

```
.
├── mergebam.log
├── mergebams.config
├── mergebams.sh
└── my_output
    ├── control_merged.bam
    ├── control_merged.bam.bai
    ├── control_merged.bw
    ├── treated_merged.bam
    ├── treated_merged.bam.bai
    └── treated_merged.bw
```

# Method

1. Input config file associates the target file (merged_bam column) with the starting files (bam column) with a dictionary
2. BAM files are merged in parallel locally or on a SLURM system using samtools merge. BAM files are then indexed once all processes are finished.
3. When BAM files are made, bamCoverage with the `--normalizeUsing CPM`, `--binSize 10`, `--smoothLength 50` flags are used to generate the bigwig files.
4. For all the steps, the exact commands are logged into the `mergebams.log` output file.

# Notes

Make sure you have samtools and deeptools in your PATH / same conda environment. 

Edit the actual script to add custom settings at each step.