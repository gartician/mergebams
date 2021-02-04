import os
import sys
import pandas as pd
import glob
from collections import defaultdict
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help = "BAM file merging scheme")
parser.add_argument("-p", "--cores", default = 4, help = "Number of cores to use (default = 4)")
parser.add_argument("-o", "--outdir", default = ".", help = "Automatically creates and export merged bams + bigwigs to this directory (default is current working directory). Do not end this option with a '/'")
args = parser.parse_args()

if not os.path.isdir(args.outdir):
	os.mkdir(args.outdir)

# read in config file and associate target file with starting files -------------------------------
config = pd.read_csv(args.config, sep = "\t", header = 0)

merge_dict = defaultdict(list)

for i, j in zip(config.bam, config.merged_bam):
	merge_dict[j].append(i)

print("Target file\tStarting files")
for target in merge_dict:
	print(target, "\t", " ".join(merge_dict[target]))

# build and execute commands to merge bams --------------------------------------------------------
print("\n")
print("Merging then indexing bam files...")
log_file = open("mergebams.log", "w")
cmd_list = [] # list of commands. each element is space-delimited list of actual cmd to run.
for target in merge_dict:
	in_files = " ".join(merge_dict[target])
	out_file = "{}/{}".format(args.outdir, target)
	cmd = "samtools merge -f -@ {p} - {i} | samtools -@ {p} -o {o}".format(p = args.cores, o = out_file, i = in_files)
	print("Making target file: ", out_file)
	print(cmd, file = log_file)
	cmd_list.append(cmd.split(" "))

# run all process and wait for them all to finish.
procs = [ subprocess.Popen(i) for i in cmd_list ]
for p in procs:
	p.wait()

log_file = open("mergebams.log", "a")
# index merged bam files serially.
print("BAM files merged, now indexing outputs")
for target in merge_dict:
	in_file = "{}/{}".format(args.outdir, target)
	cmd = "samtools index -@ {p} {i}".format(p = args.cores, i = in_file)
	print(cmd, file = log_file)
	os.system(cmd)

print("BAM files index, now making bigwigs")

# bams to bigwigs ---------------------------------------------------------------------------------
cmd_list = []
for target in merge_dict:
	in_file = "{}/{}".format(args.outdir, target)
	out_file = "{}/{}".format(args.outdir, target.replace(".bam", ".bw", 1))
	cmd = "bamCoverage -b {b} -o {o} -p {p} --normalizeUsing CPM --binSize 10 --smoothLength 50 -v".format(b = in_file, o = out_file, p = args.cores)
	print("Making target bigwig: ", out_file)
	print(cmd, file = log_file)
	cmd_list.append(cmd.split(" "))

procs = [ subprocess.Popen(i) for i in cmd_list ]
for p in procs:
	p.wait()

log_file.close()
print("Bigwig files made. Program is finished")