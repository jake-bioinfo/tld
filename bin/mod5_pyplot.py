#!/usr/bin/env python

# Imports
import sys
import argparse
import pathlib
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from math import log, sqrt

# Example argparse
parser = argparse.ArgumentParser(description='Create figure from end bam csv')
parser.add_argument("-i", "--input_csv", type=pathlib.Path, required=True,
                    help="This is the file location for the end bam csv, full path required")
parser.add_argument("-o", "--output_dir", type=pathlib.Path, required=True,
                    help="This is the output file path, full path required")
parser.add_argument("-p", "--prefix", type=str, required=True,
                    help="Prefix for output file name")
args = parser.parse_args()

print("\nStarting to produce figure of telomere lengths by Chromosome and Chromosome End\n")
print("\nThese are the options for the visualization script",
      "\n\tThis is the input file:", args.input_csv,
      "\n\tThis is the output directory:", args.output_dir,
      "\n\tThis is the prefix for the output file:", args.prefix,
      "\n")

# Check if input file and output directory exist
if Path(args.input_csv).is_file():
    print("\n{input_csv} is a file, continuing.".format(input_csv=args.input_csv))
else:
    sys.exit("\n{input_csv} does not exist, exiting.".format(input_csv=args.input_csv))

if Path(args.output_dir).is_dir():
    print("\n{out_dir} is a valid directory, continuing.\n".format(out_dir=args.output_dir))
else:
    sys.exit("\n{out_dir} is not a valid directory, exiting.\n".format(out_dir=args.output_dir))


# Setting up some interactive options
plt.ioff()

# Cleaning up data some
in_fn = args.input_csv
out_dir = str(args.output_dir)

data = pd.read_csv(in_fn)

data.rename(columns={'s.name': 'Sample Name', 's.win': 'Sliding Window (bp)',
                     'r.type': 'Read Type', 'tel.length': 'Telomere Length (kb)',
                     'norm': 'Normalization', 'seqnames': 'Chromosome', 'chr.end': 'Chromosome End'}, inplace=True)

# Process data and split into left and right ends
data_all = data[data['Read Type'].str.match('nl')]
data_all = data_all.sort_values(by='Sample Name')

# Process data into 5' End data only for plotting
l_data = data_all[data_all['Chromosome End'].str.match('5')]
l_data = l_data.rename(columns={"Telomere Length (kb)": "5' End"})

# Process data into 3' End data only for plotting
r_data = data_all[data_all['Chromosome End'].str.match('3')]
r_data = r_data.rename(columns={"Telomere Length (kb)": "3' End"})

# Set up plots
sns.set(context='paper', style='darkgrid')
sns.set_palette('Set1', color_codes=True)

# Increase fonts sizes for Figure 3 & Figure 4

plt.rc('font', size=20)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)

fig_all, axes = plt.subplots(nrows=2, sharex='all', sharey='all', figsize=(16, 12))

# Plot combined data, individual chromosome data
x = "Chromosome"
y = "5' End"

# Order
order = sorted(data_all['Chromosome'].unique())
hue = "Sample Name"

lchr = sns.violinplot(x=x, y=y, data=l_data,
                      hue=hue,
                      split=True,
                      inner="quartile",
                      height=4, aspect=1.0,
                      legend=False, cut=0,
                      palette=['r', 'b'],
                      order=order,
                      bw=0.5,
                      ax=axes[0])

lchr.legend_.remove()
lchr.set(xlabel=None)
lchr.set_xticklabels(lchr.get_xticklabels(), rotation=90)

# Reset Y for rchr
y = "3' End"

rchr = sns.violinplot(x=x, y=y, data=r_data,
                      hue=hue,
                      split=True,
                      inner="quartile",
                      height=4, aspect=1.0,
                      palette=['r', 'b'],
                      legend=False, cut=0,
                      order=order,
                      bw=0.5,
                      ax=axes[1])

rchr.legend_.remove()
rchr.set(xlabel=None)
rchr.set_xticklabels(rchr.get_xticklabels(), rotation=90)

# Add titles and formatting to chromosome end plot, figure 2
# Titles
fig_all.suptitle("Telomere Length by Chromosome End", fontsize=42,
                 y=0.95, x=0.11, ha='left')
fig_all.text(0.02, 0.5, 'Telomere Length (kb)', va='center',
             rotation='vertical', fontdict={'fontsize': 30})
fig_all.text(0.40, 0.03, 'Chromosome',
             fontdict={'fontsize': 30})

# Legends
handles, labels = axes[1].get_legend_handles_labels()
irr_status = fig_all.legend(handles, labels,
                            bbox_to_anchor=(0.988, 0.55),
                            prop={'size': 26}, title='Sample')
axes = plt.gca().add_artist(irr_status)

# Layout
plt.subplots_adjust(left=0.15, bottom=0.28, right=0.865, top=0.88)

# Save Figure 2
fil1_name = 'telo_length_by_chromosome'
plt1_f = out_dir + "/" + args.prefix + "_" + fil1_name + ".png"
plt.savefig(plt1_f, dpi=600)
plt.close('all')