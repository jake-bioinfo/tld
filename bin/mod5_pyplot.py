# Imports
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from natsort import natsorted as ns
from statannot import add_stat_annotation
from scipy import stats
from math import log, sqrt


# Setting up some interactive options
plt.ioff()
#plt.ion()
desired_width = 200
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 50)

# Cleaning up data some
in_fn = '/home/jake/tmp/docker/test_data/tld/data/o_dir/dockTest.end.bam.csv'
data = pd.read_csv(in_fn)
out_dir = '/home/jake/Desktop'

data.rename(columns={'s.name': 'Sample Name', 's.win': 'Sliding Window (bp)',
                     'r.type': 'Read Type', 'tel.length': 'Telomere Length (kb)',
                     'norm': 'Normalization', 'seqnames': 'Chromosome', 'chr.end': 'Chromosome End'}, inplace=True)

# Process data and split into left and right ends
data_all = data[data['Read Type'].str.match('nl')]
data_all = data_all.sort_values(by='Sample Name')
data_all = data_all.replace("mm", "mouse")
data_all = data_all.replace("rm", "frog")

l_data = data_all[data_all['Chromosome End'].str.match('5')]
chr_str = l_data['Chromosome'].apply(str)
l_data = l_data.assign(Chromosome=chr_str.values)
l_data.rename(columns={"Telomere Length (kb)": "5' End"}, inplace=True)
r_data = data_all[data_all['Chromosome End'].str.match('3')]
chr_str = r_data['Chromosome'].apply(str)
r_data = r_data.assign(Chromosome=chr_str.values)
r_data.rename(columns={"Telomere Length (kb)": "3' End"}, inplace=True)

# Basic stats for individual chromosome data

l_sub = l_data[['Sample Name', 'Chromosome', '5\' End']]
#chr_int = l_sub['Chromosome'].apply(int)
#l_sub = l_sub.assign(Chromosome=chr_int.values)
l_sub.sort_values(by=['Chromosome'])
l_stat = l_sub.groupby(['Sample Name', 'Chromosome']).describe()
#l_lat = l_stat.to_latex(index=True)

r_sub = r_data[['Sample Name', 'Chromosome', '3\' End']]
#chr_int = r_sub['Chromosome'].apply(int)
#r_sub = r_sub.assign(Chromosome=chr_int.values)
r_sub.sort_values(by=['Chromosome'])
r_stat = r_sub.groupby(['Sample Name', 'Chromosome']).describe()
#r_lat = r_stat.to_latex(index=True)

# Set up plots

sns.set(context='paper', style='darkgrid', font='Times New Roman')
sns.set_palette('Set1', color_codes=True)

# Increase fonts sizes for Figure 3 & Figure 4

plt.rc('font', size=20)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)

fig_all, axes = plt.subplots(nrows=2, sharex='all', sharey='all', figsize=(14, 10))

# Plot combined data, individual chromosome data
x = "Chromosome"
y = "5' End"

# Order

order = sorted(data['Chromosome'].unique())

order = ns(order)

hue = "Sample Name"
hue_order = ['mouse', 'frog']

lchr = sns.violinplot(x=x, y=y, data=l_data,
                      hue=hue,
                      split=True,
                      inner="quartile",
                      height=4, aspect=1.0,
                      legend=False, cut=1,
                      palette=['r', 'b'],
                      order=order,
                      bw=0.5,
                      ax=axes[0])

lchr.legend_.remove()
lchr.set(xlabel=None)
lchr.set_xticklabels(lchr.get_xticklabels(),rotation = 90)

# Reset Y for rchr
y = "3' End"

rchr = sns.violinplot(x=x, y=y, data=r_data,
                      hue=hue,
                      split=True,
                      inner="quartile",
                      height=4, aspect=1.0,
                      palette=['r', 'b'],
                      legend=False, cut=2,
                      order=order,
                      bw=0.5,
                      ax=axes[1])

rchr.legend_.remove()
rchr.set(xlabel=None)
rchr.set_xticklabels(rchr.get_xticklabels(),rotation = 90)

# Add titles and formatting to chromosome end plot, figure 2
# Titles
fig_all.suptitle("Telomere Length by Chromosome End", fontsize=42,
                 y=0.95, x=0.11, ha='left')
fig_all.text(0.02, 0.5, 'Telomere Length (kb)', va='center',
             rotation='vertical', fontdict={'fontsize': 30})
fig_all.text(0.45, 0.03, 'Chromosome',
             fontdict={'fontsize': 30})

# Legends
handles, labels = axes[1].get_legend_handles_labels()
irr_status = fig_all.legend(handles, labels,
                            bbox_to_anchor=(0.988, 0.65),
                            prop={'size': 26}, title='Sample')
axes = plt.gca().add_artist(irr_status)

# Layout
fig_all.tight_layout()
plt.subplots_adjust(left=0.11, bottom=0.12, right=0.865, top=0.88)

# Save Figure 2
fil1_name = 'test'
plt1_f = out_dir + "/" + fil1_name + ".png"
plt.savefig(plt1_f, optomize=True, progressive=True, dpi=600)
plt.close('all')