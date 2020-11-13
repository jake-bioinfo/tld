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

order = ['1', '2', '3', '4', '5', '6', '7',
         '8', '9', '10', '11', '12', '13', '14',
         '15', '16', '17', '18', '19', 'X', 'Y',
         'MT']
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

# Create data frames before concatenation
conditions = [

    # Irradiated Healed
    (data_all['Chromosome'].isin([1, 3])) & (data_all['Chromosome End'] == "5'") &
     (data_all['Sample Name'] == "irr"),
    (data_all['Chromosome'].isin([2])) & (data_all['Chromosome End'] == "3'") &
     (data_all['Sample Name'] == "irr"),

    # Recently healed non-irradiated
    (data_all['Chromosome'].isin([2, 6, 14]) & (data_all['Chromosome End'] == "5'") &
     (data_all['Sample Name'] == "wt")),
    (data_all['Chromosome'].isin([5, 11, 14]) & (data_all['Chromosome End'] == "3'") &
     (data_all['Sample Name'] == "wt")),

    # Recently healed irradiated
    (data_all['Chromosome'].isin([2, 6, 14]) & (data_all['Chromosome End'] == "5'") &
     (data_all['Sample Name'] == "irr")),
    (data_all['Chromosome'].isin([5, 11, 14]) & (data_all['Chromosome End'] == "3'") &
     (data_all['Sample Name'] == "irr")),

    # Non-healed irradiated
    (~data_all['Chromosome'].isin([1, 3, 2, 6, 14]) & (data_all['Chromosome End'] == "5'") &
     (data_all['Sample Name'] == "irr")),
    (~data_all['Chromosome'].isin([2, 5, 11, 14]) & (data_all['Chromosome End'] == "3'") &
     (data_all['Sample Name'] == "irr")),

    # Non-healed non-irradiated
    (~data_all['Chromosome'].isin([2, 6, 14]) & (data_all['Chromosome End'] == "5'") &
     (data_all['Sample Name'] == "wt")),
    (~data_all['Chromosome'].isin([5, 11, 14]) & (data_all['Chromosome End'] == "3'") &
     (data_all['Sample Name'] == "wt")),

]

outputs = [
    'irr-H', 'irr-H', 'wt-PH', 'wt-PH', 'irr-PH', 'irr-PH', 'irr-NH', 'irr-NH', 'wt-NH', 'wt-NH'
]


# Order
order = ['irr-H', 'irr-PH', 'irr-NH', 'wt-PH', 'wt-NH']
data_all['Status'] = np.select(conditions, outputs, 0)

# Assign group dfs and Randomly downsample
irr_H = data_all[data_all['Status'] == "irr-H"]
irr_PH = data_all[data_all['Status'] == "irr-PH"]
irr_NH = data_all[data_all['Status'] == "irr-NH"]
wt_PH = data_all[data_all['Status'] == "wt-PH"]
wt_NH = data_all[data_all['Status'] == "wt-NH"]

# Perform t-tests
min_count = min(len(irr_H), len(irr_PH), len(irr_NH), len(wt_PH), len(wt_NH))

# Perform while loop to assess all p-values and statistics
# Setup dataframes and initial values before loop
stat_df = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
sp = 5

# Loop over different samples sizes until minimum sample number adding stats
while sp <= min_count:

    # Set Sampling
    irr_H_s = irr_H.sample(sp, random_state=1)
    irr_PH_s = irr_PH.sample(sp, random_state=1)
    irr_NH_s = irr_NH.sample(sp, random_state=1)
    wt_PH_s = wt_PH.sample(sp, random_state=1)
    wt_NH_s = wt_NH.sample(sp, random_state=1)

    # Create irradiated healed and wild type non-healed stat and row
    irrH_v_wtNH = stats.ttest_ind(irr_H_s['Telomere Length (kb)'], wt_NH_s['Telomere Length (kb)'], equal_var=False)
    add_irrH_v_wtNH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_irrH_v_wtNH._set_value(sp, 'Sample Comparison', 'irr-H_wt-NH')
    add_irrH_v_wtNH._set_value(sp, 'Sample Number', sp)
    add_irrH_v_wtNH._set_value(sp, 'statistic', irrH_v_wtNH[0])
    add_irrH_v_wtNH._set_value(sp, 'pvalue', -log(irrH_v_wtNH[1]))

    # Create irradiated healed vs wild type non-healed stat and row
    irrNH_v_wtNH = stats.ttest_ind(irr_NH_s['Telomere Length (kb)'], wt_NH_s['Telomere Length (kb)'], equal_var=False)
    add_iNH_v_wtNH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_iNH_v_wtNH._set_value(sp, 'Sample Comparison', 'irr-NH_wt-NH')
    add_iNH_v_wtNH._set_value(sp, 'Sample Number', sp)
    add_iNH_v_wtNH._set_value(sp, 'statistic', irrNH_v_wtNH[0])
    add_iNH_v_wtNH._set_value(sp, 'pvalue', -log(irrNH_v_wtNH[1]))

    # Create irradiated recently healed and wild type recently healed stat and row
    irrPH_v_wtPH = stats.ttest_ind(irr_PH_s['Telomere Length (kb)'], wt_PH_s['Telomere Length (kb)'], equal_var=False)
    add_irrPH_v_wtPH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_irrPH_v_wtPH._set_value(sp, 'Sample Comparison', 'irr-PH_wt-PH')
    add_irrPH_v_wtPH._set_value(sp, 'Sample Number', sp)
    add_irrPH_v_wtPH._set_value(sp, 'statistic', irrPH_v_wtPH[0])
    add_irrPH_v_wtPH._set_value(sp, 'pvalue', -log(irrPH_v_wtPH[1]))

    # Create wild type recently healed vs wild type non-healed stat and row
    wtPH_v_wtNH = stats.ttest_ind(wt_PH_s['Telomere Length (kb)'], wt_NH_s['Telomere Length (kb)'], equal_var=False)
    add_wtPH_v_wtNH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_wtPH_v_wtNH._set_value(sp, 'Sample Comparison', 'wt-PH_wt-NH')
    add_wtPH_v_wtNH._set_value(sp, 'Sample Number', sp)
    add_wtPH_v_wtNH._set_value(sp, 'statistic', (wtPH_v_wtNH[0]))
    add_wtPH_v_wtNH._set_value(sp, 'pvalue', -log(wtPH_v_wtNH[1]))

    # Create irradiated recently healed and irradiated non-healed stat and row
    irrPH_v_irrNH = stats.ttest_ind(irr_PH_s['Telomere Length (kb)'], irr_NH_s['Telomere Length (kb)'], equal_var=False)
    add_irrPH_v_irrNH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_irrPH_v_irrNH._set_value(sp, 'Sample Comparison', 'irr-PH_irr-NH')
    add_irrPH_v_irrNH._set_value(sp, 'Sample Number', sp)
    add_irrPH_v_irrNH._set_value(sp, 'statistic', irrPH_v_irrNH[0])
    add_irrPH_v_irrNH._set_value(sp, 'pvalue', -log(irrPH_v_irrNH[1]))

    # Create irradiated healed and irradiated non-healed stat and row
    irrH_v_irrNH = stats.ttest_ind(irr_H_s['Telomere Length (kb)'], irr_NH_s['Telomere Length (kb)'], equal_var=False)
    add_irrH_v_irrNH = pd.DataFrame(data=None, columns=['Sample Comparison', 'Sample Number', 'statistic', 'pvalue'])
    add_irrH_v_irrNH._set_value(sp, 'Sample Comparison', 'irr-H_irr-NH')
    add_irrH_v_irrNH._set_value(sp, 'Sample Number', sp)
    add_irrH_v_irrNH._set_value(sp, 'statistic', irrH_v_irrNH[0])
    add_irrH_v_irrNH._set_value(sp, 'pvalue', -log(irrH_v_irrNH[1]))

    # Concatenate all add df
    add = pd.concat([add_irrH_v_wtNH,
                     add_iNH_v_wtNH,
                     #add_irrPH_v_wtPH,
                     add_wtPH_v_wtNH,
                     add_irrH_v_irrNH,
                     add_irrPH_v_irrNH])

    # Append to result dataframe
    stat_df = stat_df.append(add)

    # Increase sample size
    sp = sp + 35

# Coerce values to float
sp_num = stat_df['Sample Number'].astype(float)
stat_df = stat_df.assign(n=sp_num.values)
p_float = stat_df['pvalue'].astype(float)
stat_df = stat_df.assign(log_p=p_float.values)

# Find minimum sample count
sp_count = 35
irr_H_s = irr_H.sample(sp_count, random_state=1)
irr_PH_s = irr_PH.sample(sp_count, random_state=1)
irr_NH_s = irr_NH.sample(sp_count, random_state=1)
wt_PH_s = wt_PH.sample(sp_count, random_state=1)
wt_NH_s = wt_NH.sample(sp_count, random_state=1)

# Create basic stats text file for n=35
s_ls_df = list([irr_H_s, irr_PH_s, irr_NH_s, wt_PH_s, wt_NH_s])
s_ls_nm = list(["irr_H_s", "irr_PH_s", "irr_NH_s", "wt_PH_s", "wt_NH_s"])

len_txt_s_f = out_dir + "/" + "Telomere_grouped_lengths_35" + ".txt"
len_txt_s = open(len_txt_s_f, "w")
i = 0

# For loop over list
# while i <= (len(s_ls_df)-1):
#     print(f"\nStats that follow are n=35 for sample: {s_ls_nm[i]}", file=len_txt_s)
#     add_text = s_ls_df[i]['Telomere Length (kb)'].describe()
#     se = s_ls_df[i]['Telomere Length (kb)'].describe()['std'] / (sqrt(35))
#     print(f"{add_text}", file=len_txt_s)
#     print(f"SE:\t{se}\n", file=len_txt_s)
#     i += 1

len_txt_s.close()

# Create basic stats text file for n=1075
sp_count_1 = 1075
irr_H_1 = irr_H.sample(sp_count_1, random_state=1)
irr_PH_1 = irr_PH.sample(sp_count_1, random_state=1)
irr_NH_1 = irr_NH.sample(sp_count_1, random_state=1)
wt_PH_1 = wt_PH.sample(sp_count_1, random_state=1)
wt_NH_1 = wt_NH.sample(sp_count_1, random_state=1)

ls_df = list([irr_H_1, irr_PH_1, irr_NH_1, wt_PH_1, wt_NH_1])
ls_nm = list(["irr_H", "irr_PH", "irr_NH", "wt_PH", "wt_NH"])

len_txt_f = out_dir + "/" + "Telomere_grouped_lengths_1075" + ".txt"
len_txt = open(len_txt_f, "w")
i = 0

# For loop over list
while i <= (len(ls_df)-1):
    print(f"\nStats that follow are n=1075 for sample: {ls_nm[i]}", file=len_txt)
    add_text = ls_df[i]['Telomere Length (kb)'].describe()
    se = ls_df[i]['Telomere Length (kb)'].describe()['std'] / (sqrt(1075))
    print(f"{add_text}", file=len_txt)
    print(f"SE:\t{se}\n", file=len_txt)
    i += 1

len_txt.close()

# Concatenate into larger group for graphing
d_graph = pd.concat([irr_H_s, irr_PH_s, irr_NH_s, wt_PH_s, wt_NH_s], sort=True)

# Set Context for presentations
sns.set(context='talk', style='white', font='Times New Roman')
sns.set_palette('Set1', color_codes=True)

# Setup grouped figure
fig_grouped, ax = plt.subplots(sharex='all', sharey='all', figsize=(14, 10))

# Plot grouped by telomere healed and non telomere healed
x = "Status"
y = "Telomere Length (kb)"

# Set up statistics comparisons
box_pairs = [('irr-NH', 'wt-NH'), ('wt-PH', 'wt-NH'),
             #('irr-PH', 'wt-PH'),
             ('irr-H', 'wt-NH'),
             ('irr-PH', 'irr-NH'), ('irr-H', 'irr-NH')]

ax = sns.violinplot(x=x, y=y, data=d_graph,
                    inner="box",
                    linewidth=2,
                    legend=True,
                    cut=0, bw=0.5,
                    order=order)

# Add statistics annotations
add_stat_annotation(ax, data=d_graph, x=x, y=y,
                    box_pairs=box_pairs,
                    test='t-test_ind', text_format='star',
                    loc='inside', verbose=2, order=order)

# Add figure title
plt.title("Grouped Length Comparison", fontsize=42,
             loc='left', pad=30)
plt.xlabel('Status', fontsize=32, labelpad=15)
plt.ylabel('Telomere Length (kb)', fontsize=32, labelpad=15)

# Save figure
fil_grouped_name = 'MPM.Figure_3'
plt_grouped_f = out_dir + "/" + fil_grouped_name + ".png"
plt.savefig(plt_grouped_f, optomize=True, progressive=True, dpi=600)
plt.close('all')

# Plot p-value changes dependent on sample size
g = sns.lmplot(x='n', y='log_p', data=stat_df,
           hue='Sample Comparison',
           order=2, legend=False)

plt.legend(loc='best', prop={'size':20}, markerscale=2)
plt.xlim(0,1075)
plt.ylim(0,500)
plt.title('Grouped Telomere P-values', loc='left', pad=30, fontsize=46)
plt.xlabel('Sample Size', fontsize=32, labelpad=15)
plt.ylabel('-log(p)', fontsize=32, labelpad=15)
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.88, top=0.88)
fig = plt.gcf()
fig.set_size_inches(14, 10)

# Save figure
fil_pvalue_name = 'MPM_Figure_4'
plt_pvalue_f = out_dir + "/" + fil_pvalue_name + ".png"
plt.savefig(plt_pvalue_f, optomize=True, progressive=True, dpi=600)
plt.close('all')