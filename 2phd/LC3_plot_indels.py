#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:50:38 2023

@author: conrad
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

plotting_df = pd.read_csv("/home/conrad/les_complete/structural_variant_analysis/LC3_svimasm/LiP13-test.coords.txt",
                          sep="\t", header=0)
plotting_df['Start'] = plotting_df['Start'].astype(int)
plotting_df['End'] = plotting_df['End'].astype(int)

del_df = pd.read_csv("/home/conrad/les_complete/structural_variant_analysis/LC3_svimasm/LiP13.dels_only.txt",
                     sep="\t", header=0)
del_df['Start'] = del_df['Start'].astype(int)
del_df['End'] = del_df['End'].astype(int)

ins_df = pd.read_csv("/home/conrad/les_complete/structural_variant_analysis/LC3_svimasm/LiP13.ins_only.txt",
                     sep="\t", header=0)
ins_df['Start'] = ins_df['Start'].astype(int)
ins_df['End'] = ins_df['End'].astype(int)

# https://stackoverflow.com/questions/60009441/making-a-timeline-graph-with-a-dataframe-with-grouped-values-needing-a-for-loop
# plt.figure(figsize=(10,50))
# isolates_set = set(plotting_df['Isolate'])
# isolates = {isolate: i+1 for i, isolate in enumerate(sorted(isolates_set))}
# for isolate in isolates:
#     periods = []
#     for sample, start, end in zip(plotting_df['Isolate'], plotting_df['Start'], plotting_df['End']):
#         if isolate == sample:
#             periods.append((start, end-start))
#     plt.broken_barh(periods, (isolates[isolate] -0.45, 0.9))
# plt.yticks(ticks=range(1,206,1), labels=sorted(isolates_set))
# plt.ylim(0,206)

# plt.show()

les_prophages_gis_corrected = [
    ["P1", 706746, 721570],
    ["P2", 944892, 987088],
    ["P3", 1515027, 1558054],
    ["P4", 1829272, 1866077],
    ["P5", 3001435, 3051335],
    ["P6", 4901102, 4908700],
    ["G1", 2814147, 2861856],
    ["G2", 3062785, 3094546],
    ["G3", 3107882, 3218515],
    ["G4", 3704160, 3743729],
    ["G5", 5287607, 5318260],
    ["G6", 281883, 294308],
    ["G7", 2199910, 2216507],
    ["G8", 2920482, 2932053],
    ["G9", 2950325, 2958647],
    ["G10", 3458284, 3504146],
    ["G11", 4054701, 4109708],
    ["G12", 4299556, 4307949],
    ["G13", 5124686, 5146312],
    ["G14", 5323111, 5330434],
    ["G15", 6023041, 6047970],
    ["G16", 6123044, 6150293],
    ["G17",6595774, 6632768]]

les_prophages_gis = [
     ["P1", 665561, 680385],
     ["P2", 863875, 906018],
     ["P3", 1433825, 1476616],
     ["P4", 1684114, 1720919],
     ["P5", 2691759, 2741659],
     ["P6", 4546499, 4554097],
     ["G1", 2504700, 2552409],
     ["G2", 2753109, 2784809],
     ["G3", 2798145, 2908715],
     ["G4", 3394109, 3433537],
     ["G5", 4932837, 4962250],
     ["G6", 280017, 292442],
     ["G7", 2054752, 2071349],
     ["G8", 2611035, 2622606],
     ["G9", 2640878, 2649200],
     ["G10", 3148400, 3194166],
     ["G11", 3744468, 3798235],
     ["G12", 3947016, 3955409],
     ["G13", 4770083, 4791709],
     ["G14", 4967101, 4974424],
     ["G15", 5574296, 5585225],
     ["G16", 5658978, 5686227],
     ["G17", 6129698, 6166692]]

isolates_set = set(plotting_df['Isolate'])
isolates = {isolate: i+1 for i, isolate in enumerate(sorted(isolates_set))}

fig, ax = plt.subplots(figsize=(50,50))
for element in les_prophages_gis_corrected:
    ax.axvspan(element[1], element[2], color="lightgrey")
for isolate in isolates:
    periods = []
    # ins = []
    # dels = []
    for sample, start, end in zip(plotting_df['Isolate'], plotting_df['Start'], plotting_df['End']):
        if isolate == sample:
            periods.append((start, end-start))
    # for sample, start, end in zip(ins_df['Isolate'], ins_df['Start'], ins_df['End']):
    #     if isolate == sample:
    #         ins.append((start, end-start))
    # for sample, start, end in zip(del_df['Isolate'], del_df['Start'], del_df['End']):
    #     if isolate == sample:
    #         dels.append((start, end-start))
    ax.broken_barh(periods, (isolates[isolate] -0.45, 0.9))
    # ax.broken_barh(ins, (isolates[isolate] -0.45, 0.9), color='red')
    # ax.broken_barh(dels, (isolates[isolate] -0.45, 0.9), color='red')
    
ax.set_ylim(-1,209)
ax.set_xlim(-100000,7250000)
ax.set_yticks(range(1,207,1),sorted(isolates_set), fontsize=15.0)
ax.set_xticks(range(0,7500000, 500000))
ax.tick_params(axis='x', labelrotation=45, labelsize=25.0)
ax.set_xlabel("Position in Reference (Mbp)", fontsize=30.0)
ax.set_ylabel("Isolate", fontsize=30.0)

plt.savefig("/home/conrad/les_complete/structural_variant_analysis/LC3_svimasm/sv_plot_indels.pdf", dpi=300)

# Dels only:
fig2, ax2 = plt.subplots(figsize=(50,50))
for element in les_prophages_gis:
    ax2.axvspan(element[1], element[2], color="lightgrey")
for isolate in isolates:
    del_periods = []
    for sample, start, end in zip(del_df['Isolate'], del_df['Start'], del_df['End']):
        if isolate == sample:
            del_periods.append((start, end-start))
    ax2.broken_barh(del_periods, (isolates[isolate] -0.45, 0.9))
ax2.set_ylim(-1,209)
ax2.set_xlim(-100000, 6700000)
ax2.set_yticks(range(1,207,1), sorted(isolates_set), fontsize=15)
ax2.set_xticks(range(0,6500000, 500000))
ax2.tick_params(axis='x', labelrotation=45, labelsize=25)
ax2.set_xlabel("Position in Reference (Mbp)", fontsize=30)
ax2.set_ylabel("Isolate", fontsize=30)

plt.savefig("/home/conrad/les_complete/structural_variant_analysis/LC3_svimasm/sv_plot_dels_only.pdf", dpi=300)

# ax.set_xticks([500000, 1500000, 2500000, 3500000, 4500000, 5500000, 6500000], minor=True)
# ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(('%.1E'))) # https://stackoverflow.com/questions/25750170/show-decimal-places-and-scientific-notation-on-the-axis






