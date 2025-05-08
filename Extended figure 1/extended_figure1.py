#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:19:30 2024

@author: aliciapageau
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import zscore




# %% Import files
working_directory = '/home/alicia/Documents/antifungal_project/mardy2.0/clean_code/update_jan_2025/'
os.chdir(working_directory)

data = pd.read_csv(f"{working_directory}FungAMR_070425.csv", index_col=0)
gene_class = pd.read_csv("genes_class.csv")
gemme = pd.read_csv(f"{working_directory}gemme_scores.csv", index_col=0)
mutatex = pd.read_csv(f"{working_directory}mutatex_results.csv", index_col=0)
mutatex['ddG'] = mutatex['ddG']*-1

# Drop rows without ref_seq accession in data
data['ref_seq_uniprot_accession'] = data['ref_seq_uniprot_accession'].astype(str)
data = data[~data['ref_seq_uniprot_accession'].str.contains('nan')]

# Get database to long format (one mutation per row)
df_long = data.copy()
df_long['gene or protein name'] = df_long['gene or protein name'].apply(lambda x: x.split(','))
df_long['ortho_homolog'] = df_long['ortho_homolog'].apply(lambda x: x.split(','))
df_long['ref_seq_uniprot_accession'] = df_long['ref_seq_uniprot_accession'].apply(lambda x: x.split(','))

# Convert mutations to list of list
df_long['mutation'] = df_long.apply(lambda row: 
    [row['mutation'].split(',')] if len(row['gene or protein name']) <= 1 else 
    [item.strip().split(',') for item in row['mutation'].split('|')],
    axis=1)

# Concatenate and exploded dataframe
df_long = df_long.explode(['gene or protein name', 'mutation', 'ortho_homolog','ref_seq_uniprot_accession'], ignore_index=True)
df_long = df_long.explode('mutation').reset_index(drop=True)
df_long['gene or protein name'] = df_long['gene or protein name'].str.strip()
df_long['mutation'] = df_long['mutation'].str.strip()
df_long = df_long.drop_duplicates()

# Extract position from mutation
df_long[['wt_AA', 'position', 'alt_AA']] = df_long['mutation'].str.extract(r'(\D+)(\d+)(\D+)')

final = df_long.copy()

# %% Extended figure 1 - Ridgeline plot
# Merging gemme and mutatex datasets, removing duplicates
gemme_mutatex = pd.merge(gemme, mutatex, how='inner').drop_duplicates()

# Clipping ddG values at -10 and calculating z-scores
gemme_mutatex.loc[gemme_mutatex['ddG'] < -10, 'ddG'] = -10
gemme_mutatex['zscore_gemme'] = zscore(gemme_mutatex['gemme'])
gemme_mutatex['zscore_ddG'] = zscore(gemme_mutatex['ddG'])

# Filtering resistant mutations and setting DMS for specific evidence degrees
resistant_mutations = final[final['confidence score'] > 0].copy()
resistant_mutations['position'] = pd.to_numeric(resistant_mutations['position'], errors='coerce')
resistant_mutations.loc[resistant_mutations['confidence score'].isin([3, -3]), 'DMS'] = ' DMS'

# Merging resistant mutations with gemme and mutatex data
merged_data = pd.merge(resistant_mutations, gemme_mutatex, how='inner',
                       left_on=['species', 'ref_seq_uniprot_accession', 'ortho_homolog', 'wt_AA', 'position', 'alt_AA', 'mutation'],
                       right_on=['species', 'uniprot_accession', 'ortho_homolog', 'wt_AA', 'position', 'alt_AA', 'mutation'])

# List of species to match and rename ortho_homolog to Erg11
species_list = [
    'Candida albicans', 'Candida auris', 'Candida dubliniensis', 'Candida krusei',
    'Candida orthopsilosis', 'Candida parapsilosis', 'Candida tropicalis',
    'Cryptococcus neoformans', 'Cryptococcus neoformans/Cryptococcus deneoformans',
    'Nakaseomyces glabratus', 'Saccharomyces cerevisiae'
]

merged_data.loc[(merged_data['species'].isin(species_list)) &
                (merged_data['ortho_homolog'] == 'Cyp51'), 'ortho_homolog'] = 'Erg11'

# Creating unique identifier ID and selecting relevant columns
merged_data['ID'] = merged_data['ortho_homolog'].str.cat(merged_data['DMS'].fillna(''), sep='')
plot_data = merged_data[['ID', 'species', 'ortho_homolog', 'uniprot', 'mutation', 'gemme', 'ddG', 'zscore_gemme', 'zscore_ddG']]
plot_data = plot_data.drop_duplicates()

# Filtering genes with at least 10 unique mutations
gene_mutation_counts = plot_data.groupby('ID').mutation.nunique().reset_index(name='mut_count')
filtered_genes = gene_mutation_counts.query('mut_count >= 10').drop('mut_count', axis=1).drop_duplicates()
plot_filtered = plot_data[plot_data['ID'].isin(filtered_genes['ID'])].copy()

# Ensuring DataFrame is sorted by IDs with the highest number of resistant mutations
plot_filtered['ID'] = pd.Categorical(plot_filtered['ID'], categories=filtered_genes['ID'], ordered=True)
plot_filtered = plot_filtered.sort_values('ID')

# Calculating mean, standard deviation, and standard error of z-scores
plot_filtered['count'] = plot_filtered.groupby('ID')['mutation'].transform('nunique')
plot_filtered['mean_zscore_gemme'] = plot_filtered.groupby('ID')['zscore_gemme'].transform('mean')
plot_filtered['mean_zscore_ddG'] = plot_filtered.groupby('ID')['zscore_ddG'].transform('mean')
plot_filtered['std_zscore_gemme'] = plot_filtered.groupby('ID')['zscore_gemme'].transform('std')
plot_filtered['std_zscore_ddG'] = plot_filtered.groupby('ID')['zscore_ddG'].transform('std')
plot_filtered['se_zscore_gemme'] = plot_filtered['std_zscore_gemme'] / np.sqrt(plot_filtered['count'])
plot_filtered['se_zscore_ddG'] = plot_filtered['std_zscore_ddG'] / np.sqrt(plot_filtered['count'])

# 'plot_filtered' now contains the cleaned and processed data ready for analysis or plotting.
# Define the bin edges and corresponding labels for 'count'
bins = [0, 30, 100, 500, 1500]
labels = [0, 250, 500, 750]

# Bin 'count' into categories and convert to float
plot_filtered['size'] = pd.cut(plot_filtered['count'], bins=bins, labels=labels, include_lowest=True)
plot_filtered['size'] = plot_filtered['size'].astype(float)

# Select relevant columns and remove duplicates
plot_filtered = plot_filtered[['ID', 'ortho_homolog', 'count', 'size',
                               'zscore_gemme', 'zscore_ddG',
                               'mean_zscore_gemme', 'mean_zscore_ddG',
                               'std_zscore_gemme', 'std_zscore_ddG',
                               'se_zscore_gemme', 'se_zscore_ddG']].drop_duplicates()

# Ridgeline plot
plot = plot_filtered.copy()

gene_class['ortho_homolog'].replace('CytB','Cytochrome b',  regex=True, inplace=True)
gene_class['ortho_homolog'].replace('β-tubulin', 'Beta-tubulin', regex=True, inplace=True)
gene_class['ortho_homolog'].replace('Sqle', 'Squalene epoxidase', regex=True, inplace=True)

plot = pd.merge(plot, gene_class)

#ids and class order
ids = ['Squalene epoxidase', 'Fur1', 'Fcy1 DMS', 'Fks DMS', 'Fks', 'Dhfr DMS', 'Dhfr',
       'Erg11 DMS', 'Erg11', 'Cyp51','Cytochrome b', 'Sdh', 'Cdr1', 'Pdr1',
       'Tac1', 'Mrr1', 'Hmg1','Sur1', 'Csg2', 'Erg3', 'Erg25','Msh2']

plot = plot[plot['ID'].isin(ids)]
plot['ID'] = pd.Categorical(plot['ID'], categories=ids, ordered=True)
plot = plot.sort_values('ID')


## Ridge line plot - GEMME
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
pal = {"drug target" : "#8C0250",
       "multidrug transporter" : "#238CCD",
       "transcription factor" : "#C7B3FF",
       "other" : "#616161"}


# in the sns.FacetGrid class, the 'hue' argument is the one that is the one that will be represented by colors with 'palette'
g = sns.FacetGrid(plot, row='ID', hue='class', palette=pal, height=0.75, aspect=10)

# then we add the densities kdeplots for each month
g.map(sns.kdeplot, 'zscore_gemme',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=1, linewidth=1)

# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'zscore_gemme', 
      bw_adjust=1, clip_on=False, 
      color="w", lw=1)

# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=2, clip_on=False)

# we loop over the FacetGrid figure axes (g.axes.flat) and add the month as text with the right color
# notice how ax.lines[-1].get_color() enables you to access the last line's color in each matplotlib.Axes
for i, ax in enumerate(g.axes.flat):
    if len(ax.lines) > 0:
        id_label = plot['ID'].unique()[i]
        ax.text(-4.5, 0.1, f'{id_label}', fontweight='bold', fontsize=18, color="black")
        ax.set_xlim(-5,3)

# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.5)

# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=[])
g.set(ylabel="")
g.despine(bottom=True, left=True)
g.fig.set_dpi(600)
plt.setp(ax.get_xticklabels(), fontsize=12, fontweight='bold')
plt.xlabel('z-score gemme', fontweight='bold', fontsize=18)
plt.show()
#plt.rcParams['svg.fonttype'] = 'none'
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_gemme.png", dpi=600, bbox_inches='tight')
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_gemme.svg", bbox_inches='tight')


## Ridge line plot - MutateX
# in the sns.FacetGrid class, the 'hue' argument is the one that is the one that will be represented by colors with 'palette'
g = sns.FacetGrid(plot, row='ID', hue='class', palette=pal, height=0.75, aspect=10)

# then we add the densities kdeplots for each month
g.map(sns.kdeplot, 'zscore_ddG',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=1, linewidth=1)

# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'zscore_ddG', 
      bw_adjust=1, clip_on=False, 
      color="w", lw=1)

# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=2, clip_on=False)

# we loop over the FacetGrid figure axes (g.axes.flat) and add the month as text with the right color
# notice how ax.lines[-1].get_color() enables you to access the last line's color in each matplotlib.Axes
for i, ax in enumerate(g.axes.flat):
    if len(ax.lines) > 0:
        id_label = plot['ID'].unique()[i]
        ax.text(-4.5, 0.1, f'{id_label}', fontweight='bold', fontsize=18, color="black") 
        ax.set_xlim(-5,3)

# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.5)

# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=[])
g.set(ylabel="")
g.despine(bottom=True, left=True)
g.fig.set_dpi(600)
plt.setp(ax.get_xticklabels(), fontsize=12, fontweight='bold')
plt.xlabel('z-score ddG * -1', fontweight='bold', fontsize=18)
plt.show()
#plt.rcParams['svg.fonttype'] = 'none'
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_foldx.png", dpi=600, bbox_inches='tight')
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_foldx.svg", bbox_inches='tight')

## Grayscale bar with the number of resistant mutations
# Calculate the number of observations for grayscale
plot['count'] = plot.groupby('ID')['ID'].transform('count')
plot['norm_count'] = (plot['count'] - plot['count'].min()) / (plot['count'].max() - plot['count'].min())
bins = [0, 30, 100, 500, 1500]  # Define the bin edges
labels = ['[0, 30[', '[30, 100[', '[100, 500[', '500+']  # Define the bin labels
plot['count_cat'] = pd.cut(plot['count'], bins=bins, labels=labels, include_lowest=True)
grayscale_color = {'[0, 30[':'#DDDDDD', '[30, 100[':'#AAAAAA', '[100, 500[':'#555555', '500+':'#000000'}
plot['color'] = plot['count_cat'].map(grayscale_color)

# Plot
fig, ax = plt.subplots(figsize=(1, 22))

for i, id_ in enumerate(plot['ID'].unique()):
    color = plot.loc[plot['ID'] == id_, 'color'].values[0]
    ax.add_patch(plt.Rectangle((0, len(plot['ID'].unique()) - i - 1), 1, 1, color=color))

# Set the limits and remove axes
ax.set_xlim(0, 1)
ax.set_ylim(0, len(plot['ID'].unique()))
ax.axis('off')
 
legend_patches = [Patch(color=color, label=label) for label, color in grayscale_color.items()]

# Add legend to the plot
plt.legend(handles=legend_patches, title='Resistant mutations', bbox_to_anchor=(1.05, 1), loc='upper left', frameon = False, fontsize=12)

plt.show()
#plt.rcParams['svg.fonttype'] = 'none'
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_grayscale.png", dpi=600, bbox_inches='tight')
#plt.savefig(f"{working_directory}figures/figure_supp1_ridgeline_grayscale.svg", bbox_inches='tight')