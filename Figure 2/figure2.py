#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:32:07 2024

@author: aliciapageau
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy
from scipy.stats import zscore
from matplotlib.colors import LinearSegmentedColormap

print('pandas', pd.__version__)
print('matplotlib', plt.matplotlib.__version__)
print('numpy', np.__version__)
print('seaborn', sns.__version__)
print('scipy', scipy.__version__)


gene_colors = {'Erg11 DMS': '#8C0250',
               'Fcy1 DMS': '#71AED4',
               'Dhfr DMS': '#C7B3FF',
               'Fks DMS': '#519497',
               'Pdr1': '#616161',
               'Erg11': '#8C0250',
               'Fur1': '#DAACC6',
               'Cyp51': '#8C0250',
               'Fks': '#519497',
               'Erg25': '#83BBBD',
               'Tac1': '#C8C8C8',
               'Cytochrome b': '#D89747',
               'Mrr1': '#A67490',
               'Sur1': '#E6ADBC',
               'Dhfr': '#C7B3FF',
               'Erg3': '#0F767A',
               'Beta-tubulin': '#959595',
               'Squalene epoxidase': '#E3C398',
               'Cdr1': '#420DA5', 
               'Csg2': '#B0C4DE', 
               'Hmg1': '#663399',
               'Msh2': '#DA8632',
               'Mrr2': '#DAA540', 
               'Hap1':'#9B0E6B', 
               'Upc2': '#238CCD', 
               'Erg29':'#B0C4DE', 
               'Pdr3': '#8361CB', 
               'Dhh1':'#20B2AA', 
               'Yrr1':'#404040',
               'Sdh':'#91CEF4', 
               'Fcy1':'#404040'}

# %% Import files
working_directory = '/home/alicia/Documents/antifungal_project/mardy2.0/clean_code/update_jan_2025/Figure 2/'
os.chdir(working_directory)

data = pd.read_csv(f"{working_directory}../FungAMR_04_30_2025.csv", index_col=0)
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

data = df_long.copy()

# %% Figure 2 - Gemme and MutateX
# Merging gemme and mutatex datasets, removing duplicates
gemme_mutatex = pd.merge(gemme, mutatex, how='inner').drop_duplicates()

# Clipping ddG values at -10 and calculating z-scores
gemme_mutatex.loc[gemme_mutatex['ddG'] < -10, 'ddG'] = -10
gemme_mutatex['zscore_gemme'] = zscore(gemme_mutatex['gemme'])
gemme_mutatex['zscore_ddG'] = zscore(gemme_mutatex['ddG'])
#gemme_mutatex.loc[gemme_mutatex['zscore_ddG'] < -4, 'zscore_ddG'] = -4

gemme_mutatex['species'] = gemme_mutatex['species'].str.strip()
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Nakaseomyces glabrata", "Nakaseomyces glabratus", regex=True)
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Candida famata", "Debaryomyces hansenii", regex=True)
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Candida krusei", "Pichia kudriavzevii", regex=True)
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Candida kefyr", "Kluyveromyces marxianus", regex=True)
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Candida lusitaniae", "Clavispora lusitaniae", regex=True)
gemme_mutatex["species"] = gemme_mutatex["species"].replace("Candida auris", "Candidozyma auris", regex=True)

# Filtering resistant mutations and setting DMS for specific evidence degrees
resistant_mutations = data[data['confidence score'] > 0].copy()
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
plot_data = merged_data[['ID', 'species', 'ortho_homolog', 'ref_seq_uniprot_accession', 'mutation', 'gemme', 'ddG', 'zscore_gemme', 'zscore_ddG']]
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

# Create gradient background
x = np.linspace(1, 0, 250)
y = np.linspace(1, 0, 250)
xArray, yArray = np.meshgrid(x, y)
plotArray = np.sqrt(xArray**2 + yArray**2)**2.5

# Create the figure and axis
sns.reset_defaults()
fig = plt.figure(figsize=(12, 6), dpi=600)
ax = fig.add_subplot(111)

# Display the gradient background using imshow
colors = ["white", "silver"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
image = ax.imshow(plotArray, cmap=custom_cmap, extent=[-1.6, 1, -1.6, 0.8], origin='lower', vmin=0, vmax=1)

# Plot error bars for mean_zscore_gemme and mean_zscore_ddG
ax.errorbar(plot_filtered['mean_zscore_gemme'], plot_filtered['mean_zscore_ddG'],
            xerr=plot_filtered['se_zscore_gemme'], yerr=plot_filtered['se_zscore_ddG'],
            fmt='none', ecolor='gray', elinewidth=1, capsize=2, zorder=1)

# Scatter plot with points sized by 'size' and colored by 'ID'
scatter = sns.scatterplot(data=plot_filtered, x='mean_zscore_gemme', y='mean_zscore_ddG', hue='ID',
                          palette=gene_colors, size='size', sizes=(20, 80), ax=ax, zorder=2)

# Add ID labels to the points
for i in range(plot_filtered.shape[0]):
    x = plot_filtered['mean_zscore_gemme'].values[i]
    y = plot_filtered['mean_zscore_ddG'].values[i]
    current_id = plot_filtered['ID'].values[i]
    
    # Determine label alignment based on specific conditions
    if current_id in ['Pdr1', 'Erg11 DMS', 'Erg3', 'Fks', 'SdhB']:
        plt.text(x - 0.015, y - 0.065, current_id, fontsize=9, ha='right')
    elif current_id == 'Squalene epoxidase':
        plt.text(x - 0.015, y + 0.015, current_id, fontsize=9, ha='right')
    elif current_id == 'Erg3':
        plt.text(x + 0.015, y - 0.065, current_id, fontsize=9, ha='left')
    else:
        plt.text(x + 0.015, y + 0.015, current_id, fontsize=9, ha='left')

# Define custom size labels for the legend
size_labels = ['[0, 30[', '[30, 100[', '[100, 500[', '500+']
handles, labels = scatter.get_legend_handles_labels()
scatter.legend(handles[24:], size_labels, title='Number of resistance\nmutations',
               title_fontsize='10', loc='lower right', fontsize='9')

# Customize axis labels and ticks
plt.xlabel('Predicted impact from conservation', fontsize=12)
plt.ylabel('Predicted impact from destabilization', fontsize=12)
sns.despine()

# Set custom tick labels
ax.set_yticks([-1.5, -1, -0.5, 0, 0.5])
ax.set_yticklabels(["Destabilizes\nprotein", -1, -0.5, 0, "No effect\non protein\nstability"])
ax.set_xticks([-1.5, -1, -0.5, 0, 0.5])
ax.set_xticklabels(["Disrupts\nconserved\nresidues", -1, -0.5, 0, "Disrupts\nnon-conserved\nresidues"])

# Add and customize the colorbar
cbar = fig.colorbar(image, ax=ax, orientation='vertical')
cbar.ax.invert_yaxis()  # Reverse the colorbar
cbar.set_ticks([0.05, 0.95])  # Custom tick locations
cbar.set_ticklabels(['Low \nimpact', 'High \nimpact'])

# Set aspect ratio and adjust layout
ax.set_aspect('auto')
plt.tight_layout()
plt.show()
#plt.rcParams['svg.fonttype'] = 'none'
#plt.savefig(f"{working_directory}../figures/figure2_scatter_gemme_ddg.png", dpi=600, bbox_inches='tight')
#plt.savefig(f"{working_directory}../figures/figure2_scatter_gemme_ddg.svg", bbox_inches='tight')
