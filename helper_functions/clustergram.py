import pandas as pd
import numpy as np
import math
import re

import streamlit as st

import seaborn as sns
import matplotlib.pyplot as plt

from helper_functions.session_state import ss

class GeneHandler():
    def genes_used(self, degs, cluster_useDEG= None, cluster_textgene=None):
        if cluster_useDEG is not None:
            get_dict_values = [degs[s].index.to_list() for s in cluster_useDEG] # Get the deg list from the selected keys
            flatten_dict_values = list(set([item for sublist in get_dict_values for item in sublist]))
            gene_final = flatten_dict_values
        
        if cluster_textgene is not None:
            genes = cluster_textgene.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
            # where the user may use multiple delimiters, convert the other delimiters to commas, and then split by comma
            gene_final = [x.upper() for x in set(genes) if x != ""] # remove spaces in between and make sure to capitalise genes

        return gene_final
    
    def get_gene_vals(self, log_dict_ready, genes_used=None):
        compiled_logFC = pd.DataFrame()
        for k,v in log_dict_ready.items():
            selected_gene_vals = v.loc[~v.index.duplicated(keep='first')]
            selected_gene_vals = selected_gene_vals[selected_gene_vals.index.isin(genes_used)]
            logFC_selected_gene_vals = selected_gene_vals.filter(regex="log2FC", axis=1)
            logFC_selected_gene_vals = logFC_selected_gene_vals.add_prefix(f"{k} ")
            compiled_logFC = pd.concat([compiled_logFC, logFC_selected_gene_vals], axis=1)
        return compiled_logFC


genePP = GeneHandler()

def deg_cluster(
    proportions, log_dict, select_deg_dicts,
    resetter=False, fc_cutoff = (-1.2, 1.2),
    vminmax = (-2.0, 2.0),
    cbar_left = 0.96, cbar_bottom = 0.02, cbar_width = 0.15, cbar_height = 0.02,
    width = 10, height = 10,
    dendrogram_r = 0.2, dendrogram_c = 0.12
    ):
    '''
    Arguments
    ---------
    proportions: dict | values contain entries of up and downregulated DEG df
    log_dict: dict | values contain log-transformed ratios with pval cols attached (pre-DEG filtering)
    select_deg_dicts: list | keys from proportions dictionary to filter the log_dict by
    resetter: bool | whether to reset logFC cutoff slider (no filtering)
    fc_cutoff: tuple | lower and upper bounds of log2FC cutoffs
    vminmax: tuple | lower and upper bounds of colour scale
    cbar_left, cbar_bottom, cbar_width, cbar_height: float | cbar position and size
    width: int | clustergram width
    height: int | clustergram height
    dendrogram_r, dendrogram_c: float | relative row and column dendrogram length/height

    Returns
    -------
    Clustergram plot
    '''

#     proportions_nocounts = {k:v for k,v in proportions.items() if not type(v) is int}
#     deglist = []  # first list to add the selected DEGs

#     for l in select_deg_dicts:
        # if not resetter:
        #     fc_filter = proportions_nocounts[l][(proportions_nocounts[l].iloc[:,1].between(fc_cutoff[0], fc_cutoff[1],inclusive='both'))]
        # else:
        #     fc_filter = proportions_nocounts[l]

#         degs = fc_filter.index.tolist()
#         deglist.append(degs)

#     flattened = list(set([val for sublist in deglist for val in sublist]))
#     logFC_dfs = pd.DataFrame(index = flattened)

#     for k, v in log_dict.items():
#         log_filter = v.filter(regex=re.compile("log2FC", flags=re.I), axis=1)
#         log_filter = log_filter.add_prefix(f"{k}_")
#         log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
#         logFC_dfs = pd.concat([logFC_dfs, log_filter], axis=1)

#     if len(flattened) > 1:
#         specific_cluster = logFC_dfs.loc[flattened]

#         # plot
#         g = sns.clustermap(specific_cluster, cmap="vlag",
#                         method='average',
#                         cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
#                         center=0, 
#                         vmin = vminmax[0], vmax = vminmax[1],
#                         z_score=None,
#                         col_cluster=True,
#                         yticklabels=True,
#                         figsize=(width, height),
#                         dendrogram_ratio=(dendrogram_r, dendrogram_c),
#                         linewidths=1, linecolor='white',
#                         cbar_kws = {"label": "log2FC", 'orientation':'horizontal'})

#         g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
#         g.ax_heatmap.set_ylabel("")
#         g.fig.suptitle("Clustergram from DEGs", x=0.5, y=1.02, fontsize=14, fontweight='bold')
#         for _, spine in g.ax_heatmap.spines.items():
#             spine.set_visible(True)
#             spine.set_edgecolor("black")

#     else:
#         st.warning("Please choose more than 1 DEG or increase the log2 fold-change limit.")
#         st.stop()
    
#     return g

# ####################################################### Clustergram #################################################
# def clustergram(
#     log_dict,
#     all_df = False, select_df = None,
#     gene_list = None,
#     resetter=False, fc_cutoff = (-1.2, 1.2),
#     vminmax = (-2.0, 2.0),
#     cbar_left = 0.96, cbar_bottom = 0.02, cbar_width = 0.15, cbar_height = 0.02,
#     width = 10, height = 10,
#     dendrogram_r = 0.2, dendrogram_c = 0.12
#     ):
#     '''
#     Arguments
#     ---------
#     log_dict: dict | values contain log-transformed ratios with pval cols attached (pre-DEG filtering)
#     all_df: bool | whether to use all dfs or not (if not, select_df should be > 0)
#     select_df: list | keys from log_dict to choose which df (if multiple uploaded) to use
#     gene_list: str | text containing delimited gene names to be used in clustergram
#     resetter: bool | whether to reset logFC cutoff slider (no filtering)
#     fc_cutoff: tuple | lower and upper bounds of log2FC cutoffs
#     vminmax: tuple | lower and upper bounds of colour scale
#     cbar_left, cbar_bottom, cbar_width, cbar_height: float | cbar position and size
#     width: int | clustergram width
#     height: int | clustergram height
#     dendrogram_r, dendrogram_c: float | relative row and column dendrogram length/height

#     Returns
#     -------
#     Clustergram plot
#     '''

#     filter_on = pd.DataFrame()

#     if not all_df:
#         if len(select_df) != 0:
#             dfx = {k:dfx[k] for k in select_df}
#         else:
#             st.stop()

#     for k, v in log_dict.items():
#         log_filter = v.filter(regex=re.compile('log2FC', flags=re.I), axis=1)
#         log_filter = log_filter.add_prefix(f"{k}_")
#         log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
#         filter_on = pd.concat([filter_on, log_filter], axis=1)


#     if len(gene_list) == 1:
#         clust_expand.warning(
#             "Please enter more than one gene in the list and/or select more than one column to plot the clustergram.")

#     elif len(gene_list) > 1:
#         genes = gene_list.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
#         # where the user may use multiple delimiters, convert the other delimiters to commas, and then split by comma
#         gene_final = [x.upper() for x in set(genes) if x != ""] # remove spaces in between and make sure to capitalise genes
#         gene_cluster = filter_on.filter(items=gene_final, axis=0)

#         if not resetter:
#             specific_cluster = gene_cluster[gene_cluster.iloc[:,1].between(fc_slider[0], fc_slider[1],inclusive='both')]
#         else:
#             specific_cluster = gene_cluster

#         try:
#             g = sns.clustermap(specific_cluster, cmap="vlag",
#                         method='average',
#                         cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
#                         center=0, 
#                         vmin = vminmax[0], vmax = vminmax[1],
#                         z_score=None,
#                         col_cluster=True,
#                         yticklabels=True,
#                         figsize=(width, height),
#                         dendrogram_ratio=(dendrogram_r, dendrogram_c),
#                         linewidths=1, linecolor='white',
#                         cbar_kws = {"label": "log2FC", 'orientation':'horizontal'})

#             g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
#             g.ax_heatmap.set_ylabel("")
#             g.fig.suptitle("User-input pathway clustergram", x=0.5, y=1.02, fontsize=14, fontweight='bold')
#             for _, spine in g.ax_heatmap.spines.items():
#                 spine.set_visible(True)
#                 spine.set_edgecolor("black")
#         except ValueError:
#             st.warning("Increase the log2FC limit or add more genes.")
#             st.stop()
#     return g