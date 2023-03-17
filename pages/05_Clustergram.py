import pandas as pd
import numpy as np
import math
import re

# from scipy.spatial.distance import pdist, squareform
# from scipy.cluster.hierarchy import dendrogram, linkage

import streamlit as st

import plotly.colors as pc
import matplotlib.pyplot as plt
import seaborn as sns

from helper_functions.session_state import ss
from helper_functions.downloads import file_downloads
from helper_functions.clustergram import genePP

ss.initialise_state({'cluster_useDEG': None,
                     'cluster_textgene':'IFIT3;IL15;DDX60;ILK;IGFLR1',
                     'clust_reset':True,
                     'clust_fc':0.0})

clust_opts = st.sidebar.expander("Clustergram options", expanded=True)
degs = st.session_state['degs']
if degs is not None: # If there were DEGs already
    cluster_useDEG = clust_opts.multiselect("Select DEGs to use in clustergram",
                                            help="Leave this field blank if you wish to input a custom set of gene names",
                                            options = list(degs.keys()),
                                            default=st.session_state['cluster_useDEG'])


    if len(cluster_useDEG) == 0: # if no DEGs selected, provide text area for input
        cluster_textgene = clust_opts.text_area("Enter custom list of genes here (if not using DEGs)",
                                                value = st.session_state['cluster_textgene'],
                                                placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
        ss.save_state({'cluster_textgene':cluster_textgene,
                       'cluster_useDEG':None})
    else: # If DEGs selected, use the cluster_useDEG and revert textgene to None
        ss.save_state({'cluster_textgene': None,
                       'cluster_useDEG':cluster_useDEG})
    
else:
    cluster_textgene = clust_opts.text_area("Enter custom list of genes here",
                                            value = st.session_state['cluster_gene'],
                                            placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
    cluster_textgene = ss.save_state({'cluster_textgene':cluster_textgene})

resetter = clust_opts.checkbox("Reset clustergram to default settings", on_change=ss.binaryswitch, args=("clust_reset", ), value=st.session_state['clust_reset'])

if resetter:
    ss.save_state({'clust_fc':0.0})
clust_fc = clust_opts.number_input("Adjust absolute fold-change here",
                                   help="The app will plot the values between the user-set range",
                                   min_value=0.0, max_value=100.0, step=0.1, value= st.session_state['clust_fc'])
get_genes = genePP.genes_used(degs=degs, cluster_useDEG=st.session_state['cluster_useDEG'], cluster_textgene=st.session_state['cluster_textgene'])
gene_vals = genePP.get_gene_vals(st.session_state['log_dict_ready'], get_genes, resetter=resetter, fc_cutoff = clust_fc)
st.write(gene_vals)