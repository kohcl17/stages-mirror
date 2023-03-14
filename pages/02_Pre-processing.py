import pandas as pd
import numpy as np
import time
import math
import regex as re

# Stats modules
from statsmodels.stats import multitest
from anndata import AnnData
import decoupler as dc
import scanpy as sc
import pingouin as pg

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc
import matplotlib.pyplot as plt 

from helper_functions.preprocessing import tested, counts_pp
from helper_functions.session_state import ss

st.session_state.update(st.session_state)

ss.initialise_state({'test_fdr':'None',
                     'corrected_anova':None,
                     'to_log':False,
                     'vthresh':0,
                     'adata':None,
                     'violin_submit':False,
                     'violin1':None,
                     'violin2':None,
                     'comp_var': 0})
exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']

adjusted_dfs = {}
prep_exp = st.sidebar.expander("Pre-processing Options")

# Conditions here should mainly be
## 1. Process ANOVA data (not None)
## 2. Process Counts + metadata
## NOTE: there should not be both anova AND count data

if anovadict is not None:
    comps = tested.comparison_finder(anovadict)
    padj_mtds = {"None":None, "Bonferroni":"bonferroni", "Sidak":'sidak', 'Holm-Sidak':'holm-sidak', 'Holm':'holm', 'Simes-Hochberg':'simes-hochberg', 'Hommel':'hommel',
                    'Benjamini-Hochberg FDR':'fdr_bh', 'Benjamini-Yekutieli FDR':'fdr_by',
                    'Two-stage Benjamini-Hochberg FDR':'fdr_tsbh', 'Two-stage Benjamini-Yekutieli FDR':'fdr_tsbky'}
    test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
    ss.save_state({'test_fdr':test_fdr})
    test_fdr_match = padj_mtds[test_fdr]
    
    if test_fdr_match is not None:
        for k,v in comps.items():
            anova_file = anovadict[k]
            adj_df_per_k = pd.DataFrame()
            for comp in v:
                comp_df = anova_file.filter(regex = comp, axis=1)
                pval_col = comp_df.filter(regex = "^pval", axis=1)
                pval_col = pval_col.sort_values(by=pval_col.columns[0], ascending = True)
                pval_array = pval_col.iloc[:,0].to_numpy()
                rej, corrected, alphacSidak, alphacBonf = multitest.multipletests(pval_array,
                                                                                method=test_fdr_match,
                                                                                is_sorted = True)
                corrected_vals = pd.DataFrame(data = {f"adj_pval_{comp}":corrected}, index = pval_col.index)
                adj_df_per_k = pd.concat([adj_df_per_k, comp_df, corrected_vals], axis=1)

            adjusted_dfs[k] = adj_df_per_k
        ss.save_state({'test_fdr':test_fdr, 'comparisons':comps, 'corrected_anova':adjusted_dfs})
    
    else:
        adjusted_dfs = anovadict
        ss.save_state({'corrected_anova': adjusted_dfs, 'test_fdr':"None"})

elif exprdict is not None and metadatadict is not None: # RNAseq or microarray data
    exprdict = {k:v.groupby(v.index).mean() for k,v in exprdict.items()}
    expr_obj = exprdict[list(exprdict.keys())[0]]
    meta_obj = metadatadict[list(metadatadict.keys())[0]]

    adata = AnnData(expr_obj.T, dtype='float32') # expr data should be genes in cols, subjects in rows
    adata.obs = meta_obj # metadata should be subjects in rows, everything else in cols
    ss.save_state({"adata":adata})

    if st.session_state['file_type'] == "RNAseq Counts": # specifically RNAseq data
        split_long_violins = counts_pp.chunks(list(adata.obs_names), chunk_size=12)

        adata = st.session_state['adata']
        # Provide info here
        st.info("The violin plots From the pre-processing options in the sidebar, select a threshold value that separates the violin plots into two (narrowest point).")
        bef, aft = st.tabs(['Before pre-processing', 'After pre-processing'])
        # Create violin plot
        max_y = counts_pp.violin_maxy(adata)
        violin_thresh = prep_exp.slider(label = "Select violin plot threshold", min_value = 0, max_value= max_y,
                                  help = "This function will retain genes whose counts are above the threshold.",
                                  value=st.session_state['vthresh'])
        ss.save_state({'vthresh':violin_thresh})
        to_log = prep_exp.checkbox("Log2 transform my counts", value = st.session_state['to_log'], on_change=ss.binaryswitch, args = ('to_log', ))
        submit = prep_exp.checkbox(label = "Submit for pre-processing", value = st.session_state["violin_submit"], on_change=ss.binaryswitch, args = ('violin_submit', ))


        violin1, vaxes = counts_pp.multiviolin(adata, split_long_violins=split_long_violins)
        _ = [ax.axhline(y = st.session_state['vthresh'], color = 'red', linewidth = 1, linestyle = '--') for ax in vaxes]
        ss.save_state({'violin1':violin1})
        bef.pyplot(st.session_state['violin1'])


        if submit:
            dc.mask_features(adata, log=True, thr=st.session_state['vthresh'])
            sc.pp.filter_genes(adata, min_cells=adata.shape[0])
            violin2, vaxes2 = counts_pp.multiviolin(adata, split_long_violins=split_long_violins)
            ss.save_state({'violin2':violin2, 'adata':adata})
            aft.pyplot(st.session_state['violin2'])

            if to_log:
                adata.layers['log2'] = np.log2(st.session_state['adata'].to_df() + 1)
                logct = adata.to_df(layer = "log2").T
                ss.save_state({'adata':adata})
            else:
                unloggedct = adata.to_df().T
                ss.save_state({'adata': st.session_state['adata']})


    adata = st.session_state['adata']
    adata_obskeys = adata.obs_keys()
    comp_var = prep_exp.selectbox(label="Select variable to use for comparison", options = adata_obskeys(), index=0)
    
    # T-test requires: x, y, paired or not
