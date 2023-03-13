import pandas as pd
import numpy as np
import time
import math
import regex as re
from statsmodels.stats import multitest

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc

from helper_functions.preprocessing import tested
from helper_functions.session_state import ss

st.session_state.update(st.session_state)

ss.initialise_state({'test_fdr':'None', 'corrected_anova':None, 'to_log':False, 'preprocessed_counts':None})
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

elif exprdict is not None and metadatadict is not None:
    # Insert form here instead???
    to_log = prep_exp.checkbox("Log2 transform my counts", value = st.session_state['to_log'], on_change=ss.binaryswitch, args = ('to_log', ))
    
    if to_log is True:
        logged = {k: np.log2(v) for k,v in exprdict.items()}
        ss.save_state({'preprocessed_counts':logged})
    else:
        ss.save_state({'preprocessed_counts':st.session_state['preprocessed_counts']})

st.write(st.session_state)