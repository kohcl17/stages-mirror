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

ss.initialise_state({'test_fdr':'None', 'corrected_anova':None})
exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']

prep_exp = st.sidebar.expander("Pre-processing Options")
if anovadict is not None and st.session_state['corrected_anova'] is None:
    comps = tested.comparison_finder(anovadict)
    padj_mtds = {"None":None, "Bonferroni":"bonferroni", "Sidak":'sidak', 'Holm-Sidak':'holm-sidak', 'Holm':'holm', 'Simes-Hochberg':'simes-hochberg', 'Hommel':'hommel',
                    'Benjamini-Hochberg FDR':'fdr_bh', 'Benjamini-Yekutieli FDR':'fdr_by',
                    'Two-stage Benjamini-Hochberg FDR':'fdr_tsbh', 'Two-stage Benjamini-Yekutieli FDR':'fdr_tsbky'}
    test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
    test_fdr_match = padj_mtds[test_fdr]
    
    adjusted_dfs = {}
    for k,v in comps.items():
        anova_file = anovadict[k]
        adj_df_per_k = pd.DataFrame()
        for comp in v:
            comp_df = anova_file.filter(regex = comp, axis=1)
            pval_array = comp_df.filter(regex = "^pval", axis=1)
            pval_array = pval_array.sort_values(by=pval_array.columns[0], ascending = True)
            rej, corrected, alphacSidak, alphacBonf = multitest.multipletests(pval_array.iloc[:,0].to_numpy(),
                                                                              method=test_fdr_match,
                                                                              is_sorted = True)
            corrected_vals = pd.DataFrame(data = {f"adj_pval_{comp}":corrected}, index = pval_array.index)
            adj_df_per_k = pd.concat([adj_df_per_k, comp_df, corrected_vals], axis=1)

        adjusted_dfs[k] = adj_df_per_k

    ss.save_state({'test_fdr':test_fdr, 'comparisons':comps, 'corrected_anova':adjusted_dfs})
else:
    ss.save_state({'anova_dict':st.session_state['anova_dict'], 'corrected_anova':None})
