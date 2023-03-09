import pandas as pd
import numpy as np
import time
import math
import regex as re

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc

from helper_functions.preprocessing import tested
from helper_functions.session_state import ss

st.session_state.update(st.session_state)

ss.initialise_state({'test_fdr':'None'})
exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']

prep_exp = st.sidebar.expander("Pre-processing Options")
if st.session_state['file_type'] == 'Fold-Changes and P-values':
    comps = tested.comparison_finder(anovadict)
    padj_mtds = {"None":None, "Bonferroni":"bonferroni", "Sidak":'sidak', 'Holm-Sidak':'holm-sidak', 'Holm':'holm', 'Simes-Hochberg':'simes-hochberg', 'Hommel':'hommel',
                    'Benjamini-Hochberg FDR':'fdr_bh', 'Benjamini-Yekutieli FDR':'fdr_by'}
    test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
    


    ss.save_state({'test_fdr':test_fdr, 'comparisons':comps})

