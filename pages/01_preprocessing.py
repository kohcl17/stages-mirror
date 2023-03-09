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
ss.initialise_state({'file_type':0, 'demo':False, 'cleandict':None, 'view_df':True, 'meta_dict':None, 'test_fdr':0})

cleandict = st.session_state['cleandict']

prep_exp = st.sidebar.expander("Pre-processing Options")
if st.session_state['file_type'] == "Fold-Changes and P-values":
    comps = tested.comparison_finder(cleandict)
    padj_mtds = {"None":None, "Bonferroni":"bonferroni", "Sidak":'sidak', 'Holm-Sidak':'holm-sidak', 'Holm':'holm', 'Simes-Hochberg':'simes-hochberg', 'Hommel':'hommel',
                    'Benjamini-Hochberg FDR':'fdr_bh', 'Benjamini-Yekutieli FDR':'fdr_by'}
    test_fdr = st.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = st.session_state['test_fdr'])
    ss.save_state({'test_fdr':list(padj_mtds.keys()).index(test_fdr)})

