import pandas as pd
import numpy as np
import time
import math
import regex as re

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc


class FC_class():
    def comparison_finder(self, cleandict):
        comparison_regex = r"(ratio|p[\.\-value]*)[_\-\s\.](.*[_\-\s\.]vs[_\-\s\.].*)"
        comparison_dict = {}
        for k,v in cleandict.items():
            comparison = list(dict.fromkeys([re.match(pattern=comparison_regex, string=i, flags=re.I).group(2) for i in v.columns if re.match(pattern=comparison_regex, string=i, flags=re.I) is not None]))
            comparison_dict[k] = comparison
        return comparison_dict


    def nclrs(self, comparison_dict):
        combined_comparisons = set([s for v in comparison_dict.values() for s in v])
        n_comps = len(combined_comparisons)
        plotly_clrs = pc.qualitative.Plotly
        if n_comps > 10:
            colors = pc.sample_colorscale(plotly_clrs, [n/(n_comps -1) for n in range(n_comps)], colortype='tuple')
        else:
            colors = plotly_clrs[0:n_comps]
        return colors

    def log_transform(self, cleandict, comparison_dict):
        log_dict = {}
        for k,v in cleandict.items():
            new_combined = pd.DataFrame()
            comps_per_df = comparison_dict[k]
            for comp in comps_per_df:
                ratio = v.filter(regex=re.compile(f"ratio[_\-\s\.]{comp}", flags=re.I), axis=1)
                pval = v.filter(regex=re.compile(f"^(p[\.\-valuedj]*)[_\-\s\.]{comp}", flags=re.I), axis=1)
                comp_df = pd.concat([np.log2(ratio), np.log10(pval)*(-1)], axis=1)
                comp_df.columns = [f"log2FC_{comp}", f"negative_log_pval_{comp}"]
                new_combined = pd.concat([new_combined, comp_df], axis=1)
            log_dict[k] = new_combined
        return log_dict


# class RNAseq():
#     def get_ttest_stats(self, cleandict):
        

tested = FC_class()