import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import textwrap

import numpy as np
import pandas as pd
import gseapy as gp

import streamlit as st

class Enrichr_STAGES():
    '''
    Class to run enrichr functions for STAGES
    '''
    def execute_enrichr(self, gene_dict, select_dataset, enr_pthresh=0.05, enr_showX=10):
        '''
        Parameters
        ----------
        gene_dict: dict | keys containing deg keys or user_genes and values containing list of genes to use
        select_dataset: str | one of the gene sets for enrichr
        enr_pthresh: float | pvalue to filter pathways by
        enr_showX: int | number of pathways to display from filtered dataset
        '''
        enr_significant, enr_all = {}, {}
        non_zero = {k:gene_dict[k] for k in gene_dict.keys() if len(gene_dict[k]) !=0} # in case deg sets with 0 genes were selected
        for k,v in non_zero.items():
            enr = gp.enrichr(
                gene_list=v,
                gene_sets=select_dataset,
                outdir=None,
                no_plot=True,
                cutoff=0.5
            )

            # Sort values by adjusted p-value
            data = enr.results
            data.set_index("Term", inplace=True)
            data = data.sort_values(by=['Adjusted P-value'])
            # Drop the unimportant variables
            data_truncated = data.drop("Gene_set", axis=1)
            # Filter data by adjusted p-value and rank
            data_sig = data_truncated[(data_truncated['Adjusted P-value'] < enr_pthresh)]
            # Calculate -logP
            data_sig['-logadjP'] = np.log10(data_sig['Adjusted P-value']) * (-1)

            # Sort and return
            enr_significant[k] = data_sig.sort_values(by = "-logadjP", ascending = True).tail(enr_showX)
            enr_all[k] = data
        return enr_all, enr_significant
            
    def enr_barplot(self, enr_significant, enr_useDEG=None, enr_pthresh=0.05, enr_showX=10, enr_ht=500):
        if enr_useDEG is not None: # which implies it will require DEGs
            fig = make_subplots(rows=len(enr_significant), cols=1, subplot_titles=list(enr_significant.keys()),
                                x_title="-log10 (adjusted p-value)", shared_xaxes=True,
                                vertical_spacing=0.05) # make a subplot regardless
                
            i = 1
            for k,v in enr_significant.items(): # need to rethink this part as it creates a subplot for every comparison
                wrap_enr_sig = ["<br>".join(textwrap.wrap(a, 30)) for a in v.index]
                marker_clr = "#EF553B" if "UP" in k else "#636EFA"
                fig.add_trace(go.Bar(x=v['-logadjP'], y=wrap_enr_sig, 
                                     orientation='h', marker_color=marker_clr,
                                     name="",
                                     hovertemplate="<b>%{y}</b><br>-logadjP: %{x}<br>Enriched genes: %{customdata}"),
                                     row = i, col = 1)
                i += 1
            fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
            fig.update_layout(title=f"Enriched Pathways (Top {enr_showX}), adjusted p-value < {enr_pthresh}", title_x=0.5,
                              showlegend=False,
                              yaxis={'tickmode': 'linear'},
                              font=dict(family='Arial', size=14),
                              width = 750, height = enr_ht)
            

        else:
            fig = go.Figure()
            user_df = enr_significant['user_genes']
            wrap_usergenes = ["<br>".join(textwrap.wrap(a, 50)) for a in user_df.index]
            fig.add_trace(go.Bar(x=user_df['-logadjP'], y=wrap_usergenes,
                                 customdata=user_df['Genes'],
                                 orientation='h',
                                 marker_color="#4FC04F",
                                 name="",
                                 hovertemplate="<b>%{y}</b><br>-logadjP: %{x}<br>Enriched genes: %{customdata}")
                                 )
            fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
            fig.update_layout(title=f"Enriched Pathways (Top {enr_showX}) for user-input genes<br>adjusted p-value < {enr_pthresh}", title_x=0.5,
                            showlegend=False,
                            yaxis={'tickmode': 'linear'},
                            font=dict(family='Arial', size=14),
                            width = 750, height = enr_ht)
        return fig

enr = Enrichr_STAGES()