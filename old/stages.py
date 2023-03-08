#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import time
import math
import re
import base64
from PIL import Image
from io import BytesIO

import requests

import gseapy as gp
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import stats
import phik

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import plotly.express as px
import plotly.colors as pc
import matplotlib.pyplot as plt
import seaborn as sns

from date_gene import *

# What's new: recoded the data cleaning codes, cleaned up code flow to only run data cleaning once per session

st.set_page_config(page_title="STAGEs", page_icon="📊")
################################################ for df download #######################################################
def convert_df(df):
    return df.to_csv().encode('utf-8')


def to_excel(df, sheetnames = None):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    if sheetnames is not None:
        for d, i in zip(df, sheetnames):
            d.to_excel(writer, sheet_name = i)
    else:
        for d, i in zip(df, range(len(df))):
            d.to_excel(writer, sheet_name=f'Sheet {i + 1}')
    writer.save()
    processed_data = output.getvalue()
    return processed_data


def get_table_download_link(df, purpose, sheetnames = None):  # downloads without needing to reset the whole scripts
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    val = to_excel(df, sheetnames=sheetnames)
    b64 = base64.b64encode(val)
    return f'<a style = "border:1px solid #31333f33; border-radius:0.25rem; background-color:#f9f9fb; text-decoration:none; color:black; padding:0.50rem 0.75rem" href="data:application/octet-stream;base64,{b64.decode()}" download="{purpose}.xlsx">' \
           f'📥 Download {purpose} as Excel file 📥</a>'  # decode b'abc' => abc

# def zip_file(dfs, keys):
#     with zipfile.ZipFile("cleanedfiles.zip", 'w') as compress:
#         for df, k in zip(dfs, keys):
#             file = df.to_csv().encode('utf-8')
#             compress.writestr(f"cleaned_{k}.csv", file)

#     with open("cleanedfiles.zip", "rb") as fp:
#           btn = st.download_button(
#               label="Download ZIP",
#               data=fp,
#               file_name="cleanedfiles.zip",
#               mime="application/octet-stream"
#               )

def create_pdf(fig, fn, graph_module = "pyplot"):
    buf = BytesIO()
    if graph_module == 'pyplot':
        fig.savefig(buf, format = 'pdf', bbox_inches = 'tight')
    elif graph_module == 'plotly':
        fig.write_image(file = buf, format = 'pdf')
    st.download_button(
        label = f"Download {fn.replace('_', ' ')} as pdf",
        data = buf,
        file_name = f"{fn}.pdf",
        mime='application/pdf'
        )

st.title("STAGEs Dashboard \U0001F4CA")

################################################# Documentation ########################################################
documentation = st.sidebar.checkbox("Read the Docs", value=False, key='documentation')
################################################# Read the Docs #######################################################
def read_docs():
    st.subheader("Static and Temporal Analysis of Gene Expression Studies (STAGES) documentation")
    st.image("images/STAGESgraphicalabstract_v4.png",width=None) #changed
    st.markdown(
    '''
    <div style="text-align: justify">
    STAGES is an easy-to-use web tool that integrates data visualisation and pathway enrichment analysis for both static and temporal gene expression studies. STAGES is free and open to all users and there is no login requirement. The web tool works by running the Python programming language at backend to perform the data analysis and graph plotting, while the Streamlit framework is used to display the output data tables and graphs at frontend. Overall, STAGEs allow users to perform the following:
    

    1.	Plot interactive volcano plots
    2.	Filter data to identify and quantify number of differentially expressed genes based on users’ pre-defined fold change and p-value cut-offs
    3.	Pathway analysis by Enrichr against Gene Ontology (GO) Biological Processes, GO Molecular Function. GO Cellular Component, Reactome databases. Also allows pathway enrichment analysis against customised gene sets such as the Blood Transcriptomic Modules (BTM) and our in-house curated database (Vaccinomics)
    4.	GSEA pre-ranked analysis against the Reactome database, BTM and Vaccinomics
    5.	Plot clustergrams based on DEGs, genes from Enrichr pathway enrichment analysis or leading edge genes from GSEA
    6.  STRING query based on DEGs or user-input gene list.
    7.	Correlation matrix comparing transcriptomics responses between different experimental conditions
    
    All data uploaded and analysed are safe and is never stored anywhere.
    
    ### Getting Started
    To use the web tool, you will need at least one comparison file which contain:

    1.	**Gene names** on the first column
    2.	**Ratio values** (relative transcript expression versus control or baseline)
    3.	**Adjusted p-values** (or p-values)

    For the tool to be able to recognise your ratio and p-value columns, please format the ratio and p-value header columns to be parsed by underscores ( _ ):

    1.	Ratio as ratio_X_vs_Y
    2.	Adjusted p-values (or p-values) as pval_X_vs_Y,

    where X and Y are comparison variables parsed by underscores ( _ ). The X and Y variables can be time-point comparisons (e.g.ratio_hr6_vs_0) or experimental-control comparisons (e.g. ratio_drugA_vs_placebo, ratio_virus_vs_ctrl).
    
    For multiple comparisons to be analysed simultaneously, users can add more ratio and p-value columns (e.g. ratio_A_vs_Y, pval_A_vs_Y, ratio_B_vs_Y, pval_B_vs_Y). Users do not need to manually remove any other column labels and values that are present within the file.
    
    To analyse different experimental conditions across various time-points, users can upload multiple .csv or .xlsx files. However, the time-points sampled and the naming of header columns must be matched across the different experimental conditions for valid comparisons to be made.
   
    ### Browser Compatibility
    Streamlit has been demonstrated to work on Google Chrome, Firefox, Safari and Microsoft Edge (https://docs.streamlit.io/knowledge-base/using-streamlit/supported-browsers), and on MacOS, Windows and Linux (https://docs.streamlit.io/library/get-started/installation). While we have not tested STAGEs in all these browsers and OS, the web tool has been tested and evaluated in these below settings:
    |     OS         |     Version    |     Chrome          |     Firefox          |     Microsoft Edge    |     Safari    |
    |----------------|----------------|---------------------|----------------------|-----------------------|---------------|
    |     MacOS      |     Big Sur    |     96.0.4664.93    |     95.0 (64 bit)    |     not tested        |     14.1.2    |
    |     Windows    |     10         |     96.0.4664.93    |     95.0 (64 bit)    |     not tested        |     14.1.2    |


    <br> </br>
    ### Data Output from STAGEs
    #### STAGEs Dashboard
    The STAGEs dashboard displays the output data tables and charts, while the sidebar allows users to upload data files, select the data output files to render and customise graph settings. If no dataset is uploaded, an option to use a demo dataset showing the gene transcript expression levels in seronegative subjects after MERCK Ad5/HIV vaccination by Zak et al., 2012 will be given.
    Users can click on the checkbox to show uploaded demo dataset in the dashboard. Clicking the header columns will sort the numerical values in ascending/descending order.
    </div>
    ''', unsafe_allow_html=True)

    st.image("images/dashboard1.png", width=None) #keep
    st.markdown('''
    <div style="text-align: justify">

    ##### Gene Updater
    Gene Updater is incorporated into STAGEs at backend, and auto-converts old gene names and dates back into the new gene names as recommended by the HUGO Gene Nomenclature Committee (HGNC).
    For more details, please see [article here](https://doi.org/10.1038/s41598-022-17104-3).

    ##### Volcano Plot Output
    The first graph that can be rendered is the volcano plot, where users can visualise the distribution of fold change and p-values of all data points. If multiple comparisons are indicated in the uploaded dataset, the volcano plots of each comparison will be overlaid, allowing users to compare distribution of volcano plots between the different experimental conditions. <br> 
    The top 10 most upregulated and downregulated are also annotated in the volcano plots. The default settings allow all data points to be plotted. An example of the output from the demo dataset is as shown below. <br>
    To define the range of x and y-axis to be plotted, users can click on the expander for volcano plots at the side bar. In this case, we have changed the settings of the log2-transformed fold change values from -3.00 to 6.50, and the dashboard will immediately update the dashboard. <br>
    To visualise the contents of every single data-point, users can click on the checkbox on “Show interactive volcano plot” and the dashboard will show a new volcano plot output. This graph will allow users to show data characteristics upon mouseover as seen in b). <br>
    To reset the settings to default settings, users just have to check on the “Reset to default settings” checkbox.
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/Fig2.png", width=None) #changed
    st.markdown('''
    <div style="text-align: justify">

    ##### DEGs Stacked Bar Output
    For the DEGs analysis, users can define the fold-change and p-value cutoffs at the side-bar, and the corresponding stacked bar chart showing the number of upregulated and downregulated DEGs will be updated in real-time on the STAGEs dashboard. The default cutoffs are fold-change ≥ 2 and p-values < 0.05, but if required, users can adjust these parameters located at the sidebar. Users can also hover the mouse cursor over the bar charts to display the exact number of up- or down-regulated DEGs. An example of the DEGs stacked bar output from the demo dataset is as shown below:
    </div>
    ''', unsafe_allow_html=True)

    st.image("images/Fig3A.png", width=None) #changed-split Fig3
    st.markdown('''
    <div style="text-align: justify">
    
    To determine the identity of the DEGs and their respective log2 fold-change values and p-values, users can click on the **Data** tab and the data table showing the respective values will be displayed. At the bottom of the page, users can download the data as an Excel file to easy visualisation of tables. An example of the output is as displayed:
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/dashboard6.png", width=None) #keep
    st.markdown('''
    <div style="text-align: justify">
    To plot clustergrams from DEGs, users can click on the expander at the sidebar and check on the checkbox to plot DEGs in the clustergram. Next, users can specify the DEGs to plot on the clustergram with the widget: “Select DEGs to plot.” The default clustergram will then be plotted. Users can uncheck the default settings and adjust the sliders for the range of log2-fold change values to be displayed if a customised clustergram is preferred. With the example below, we have used the upregulated DEGs from day 1 of the demo set to plot the clustergram to examine the relative expression between the different time-points.
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/Fig3B.png", width=None) #changed-split Fig3
    st.markdown('''
    <div style="text-align: justify">

    ##### Enrichr Output

    DEGs can be queried against curated pathway databases such as Gene Ontology (GO), Reactome, KEGG, and HumanCyc to understand the role of DEGs in biological processes, functions and their cell localisation. Once the users selected the enrichr app to display on dashboard, the databases available for database will be displayed at the sidebar after clicking on the “Expand for Enrichr pathway analysis” expander. <br>

    First, users select a geneset as the database to query against. <br>

    Next, using the widget, users can select the DEGs to be use for analysis. Alternatively, users can also manually input the genes by selecting the “Add manually” option. <br>
    
    Finally, users can click on the checkbox on “Run Enrichr” to peform the enrichr pathway analysis. <br>

    Besides the established databases, we have also included the Blood Transcriptomic Modules (BTM) annotated by Li et al., 2014, and also a customised in-house dataset, named as Vaccinomics which curates the different vaccine signatures that are published to date. In the example below, we have selected the BTM database, and used upregulated DEGs from day 1 as the gene list for database query. The top 10 enriched pathways will then be plotted in the **Bar Chart** tab. <br>

    To obtain the data table showing the full list of pathways, p-values and the identity of DEGs that are captured in the respective pathways, users can click on the **Data** tab to display the enrichr dataframe as shown below.  The Excel file capturing the data table information can also be exported by clicking on the download hyperlink.
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/Fig4.png", width=None) #changed
    
    st.markdown('''
    <div style="text-align: justify">

    ##### Prerank Analysis Output
    Another method for pathway analysis is to use the Gene Set Enrichment Analysis, which relies on ranking of ratio values to determine leading edge genes responsible for pathway enrichment. As the pathway analysis utilises on the full list of genes in your dataset, the time taken for data analysis will be much longer than in Enrichr. <br>
    
    At this point, we also recommend users to save all the graph plots and tables in your computer and remove (or uncheck) the app rendering for DEGs and Enrichr before proceeding with prerank analysis, as this will drastically improve the analysis speed. Similar to enrichr analysis, users can select the experimental condition to perform the prerank analysis, and select the geneset database for pathway analysis. Finally, users can apply the FDR<0.05 cutoff and run the GSEA prerank to render the bar charts. As shown below, we have unchecked the apps for DEGs and enrichr, and proceeded with prerank analysis, where we analysed the day 1 for prerank analysis. We also selected to apply the FDR<0.05 to display on the charts. <br>
    
    The top 10 positively and negatively enriched pathways are presented below under the **Bar Chart** tab.

    To obtain the data table showing the full list of pathways, leading edge genes, normalised enrichment scores and FDR values users can click on the **Data** tab to display the prerank dataframe as shown below. The Excel file capturing the data table information can also be exported by clicking on the download hyperlink.
    </div>
    ''', unsafe_allow_html=True)

    st.image("images/Fig5_and_dashboard11.png", width=None) #changed

    st.markdown('''
    <div style="text-align: justify">
    <br></br>

    ##### Pathway Clustergram Output
    After pathway enrichment analysis, users can render the pathway clustergram to highlight the gene expression levels of genes involved in the respective pathways. Users can select the dataframe, and perform a copy-and-paste from the data table of the genes obtained from the enrichr or prerank analysis. The default clustergrams will then the plotted on the dashboard. <br>

    In the example below, we used the demo dataset and used the leading edge genes from the DC surface signature enrichment analysis to render the pathway clustergram. The default settings are checked. To customise the clustergram, users can simply uncheck the default settings and use the sliders to select the data range to be plotted on the clustergram.
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/dashboard12.png", width=None) #keep

    st.markdown('''
    <div style="text-align: justify">
    <br></br>

    ##### STRING Query Output
    Users may also use the STRING query function to render a STRING protein-protein interaction network with their own gene list or from selected DEGs.
    </div>
    ''', unsafe_allow_html=True)
    st.image("images/Fig7.png", width=None)
    
    st.markdown('''
    <div style="text-align: justify">
    <br></br>

    ##### Correlation Matrix Output
    The correlation matrix can be used to compare the similarities in host transcriptomics responses between different experimental conditions. The ratio values are converted to log2-transformed fold change values at backend, and the correlation matrices are generated through pairwaise correlations of the log2-transformed fold changes between the different experimental conditions. Users may also select their preferred correlation coefficient from the sidebar. <br>
    
    If the transcriptomics responses are similar, then the correlation coefficient will be close to 1 (positive correlation) or -1 (Negative correlation). In the demo example, we used STAGEs to correlated the transcriptional responses between the different time-points.
    </div>
    ''', unsafe_allow_html=True)

    st.image("images/Fig8.png", width=None) #changed
    st.markdown('''
    <div style="text-align: justify">

    ### Contributors
    These apps are jointly made by myself (Kuan Rong Chan), Clara Koh and Justin Ooi from Duke-NUS, Department of Emerging Infectious Diseases. I am also thankful for Eugenia Ong and Ayesa Syenina from VIREMICS for their constructive feedback. For more details on what we do, feel free to visit us at [kuanrongchan.com](kuanrongchan.com).
    </div>
    ''', unsafe_allow_html=True)


if documentation:
    read_docs()

################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
df_query = file_opts.file_uploader('Upload your .csv/.xlsx file', type = ['csv', 'xlsx'], accept_multiple_files=True)

df_dict = {} # initial, but should use cleaned_dict after uploading and qc checks
df_names = []

if len(df_query) != 0:
    for d in df_query:
        head, sep, tail = str(d.name).partition(".")
        if tail == 'csv':
            data = st.experimental_memo(pd.read_csv)(d, index_col=0)
            df_dict[head] = data
            df_names.append(head)

        elif tail == 'xlsx':
            x = st.experimental_memo(pd.read_excel)(d, index_col=0, sheet_name=None, engine='openpyxl')
            selected_sheet = file_opts.multiselect(label="* Select which sheet to read in", options=x.keys())
            if len(selected_sheet) != 0:
                for i in selected_sheet:
                    df_dict[i] = x[i]
                    df_names.append(i)
            else:
                st.stop()

else:
    if file_opts.checkbox("Use demo dataset", value=False):
        testdata = st.experimental_memo(pd.read_csv)("demo_dataframe_corrected.csv", index_col=0)
        testname = "Demo"
        df_dict[testname] = testdata
        df_names.append(testname)
    else:
        st.stop()

for df in df_dict.values():
    df.index = df.index.astype(str, copy=False) # expand to format actual dates from excel sheets as text
    df.index = df.index.str.upper()


def comparison_finder(cleandict):
    comparison_regex = r"(ratio|p[\.\-value]*)[_\-\s\.](.*[_\-\s\.]vs[_\-\s\.].*)"
    comparison_dict = {}
    for k,v in cleandict.items():
        comparison = list(dict.fromkeys([re.match(pattern=comparison_regex, string=i, flags=re.I).group(2) for i in v.columns if re.match(pattern=comparison_regex, string=i, flags=re.I) is not None]))
        comparison_dict[k] = comparison
    return comparison_dict

####################################################################### Volcano Plot #######################################################################
def nclrs(comparison_dict):
    combined_comparisons = set([s for v in comparison_dict.values() for s in v])
    n_comps = len(combined_comparisons)
    plotly_clrs = pc.qualitative.Plotly
    if n_comps > 10:
        colors = pc.sample_colorscale(plotly_clrs, [n/(n_comps -1) for n in range(n_comps)], colortype='tuple')
    else:
        colors = plotly_clrs[0:n_comps]
    return colors

def log_transform(cleandict, comparison_dict):
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

def volcano(
    user_log, comparison_dict,
    reset=False, xaxes = (0.0, 0.0), yaxes = (0.0),
    interactive_volcano = False
    ):
    
    xaxes = (0.0, 0.0) if reset else xaxes
    yaxes = 0.0 if reset else yaxes

    plt.style.use("ggplot")
    top10annotation, bottom10annotation = [], []
    
    unlist_comparisons = sorted(list(set([item for sublist in comparison_dict.values() for item in sublist])), reverse=True)
    colorlist = nclrs(comparison_dict = comparison_dict)
    legend_dict = {}
    for a, c in zip(unlist_comparisons, colorlist):
        legend_dict[a.replace("_"," ").replace("-", " ")] = c

    if len(user_log) == 1:
        volcano1 = go.Figure()
        fig, ax = plt.subplots()
        for k, df in user_log.items():
            comps = comparison_dict[k] # a list of comparisons made for each dataframe that the user uploads
            for i, tp in enumerate(comps):
                complabels = tp.replace("_", " ").replace("-", " ")
                hex_clr = legend_dict[complabels]

                #### selecting the required FC and pval for plotting
                pval_name = f'negative_log_pval_{tp}'
                fc_name = f'log2FC_{tp}'
                mini_df = df[[fc_name, pval_name]]
                max_y = math.ceil(max(mini_df[pval_name]))
                
                if xaxes != (0.0, 0.0) and yaxes != (0.0):
                    user_filter = mini_df[(mini_df[pval_name] <= yaxes) & (mini_df[fc_name].between(xaxes[0], xaxes[1],inclusive='both'))]
                elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                    user_filter = mini_df[(mini_df[pval_name] <= yaxes)]

                elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                    user_filter = mini_df[(mini_df[fc_name].between(xaxes[0], xaxes[1], inclusive='both'))]
                else:
                    user_filter = mini_df
                
                top_10 = user_filter.sort_values(by=fc_name, ascending=False).head(10)
                bottom_10 = user_filter.sort_values(by=fc_name, ascending=False).tail(10)
                
                bottom10annotation.append(
                    bottom_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))
                top10annotation.append(
                    top_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))

                ax.grid(visible=True, which="major", axis="both", alpha=0.5)
                plt.scatter(user_filter[fc_name], user_filter[pval_name], alpha=0.8, label=complabels, c = [hex_clr])
                plt.title("Volcano plot across comparisons", loc='center')
                plt.xlabel('log2(Fold-change)')
                plt.ylabel('-log10(p-value)')
                plt.axhline(y=0, color='r', linestyle='dashed')
                plt.axvline(x=0, linestyle='dashed')

                leg = ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, borderaxespad=0.3)

                if xaxes != (0.0,0.0):
                    plt.xlim([xaxes[0], xaxes[1]])
                else:
                    pass

                if yaxes != (0.0):
                    plt.ylim(0.0, yaxes)
                else:
                    plt.ylim(0.0, max_y)

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=user_filter[fc_name], y=user_filter[pval_name],
                                                  mode='markers',
                                                  name = complabels, hovertext=list(user_filter.index),
                                                  marker=dict(color=hex_clr, size=8), opacity=0.9,
                                                  legendgroup=tp
                                                  )
                                       )
            annotationconcat_top = pd.concat(top10annotation, axis=0)
            annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

            annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
            annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

            # in case I forget for annotation:
            # 1. get the top 10 of each comparison and rename the annotation cols to log2FC and negative_log_pval (generic)
            # 2. concat all the top 10s for each comparison along axis=0 and sort by highest log2FC
            # 3. the result is a long df of 2 cols (log2FC, neglogpval), where my xy annotation coordinates is x = log2FC (1st col), y = neglogpval (2nd col)

            for i in range(len(annotationconcat_top)):
                plt.annotate(text=annotationconcat_top.index[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # top 10
            for i in range(len(annotationconcat_bottom)):
                plt.annotate(text=annotationconcat_bottom.index[i],
                             xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # bottom 10
        names = set()
        volcano1.for_each_trace(lambda trace:trace.update(showlegend=False)if trace.name in names else names.add(trace.name))

        volcano1.update_layout(showlegend=True,
                               title="Interactive volcano across comparisons",
                               legend_title_text="Comparisons",
                               font=dict(family='Arial', size=14),
                               xaxis_title="log2(Fold-change)",yaxis_title="-log10(p-value)")

        if xaxes != (0.0,0.0):
            volcano1.update_xaxes(range=[xaxes[0], xaxes[1]])
        else:
            pass
        
    else:
        i = 1
        if len(user_log) % 2 == 0:
            nrows = math.ceil(len(user_log) / 2)
            extras = nrows*2 - len(user_log)
            volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(user_log.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)", shared_xaxes=True, shared_yaxes=True)
            v_row, v_col = 1, 1
            j = 1
            fig, axs = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey = True, figsize=(8, 7))

        else:
            nrows = math.ceil(len(user_log) / 3)
            extras = nrows*3 - len(user_log)
            volcano1 = make_subplots(rows=nrows, cols=3, subplot_titles=(list(user_log.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)", shared_xaxes=True, shared_yaxes=True)
            v_row, v_col = 1, 1
            j = 1
            fig, axs = plt.subplots(nrows=nrows, ncols=3, sharex=True, sharey=True, figsize=(8,7))

        min_x, max_x, max_y = 0,0,0
        for k, df in user_log.items():
            comps = comparison_dict[k] # a list of comparisons made for each dataframe that the user uploads
            for i, tp in enumerate(comps):
                complabels = tp.replace("_", " ").replace("-", " ")
                hex_clr = legend_dict[complabels]
                #### selecting the required FC and pval for plotting
                pval_name = f'negative_log_pval_{tp}'
                fc_name = f'log2FC_{tp}'
                mini_df = df[[fc_name, pval_name]]

                min_x = math.floor(min(mini_df[fc_name])) if math.floor(min(mini_df[fc_name])) < min_x else min_x
                max_x = math.ceil(max(mini_df[fc_name])) if math.ceil(max(mini_df[fc_name])) > max_x else max_x
                max_y = math.ceil(max(mini_df[pval_name])) if math.ceil(max(mini_df[pval_name])) > max_y else max_y

                if xaxes != (0.0, 0.0) and yaxes != (0.0):
                    user_filter = mini_df[(mini_df[pval_name] <= yaxes) & (mini_df[fc_name].between(xaxes[0], xaxes[1],inclusive='both'))]
                elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                    user_filter = mini_df[(mini_df[pval_name] <= yaxes)]

                elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                    user_filter = mini_df[(mini_df[fc_name].between(xaxes[0], xaxes[1], inclusive='both'))]
                else:
                    user_filter = mini_df

                top_10 = user_filter.sort_values(by=fc_name, ascending=False).head(10)
                bottom_10 = user_filter.sort_values(by=fc_name, ascending=True).head(10)

                bottom10annotation.append(
                    bottom_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))
                top10annotation.append(
                    top_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))

                ax = plt.subplot(nrows, 2, j) if len(log_dict) % 2 == 0 else plt.subplot(nrows, 3, j)
                ax.grid(visible=True, which="major", axis="both", alpha=0.5)
                ax.scatter(user_filter[fc_name], user_filter[pval_name], alpha=0.9, label = complabels, c = [hex_clr])
                ax.axhline(y=0, color='r', linestyle='dashed')
                ax.axvline(x=0, linestyle='dashed')
                ax.set_title(f"{k}", fontdict={'fontsize':10})

                if xaxes != (0.0,0.0):
                    ax.set_xlim([xaxes[0], xaxes[1]])
                else:
                    ax.set_xlim([min_x, max_x])

                if yaxes != (0.0):
                    ax.set_ylim([-1.0, yaxes])
                else:
                    ax.set_ylim([-1.0, max_y])

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=user_filter[fc_name], y=user_filter[pval_name],
                                                  mode='markers',
                                                  name=complabels, hovertext=list(user_filter.index),
                                                  marker=dict(color=hex_clr, size=8, opacity=0.9), legendgroup=str(i)),
                                       row=v_row, col=v_col)
                    i += 1

            annotationconcat_top = pd.concat(top10annotation, axis=0)
            annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

            annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
            annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

            for i in range(len(annotationconcat_top)):
                plt.annotate(text=annotationconcat_top.index[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # top 10
            for i in range(len(annotationconcat_bottom)):
                plt.annotate(text=annotationconcat_bottom.index[i],
                             xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # bottom 10

            top10annotation.clear()
            bottom10annotation.clear()

            i = 1
            j += 1
            v_col += 1

            if (len(user_log) % 2 == 0) and v_col > 2:
                v_col = 1
                v_row += 1
            if (len(user_log) % 2 != 0) and v_col > 3:
                v_col = 1
                v_row += 1
        
        handles = {}
        for a in fig.axes:
            hn, lb = a.get_legend_handles_labels()
            for h,l in zip(hn, lb):
                handles[l] = h
 
        leg = plt.legend(handles = handles.values(), labels = handles.keys(), bbox_to_anchor=(1.1, 1), loc='upper left', frameon=False, borderaxespad=0.3)

        if extras == 1:
            axs[nrows-1, 2].remove()
        elif extras == 2:
            axs[nrows-1, 2].remove()
            axs[nrows-1, 1].remove()
        else:
            pass
        
        fig.add_subplot(111, frame_on=False)
        plt.grid(visible=False)
        plt.tick_params(labelcolor="none", bottom=False, left=False)
        fig.suptitle("Volcano plot across comparisons", fontsize=14)
        plt.xlabel("log2(Fold-change)")
        plt.ylabel("-log10(p-value)")

        if xaxes != (0.0,0.0):
            volcano1.update_xaxes(range=[xaxes[0], xaxes[1]])
        else:
            pass
        ## removing duplicate legends
        names = set()
        volcano1.for_each_trace(lambda trace: trace.update(showlegend=False) if (trace.name in names) else names.add(trace.name))

        volcano1.update_layout(showlegend=True,
                               title="Interactive volcano across comparisons", title_x=0.5,
                               legend_title_text="Comparisons",
                               font=dict(family='Arial', size=14)
                               )
    return fig, volcano1, interactive_volcano


############################################# Stacked DEGs ############################################################
def degs(cleaned_dict, comparison_dict, pval_cutoff=0.0, fc_cutoff=0.0, u_width = 800, u_height=600):
    log2fc_cutoff = np.log2(fc_cutoff)
    ####################################### Filter DF by Pvals and FC #################################################

    def organised_df():
        log_dict = {}
        for k,v in cleaned_dict.items():
            new_combined = pd.DataFrame()
            comps_per_df = comparison_dict[k]
            for comp in comps_per_df:
                ratio = v.filter(regex=re.compile(f"ratio[_\-\s\.]{comp}", flags=re.I), axis=1)
                pval = v.filter(regex=re.compile(f"^(p[\.\-valuedj]*)[_\-\s\.]{comp}", flags=re.I), axis=1)
                comp_df = pd.concat([np.log2(ratio), pval], axis=1)
                comp_df.columns = [f"log2FC_{comp}", f"pval_{comp}"]
                new_combined = pd.concat([new_combined, comp_df], axis=1)
            log_dict[k] = new_combined
        return log_dict

    log_dict = organised_df()

    proportions, deg_dict = {}, {}
    for k,v in log_dict.items():
        comps = comparison_dict[k]
        up, down = [], []
        for cmp in comps:
            pval_name = f"pval_{cmp}"
            logfc_name = f"log2FC_{cmp}"
            filtered = v[[pval_name, logfc_name]]
            degs = filtered[(filtered[pval_name] < pval_cutoff) & ((filtered[logfc_name] > log2fc_cutoff)|(filtered[logfc_name] < -log2fc_cutoff))]
            deg_dict[f"{k}_{cmp}"] = degs
            # calculating proportion of DEGs for pie chart
            upreg_deg = degs[degs.loc[:, logfc_name] > 0]
            up.append(len(upreg_deg))
            proportions[f"{k}_upcount"] = up
            proportions[f"UP_{k}_{cmp}"] = upreg_deg

            downreg_deg = degs[degs.loc[:, logfc_name] < 0]
            down.append(len(downreg_deg))
            proportions[f"{k}_downcount"] = down
            proportions[f"DOWN_{k}_{cmp}"] = downreg_deg
            
    if len(log_dict) == 1:
        stacked1 = go.Figure()
    else:
        stacked_row = 1
        stacked_col = 1
        if len(log_dict) % 2 == 0:
            nrows = math.ceil(len(log_dict) / 2)
            stacked1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(log_dict.keys())),
                                    y_title='Number of DEGs',
                                    vertical_spacing = 0.5, shared_yaxes=True)
        else:
            nrows = math.ceil(len(log_dict) / 3)
            stacked1 = make_subplots(rows=nrows, cols=3, subplot_titles=(list(log_dict.keys())),
                                    y_title='Number of DEGs', 
                                    vertical_spacing=0.5, horizontal_spacing=0.02,
                                    shared_yaxes=True)
        
    for k,v in log_dict.items():
        if len(log_dict) == 1:
            # Stacked Bar
            stacked1.add_trace(
                go.Bar(x=comps, y=proportions[f'{k}_downcount'], name="Downregulated", marker_color="#636EFA"))
            stacked1.add_trace(
                go.Bar(x=comps, y=proportions[f'{k}_upcount'], name="Upregulated", marker_color="#EF553B"))
        else:
            # Stacked Bar
            stacked1.add_trace(
                go.Bar(x=comps, y=proportions[f'{k}_downcount'], name="Downregulated", marker_color="#636EFA",
                        legendgroup="A"),
                row=stacked_row, col=stacked_col)
            stacked1.add_trace(
                go.Bar(x=comps, y=proportions[f'{k}_upcount'], name="Upregulated", marker_color="#EF553B",
                        legendgroup="B"),
                row=stacked_row, col=stacked_col)

            stacked_col += 1
            if len(log_dict) % 2 == 0 and stacked_col > 2:
                stacked_col = 1
                stacked_row += 1
            elif len(log_dict) % 2 != 0 and stacked_col > 3:
                stacked_col = 1
                stacked_row += 1

            ## removing duplicate legends
            names = set()
            stacked1.for_each_trace(lambda trace:trace.update(showlegend=False) if (trace.name in names) else names.add(trace.name))
    
    stacked1.update_layout(showlegend=True, barmode='stack',
                           title="Number of DEGs across comparisons (based on selected cutoffs)",
                           title_x=0.5,
                           legend_title_text='DEGs:',
                           font=dict(family='Arial', size=14), width=u_width, height=u_height)
    stacked1.update_xaxes(automargin=True)


    return stacked1, proportions, deg_dict, log_dict

############################################ log2FC clustergram based on DEGs #############################################
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

    proportions_nocounts = {k:v for k,v in proportions.items() if not type(v) is int}
    deglist = []  # first list to add the selected DEGs

    for l in select_deg_dicts:
        if not resetter:
            fc_filter = proportions_nocounts[l][(proportions_nocounts[l].iloc[:,1].between(fc_cutoff[0], fc_cutoff[1],inclusive='both'))]
        else:
            fc_filter = proportions_nocounts[l]

        degs = fc_filter.index.tolist()
        deglist.append(degs)

    flattened = list(set([val for sublist in deglist for val in sublist]))
    logFC_dfs = pd.DataFrame(index = flattened)

    for k, v in log_dict.items():
        log_filter = v.filter(regex=re.compile("log2FC", flags=re.I), axis=1)
        log_filter = log_filter.add_prefix(f"{k}_")
        log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
        logFC_dfs = pd.concat([logFC_dfs, log_filter], axis=1)

    if len(flattened) > 1:
        specific_cluster = logFC_dfs.loc[flattened]

        # plot
        g = sns.clustermap(specific_cluster, cmap="vlag",
                        method='average',
                        cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
                        center=0, 
                        vmin = vminmax[0], vmax = vminmax[1],
                        z_score=None,
                        col_cluster=True,
                        yticklabels=True,
                        figsize=(width, height),
                        dendrogram_ratio=(dendrogram_r, dendrogram_c),
                        linewidths=1, linecolor='white',
                        cbar_kws = {"label": "log2FC", 'orientation':'horizontal'})

        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
        g.ax_heatmap.set_ylabel("")
        g.fig.suptitle("Clustergram from DEGs", x=0.5, y=1.02, fontsize=14, fontweight='bold')
        for _, spine in g.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_edgecolor("black")

    else:
        st.warning("Please choose more than 1 DEG or increase the log2 fold-change limit.")
        st.stop()
    
    return g


####################################################### Clustergram #################################################
def clustergram(
    log_dict,
    all_df = False, select_df = None,
    gene_list = None,
    resetter=False, fc_cutoff = (-1.2, 1.2),
    vminmax = (-2.0, 2.0),
    cbar_left = 0.96, cbar_bottom = 0.02, cbar_width = 0.15, cbar_height = 0.02,
    width = 10, height = 10,
    dendrogram_r = 0.2, dendrogram_c = 0.12
    ):
    '''
    Arguments
    ---------
    log_dict: dict | values contain log-transformed ratios with pval cols attached (pre-DEG filtering)
    all_df: bool | whether to use all dfs or not (if not, select_df should be > 0)
    select_df: list | keys from log_dict to choose which df (if multiple uploaded) to use
    gene_list: str | text containing delimited gene names to be used in clustergram
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

    filter_on = pd.DataFrame()

    if not all_df:
        if len(select_df) != 0:
            dfx = {k:dfx[k] for k in select_df}
        else:
            st.stop()

    for k, v in log_dict.items():
        log_filter = v.filter(regex=re.compile('log2FC', flags=re.I), axis=1)
        log_filter = log_filter.add_prefix(f"{k}_")
        log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
        filter_on = pd.concat([filter_on, log_filter], axis=1)


    if len(gene_list) == 1:
        clust_expand.warning(
            "Please enter more than one gene in the list and/or select more than one column to plot the clustergram.")

    elif len(gene_list) > 1:
        genes = gene_list.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
        # where the user may use multiple delimiters, convert the other delimiters to commas, and then split by comma
        gene_final = [x.upper() for x in set(genes) if x != ""] # remove spaces in between and make sure to capitalise genes
        gene_cluster = filter_on.filter(items=gene_final, axis=0)

        if not resetter:
            specific_cluster = gene_cluster[gene_cluster.iloc[:,1].between(fc_slider[0], fc_slider[1],inclusive='both')]
        else:
            specific_cluster = gene_cluster

        try:
            g = sns.clustermap(specific_cluster, cmap="vlag",
                        method='average',
                        cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
                        center=0, 
                        vmin = vminmax[0], vmax = vminmax[1],
                        z_score=None,
                        col_cluster=True,
                        yticklabels=True,
                        figsize=(width, height),
                        dendrogram_ratio=(dendrogram_r, dendrogram_c),
                        linewidths=1, linecolor='white',
                        cbar_kws = {"label": "log2FC", 'orientation':'horizontal'})

            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
            g.ax_heatmap.set_ylabel("")
            g.fig.suptitle("User-input pathway clustergram", x=0.5, y=1.02, fontsize=14, fontweight='bold')
            for _, spine in g.ax_heatmap.spines.items():
                spine.set_visible(True)
                spine.set_edgecolor("black")
        except ValueError:
            st.warning("Increase the log2FC limit or add more genes.")
            st.stop()
    return g

################################################# Enrichr ##############################################################
def select_enrichr_dataset():
    geneset_dict = {
        "Blood Transcriptomic Modules (BTM)": "BTM.gmt",
        "Reactome 2021": "Reactome.gmt",
        "Vaccinomics (In-house)": "Vaccinomics.gmt", 
        "GO Biological Process 2021": "GO_Biological_Process_2021",
        "GO Molecular Function 2021": "GO_Molecular_Function_2021",
        "GO Cellular Component 2021": "GO_Cellular_Component_2021",
        "KEGG 2021 Human": "KEGG_2021_Human",
        "KEGG 2019 Mouse":"KEGG_2019_Mouse",
        "HumanCyc 2016": "HumanCyc_2016"
    }

    # Selecting genesets (BTM or reactome) to plot from a list
    geneset = enrichr_exp.radio(label='Select a geneset', options=geneset_dict.keys(), key='enrichr')
    return geneset_dict[geneset]


def genes_used(premade_dict=None, key = None):
    gene_dict = {}

    if premade_dict is not None:
        choose_genetype = st.selectbox("Select whether to use DEGs or manually add gene list",
                                                options=["DEGs", "Add manually"], key = f"{key}_choice")
        if choose_genetype == "DEGs":
            proportion_keys = list(premade_dict.keys())
            proportion_keys = [i for i in proportion_keys if re.search("^UP|DOWN", i)]

            select_DEGs = st.multiselect("Select DEGs to use",
                                                  options=sorted(proportion_keys, key=str.casefold), key = f"{key}_deg")

            upkeys = [s for s in select_DEGs if re.search("^UP", s)]
            downkeys = [s for s in select_DEGs if re.search("^DOWN", s)]

            uplist = [premade_dict[s].index.to_list() for s in upkeys]
            downlist = [premade_dict[s].index.to_list() for s in downkeys]

            flattenedup = [val for sublist in uplist for val in sublist]  # all the upDEGs in 1 list
            flatteneddown = [val for sublist in downlist for val in sublist]  # all the downDEGs in 1 list

            gene_dict['UP'], gene_dict['DOWN'] = flattenedup, flatteneddown

        elif choose_genetype == "Add manually":
            gene_in = st.text_area(label="Input list of at least 3 genes here",
                                            help="Please use one of the following delimiters:line breaks, commas, or semicolons", key = f"{key}_manual")
            genes = gene_in.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
            # where the user may use multiple delimiters
            gene_final = [x.upper() for x in set(genes) if x != ""]

            gene_dict['User'] = gene_final
    else:
        gene_in = enrichr_exp.text_area(label="Input list of at least 3 genes here",
                                        help="Please use one of the following delimiters:line breaks, commas, or semicolons", key = f"{key}_manual")

        genes = gene_in.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
        # where the user may use multiple delimiters
        gene_final = [x.upper() for x in set(genes) if x != ""]
        gene_dict['User'] = gene_final


    return gene_dict


def execute_enrichr(gene_dict, select_dataset, w_enr = 850, h_enr = 600):
    '''
    Parameters
    ----------
    gene_dict: dict | keys containing any of UP, DOWN, USER and values containing consolidated list of genes to use
    select_dataset: str | one of the gene sets for enrichr
    h_enr = int | height of plots
    '''
    enr_significant, enr_all = {}, {}
    non_zero = {k:gene_dict[k] for k in gene_dict.keys() if len(gene_dict[k]) !=0}
    for k,v in non_zero.items():
        enr = st.experimental_memo(gp.enrichr)(
            gene_list=sorted(v),
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
        data_sig = data_truncated[(data_truncated['Adjusted P-value'] < 0.05)]
        # Calculate -logP
        data_sig['-logP'] = np.log10(data_sig['Adjusted P-value']) * (-1)

        enr_significant[k] = data_sig
        enr_all[k] = data

    if "UP" in enr_significant.keys() and "DOWN" in enr_significant.keys():
        toplot_up = enr_significant['UP'].sort_values(by = "-logP", ascending = True).tail(10)
        toplot_down = enr_significant['DOWN'].sort_values(by = "-logP", ascending = True).tail(10)
        fig = make_subplots(rows=2, cols=1, subplot_titles=["Upregulated DEGs", "Downregulated DEGs"],
                                x_title="-logP", shared_xaxes=True)
        
        fig.add_trace(go.Bar(x=toplot_up['-logP'], y=toplot_up.index,
                                orientation='h', marker_color="#EF553B"),
                        row=1, col=1)
        
        fig.add_trace(go.Bar(x=toplot_down['-logP'], y=toplot_down.index,
                                orientation='h', marker_color="#636EFA"),
                        row=2, col=1)

        
    elif "UP" in enr_significant.keys():
        toplot_up = enr_significant['UP'].sort_values(by = "-logP", ascending = True).tail(10)
        fig = go.Figure()
        fig.add_trace(go.Bar(x=toplot_up['-logP'], y=toplot_up.index, orientation='h', marker_color="#EF553B"))

        
    elif "DOWN" in enr_significant.keys():
        toplot_down = enr_significant['DOWN'].sort_values(by = "-logP", ascending = True).tail(10)
        fig = go.Figure()
        fig.add_trace(go.Bar(x=toplot_down['-logP'], y=toplot_down.index, orientation='h', marker_color="#636EFA"))

    elif "USER" in enr_significant.keys():
        toplot = data_sig.sort_values("-logP", ascending=True).tail(10)
        fig = go.Figure()
        fig.add_trace(data=go.Bar(x=toplot['-logP'], y=toplot.index, orientation='h', marker_color="#4FC04F"))
        
    fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
    fig.update_layout(title="Enriched Pathways (Top 10), adjusted p-value < 0.05", title_x=0.5,
                      showlegend=False,
                      yaxis={'tickmode': 'linear'},
                      font=dict(family='Arial', size=16),
                      width = w_enr, height = h_enr)
    
    return fig, enr_all


############################################ Prerank and Visualisation #################################################
########## Choose the prerank geneset to use #############
def select_prerank_dataset():
    geneset_dict = {
        "Blood Transcriptomic Modules (BTM)": "BTM.gmt",
        "Reactome 2021": "Reactome.gmt",
        "Vaccinomics (In-house)": "Vaccinomics.gmt",
        "GO Cellular Component 2021": "GO_Cellular_Component_2021",
        "GO Biological Process 2021": "GO_Biological_Process_2021",
        "GO Molecular Function 2021": "GO_Molecular_Function_2021",
        "KEGG 2021 Human": "KEGG_2021_Human",
        "KEGG 2019 Mouse":"KEGG_2019_Mouse",
        "HumanCyc 2016": "HumanCyc_2016"
    }

    # Selecting genesets (BTM or reactome) to plot from a list
    geneset = prnk_opts.radio(label='Select a geneset', options=geneset_dict.keys(), key="prerank")
    return geneset_dict[geneset]


######### Obtain columns for ranking ###################
def find_cols(df, comparisons):
    '''
    Parameters
    ----------
    df: dataframe | log-transformed values (log2FC, -logP)
    comparisons: list | contains comparisons of the selected df
    '''
    
    col_storage = {}
    for comp in comparisons:
        logfc_comp = df.loc[:,f"log2FC_{comp}"]
        logfc_comp = logfc_comp.reset_index()
        logfc_comp.columns = [0,1]
        logfc_comp = logfc_comp.sort_values(by=1, ascending=False)
        col_storage[comp] = logfc_comp
    return col_storage


########### Run Prerank #############################
def execute_prerank(col_dict, geneset):
    prerank_results_dict = {}
    for key, data in col_dict.items():
        running = st.experimental_memo(gp.prerank)(rnk=data,
                             gene_sets=geneset,
                             permutation_num=200,  # reduce number to speed up testing
                             outdir=None,
                             seed=123,
                             no_plot=True)
        prerank_results_dict[key] = running
    plots = []
    for key, result in prerank_results_dict.items():
        terms = result.res2d.Term
        ledge = result.res2d['Lead_genes']
        nes = result.res2d['NES']
        fdr = result.res2d['FDR q-val']

        ranked = pd.concat([terms, ledge, nes, fdr], axis=1)
        ranked = ranked.set_index("Term")
        ranked.index = [i.replace('"', "") for i in ranked.index]

        if applyfdr:
            pos_nes = ranked[(ranked["NES"] > 0) & (ranked["FDR q-val"] < 0.05)]
            neg_nes = ranked[(ranked["NES"] < 0) & (ranked["FDR q-val"] < 0.05)]
            neg_nes["negative NES"] = neg_nes["NES"] * -1
        else:
            pos_nes = ranked[ranked["NES"] > 0]
            neg_nes = ranked[ranked["NES"] < 0]
            neg_nes["negative NES"] = neg_nes["NES"] * -1

        pos_nes_sort = pos_nes.sort_values(by=['NES']).tail(10)
        pos_nes_sort.reset_index(inplace=True, names='Term')
        pos_nes_sort['direction'] = "positive"
        top_pos = len(pos_nes_sort)  # For formatting plot title

        neg_nes_sort = neg_nes.sort_values(by=['negative NES']).tail(10)
        neg_nes_sort.reset_index(inplace=True, names='Term')
        neg_nes_sort['direction'] = "negative"
        top_neg = len(neg_nes_sort)  # For formatting plot title

        pos = px.bar(pos_nes_sort, x="NES", y="Term", orientation="h",
                     title=f"Positive enrichment ({key})")
        neg = px.bar(neg_nes_sort, x="negative NES", y="Term", orientation="h",
                     title=f"Negative enrichment ({key})")

        pos.update_traces(marker_color="#EF553B")
        neg.update_traces(marker_color="#636EFA")

        pos.update_layout(title=f"Top {top_pos} positively enriched pathways",
                          xaxis_title='NES',
                          font=dict(family='Arial', size=16),
                          yaxis_automargin = True
                          )
        neg.update_layout(title=f"Top {top_neg} negatively enriched pathways",
                          xaxis_title='NES',
                          font=dict(family='Arial', size=16),
                          yaxis_automargin=True
                          )
        
        plots.append(pos)
        plots.append(neg)

    return plots, prerank_results_dict


##################################### Correlation Matrix (using original df) #########################################
def corr_matrix(log_dict, method):
    updated_df_list = []
    for k, v in log_dict.items():
        filtered = v.filter(regex="log2FC", axis=1)
        df_new = filtered.add_prefix(f"{k}_")
        df_new = df_new.loc[~df_new.index.duplicated(keep='first')]
        updated_df_list.append(df_new)

    concat_fc = pd.concat(updated_df_list, axis=1)
    concat_fc = concat_fc.dropna()
    if method != 'phik':
        concat_corr = concat_fc.corr(method=method)
    else:
        concat_corr = concat_fc.phik_matrix()

    # plot
    mask = np.triu(np.ones_like(concat_corr, dtype=bool))
    df_mask = concat_corr.mask(mask)
    z_raw = df_mask.to_numpy().round(3)
    corr_matrix = px.imshow(z_raw, x = df_mask.columns.to_list(),
                            y = df_mask.columns.to_list(),
                            color_continuous_scale = px.colors.diverging.RdBu_r,
                            color_continuous_midpoint=0.0,
                            text_auto = True,
                            labels = {"color":f"{method} correlation"})
    
    corr_matrix.update_layout(
        title_text="Comparison corrlation matrix (based on log2FC values)",
        title_x=0.5,
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        xaxis_zeroline=False,
        yaxis_zeroline=False,
        yaxis_autorange='reversed',
        template='plotly_white',
        font=dict(
            family='Arial', size=14)
    )

    # NaN values are not handled automatically and are displayed in the figure
    # So we need to get rid of the text manually
    # for i in range(len(corr_matrix.layout.annotations)):
        # corr_matrix.layout.annotations[i].font.size = 14
        # if corr_matrix.layout.annotations[i].text == '0.000':
        #     corr_matrix.layout.annotations[i].text = ""

    return corr_matrix

##################################### STRINGdb from DEGs or input list ###############################################
def string_query(gene_dict):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "highres_image"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])

    gene_choice = [item for sublist in gene_dict.values() for item in sublist]
                
    if string_exp.checkbox("Run STRING query"):
        gene_ready  = "%0d".join(gene_choice)

        params = {
            "identifiers" : gene_ready,
            "species" : 9606, # species NCBI identifier 
            "network_flavor": "confidence", # show confidence links
            "caller_identity" : "stages" # your app name
            }
        if len(gene_ready) != 0:
            response = st.experimental_memo(requests.post)(request_url, data=params)
            in_memory_file = BytesIO(response.content)
            im = Image.open(in_memory_file)

            st.image(im, use_column_width=True)
    else:
        im = None
    return im
    
##################################### Choose app to build dashboard ##################################################
choose_app = st.sidebar.multiselect("Choose an app to render in the main page 👉",
                                    options=["volcano plot", "DEGs", "enrichr", "GSEA prerank", "pathway clustergram",
                                             "STRING query", "correlation matrix"])

# data processing and prepping for downstream analysis
cleaned_dict = qc_df(df_dict)
comparison_dict = comparison_finder(cleaned_dict)
colours = nclrs(comparison_dict)
log_dict = log_transform(cleaned_dict, comparison_dict)
#####################################################

if file_opts.checkbox("Show uploaded/demo dataframe"):
    for k, v in cleaned_dict.items():
        st.markdown(f"**{k}**")
        st.info("Click anywhere within the dataframe and use Cmd+F or Ctrl+F to search the dataframe.")
        st.dataframe(v)
        # downloadzip = zip_file(list(cleaned_dict.values()), cleaned_dict.keys())

for c in choose_app:
    with st.spinner("🔨Building your dashboard 🔨"):
        time.sleep(0.25)
        if c == "volcano plot":
            st.subheader("Volcano plot analysis")
            vol_opts = st.sidebar.expander("Volcano plot options", expanded=False)
            reset = vol_opts.checkbox("Reset to default settings", value=False)
            xaxes = vol_opts.slider("Choose log2 fold-change boundaries for volcano plot",
                                        help="The app will plot the values between the user-set range",
                                        min_value=-10.0, max_value=10.0, step=0.1, value=(0.0, 0.0))
            if reset:
                xaxes = (0.0, 0.0)
            yaxes = vol_opts.slider("Choose negative log10 p-value boundaries for volcano plot",
                                        help="The app will plot the values less than or equal to the user-set value",
                                        min_value=0.0, max_value=50.0, step=0.1, value=0.0)
            interactive_volcano = vol_opts.checkbox(label="Show interactive volcano plot", value=False,
                                                        help="Facilitates gene name display on hover. This may cause lag")
            vol_plot, iplot, interactive = volcano(
                user_log=log_dict,
                comparison_dict=comparison_dict,
                reset = reset, xaxes=xaxes, yaxes=yaxes,
                interactive_volcano=interactive_volcano)
            if interactive:
                ivolcano = st.success("Plot complete!")
                time.sleep(0.25)
                ivolcano.empty()

                st.pyplot(vol_plot)
                create_pdf(vol_plot, "volcano_plot", graph_module='pyplot')

                st.plotly_chart(iplot, theme = None, use_container_width=True)
                create_pdf(iplot, "interactive_volcano_plot", graph_module="plotly")

            else:
                svolcano = st.success("Plot complete!")
                time.sleep(0.25)
                svolcano.empty()
                st.pyplot(vol_plot)
                create_pdf(vol_plot, "volcano_plot", graph_module='pyplot')

        elif c == "DEGs":
            st.subheader('Differential expression analysis')
            deg_opts = st.sidebar.expander("Expand for stacked bar chart", expanded=False)
            pval_cutoff = deg_opts.slider(label="Adjust p-value cutoff here", value=0.05, min_value=0.0, max_value=1.0, step=0.01)
            fc_cutoff = deg_opts.slider(label="Adjust fold-change cutoff here ", min_value=0.0, max_value=5.0, value=2.0, step=0.1)
            u_width = deg_opts.slider(label="Adjust plot width (in px)", min_value=300, max_value=1200, value=800, step=50)
            u_height = deg_opts.slider(label="Adjust plot height (in px)", min_value=300, max_value=1200, value=600, step=50)
            stacked1, proportions, deg_dict, log_dict = degs(cleaned_dict, comparison_dict, pval_cutoff=pval_cutoff, fc_cutoff=fc_cutoff, u_width=u_width, u_height=u_height)
            stack_success = st.success("Plot complete!")
            time.sleep(0.25)
            stack_success.empty()
            deg_graph, deg_data  = st.tabs(["Bar Chart", "Data"])
            deg_graph.plotly_chart(stacked1, theme=None, use_container_width=False)
            create_pdf(stacked1, "stacked_DEG_plot", graph_module='plotly')
            
            with deg_data:
                st.info("Users may select the gene names within this dataframe and copy them for downstream analysis.")
                for k, v in deg_dict.items():
                    st.write(f"**{k}**", v)
                st.download_button(label="Download DEGs", data=to_excel(deg_dict.values()), file_name="DEGs.xlsx")
                # st.markdown(get_table_download_link(deg_dict.values(), "DEGs"), unsafe_allow_html=True)


            postdeg = st.sidebar.expander("Expand to plot clustergram based on DEGs", expanded=False)
            with postdeg:
                deg_keys = [i for i in proportions.keys() if re.search("^UP|DOWN", i)]
                select_deg_dicts = st.multiselect("Select DEGs to plot", options=sorted(deg_keys, key=str.casefold))
                st.info("Note that you should deselect the default settings checkbox before setting your log2 fold-change to see the changes.")
                resetter = st.checkbox("Default settings", help="Do not filter by log2 fold-change cutoff", value=True, key='degbased')
                fc_slider = st.slider("Adjust log2 fold-change here", help="The app will plot the values between the user-set range",
                                        min_value=-10.0, max_value=10.0, step=0.1, value=(-2.0,2.0), key='degbasedfc')

                f_width = st.slider("Clustergram width (in inches)", min_value=4, max_value=20,step=1, value=5, key = 'degcl_w')
                f_height = st.slider("Clustergram height (in inches)", min_value=4, max_value=50, step=1, value=10, key = 'degcl_h')
                vminmax = st.slider("Adjust minimum and maximum values of the colour bar", min_value=-10.0, max_value=10.0, step = 0.5, value=(-2.0, 2.0), key='degcl_vminmax')
                cbar_left = st.slider("Adjust colourbar position from left", min_value=0.0, max_value=1.0, step=0.01, value = 0.78, key='degcl_cbarl')
                cbar_bottom = st.slider("Adjust colourbar position from bottom", min_value=0.0, max_value=1.0, step=0.01, value = 0.05, key='degcl_cbarb')
                cbar_width = st.slider("Adjust colourbar width", min_value=0.0, max_value=1.0, step=0.01, value = 0.15, key='degcl_cbarw')
                cbar_height = st.slider("Adjust colourbar height", min_value=0.0, max_value=1.0, step=0.01, value = 0.02, key='degcl_cbarh')
                dendrogram_r = st.slider("Adjust relative row dendrogram length", min_value=0.01, max_value=1.00, step=0.01, value=0.12, key='degcl_dendror')
                dendrogram_c = st.slider("Adjust relative column dendrogram height", min_value=0.01, max_value=1.00, step=0.01, value=0.08, key='degcl_dendroc')

                plot_deg_clust = st.checkbox("Plot DEGs in pathway clustergram", value=False)
            if plot_deg_clust:
                st.subheader("Clustergram from DEGs")
                deg_cl = deg_cluster(
                    proportions=proportions, log_dict=log_dict, select_deg_dicts=select_deg_dicts,
                    resetter = resetter, fc_cutoff=fc_slider,
                    vminmax = vminmax,
                    cbar_left=cbar_left, cbar_bottom=cbar_bottom, cbar_width=cbar_width, cbar_height=cbar_height,
                    dendrogram_c = dendrogram_c, dendrogram_r = dendrogram_r,
                    width = f_width, height = f_height)
                st.pyplot(deg_cl)
                create_pdf(fig = deg_cl, fn="DEG_clustergram", graph_module='pyplot')

        elif c == 'enrichr':
            enrichr_exp = st.sidebar.expander("Expand for Enrichr pathway analysis", expanded=False)
            st.subheader("Enrichr Analysis")
            enr_graph, enr_data = st.tabs(['Bar Chart', 'Data'])
            if "DEGs" in choose_app:
                select_dataset = select_enrichr_dataset()
                with enrichr_exp:
                    gene_dict = genes_used(premade_dict=proportions, key = 'enrichr')
                
            else:
                select_dataset = select_enrichr_dataset()
                with enrichr_exp:
                    gene_dict = genes_used(premade_dict=None)
                
            enr_w = enrichr_exp.slider(label="Adjust plot width (in px)", min_value=300, max_value=1200, value=750, step=50, key = 'enrichrW')
            enr_ht = enrichr_exp.slider(label="Adjust plot height (in px)", min_value=300, max_value=1200, value=600, step=50, key = 'enrichrH')
            run_enrichr = enrichr_exp.checkbox("Run Enrichr", value=False)
            if run_enrichr:
                check_DEGs = sum([len(v) for v in gene_dict.values()])
                if check_DEGs == 0:
                    st.error("No genes to conduct analysis!")
                else:
                    enr_fig, enr_res = execute_enrichr(gene_dict=gene_dict, select_dataset=select_dataset, h_enr = enr_ht, w_enr = enr_w)
                    enr_graph.info("Expand the plot to view all of the terms.")
                    enr_graph.plotly_chart(enr_fig, theme='streamlit', use_container_width=True)
                    enr_data.markdown(f"**Enrichr analysis using {select_dataset.replace('.gmt', '')}**")
                    for k,v in enr_res.items():
                        if k == "USER":
                            enr_data.write("Enrichr with user-input genes")
                        else:
                            enr_data.write(f"Enrichr with {k.lower()}regulated DEGs")
                        enr_data.dataframe(v)
                    enr_data.download_button(label="Download Enrichr results", data=to_excel(list(enr_res.values())),
                                            file_name=f"Enrichr_{select_dataset.replace('.gmt', '')}.xlsx")

        elif c == 'GSEA prerank':
            prnk_opts = st.sidebar.expander("Expand for Prerank Analysis", expanded=False)
            select_df = prnk_opts.selectbox("Select your dataset to use", options=cleaned_dict.keys())

            use_tps = prnk_opts.multiselect("Select comparisons for prerank", options=comparison_dict[select_df])

            prerank_dataset = select_prerank_dataset()  # same datasets as enrichr
            prerank_cols = find_cols(log_dict[select_df], use_tps)
            applyfdr = prnk_opts.checkbox("Apply FDR < 0.05 cutoff to the bar plots", value=False)
            run_prerank = prnk_opts.checkbox("Run GSEA prerank", value=False)
            if run_prerank:
                st.subheader("GSEA Prerank Analysis")
                prerank_graph, prerank_data = st.tabs(["Bar Chart", "Data"])
                prnk_plots, prnk_results = execute_prerank(prerank_cols, prerank_dataset)
                _ = [prerank_graph.plotly_chart(i, theme = 'streamlit', use_container_width = True) for i in prnk_plots]
                for k,v in prnk_results.items():
                    prerank_data.markdown(f"**{k}**")
                    prerank_data.dataframe(v.res2d)
                
                prerank_data.download_button("Download prerank results", 
                                             data = to_excel([v.res2d for v in prnk_results.values()], sheetnames=prnk_results.keys()),
                                             mime='application/octet-stream', 
                                             file_name = f"{select_df}_prerank.xlsx")


        elif c == "pathway clustergram":
            st.subheader("Pathway Clustergram")
            clust_expand = st.sidebar.expander("Expand for user-input pathway clustergram", expanded=False)
            with clust_expand:
                all_df = st.checkbox("Use all dataframes", value = False)
                select_df = st.multiselect("Select dataframes to use", options = sorted(log_dict.keys()))
                gene_input = st.text_area(label="Input list of genes here", help="Please use one of the following delimiters: line breaks, commas, or semicolons")
                st.info("Note that you should deselect the default settings checkbox before setting your log2 fold-change to see the changes.")
                resetter = st.checkbox("Default settings", help="Do not filter by log2 fold-change cutoff", value=True, key='normalcl')
                fc_slider = st.slider("Adjust log2 fold-change here", help="The app will plot the values between the user-set range",
                                        min_value=-10.0, max_value=10.0, step=0.1, value=(-2.0,2.0), key='normalfc')

                f_width = st.slider("Clustergram width (in inches)", min_value=4, max_value=20,step=1, value=5, key = 'cl_w')
                f_height = st.slider("Clustergram height (in inches)", min_value=4, max_value=50, step=1, value=10, key = 'cl_h')
                vminmax = st.slider("Adjust minimum and maximum values of the colour bar", min_value=-10.0, max_value=10.0, step = 0.5, value=(-2.0, 2.0), key='cl_vminmax')
                cbar_left = st.slider("Adjust colourbar position from left", min_value=0.0, max_value=1.0, step=0.01, value = 0.78, key='cl_cbarl')
                cbar_bottom = st.slider("Adjust colourbar position from bottom", min_value=0.0, max_value=1.0, step=0.01, value = 0.05, key='cl_cbarb')
                cbar_width = st.slider("Adjust colourbar width", min_value=0.0, max_value=1.0, step=0.01, value = 0.15, key='cl_cbarw')
                cbar_height = st.slider("Adjust colourbar height", min_value=0.0, max_value=1.0, step=0.01, value = 0.02, key='cl_cbarh')
                dendrogram_r = st.slider("Adjust relative row dendrogram length", min_value=0.01, max_value=1.00, step=0.01, value=0.12, key='cl_dendror')
                dendrogram_c = st.slider("Adjust relative column dendrogram height", min_value=0.01, max_value=1.00, step=0.01, value=0.08, key='cl_dendroc')
                plot_clust = st.checkbox("Plot pathway clustergram", value=False, key='normalcl_plot')

            if plot_clust:
                normal_cl = clustergram(
                    log_dict = log_dict, all_df = all_df, select_df=select_df, gene_list = gene_input,
                    resetter = resetter, fc_cutoff=fc_slider,
                    vminmax = vminmax,
                    cbar_left=cbar_left, cbar_bottom=cbar_bottom, cbar_width=cbar_width, cbar_height=cbar_height,
                    dendrogram_c = dendrogram_c, dendrogram_r = dendrogram_r,
                    width = f_width, height = f_height)
                st.pyplot(normal_cl)

        elif c == "correlation matrix":
            corr_exp = st.sidebar.expander("Expand for correlation matrix", expanded=False)
            method = corr_exp.selectbox("Choose the correlation coefficient to use", options=['pearson', 'kendall', 'spearman', 'phik'], format_func = lambda x: x.title())
            st.subheader("Correlation matrix")
            mtx = corr_matrix(log_dict, method = method)
            st.info("If the plot is too small, please hover over the plot and click the expand button on the top right corner of the plot.")
            st.plotly_chart(mtx, theme = None, use_container_width=False)

        elif c == "STRING query":
            string_exp = st.sidebar.expander("Expand for STRING query", expanded=False)
            st.subheader("STRING Interaction Network")
            if "DEGs" in choose_app:
                with string_exp:
                    str_genes = genes_used(premade_dict=proportions, key = "string")
                string_fig = string_query(str_genes)
            else:
                with string_exp:
                    str_genes = genes_used(premade_dict=None)
                    string_fig = string_query(str_genes)

