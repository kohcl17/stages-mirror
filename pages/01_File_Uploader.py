import streamlit as st
import pandas as pd
from io import StringIO

from helper_functions.downloads import file_downloads
from helper_functions.preprocessing import tested
from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads

###################################################### SESSION STATES ##################################################
st.session_state.update(st.session_state)
ss.initialise_state(state_dict = {'file_type':'Fold-Changes and P-values', 'demo':True, 'view_df':True})
########################################################################################################################


################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
use_demo = file_opts.checkbox("Use demo dataset", value=st.session_state['demo'])
file_type = file_opts.radio(label="Select data type for upload", options = ['RNAseq Counts', 'Log2-Normalised Data', 'Fold-Changes and P-values'],
                             index = ['RNAseq Counts', 'Log2-Normalised Data', 'Fold-Changes and P-values'].index(st.session_state['file_type']))
ss.save_state(dict(file_type = file_type))
df_query = file_opts.file_uploader('Upload your file', type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True, help='Note that excel files take longer to upload')


if use_demo:
    testdata = pd.read_csv("accessory_files/demo_dataframe_corrected.csv", index_col=0)
    cleandict = {"Demo":testdata}
    ss.save_state({'demo':True, 'anova_dict':cleandict, 'expr_dict':None, 'meta_dict':None})
else:
    ss.save_state({'demo':False, 'anova_dict':None, 'expr_dict':None, 'meta_dict':None})

if len(df_query) != 0:
    with file_opts:
        cleandict = st.cache_data(experimental_allow_widgets = True)(fileuploads.read_xfile)(df_query)
    cleandict = fileuploads.capslock_genes(cleandict)

    if file_type == "Fold-Changes and P-values":
        ss.save_state({'anova_dict':cleandict, 'meta_dict':None, 'expr_dict':None})

    else:
        with file_opts:
            metadata = st.file_uploader(label="Upload gene expression's metadata here", type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True)
            if len(metadata) != 0:
                meta_dict = fileuploads.read_xfile(metadata)
                ss.save_state({'expr_dict':cleandict,'meta_dict':meta_dict, 'anova_dict':None})
            else:
                ss.save_state({'expr_dict':cleandict, 'meta_dict':None, 'anova_dict':None})

view_df = file_opts.checkbox("View demo/uploaded gene expression dataset", value = st.session_state['view_df'])

if not view_df:
    ss.save_state(dict(view_df = False))

else:
    main_expr, meta_expr, anova_expr = st.tabs(['Gene Expression Data', 'Metadata', 'Fold-Change and P-value Data'])
    exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']
    _ = {main_expr.subheader(k):main_expr.dataframe(v) for k,v in exprdict.items()} if exprdict is not None else main_expr.info("No gene expression counts uploaded")
    _ = {meta_expr.subheader(k):meta_expr.dataframe(v) for k,v in metadatadict.items()} if metadatadict is not None else meta_expr.info("No metadata uploaded")
    _ = {anova_expr.subheader(k):anova_expr.dataframe(v) for k,v in anovadict.items()} if anovadict is not None else anova_expr.info("No fold-changes and p-values uploaded")
