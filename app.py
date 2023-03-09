import streamlit as st
import pandas as pd

from helper_functions.downloads import file_downloads
from helper_functions.preprocessing import tested
from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads

###################################################### SESSION STATES ##################################################
st.session_state.update(st.session_state)
ss.initialise_state(state_dict = {'file_type':2, 'demo':False, 'cleandict':None, 'view_df':True, 'meta_dict':None})
########################################################################################################################


st.title("STAGEs Dashboard \U0001F4CA")
################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
file_type = file_opts.radio(label="Select data type for upload", options = ['RNAseq Counts', 'Log2-Normalised Data', 'Fold-Changes and P-values'], index = st.session_state['file_type'])
ss.save_state(dict(file_type = ['RNAseq Counts', 'Log2-Normalised Data', 'Fold-Changes and P-values'].index(file_type)))
df_query = file_opts.file_uploader('Upload your file', type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True, help='Note that excel files take longer to upload')
use_demo = file_opts.checkbox("Use demo dataset", value=st.session_state['demo'])
view_df = file_opts.checkbox("View demo/uploaded gene expression dataset", value = st.session_state['view_df'])

df_dict = {} # initial, but should use cleaned_dict after uploading and qc checks
if len(df_query) != 0:
    cleandict = fileuploads.read_xfile(df_query)
    cleandict = fileuploads.capslock_genes(cleandict)
    ss.save_state({'cleandict':cleandict})
else:
    if use_demo:
        testdata = pd.read_csv("accessory_files/demo_dataframe_corrected.csv", index_col=0)
        cleandict = {"Demo":testdata}
        ss.save_state({'demo':True, 'cleandict':cleandict})
    else:
        ss.save_state({'demo':False})
        st.stop()

if view_df:
    for k,v in cleandict.items():
        st.subheader(k)
        st.dataframe(v)
    ss.save_state(dict(view_df = True))
else:
    ss.save_state(dict(view_df = False))


if file_type != "Fold-Changes and P-values":
    metadata = file_opts.file_uploader(label="Upload gene expression's metadata here", type = ['csv', 'txt', 'xlsx'], accept_multiple_files=False)
    meta_dict = fileuploads.read_xfile(metadata)
    ss.save_state({'meta_dict':meta_dict})
