import streamlit as st
import pandas as pd

from helper_functions.downloads import file_downloads
from helper_functions.preprocessing import tested
from helper_functions.date_gene import qc_df

################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
file_type = file_opts.radio(label="Select data type for upload", options = ['RNAseq Counts', 'Log2-Normalised Data', 'Fold-Changes and P-values'])
df_query = file_opts.file_uploader('Upload your file', type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True, help='Note that excel files take longer to upload')

df_dict = {} # initial, but should use cleaned_dict after uploading and qc checks
if len(df_query) != 0:
    for d in df_query:
        head, sep, tail = str(d.name).partition(".")
        if tail == 'csv':
            data = pd.read_csv(d, index_col=0)
            df_dict[head] = data
        
        elif tail == 'txt':
            data = pd.read_csv(d, sep='\t', index_col=0)
            df_dict[head] = data

        elif tail == 'xlsx':
            x = pd.read_excel(d, index_col=0, sheet_name=None, engine='openpyxl')
            selected_sheet = file_opts.multiselect(label="* Select which sheet to read in", options=x.keys())
            if len(selected_sheet) != 0:
                for i in selected_sheet:
                    df_dict[f"{head}_{i}"] = x[i]
            else:
                st.stop()

else:
    if file_opts.checkbox("Use demo dataset", value=False):
        testdata = pd.read_csv("accessory_files/demo_dataframe_corrected.csv", index_col=0)
        testname = "Demo"
        df_dict[testname] = testdata
    else:
        st.stop()

for df in df_dict.values():
    df.index = df.index.astype(str, copy=False) # expand to format actual dates from excel sheets as text
    df.index = df.index.str.upper()

cleandict = qc_df(df_dict)

if file_type == "Fold-Changes and P-values":
    comps = tested.comparison_finder(cleandict)
    st.write(comps)

elif file_type == "RNAseq Counts":
    metadata = file_opts.file_uploader(label="Upload gene expression's metadata here", type = ['csv', 'txt', 'xlsx'], accept_multiple_files=False)
