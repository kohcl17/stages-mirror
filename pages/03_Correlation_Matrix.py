import streamlit as st
from helper_functions.session_state import ss
from helper_functions.corr_matrix import cmatrix
from helper_functions.preprocessing import tested

st.session_state.update(st.session_state)

ss.initialise_state({'use_corrected_pval':False,
                     'corr_mtd': "pearson"})

corr_exp = st.sidebar.expander("Expand for correlation matrix", expanded=False)
corr_opts = ['pearson', 'kendall', 'spearman', 'phik']
log_dict = st.session_state['log_dict_ready']

if st.session_state['ready'] is not None:
    corr_mtd = corr_exp.selectbox("Choose the correlation coefficient to use", options=corr_opts, format_func = lambda x: x.title(), index = corr_opts.index(st.session_state['corr_mtd']))
    ss.save_state({'corr_mtd':corr_mtd})

    st.header("Correlation Matrix")
    mtx = cmatrix.corr_matrix(log_dict, method = st.session_state['corr_mtd'])
    ss.save_state({'corr_matrix_plot':mtx})
    st.info("If the plot is too small, please hover over the plot and click the expand button on the top right corner of the plot.")
    st.plotly_chart(st.session_state['corr_matrix_plot'], theme = None, use_container_width=False)

    st.subheader("Pre-processed Data")
    data_exp = st.expander("Expand to show pre-processed data", expanded=False)
    for k,v in st.session_state['ready'].items():
        data_exp.subheader(k)
        data_exp.dataframe(v)