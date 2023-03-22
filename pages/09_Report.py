import streamlit as st
import streamlit.components.v1 as components

import jinja2
import pdfkit
import base64
from io import BytesIO

import plotly
import matplotlib
from PIL import Image
from datetime import datetime

from helper_functions.downloads import file_downloads

st.header("Report Generation")
# Get time
today = datetime.now()
dt_string = today.strftime("%d %B %Y %I:%M:%S %p")

# Converting all the plots to output to buffer with BytesIO
cmatrix = file_downloads.plot_to_bytes(st.session_state['corr_matrix_plot'], graph_module="plotly", format="png")
volplot = file_downloads.plot_to_bytes(st.session_state['volcano_plots'][0], graph_module="pyplot", format="png") # As a static plot, use the matplotlib one only
cdf = file_downloads.plot_to_bytes(st.session_state['cdf_plot'], graph_module="plotly", format="png")
barplot = file_downloads.plot_to_bytes(st.session_state['barplot'], graph_module="plotly", format="png")
clustergram = file_downloads.plot_to_bytes(st.session_state['clustergram_plot'], graph_module="pyplot", format="png")
enrichr = file_downloads.plot_to_bytes(st.session_state['enrichr_plots'], graph_module="plotly", format="png")
prerank = file_downloads.plot_to_bytes(st.session_state['prerank_plots'], graph_module="plotly", format="png")
string = {}
for k,v in st.session_state['string_plots'].items():
    in_memory_file = BytesIO(v)
    im = Image.open(in_memory_file)
    im = im.resize((round(im.size[0]*0.3), round(im.size[1]*0.3)))
    b = BytesIO()
    im.save(b, 'png')
    data = base64.b64encode(b.getbuffer()).decode("ascii")
    string[k] = data
pval_fmt = "adjusted p-value" if st.session_state['use_corrected_pval'] else "p-value"

templateLoader = jinja2.FileSystemLoader(searchpath="accessory_files/")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "output_report_template.html"
template = templateEnv.get_template(TEMPLATE_FILE)


outputText = template.render(date = dt_string,
                             cmatrix = cmatrix,
                             volplot = volplot,
                             cdf= cdf,
                             barplot=barplot,
                             pval_fmt=pval_fmt,
                             bar_pval=st.session_state['bar_pval'],
                             bar_fc=st.session_state['bar_fc'],
                             clustergram=clustergram,
                             geneset_enr=st.session_state['geneset_enr'],
                             enr_genedict=st.session_state['enr_genedict'],
                             enr_pthresh=st.session_state['enr_pthresh'],
                             enr_showX=st.session_state['enr_showX'],
                             enrichr=enrichr,
                             geneset_prerank=st.session_state['geneset_prerank'],
                             prerank_choose_col=st.session_state['prerank_choose_col'],
                             prerank_pthresh=st.session_state['prerank_pthresh'],
                             prerank_showX=st.session_state['prerank_showX'],
                             prerank=prerank,
                             string_dict=string)
html_file = open("STAGES_report.html", 'w')
html_file.write(outputText)
html_file.close()

pdf_out = pdfkit.from_string(outputText, False)
st.download_button("Download STAGES report as PDF here", data=pdf_out, file_name="STAGES_report.pdf", mime='application/octet-stream')