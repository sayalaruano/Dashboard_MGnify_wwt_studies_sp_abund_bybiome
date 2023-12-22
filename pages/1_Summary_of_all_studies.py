# Web app
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

# OS and file management
from PIL import Image
import os
import glob

# General options 
im = Image.open("img/favicon.ico")
st.set_page_config(
    page_title="Mgnify waste water treatment studies and samples summary",
    page_icon=im,
    layout="wide",
)

# Attach customized ccs style
with open('style.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

# Function to load abundance table for a specific study
@st.cache_data
def load_abund_table(selected_study):
    # Set the file name 
    folder_path = f"Abundance_tables/{selected_study}/"
    filename = os.path.join(folder_path, f"{selected_study}*phylum_taxonomy*.tsv")

    print(filename)
    # Use glob to find files that match the pattern
    file = glob.glob(filename)
    print(file)

    if not file:
        print(f"No files found for the study '{selected_study}' in folder '{folder_path}'.")
        return None
    
    # Take the first file found (you can modify this if you expect multiple files)
    file_path = file[0]
    
    # Load abundance table for the study
    abund_table = pd.read_csv(file_path, sep='\t')
    
    return abund_table

# Function to preprocess abundance table for a specific study
def preprocess_abund_table(abund_table):
    # Delete kingdom and superkingdom columns
    if 'superkingdom' in abund_table.columns:
        abund_table = abund_table.drop(columns=['kingdom', 'superkingdom'])
    else:
        abund_table = abund_table.drop(columns=['kingdom'])

    # Delete rows with unassigned phylum
    abund_table = abund_table[abund_table['phylum'] != 'Unassigned']

    # Delete rows that contain 'Candidatus' or 'candidate'
    abund_table = abund_table[~abund_table['phylum'].str.contains('Candidatus')]
    abund_table = abund_table[~abund_table['phylum'].str.contains('candidate')]
    
    # Reset index
    abund_table = abund_table.reset_index(drop=True)
    
    return abund_table

# Add a title and info about the app
st.title('Summary of waste water treatment studies from Mgnify with more than 10 samples')
st.header('Plots to summarize all studies')

# Create plot the number of studies per study
st.subheader("Number of samples per study")

plot1 = px.histogram(st.session_state.studies_data, x="n_samples", nbins=10, opacity=0.8).update_layout(
                    xaxis=dict(title='Number of samples',tickfont=dict(size=18), titlefont=dict(size=20)),
                    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20))).update_xaxes(
                    showgrid=False).update_yaxes(showgrid=False)

st.plotly_chart(plot1, use_container_width=True)

# Create a pie plot for the sampling countries of the studies
st.subheader("Sampling countries")

plot2 = px.pie(values=st.session_state.studies_data["sampling_country"].value_counts(),
                names=st.session_state.studies_data["sampling_country"].value_counts().index, 
                opacity=0.8).update_traces(textposition='inside', textinfo='value', 
                                           insidetextfont=dict(size=18)).update_layout(legend=dict(font=dict(size=20)))

st.plotly_chart(plot2, use_container_width=True)

# Create bar plot for the data types of the studies
st.subheader("Data types")

plot3 = px.histogram(st.session_state.studies_data, x="experiment_type", opacity=0.8, color='experiment_type').update_layout(
                    xaxis=dict(title='Number of samples',tickfont=dict(size=18), titlefont=dict(size=20)),
                    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20)), 
                    showlegend=False).update_xaxes(showgrid=False).update_yaxes(showgrid=False)

st.plotly_chart(plot3, use_container_width=True)

# Create pie plot for the biomes of the studies
st.subheader("Biomes")

plot4 = px.pie(values=st.session_state.studies_data["biomes"].value_counts(),
                names=st.session_state.studies_data["biomes"].value_counts().index,
                opacity=0.8).update_traces(textposition='inside', textinfo='value', 
                                           insidetextfont=dict(size=18)).update_layout(legend=dict(font=dict(size=20)))

st.plotly_chart(plot4, use_container_width=True)

# Creat PCA plot for all studies
# st.subheader("PCA plot")

# Obtain list of studies
# studies_list = list(st.session_state.studies_data['study_id'].unique())

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')
    


