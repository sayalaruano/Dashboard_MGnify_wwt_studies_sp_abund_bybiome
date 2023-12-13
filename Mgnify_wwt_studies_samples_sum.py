import streamlit as st
import pandas as pd
from st_aggrid import AgGrid

# Function to load data for a specific study
@st.cache_data
def load_study_data(all_data, selected_study):
    # Filter data for the selected study
    study_info = all_data[all_data['study_id'] == selected_study]
    
    # Load sample information for the study
    sample_info = pd.read_csv(f"Samples_metadata/{selected_study}/{selected_study}_samples_metadata.csv")
    
    return study_info, sample_info

# Function to load all data
@st.cache_data
def load_all_data():
    all_data = pd.read_csv("Mgnify_studies_wwt_shot_metag_assembly_more10samples.csv")
    return all_data

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

# Load all data
all_data = load_all_data()

# Get unique study names
studies_list = all_data['study_id'].unique()

# Add a title and info about the app
st.title('Summary of waste water treatment studies from Mgnify with more than 10 samples')

with st.expander('About this app'):
    st.write('''
    This app shows a summary of the waste water treatment studies from Mgnify with more than 10 samples. It includes informations about the studies and its samples.
    
    **Credits**
    - Developed by [Sebastián Ayala Ruano](https://sayalaruano.github.io/).
    - The data was obtained using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI).
      ''')

# Select study
st.sidebar.header('How it works?')

st.sidebar.write('Select a study in the sidebar to display its information in the main panel.')

# Create a selectbox to choose the study
selected_study = st.sidebar.selectbox("Mgnify study:", studies_list)

st.sidebar.header('Data')

st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')

st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

# Load data for the selected study
study_info, sample_info = load_study_data(all_data, selected_study)

# Display study information
st.subheader(f'Description - {selected_study}')
description = study_info['study_name'].values[0]
st.text(description)

biome = study_info['biomes'].values[0]
st.subheader('Biome')
st.text(biome)

centre_name = study_info['centre_name'].values[0]
st.subheader('Centre name')
st.text(centre_name)

sampling_country = study_info['sampling_country'].values[0]
st.subheader('Sampling country')
st.text(sampling_country)

experiment_type = study_info['experiment_type'].values[0]
st.subheader('Experiment type')
st.text(experiment_type)

n_samples = study_info['n_samples'].values[0]
pipeline_version = study_info['pipeline_version'].values[0]

col1, col2 = st.columns(2)
col1.metric("Number of samples", n_samples)
col2.metric("Pipeline version", pipeline_version)

# st.table(study_info)

# Display sample information for the selected study
st.subheader("Sample Information:")
AgGrid(sample_info, editable=True)

# Button to download the data
st.subheader("Download sample data")

sample_info_csv = convert_df(sample_info)
st.download_button(
    label="Download data as CSV",
    data=sample_info_csv,
    file_name=f'sample_info_{selected_study}.csv',
    mime='text/csv',
)

