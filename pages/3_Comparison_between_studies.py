# Web app
import streamlit as st
from st_aggrid import AgGrid
import pandas as pd
import numpy as np
import glob
import os
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

# OS and file management
from PIL import Image

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

# Function to load data for a specific study
@st.cache_data
def load_study_data(all_data, selected_study):
    
    # Filter data for the selected study
    study_info = all_data[all_data['study_id'] == selected_study]

    return study_info

# Add a title and info about the app
st.title('Summary of waste water treatment studies from Mgnify with more than 10 samples')

# Create a df with study ids and biomes
filt_columns = ['study_id', 'biomes']
studies_biomes = st.session_state.studies_data[filt_columns].drop_duplicates()
studies_biomes = studies_biomes.dropna()

# Show the df
st.subheader('List of studies and their biomes')

studies_biomes.columns = ['Study ID', 'Biomes']

st.data_editor(
    studies_biomes,
    width = 1000,
    column_config={
        "widgets": st.column_config.TextColumn(
            help= f"List of studies and their biomes"
        )
    },
    hide_index=True,
)

# Make comparison between studies
st.subheader('Comparison between studies')

# Create a selectbox to choose the two studies
selected_studies = st.multiselect("Selected studies:", studies_biomes['Study ID'].unique())

if len(selected_studies) == 0:
    st.write('Please select two studies')
elif len(selected_studies) == 1:
    st.write('Please select another study')
else:
    # Display summary for the selected studies
    st.subheader(f'Studies {selected_studies[0]} and {selected_studies[1]}')

    # Load data for the selected studies
    study1_info = load_study_data(st.session_state.studies_data, selected_studies[0])
    study1_info = study1_info.T
    study1_info = study1_info.drop(index='study_id')
    study1_info.columns = [selected_studies[0]]

    study2_info = load_study_data(st.session_state.studies_data, selected_studies[1])
    study2_info = study2_info.T
    study2_info = study2_info.drop(index='study_id')
    study2_info.columns = [selected_studies[1]]

    # Merge the study info dfs
    studies_info = pd.merge(study1_info, study2_info, left_index=True, right_index=True)

    # Replace indices and column names
    studies_info.index = ['Description', 'Bioproject ID', 'Centre name', 'Number of samples',
                    'Biome','Experiment type', 'Pipeline version', 'Instrument platform',
                    'Sampling country']
    
    st.data_editor(
    studies_info,
    width = 2000,
    column_config={
        "widgets": st.column_config.TextColumn(
            f"Summary {selected_studies[0]} and {selected_studies[1]}",
            help= f"Summary of the studies {selected_studies[0]} and {selected_studies[1]}"
        )
        },
        hide_index=False,
        )
    
    # Load and preprocess abundance tables for the selected studies and store in independent variables
    abund_table1 = load_abund_table(selected_studies[0])
    abund_table1 = preprocess_abund_table(abund_table1)
    abund_table2 = load_abund_table(selected_studies[1])
    abund_table2 = preprocess_abund_table(abund_table2)

    # Merge the abundance tables, keeping only the phyla that are present in both studies
    abund_table_merged = pd.merge(abund_table1, abund_table2, on='phylum', how='inner')

    # Transpose the abundance table and keep it as a df
    abund_table_merged = abund_table_merged.set_index('phylum').T
    #abund_table_merged.T

    # Set the first row as the header
    new_header = abund_table_merged.iloc[0]
    abund_table_merged = abund_table_merged[1:]
    abund_table_merged.columns = new_header
    
    # Show the merged abundance table
    st.subheader(f"Abundance table for {selected_studies[0]} and {selected_studies[1]}")
    # AgGrid(abund_table_merged, editable=True)

    # Create a list of the samples from the merge abundance table
    samples_list = abund_table_merged.index.tolist()

    # Create a df with the samples and their study ids
    samples_df = pd.DataFrame(samples_list, columns=['sample_id'])

    # Add the study ids to the samples df
    for sample in samples_list:
        if sample in abund_table1.columns:
            samples_df.loc[samples_df['sample_id'] == sample, 'study_id'] = selected_studies[0]
        else:
            samples_df.loc[samples_df['sample_id'] == sample, 'study_id'] = selected_studies[1]

    # Create a pca plot with the merged abundance table
    st.subheader('PCA plot')
    st.write('The PCA plot shows the distribution of the samples from the two selected studies based on the abundance of the phyla present in both studies.')

    # Create a PCA object
    pca = PCA(n_components=2)

    # Standardize the data
    abund_values_std = StandardScaler().fit_transform(abund_table_merged)

    # Fit the PCA object to the standardized data
    pca.fit_transform(abund_values_std)

    # Transform the standardized data using the fitted PCA object
    pca_data = pca.transform(abund_values_std)

    # Create a df with the PCA data
    pca_df = pd.DataFrame(data = pca_data, columns = ['PC1', 'PC2'])

    # Add the study ids to the PCA df
    pca_df['study_id'] = samples_df['study_id'] 

    # Get biomes for the selected studies
    biome_study1 = studies_biomes.loc[studies_biomes['Study ID'] == selected_studies[0], 'Biomes']
    biome_study2 = studies_biomes.loc[studies_biomes['Study ID'] == selected_studies[1], 'Biomes']

    # Define conditions and values for filling the biome column
    conditions = [
        (pca_df['study_id'] == selected_studies[0]),  
        (pca_df['study_id'] == selected_studies[1]) 
    ]

    values = [biome_study1, biome_study2]

    # Create a new column based on conditions
    pca_df['biome'] = np.select(conditions, values, default='Unknown')
    
    # Add biome in the study id column 
    pca_df['study_id'] = pca_df['study_id'].str.cat(pca_df['biome'], sep=' - ')

    # Create a plotly figure
    pca_plot = px.scatter(pca_df, x='PC1', y='PC2', opacity=0.8, color='study_id', 
                      hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Plotly)

    # Add title and axis labels
    pca_plot.update_layout(
        xaxis=dict(
            title='PC1',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        yaxis=dict(
            title='PC2',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        legend_title_text ='Study ID'
    )

    # Show pca plot
    st.plotly_chart(pca_plot, use_container_width=True)

    # Clustering plot
    # Run cluster analysis with k-means
    st.subheader('Clustering analysis')
    st.write('The clustering analysis was performed using the kmeans algorithm. The number of clusters was set to 2, based on the number of studies selected.')

    # Specify the number of clusters (k)
    k = 2

    # Fit K-Means model
    kmeans = KMeans(n_clusters=k)
    clusters = kmeans.fit_predict(abund_values_std)

    # Add the cluster labels to the PCA df
    pca_df['Cluster'] = clusters

    # Set column type to string
    pca_df['Cluster'] = pca_df['Cluster'].astype(str)

    # Create a plotly figure
    clust_plot = px.scatter(pca_df, x='PC1', y='PC2', opacity=0.8, color='Cluster',
                            color_discrete_sequence=px.colors.qualitative.Plotly)
    
    # Add title and axis labels
    clust_plot.update_layout(
        xaxis=dict(
            title='PC1',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        yaxis=dict(
            title='PC2',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        legend_title_text ='Cluster'
    )

    # Show the plot
    st.plotly_chart(clust_plot, use_container_width=True)

    # Create PcoA plot
    st.subheader('PCoA plot with Bray-Curtis distance')
    st.write('The PCoA plot shows the distribution of the samples from the two selected studies based on the Bray-Curtis distance between them.')
    
    # Extract OTU names and sample names as numpy arrays
    otu_names = list(abund_table_merged.columns.values)
    sample_names = list(abund_table_merged.index.values)

    # Convert abundance table to numpy array
    abund_table_merged_mat = abund_table_merged.to_numpy()

    # Obtain bray-curtis distance matrix
    bc_mat = beta_diversity("braycurtis", abund_table_merged_mat, sample_names)

    # Replace NaN values with 0
    bc_mat = np.nan_to_num(bc_mat.data, nan=0.0)
    
    # Run PCoA
    bc_pcoa = pcoa(bc_mat)

    # Extract the data to plot the PcoA
    bc_pcoa_data = pd.DataFrame(data = bc_pcoa.samples[['PC1', 'PC2']])
    
    # Reset index
    bc_pcoa_data = bc_pcoa_data.reset_index(drop=True)

    # Add the study ids to the PCA df
    bc_pcoa_data['study_id'] = samples_df['study_id']

    # Define conditions and values for filling the biome column
    conditions = [
        (bc_pcoa_data['study_id'] == selected_studies[0]),  
        (bc_pcoa_data['study_id'] == selected_studies[1]) 
    ]

    # Create a new column based on conditions
    bc_pcoa_data['biome'] = np.select(conditions, values, default='Unknown')
    
    # Add biome in the study id column 
    bc_pcoa_data['study_id'] = bc_pcoa_data['study_id'].str.cat(bc_pcoa_data['biome'], sep=' - ')

    # Create a plotly figure
    pcoa_plot = px.scatter(bc_pcoa_data, x='PC1', y='PC2', opacity=0.8, color='study_id', 
                      hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Plotly)
    
    # Add title and axis labels
    pcoa_plot.update_layout(
        xaxis=dict(
            title='PC1',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        yaxis=dict(
            title='PC2',
            tickfont=dict(size=18),
            titlefont=dict(size=20)
        ),
        legend_title_text ='Study ID'
    )

    # Show pcoa plot
    st.plotly_chart(pcoa_plot, use_container_width=True)

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')
    


