# Web app
import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

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
st.title('Summary and EDA of waste water treatment studies from Mgnify')
st.header('Plots to summarize all studies')

# Create plot the number of studies per study
st.subheader("Number of samples per study")

hist_samples = px.histogram(st.session_state.studies_data, x="n_samples", 
                            nbins=10, opacity=0.8, color_discrete_sequence=px.colors.qualitative.Plotly)

hist_samples.update_layout(
    xaxis=dict(title='Number of samples', tickfont=dict(size=18), titlefont=dict(size=20)),
    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20))
    ).update_xaxes(
        showgrid=False
    ).update_yaxes(
        showgrid=False)

st.plotly_chart(hist_samples, use_container_width=True)

# Create a pie plot for the sampling countries of the studies
st.subheader("Sampling countries")

pie_plot_countries = px.pie(values=st.session_state.studies_data["sampling_country"].value_counts(),
                names=st.session_state.studies_data["sampling_country"].value_counts().index, 
                opacity=0.8, color_discrete_sequence=px.colors.qualitative.Plotly)

pie_plot_countries.update_traces(
    textposition='inside',
    textinfo='value',
    insidetextfont=dict(size=18)
    ).update_layout(
        legend_title=dict(text='Country', font=dict(size=19)),
        legend=dict(font=dict(size=17))
    )

st.plotly_chart(pie_plot_countries, use_container_width=True)

# Create bar plot for the data types of the studies
st.subheader("Data types")

bar_plot_dtypes = px.histogram(st.session_state.studies_data, x="experiment_type", opacity=0.8, 
                               color='experiment_type', color_discrete_sequence=px.colors.qualitative.Plotly)

bar_plot_dtypes.update_layout(
    xaxis=dict(title='Number of samples', tickfont=dict(size=18), titlefont=dict(size=20)),
    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20)),
    showlegend=False
    ).update_xaxes(
        showgrid=False
    ).update_yaxes(
        showgrid=False)

st.plotly_chart(bar_plot_dtypes, use_container_width=True)

# Create pie plot for the biomes of the studies
st.subheader("Biomes")

pie_plot_biomes = px.pie(values=st.session_state.studies_data["biomes"].value_counts(),
                names=st.session_state.studies_data["biomes"].value_counts().index,
                opacity=0.8, color_discrete_sequence=px.colors.qualitative.Plotly)

pie_plot_biomes.update_traces(
    textposition='inside',
    textinfo='value',
    insidetextfont=dict(size=18)
    ).update_layout(
        legend_title=dict(text='Biome', font=dict(size=18)),
        legend=dict(font=dict(size=16))
    )

st.plotly_chart(pie_plot_biomes, use_container_width=True)

# Creat PCA plots for all studies
# Load the merged abundance table
merged_df = pd.read_csv('Abundance_tables/merged_all_abund_tables.csv', index_col=0)

# Extract the study_id column and remove it from the DataFrame
study_id = merged_df['study_id']
merged_df = merged_df.drop(columns=['study_id'])

# Transpose the DataFrame
merged_df_transp = merged_df.T

# Create a PCA object
pca = PCA(n_components=2)

# Standardize the data
abund_values_std = StandardScaler().fit_transform(merged_df)

# Fit the PCA object to the standardized data
pca.fit_transform(abund_values_std)

# Transform the standardized data using the fitted PCA object
pca_data = pca.transform(abund_values_std)

# Create a df with the PCA data
pca_df = pd.DataFrame(data = pca_data, columns = ['PC1', 'PC2'], index=merged_df_transp.columns)

# Add study_id column
pca_df['study_id'] = study_id

# Add biome column to the PCA df 
pca_df = pca_df.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

# Add biome in the study id column 
pca_df['study_id'] = pca_df['study_id'].str.cat(pca_df['biomes'], sep=' - ')

# Plot PCA colored by biome
st.subheader("PCA plot of the analyses from all studies colored by biome")

pca_biome = px.scatter(pca_df, x='PC1', y='PC2', opacity=0.8, color='biomes',
                       hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Plotly)

# Add title and axis labels
pca_biome.update_traces(
    marker=dict(size=6)
    ).update_layout(
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
    legend_title=dict(text='Biome', font=dict(size=20)),
    legend=dict(font=dict(size=16))
)

st.plotly_chart(pca_biome, use_container_width=True)

# Create PCoA plots for all studies
# Extract analyses names as numpy arrays 
analyses_names = list(merged_df.index.values)

# Convert abundance table to numpy array
abund_table_mat = merged_df.to_numpy()

# Obtain bray-curtis distance matrix
bc_mat = beta_diversity("braycurtis", abund_table_mat, analyses_names)

# Replace NaN values with 0
bc_mat = np.nan_to_num(bc_mat.data, nan=0.0)

# Run PCoA
bc_pcoa = pcoa(bc_mat)

# Extract the data to plot the PCoA
bc_pcoa_data = pd.DataFrame(data = bc_pcoa.samples, columns = ['PC1', 'PC2'])

# Reset index
bc_pcoa_data = bc_pcoa_data.reset_index(drop=True)

# Add analyses names as index
bc_pcoa_data.index = merged_df_transp.columns

# Add study_id column to the PCoA df
bc_pcoa_data['study_id'] = study_id

# Add biome column to the PCoA df 
bc_pcoa_data = bc_pcoa_data.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

# Add biome in the study id column 
bc_pcoa_data['study_id'] = bc_pcoa_data['study_id'].str.cat(bc_pcoa_data['biomes'], sep=' - ')

# Plot PCoA colored by biome
st.subheader("PCoA plot (Bray Curtis distance) of the analyses from all studies colored by biome")

pcoa_biome = px.scatter(bc_pcoa_data, x='PC1', y='PC2', opacity=0.8, color='biomes',
                        hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Plotly)

# Add title and axis labels
pcoa_biome.update_traces(
    marker=dict(size=6)
    ).update_layout(
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
    legend_title=dict(text='Biome', font=dict(size=20)),
    legend=dict(font=dict(size=16))
)

st.plotly_chart(pcoa_biome, use_container_width=True)

# Plot PCA colored by study id and biome
st.subheader("PCA plot of the analyses from all studies colored by study id and biome")

pca_studyid_biome = px.scatter(pca_df, x='PC1', y='PC2', opacity=0.8, color='study_id', 
                       hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Dark24)

# Add title and axis labels
pca_studyid_biome.update_traces(
    marker=dict(size=6)
    ).update_layout(
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
    legend_title=dict(text='Study ID - Biome', font=dict(size=16)),
    legend=dict(font=dict(size=12))
)

st.plotly_chart(pca_studyid_biome, use_container_width=True)

# Plot PCoA colored by study id and biome
st.subheader("PCoA plot (Bray Curtis distance) of the analyses from all studies colored by study id and biome")

pcoa_studyid_biome = px.scatter(bc_pcoa_data, x='PC1', y='PC2', opacity=0.8, color='study_id',
                                hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Dark24)

# Add title and axis labels
pcoa_studyid_biome.update_traces(
    marker=dict(size=6)
    ).update_layout(
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
    legend_title=dict(text='Study ID - Biome', font=dict(size=16)),
    legend=dict(font=dict(size=12))
)

st.plotly_chart(pcoa_studyid_biome, use_container_width=True)

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')
    


