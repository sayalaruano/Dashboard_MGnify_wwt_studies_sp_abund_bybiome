# Web app
import streamlit as st
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
import pandas as pd
import numpy as np
import glob
import os
import plotly.express as px
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

@st.cache_data
# Function to load abundance table for a specific study
def load_abund_table(selected_study, phylum):
    # Set the folder name 
    folder_path = f"Abundance_tables/{selected_study}/"

    # Broad pattern to initially match files
    broad_pattern = f"{selected_study}*taxonomy*.csv"
    file_list = glob.glob(os.path.join(folder_path, broad_pattern))

    if phylum:
        # Filtering for phylum taxonomy files
        filtered_files = [f for f in file_list if 'phylum_taxonomy' in f]
    else:
        # Filtering out unwanted files (those with '_phylum_')
        filtered_files = [f for f in file_list if '_phylum_' not in f]

    # Check if the filtered list is not empty
    if filtered_files:
        filename = filtered_files[0]  # Selecting the first matching file
    else:
        print(f"No files found for the study '{selected_study}' in folder '{folder_path}'.")
        return None

    # Load abundance table for the study
    abund_table = pd.read_csv(filename, sep=',')
    
    return abund_table

# Function to preprocess abundance table for a specific study
def preprocess_abund_table(abund_table, phylum):
    # Delete NaN rows
    abund_table = abund_table.dropna(how='all')
    if phylum:
        # Delete kingdom and superkingdom columns
        if 'superkingdom' in abund_table.columns:
            abund_table = abund_table.drop(columns=['superkingdom', 'kingdom'])
        else:
            abund_table = abund_table.drop(columns=['kingdom'])
        
         # Set the phylum column as index
        abund_table = abund_table.set_index('phylum')
    
    else:
        # Delete extra taxonomic columns 
        if 'Superkingdom' in abund_table.columns:
            abund_table = abund_table.drop(columns=['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Species'])
        else:
            abund_table = abund_table.drop(columns=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Species'])
        
         # Set the genus column as index
        abund_table = abund_table.set_index('Genus')
    
    return abund_table

# Function to load data for a specific study
@st.cache_data
def load_study_data(all_data, selected_study):
    # Filter data for the selected study
    study_info = all_data[all_data['study_id'] == selected_study]

    return study_info

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

# Add info on the sidebar
# Select study
st.sidebar.header('How it works?')

st.sidebar.write('Select a taxonomic rank in the sidebar to make the comparison between the studies in the main panel.')

# Create a selectbox to choose the taxonomic rank
tax_rank = st.sidebar.selectbox(
    "Taxonomic Rank:",
    options=["Phylum", "Genus"]
)

# Add aditional info
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')

# Add a title and info about the app
st.title('Summary and EDA of waste water treatment studies from Mgnify')

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
    if tax_rank == 'Phylum':
        abund_table1 = load_abund_table(selected_studies[0], phylum=True)
        abund_table1 = preprocess_abund_table(abund_table1, phylum=True)
        abund_table2 = load_abund_table(selected_studies[1], phylum=True)
        abund_table2 = preprocess_abund_table(abund_table2, phylum=True)
    else:  # Genus
        abund_table1 = load_abund_table(selected_studies[0], phylum=False)
        abund_table1 = preprocess_abund_table(abund_table1, phylum=False)
        abund_table2 = load_abund_table(selected_studies[1], phylum=False)
        abund_table2 = preprocess_abund_table(abund_table2, phylum=False)

    # Add a row with the study ID for all samples in the current study
    study_id_row1 = pd.Series(selected_studies[0], index=abund_table1.columns)
    study_id_row2 = pd.Series(selected_studies[1], index=abund_table2.columns)

    # Add the study id row to the abundance table
    abund_table1 = pd.concat([pd.DataFrame([study_id_row1]), abund_table1])
    abund_table2 = pd.concat([pd.DataFrame([study_id_row2]), abund_table2])

    # Merge the abundance tables
    merged_df = pd.concat([abund_table1, abund_table2], axis=1)

    # Fill NaN values with 0
    merged_df = merged_df.fillna(0)

    # Transpose the DataFrame
    merged_df_transp = merged_df.T

    # Rename study_id column
    merged_df_transp = merged_df_transp.rename(columns={0: 'study_id'})

    # Extract the study_id column and remove it from the DataFrame
    study_id = merged_df_transp['study_id']
    merged_df_transp = merged_df_transp.drop(columns=['study_id'])
    merge_df = merged_df_transp.T

    # Show the merged abundance table
    st.subheader(f"Abundance table for {selected_studies[0]} and {selected_studies[1]}")
    
    # Display abundance table for the selected study
    builder = GridOptionsBuilder.from_dataframe(merged_df_transp)
    builder.configure_default_column(editable=True, groupable=True)
    builder.configure_side_bar(filters_panel = True, columns_panel = True)
    builder.configure_selection(selection_mode="multiple")
    builder.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
    go = builder.build()

    AgGrid(merged_df_transp, gridOptions=go)

    # Button to download the data
    abund_table_csv = convert_df(merged_df_transp)
    st.download_button(
        label="Download abundance data as CSV",
        data=abund_table_csv,
        file_name=f'abund_table_{selected_studies[0]}_and_{selected_studies[1]}.csv',
        mime='text/csv',
    )

    # Create PcoA plot
    st.subheader('PCoA plot with Bray-Curtis distance')
    st.write('The PCoA plot shows the distribution of the analyses from the two selected studies based on the Bray-Curtis distance between them.')
    
    # Extract analyses names as numpy arrays
    analyses_names = list(merged_df_transp.index.values)

    # Convert abundance table to numpy array
    abund_table_merged_mat = merged_df_transp.apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy()
    # print(type(abund_table_merged_mat))

    # Obtain bray-curtis distance matrix
    bc_mat = beta_diversity("braycurtis", abund_table_merged_mat, analyses_names)

    # Replace NaN values with 0
    bc_mat = np.nan_to_num(bc_mat.data, nan=0.0)
    
    # Run PCoA
    bc_pcoa = pcoa(bc_mat)

    # Extract the data to plot the PcoA
    bc_pcoa_data = pd.DataFrame(data = bc_pcoa.samples, columns = ['PC1', 'PC2'])
    
    # Reset index
    bc_pcoa_data = bc_pcoa_data.reset_index(drop=True)

    # Add analyses names as index
    bc_pcoa_data.index = merged_df.columns

    # Add study_id column to the PCoA df
    bc_pcoa_data['study_id'] = study_id

    # Add biome column to the PCoA df 
    bc_pcoa_data = bc_pcoa_data.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

    # Add biome in the study id column 
    bc_pcoa_data['study_id'] = bc_pcoa_data['study_id'].str.cat(bc_pcoa_data['biomes'], sep=' - ')

    # Get explained variance ratio
    explained_var_ratio = bc_pcoa.proportion_explained

    # Create a plotly figure
    pcoa_plot = px.scatter(bc_pcoa_data, x='PC1', y='PC2', opacity=0.8, color='study_id', 
                      hover_data=['study_id'], color_discrete_sequence=px.colors.qualitative.Plotly)

    # Add title and axis labels
    pcoa_plot.update_traces(
        marker=dict(size=7)
        ).update_layout(
        xaxis=dict(
            title=f'PCo1 ({explained_var_ratio[0]:.2%})',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        yaxis=dict(
            title=f'PCo2 ({explained_var_ratio[1]:.2%})',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        legend_title=dict(text='Study ID - Biome', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

    # Show pcoa plot
    st.plotly_chart(pcoa_plot, use_container_width=True)

    # Clustering plot
    # Run cluster analysis with k-means
    st.subheader('Clustering analysis')
    st.write('The clustering analysis was performed using the kmeans algorithm. The number of clusters was set to 2, based on the number of studies selected.')

    # Specify the number of clusters (k)
    k = 2

    # Fit K-Means model
    kmeans = KMeans(n_clusters=k)
    clusters = kmeans.fit_predict(abund_table_merged_mat)

    # Add the cluster labels to the PCA df
    bc_pcoa_data['Cluster'] = clusters

    # Set column type to string
    bc_pcoa_data['Cluster'] = bc_pcoa_data['Cluster'].astype(str)

    # Create a plotly figure
    clust_plot = px.scatter(bc_pcoa_data, x='PC1', y='PC2', opacity=0.8, color='Cluster',
                            color_discrete_sequence=px.colors.qualitative.Plotly)
    
    # Add title and axis labels
    clust_plot.update_layout(
        xaxis=dict(
            title='PC1',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        yaxis=dict(
            title='PC2',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        legend_title=dict(text='Cluster', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

    # Show the plot
    st.plotly_chart(clust_plot, use_container_width=True)
    


