# Web app
import glob
import os
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
from scipy.spatial.distance import squareform

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

# Function to load data for a specific study
@st.cache_data
def load_study_data(all_data, selected_study):
    # Filter data for the selected study
    study_info = all_data[all_data['study_id'] == selected_study]
    
    # Load sample information for the study
    sample_info = pd.read_csv(f"Samples_metadata/{selected_study}/{selected_study}_samples_metadata.csv")
    
    return study_info, sample_info

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

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
def preprocess_abund_table(abund_table,phylum):
    # Delete rows with NANs in all columns
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
        # Check available taxonomic levels and drop the corresponding columns
        taxonomic_levels = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Family_Genus']
        for level in taxonomic_levels:
            if level in abund_table.columns:
                abund_table = abund_table.drop(columns=level)
        
        # Set the genus column as index
        abund_table = abund_table.set_index('Genus')

    return abund_table

# Get unique study names
studies_list = st.session_state.studies_data['study_id'].unique()
studies_list = studies_list[~pd.isnull(studies_list)]

# Add info on the sidebar
# Select study
st.sidebar.header('How it works?')

st.sidebar.write('Select a study and taxonomic rank in the sidebar to display its information in the main panel.')

# Create a selectbox to choose the study
selected_study = st.sidebar.selectbox("Mgnify study:", studies_list)

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
st.title('Summary and EDA of waste water treatment studies from Mgnify combined by biomes')

# Load data for the selected study
study_info, sample_info = load_study_data(st.session_state.studies_data, selected_study)

# Display study information
st.subheader(f'Study {selected_study}')

# Create a df to summarize the study information
summ_df = pd.DataFrame(study_info).T

# Remove study_id row
summ_df = summ_df.drop(index='study_id')

# Replace indices and column names
summ_df.index = ['Description', 'Bioproject ID', 'Centre name', 'Number of samples', 
                 'Biome','Experiment type', 'Pipeline version', 'Instrument platform',
                 'Sampling country']

summ_df.columns = ['Summary metrics']

st.data_editor(
    summ_df,
    width = 1000,
    column_config={
        "widgets": st.column_config.TextColumn(
            f"Summary {selected_study}",
            help= f"Summary of the study {selected_study}"
        )
    },
    hide_index=False,
)

# Display sample information for the selected study
st.subheader("Sample Information")

builder = GridOptionsBuilder.from_dataframe(sample_info)
builder.configure_default_column(editable=True, groupable=True)
builder.configure_side_bar(filters_panel = True, columns_panel = True)
builder.configure_selection(selection_mode="multiple")
builder.configure_pagination(paginationAutoPageSize=False, paginationPageSize=15)
go = builder.build()

AgGrid(sample_info, gridOptions=go)

# Button to download the data
sample_info_csv = convert_df(sample_info)
st.download_button(
    label="Download sample data as CSV",
    data=sample_info_csv,
    file_name=f'sample_info_{selected_study}.csv',
    mime='text/csv',
)

# Load and preprocess abundance table for the selected study and taxonomic rank. Store it in an independent variable
if tax_rank == 'Phylum':
    abund_table = load_abund_table(selected_study, phylum=True)
    abund_table = preprocess_abund_table(abund_table, phylum=True)
else:  # Genus
    abund_table = load_abund_table(selected_study, phylum=False)
    abund_table = preprocess_abund_table(abund_table, phylum=False)

# Display abundance table for the selected study
st.subheader(f"Abundance table at {tax_rank} level")

# Capture the name of the index, or default to 'index' if it's unnamed
index_name = abund_table.index.name if abund_table.index.name else 'index'

# Reset the index to make it a column
abund_table_with_index = abund_table.reset_index()

# Ensure the index column is the first one
column_order = [index_name] + [col for col in abund_table_with_index.columns if col != index_name]
abund_table_with_index = abund_table_with_index[column_order]

builder = GridOptionsBuilder.from_dataframe(abund_table_with_index)
builder.configure_default_column(editable=True, groupable=True)
builder.configure_side_bar(filters_panel = True, columns_panel = True)
builder.configure_selection(selection_mode="multiple")
builder.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
go = builder.build()

AgGrid(abund_table_with_index, gridOptions=go)

# Button to download the data
abund_table_csv = convert_df(abund_table_with_index)
st.download_button(
    label="Download abundance data as CSV",
    data=abund_table_csv,
    file_name=f'abund_table_{tax_rank}_{selected_study}.csv',
    mime='text/csv',
)

# Add condition to check if the study has more than one sample
if len(sample_info) > 2:

    # Reshape abundance df so that the assembly_run_ids become a column, considering the phylum or genus column as the index
    if tax_rank == 'Phylum':
        abund_df_reshaped = abund_table.reset_index().melt(id_vars='phylum', var_name='assembly_run_ids', value_name='count')
    else:  # Genus
        abund_df_reshaped = abund_table.reset_index().melt(id_vars='Genus', var_name='assembly_run_ids', value_name='count')
    
    # Split the multiple IDs in the assembly_run_ids column of df1
    sample_info['assembly_run_ids'] = sample_info['assembly_run_ids'].str.split(';')

    # Explode the dataframe based on the assembly_run_ids column
    sample_info_exploded = sample_info.explode('assembly_run_ids')

    # Merge dataframes on the 'assembly_run_ids' column
    samples_df = pd.merge(sample_info_exploded, abund_df_reshaped, left_on='assembly_run_ids', right_on='assembly_run_ids')

    # Filter the merged df to keep unique rows based on the assembly_run_ids column
    samples_df = samples_df.drop_duplicates(subset=['assembly_run_ids'])

    # Filter the merged df to keep only the relevant columns
    samples_df = samples_df[['assembly_run_ids', 'sample_id', 'biome_feature', 'biome_material']].reset_index(drop=True)

    # Combine biome_feature and biome_material columns if there is info on the columns, otherwise add 'NA'
    samples_df['biome_feature'] = samples_df['biome_feature'].fillna('NA')
    samples_df['biome_material'] = samples_df['biome_material'].fillna('NA')
    samples_df['biome'] = samples_df['biome_feature'] + ' - ' + samples_df['biome_material']

    # Sort df by assembly_run_ids
    samples_df = samples_df.sort_values(by=['assembly_run_ids']).reset_index(drop=True)

    # Create PcoA plot at taxonomic rank level
    st.subheader(f'PCoA plot at {tax_rank} level with Bray-Curtis distance')
    st.write('The PCoA plot shows the distribution of the analyses from the selected study based on the Bray-Curtis distance between them.')

    # Transpose abundance table
    abund_table = abund_table.T

    # Extract analyses names as numpy arrays
    analyses_names = list(abund_table.index.values)

    # Convert abundance table to numpy array
    abund_table_mat = abund_table.to_numpy()

    # Obtain bray-curtis distance matrix
    bc_mat = beta_diversity("braycurtis", abund_table_mat, analyses_names)

    # Replace NaN values with 0
    bc_mat = np.nan_to_num(bc_mat.data, nan=0.0)

    # Run PCoA
    bc_pcoa = pcoa(bc_mat)

    # Extract the data to plot the PcoA
    bc_pcoa_data = pd.DataFrame(data = bc_pcoa.samples[['PC1', 'PC2']])

    # Reset index
    bc_pcoa_data = bc_pcoa_data.reset_index(drop=True)

    # Add the biome to the PcoA df
    bc_pcoa_data['biome'] = samples_df['biome']

    # Get explained variance ratio
    explained_var_ratio = bc_pcoa.proportion_explained

    # Create a plotly figure
    pcoa_plot = px.scatter(bc_pcoa_data, x='PC1', y='PC2', 
                        opacity=0.8, color='biome', 
                        hover_name=abund_table.index,
                        color_discrete_sequence=px.colors.qualitative.Plotly)

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
        legend_title=dict(text="Biome", font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

    # Show pcoa plot at the chosen taxonomic rank level
    st.plotly_chart(pcoa_plot, use_container_width=True)

    # Create a boxplot to show the distribution of the distances within the study
    st.subheader(f'Distribution of the Bray-Curtis distances at {tax_rank} level')

    # Assuming bc_mat is the Bray-Curtis distance matrix you have
    # Convert the distance matrix to a 1D array of distances
    dist_within = squareform(bc_mat)

    # Filter out zero distances (self-comparisons)
    dist_within = dist_within[dist_within != 0]

    # Prepare data for violin plot
    data = {'Distance': dist_within}
    df_violin_plot = pd.DataFrame(data)

    # Create the violin plot
    violin_plot = px.violin(df_violin_plot, y='Distance', box=True, color_discrete_sequence=px.colors.qualitative.Plotly)

    # Add title and axis labels
    violin_plot.update_layout(
        yaxis=dict(
            title='Bray-Curtis Distance',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        legend_title=dict(text='Distance', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

    # Show violin plot at the chosen taxonomic rank level
    st.plotly_chart(violin_plot, use_container_width=True)

    # Create a datframe with the top 10 most abundant taxa
    # Reset the index of the DataFrame, this will add the index as a new column
    abund_table_reset = abund_table.T.reset_index()

    # Remove the first row (which is the index)
    abund_table_reset = abund_table_reset.drop(abund_table_reset.index[0])

    if tax_rank == 'Phylum':
        # Now use the taxa column as id_vars in melt
        abund_table_top10 = abund_table_reset.melt(id_vars='phylum', var_name='assembly_run_ids', value_name='count')
        # Rename the column
        abund_table_top10 = abund_table_top10.rename(columns={'phylum': 'taxa'})
    else:  # Genus
        abund_table_top10 = abund_table_reset.melt(id_vars='Genus', var_name='assembly_run_ids', value_name='count')
        abund_table_top10 = abund_table_top10.rename(columns={'Genus': 'taxa'})

    # Group by taxonomic rank and sum the counts
    abund_table_top10 = abund_table_top10.groupby('taxa').sum().reset_index()

    # Sort values by count
    abund_table_top10 = abund_table_top10.sort_values(by='count', ascending=False).reset_index(drop=True)

    # Keep only the top 10 most abundant taxa
    abund_table_top10 = abund_table_top10.head(10)

    # Create a bar plot to show the top 10 most abundant taxa
    st.subheader(f'Top 10 most abundant taxa at {tax_rank} level')

    # Create a plotly figure
    top10_plot = px.bar(abund_table_top10, x="taxa", y='count', color_discrete_sequence=px.colors.qualitative.Plotly)

    # Add title and axis labels
    top10_plot.update_layout(
        xaxis=dict(
            title=f'{tax_rank}',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        yaxis=dict(
            title='Count',
            tickfont=dict(size=18),
            titlefont=dict(size=20),
            showgrid=False
        ),
        legend_title=dict(text='Count', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

    # Show top 10 most abundant taxa plot at the chosen taxonomic rank level
    st.plotly_chart(top10_plot, use_container_width=True)

else:
    st.warning('The selected study has less than two samples, so the PCoA plot and the distribution of the Bray-Curtis distances cannot be created.', icon='⚠️')

