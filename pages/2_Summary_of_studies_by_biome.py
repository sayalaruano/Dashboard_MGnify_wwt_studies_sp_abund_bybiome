# Web app
import streamlit as st
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
import plotly.express as px
import plotly.graph_objects as go
import colorcet as cc
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import squareform

# OS and file management
from PIL import Image

# General options 
im = Image.open("img/favicon.ico")
st.set_page_config(
    page_title="Summary and EDA of waste water treatment studies from Mgnify at species level and combined by biome",
    page_icon=im,
    layout="wide",
)

# Attach customized ccs style
with open('style.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

# Add a title and info about the app
st.title('Summary and EDA of waste water treatment studies from Mgnify at species level and combined by biome')
st.header('Plots to summarize studies by biome')

# Get the unique biomes removing nan values
biomes = st.session_state.studies_data["biomes"].unique()
biomes = [b for b in biomes if str(b) != 'nan']

# Sort the biomes
biomes.sort()

# Selectbox to choose a biome
st.subheader("Select a biome to create the plots")
biome = st.selectbox("Biome:", biomes) 

# Replace the ":" and " " characters by "_" in the biome name
biome = biome.replace(":", "_").replace(" ", "_")

# Load the merged sample data for the selected biome
sample_info = pd.read_csv(f"Samples_metadata/Merged_tables/{biome}_merged_samples_metadata.csv")

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
    file_name=f'sample_info_{biome}.csv',
    mime='text/csv',
)

# Create PCoA plots for all studies at species level
# Load the merged abundance and sample data
abund_df_species = pd.read_csv(f"Abundance_tables/Merged_tables/{biome}/{biome}_merged_abund_tables_species.csv", index_col=0)
tax_df_species = pd.read_csv(f"Abundance_tables/Merged_tables/{biome}/{biome}_merged_taxa_tables_species.csv", index_col=0)
tax_df_species['Genus_Species'] = tax_df_species['Genus'] + '_' + tax_df_species['Species']
study_ids = pd.read_csv(f"Abundance_tables/Merged_tables/{biome}/{biome}_studies_per_sample.csv", index_col=0)

# Transpose the DataFrame
abund_df_species_transp = abund_df_species.T
tax_df_species_transp = tax_df_species.T

# Reshape the abundance table
abund_df_reshaped = abund_df_species.reset_index().melt(id_vars='OTU', var_name='assembly_run_ids', value_name='count')

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

# Extract analyses names as numpy arrays
analyses_names = list(abund_df_species_transp.index.values)

# Convert abundance table to numpy array
abund_table_mat_species = abund_df_species_transp.to_numpy()

# Obtain bray-curtis distance matrix
bc_mat_species = beta_diversity("braycurtis", abund_table_mat_species, analyses_names)

# Replace NaN values with 0
bc_mat_species = np.nan_to_num(bc_mat_species.data, nan=0.0)

# Run PCoA
bc_pcoa_species = pcoa(bc_mat_species)

# Extract the data to plot the PCoA
bc_pcoa_species_data = pd.DataFrame(data = bc_pcoa_species.samples, columns = ['PC1', 'PC2'])

# Reset index
bc_pcoa_species_data = bc_pcoa_species_data.reset_index(drop=True)

# Add specific biome column to the PCoA df
bc_pcoa_species_data['specific_biome'] = samples_df['biome']

# Add analyses names as index
bc_pcoa_species_data.index = abund_df_species.columns

# Add study_id column to the PCoA df
bc_pcoa_species_data['study_id'] = study_ids["study_id"]

# Add general biome column to the PCoA df
bc_pcoa_species_data = bc_pcoa_species_data.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

# If specific biome is "NA - NA", replace it with the general biome
bc_pcoa_species_data['specific_biome'] = np.where(bc_pcoa_species_data['specific_biome'] == 'NA - NA', bc_pcoa_species_data['biomes'], bc_pcoa_species_data['specific_biome'])

# Add country column to the PCoA df
bc_pcoa_species_data = bc_pcoa_species_data.merge(st.session_state.studies_data[['study_id', 'sampling_country']], on='study_id')

# Add type of data column to the PCoA df
bc_pcoa_species_data = bc_pcoa_species_data.merge(st.session_state.studies_data[['study_id', 'experiment_type']], on='study_id')

# Add pipeline version column to the PCoA df
bc_pcoa_species_data = bc_pcoa_species_data.merge(st.session_state.studies_data[['study_id', 'pipeline_version']], on='study_id')

# Change the pipeline version to a string column
bc_pcoa_species_data['pipeline_version'] = bc_pcoa_species_data['pipeline_version'].astype(str)

# Add sequecing platform column to the PCoA df
bc_pcoa_species_data = bc_pcoa_species_data.merge(st.session_state.studies_data[['study_id', 'instrument_platform']], on='study_id')

# Add biome in the study id column
bc_pcoa_species_data['study_id'] = bc_pcoa_species_data['study_id'].str.cat(bc_pcoa_species_data['biomes'], sep=' - ')

# Get explained variance ratio
explained_var_ratio = bc_pcoa_species.proportion_explained

# Create a boxplot to show the distribution of the distances within and between the studies
st.subheader(f'Distribution of the Bray-Curtis distances at species level')
# Initialize an empty DataFrame to store distances and corresponding study IDs
distances_df = pd.DataFrame()

# Loop over each unique study ID
for study in np.unique(study_ids):
    # Get the indices of samples belonging to this study
    indices = np.where(study_ids == study)[0]

    # Extract the submatrix of distances for these samples
    submatrix = bc_mat_species[np.ix_(indices, indices)]

    # Convert the submatrix to a 1D array of distances
    dist_within = squareform(submatrix)

    # Filter out zero distances (self-comparisons)
    dist_within = dist_within[dist_within != 0]

    # Add these distances to the DataFrame
    temp_df = pd.DataFrame({'Distance': dist_within, 'Study': study})
    distances_df = pd.concat([distances_df, temp_df], ignore_index=True)

# Create the violin plot for each study
violin_plot_species = px.violin(distances_df, y='Distance', x='Study', box=True,
                color='Study', color_discrete_sequence=px.colors.qualitative.Dark24)

# Rotate x-axis labels to display them vertically
violin_plot_species.update_layout(
    xaxis=dict(
        tickfont=dict(size=18),
        titlefont=dict(size=20)
    ),
    yaxis=dict(
        title='Bray-Curtis distance',
        tickfont=dict(size=18),
        titlefont=dict(size=20),
        showgrid=False
    ),
    legend_title=dict(text='Study', font=dict(size=24)),
    legend=dict(font=dict(size=20))
)

st.plotly_chart(violin_plot_species, use_container_width=True)
st.warning('Some studies have just a few samples, so they are not shown in this plot', icon="⚠️")

# Create a barplot to show the top 5 genera per study
st.subheader(f'Top 5 genera per study')

# Merge the abund_df_species and tax_df_species by index
abund_df_species_merged = abund_df_species.merge(tax_df_species, left_index=True, right_index=True)

# Set "Genus_Species" column as index for the merged DataFrame
abund_df_species_merged.index = abund_df_species_merged['Genus_Species']

# Delete extra taxonomic columns
abund_df_species_merged.drop(columns=['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genus_Species'], inplace=True)

# Transpose the DataFrame
abund_df_species_merged_transp = abund_df_species_merged.T

# Add study_id back to merged_df_species
abund_df_species_merged_transp['study_id'] = study_ids

# Initialize an empty DataFrame for the top 5 genera data
top_genera_df = pd.DataFrame()

# Loop over each study to find the top 5 genera
for study in abund_df_species_merged_transp['study_id'].unique():
    study_data = abund_df_species_merged_transp[abund_df_species_merged_transp['study_id'] == study]
    top_genera = study_data.drop(columns='study_id').sum().nlargest(5)
    total_abundance = top_genera.sum()
    
    # Normalize the abundance for relative comparison
    top_genera_relative = (top_genera / total_abundance) * 100

    temp_df = pd.DataFrame({
        'Study': study,
        'Species': top_genera_relative.index,
        'Relative Abundance': top_genera_relative.values
    })
    top_genera_df = pd.concat([top_genera_df, temp_df])

# Generate a palette with many unique colors and convert the colorcet palette to HEX format
palette_hex = ['#' + ''.join([f'{int(c*255):02x}' for c in rgb]) for rgb in cc.glasbey_bw_minc_20]

# Select colors based on the unique values of the species column
unique_values_top_genera = top_genera_df["Species"].nunique()
selected_palette_top_genera = palette_hex[:unique_values_top_genera]

# Create a stacked bar chart for the top 5 genera
top_genera_plot = px.bar(top_genera_df, x='Study', y='Relative Abundance', color = 'Species',
             category_orders={"Species": top_genera_df['Species'].unique()}, opacity=0.8,
             color_discrete_sequence=selected_palette_top_genera)

# Update layout to adjust the margin, if necessary, to ensure text is not cut off
top_genera_plot.update_layout(
    margin=dict(l=40, r=40, t=40, b=40),
    xaxis=dict(
        title='Study',
        tickfont=dict(size=18),
        titlefont=dict(size=20),
        showgrid=False
    ),
    yaxis=dict(
        title='Relative abundance (%)',
        tickfont=dict(size=18),
        titlefont=dict(size=20),
        showgrid=False
    ),
    legend_title=dict(font=dict(size=24)),
    legend=dict(font=dict(size=20))
)

# Show the plot
st.plotly_chart(top_genera_plot, use_container_width=True)

# Display the top genera data
builder = GridOptionsBuilder.from_dataframe(top_genera_df)
builder.configure_default_column(editable=True, groupable=True)
builder.configure_side_bar(filters_panel = True, columns_panel = True)
builder.configure_selection(selection_mode="multiple")
builder.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
grid_opt = builder.build()

AgGrid(top_genera_df, gridOptions=grid_opt)

# Button to download the data
top_genera_csv = convert_df(top_genera_df)
st.download_button(
    label="Download top5 genera for all studies data as CSV",
    data=top_genera_csv,
    file_name=f'top5_genera_data_allstudies.csv',
    mime='text/csv',
)

# Plot PCoA colored by biome
st.subheader(f"PCoA plot (Bray Curtis distance) of the analyses from all studies at species level")

# Dropdown menu for selecting the color variable
color_option = st.selectbox("Select a variable to color by:", 
                            ('Biomes', 'Study ID and Biome', 'Sampling country', 'Experiment type',
                             'MGnify pipeline', 'Sequencing platform'))

# Create a function to update the figure
def update_figure(selected_variable):
    if selected_variable == 'Biomes':
        return 'specific_biome'
    elif selected_variable == 'Study ID and Biome':
        return 'study_id'
    elif selected_variable == 'Sampling country':
        return 'sampling_country'
    elif selected_variable == 'Experiment type':
        return 'experiment_type'
    elif selected_variable == 'MGnify pipeline':
        return 'pipeline_version'
    elif selected_variable == 'Sequencing platform':
        return 'instrument_platform'

# Select colors based on the unique values of the selected variable
color_var = update_figure(color_option)
unique_values_pcoa = bc_pcoa_species_data[color_var].nunique()
selected_palette_pcoa = palette_hex[:unique_values_pcoa]

# Make the plot
pcoa_species = px.scatter(bc_pcoa_species_data, x='PC1', y='PC2', opacity=0.8, color=color_var,
                            hover_data=['study_id'], color_discrete_sequence=selected_palette_pcoa)

# Add title and axis labels
pcoa_species.update_traces(
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
    legend_title=dict(text=color_option, font=dict(size=24)),
    legend=dict(font=dict(size=20))
)

# Show the plot
st.plotly_chart(pcoa_species, use_container_width=True)

# Display merged abundance data
biome_title = biome.replace("Wastewater_", "").replace("_", " ")
st.subheader(f'Abundance table for the {biome_title} biome')

# Reset the index to make it a column
abund_table_with_index = abund_df_species_merged.reset_index()

# Ensure the index column is the first one
column_order = ["Genus_Species"] + [col for col in abund_table_with_index.columns if col != "Genus_Species"]
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
    label=f"Download abundance data as CSV",
    data=abund_table_csv,
    file_name=f'abund_table_{biome}.csv',
    mime='text/csv',
)

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')