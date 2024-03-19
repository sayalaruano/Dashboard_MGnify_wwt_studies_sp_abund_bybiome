# Web app
import streamlit as st
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import colorcet as cc
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import squareform

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

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

# Add a title and info about the app
st.title('Summary and EDA of waste water treatment studies from Mgnify combined by biomes')
st.header('Plots to summarize all biomes')

# Create plot the number of studies per study
st.subheader("Number of samples per study")

hist_samples = px.histogram(st.session_state.studies_data, x="n_samples", 
                            nbins=40, opacity=0.8, color_discrete_sequence=px.colors.qualitative.Plotly)

hist_samples.update_layout(
    xaxis=dict(title='Number of samples', tickfont=dict(size=18), titlefont=dict(size=20), dtick=10, tick0=0),
    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20))
    ).update_xaxes(
        showgrid=False
    ).update_yaxes(
        showgrid=False)

st.plotly_chart(hist_samples, use_container_width=True)

# Create a pie plot for the sampling countries of the studies
st.subheader("Sampling countries for all studies")

pie_plot_countries = px.pie(values=st.session_state.studies_data["sampling_country"].value_counts(),
                names=st.session_state.studies_data["sampling_country"].value_counts().index, 
                opacity=0.8, color_discrete_sequence=px.colors.qualitative.Plotly)

pie_plot_countries.update_traces(
    textposition='inside',
    textinfo='value',
    insidetextfont=dict(size=18)
    ).update_layout(
        legend_title=dict(text='Country', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

st.plotly_chart(pie_plot_countries, use_container_width=True)

# Create bar plot for the data types of the studies
st.subheader("Data types for all studies")

bar_plot_dtypes = px.histogram(st.session_state.studies_data, x="experiment_type", opacity=0.8, 
                               color='biomes', color_discrete_sequence=px.colors.qualitative.Plotly)

bar_plot_dtypes.update_layout(
    xaxis=dict(title='Experiment type', tickfont=dict(size=18), titlefont=dict(size=20)),
    yaxis=dict(title='Number of studies', tickfont=dict(size=18), titlefont=dict(size=20)),
    legend_title=dict(text='Biome', font=dict(size=24)),
    legend=dict(font=dict(size=20)),
    ).update_xaxes(
        showgrid=False
    ).update_yaxes(
        showgrid=False)

st.plotly_chart(bar_plot_dtypes, use_container_width=True)

# Load the merged abundance and taxonomic data
abund_df_genus_wwt = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater/Wastewater_merged_abund_tables_genus.csv", index_col=0)
tax_df_genus_wwt = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater/Wastewater_merged_taxa_tables_genus.csv", index_col=0)
abund_df_genus_wwt_ws = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Water_and_sludge/Wastewater_Water_and_sludge_merged_abund_tables_genus.csv", index_col=0)
tax_df_genus_wwt_ws = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Water_and_sludge/Wastewater_Water_and_sludge_merged_taxa_tables_genus.csv", index_col=0)
abund_df_genus_wwt_ind = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Industrial_wastewater/Wastewater_Industrial_wastewater_merged_abund_tables_genus.csv", index_col=0)
tax_df_genus_wwt_ind = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Industrial_wastewater/Wastewater_Industrial_wastewater_merged_taxa_tables_genus.csv", index_col=0)
abund_df_genus_wwt_as = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Activated_Sludge/Wastewater_Activated_Sludge_merged_abund_tables_genus.csv", index_col=0)
tax_df_genus_wwt_as = pd.read_csv(f"Abundance_tables/Merged_tables/Wastewater_Activated_Sludge/Wastewater_Activated_Sludge_merged_taxa_tables_genus.csv", index_col=0)

# Create pie plot for the biomes of the studies
st.subheader("Number of studies and samples per biome")

# Calculate the number of studies per biome
studies_per_biome = st.session_state.studies_data["biomes"].value_counts()
biomes = studies_per_biome.index

# Calculate the number of samples per biome by counting the number of columns in the abundance table, excluding the "Genus" column
samples_per_biome = [abund_df_genus_wwt_as.shape[1], abund_df_genus_wwt.shape[1], abund_df_genus_wwt_ws.shape[1], abund_df_genus_wwt_ind.shape[1]]

# Create subplots
fig = make_subplots(1, 2, specs=[[{'type': 'domain'}, {'type': 'domain'}]],
                    subplot_titles=['Studies', 'Samples'])

# Add the studies pie chart
fig.add_trace(go.Pie(labels=biomes, values=studies_per_biome.values, name='Studies',
                     marker=dict(colors=px.colors.qualitative.Plotly)), 1, 1)

# Add the samples pie chart
fig.add_trace(go.Pie(labels=biomes, values=samples_per_biome, name='Samples',
                     marker=dict(colors=px.colors.qualitative.Plotly)), 1, 2)

# Update layout
fig.update_layout(
    legend_title=dict(text='Biomes', font=dict(size=24)),
    legend=dict(font=dict(size=20))
)
fig.update_traces(textposition='inside', textinfo='value+percent', insidetextfont=dict(size=18))

# Display the plot in Streamlit
st.plotly_chart(fig, use_container_width=True)

# Create a stacked bar plot for the top 5 genera by biome
st.subheader("Top 5 genera by biome")

# Initialize an empty DataFrame for the top 5 genera data
top_genera_df_all_biomes = pd.DataFrame()

# List of biome dataframes and their names
abund_dfs = [abund_df_genus_wwt, abund_df_genus_wwt_ws, abund_df_genus_wwt_ind, abund_df_genus_wwt_as]
tax_dfs = [tax_df_genus_wwt, tax_df_genus_wwt_ws, tax_df_genus_wwt_ind, tax_df_genus_wwt_as]
biome_names = ['Wastewater', 'Water and sludge', 'Industrial wastewater', 'Activated sludge']

# Process each biome
for abund_df, tax_df, biome_name in zip(abund_dfs, tax_dfs, biome_names):
    # Merge the abundance and taxonomy dataframes by index
    abund_tax_merged = abund_df.merge(tax_df, left_index=True, right_index=True)
    abund_tax_merged.index = abund_tax_merged['Genus']
    abund_tax_merged.drop(columns=['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'], inplace=True)
    abund_tax_merged_transp = abund_tax_merged.T

    # Calculate top 5 genera
    top_genera = abund_tax_merged_transp.sum().nlargest(5)
    total_abundance = top_genera.sum()
    top_genera_relative = (top_genera / total_abundance) * 100

    temp_df = pd.DataFrame({
        'Biome': biome_name,
        'Genus': top_genera_relative.index,
        'Relative Abundance': top_genera_relative.values
    })
    top_genera_df_all_biomes = pd.concat([top_genera_df_all_biomes, temp_df])

# Create a stacked bar chart for the top 5 genera in each biome
top_genera_plot_biome = px.bar(top_genera_df_all_biomes, x='Biome', y='Relative Abundance', color='Genus',
             category_orders={"Genus": top_genera_df_all_biomes['Genus'].unique()}, opacity=0.8,
             color_discrete_sequence=px.colors.qualitative.Dark24)

# Update layout (similar to your existing code, but adjust titles and labels accordingly)
top_genera_plot_biome.update_layout(
    margin=dict(l=40, r=40, t=40, b=40),
    xaxis=dict(
        title='Biome',
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
st.plotly_chart(top_genera_plot_biome, use_container_width=True)

# Load the merged sample data for the selected biome
sample_info = pd.read_csv(f"Samples_metadata/Merged_tables/all_biomes_merged_samples_metadata.csv")

# Create PCoA plots for all studies at genus level
# Load the merged abundance and sample data
abund_df_genus = pd.read_csv(f"Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_abund_tables_genus.csv", index_col=0)
tax_df_genus = pd.read_csv(f"Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_taxa_tables_genus.csv")
study_ids = pd.read_csv(f"Abundance_tables/Merged_tables/All_biomes/all_biomes_studies_per_sample.csv", index_col=0)

# Transpose the DataFrame
abund_df_genus_transp = abund_df_genus.T
tax_df_genus_transp = tax_df_genus.T

# Reshape the abundance table
abund_df_reshaped = abund_df_genus.reset_index().melt(id_vars='OTU', var_name='assembly_run_ids', value_name='count')

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

# Sort df by assembly_run_ids
samples_df = samples_df.sort_values(by=['assembly_run_ids']).reset_index(drop=True)

# Extract analyses names as numpy arrays
analyses_names = list(abund_df_genus_transp.index.values)

# Convert abundance table to numpy array
abund_table_mat_genus = abund_df_genus_transp.to_numpy()

# Obtain bray-curtis distance matrix
bc_mat_genus = beta_diversity("braycurtis", abund_table_mat_genus, analyses_names)

# Replace NaN values with 0
bc_mat_genus = np.nan_to_num(bc_mat_genus.data, nan=0.0)

# Run PCoA
bc_pcoa_genus = pcoa(bc_mat_genus)

# Extract the data to plot the PCoA
bc_pcoa_genus_data = pd.DataFrame(data = bc_pcoa_genus.samples, columns = ['PC1', 'PC2'])

# Reset index
bc_pcoa_genus_data = bc_pcoa_genus_data.reset_index(drop=True)

# Add analyses names as index
bc_pcoa_genus_data.index = abund_df_genus.columns

# Add study_id column to the PCoA df
bc_pcoa_genus_data['study_id'] = study_ids["study_id"]

# Add general biome column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

# Add country column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'sampling_country']], on='study_id')

# Add type of data column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'experiment_type']], on='study_id')

# Add pipeline version column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'pipeline_version']], on='study_id')

# Change the pipeline version to a string column
bc_pcoa_genus_data['pipeline_version'] = bc_pcoa_genus_data['pipeline_version'].astype(str)

# Add sequecing platform column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'instrument_platform']], on='study_id')

# Add biome in the study id column
bc_pcoa_genus_data['study_id'] = bc_pcoa_genus_data['study_id'].str.cat(bc_pcoa_genus_data['biomes'], sep=' - ')

# Get explained variance ratio
explained_var_ratio = bc_pcoa_genus.proportion_explained

# Generate a palette with many unique colors and convert the colorcet palette to HEX format
palette_hex = ['#' + ''.join([f'{int(c*255):02x}' for c in rgb]) for rgb in cc.glasbey_bw_minc_20]

# Plot PCoA colored by biome
st.subheader(f"PCoA plot (Bray Curtis distance) of the analyses from all studies at Genus level")

# Dropdown menu for selecting the color variable
color_option = st.selectbox("Select a variable to color by:", 
                            ('Biomes', 'Study ID and Biome', 'Sampling country', 'Experiment type',
                             'MGnify pipeline', 'Sequencing platform'))

# Create a function to update the figure
def update_figure(selected_variable):
    if selected_variable == 'Biomes':
        return 'biomes'
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
unique_values_pcoa = bc_pcoa_genus_data[color_var].nunique()
selected_palette_pcoa = palette_hex[:unique_values_pcoa]

# Make the plot
pcoa_genus = px.scatter(bc_pcoa_genus_data, x='PC1', y='PC2', opacity=0.8, color=color_var,
                            hover_data=['study_id'], color_discrete_sequence=selected_palette_pcoa)

# Add title and axis labels
pcoa_genus.update_traces(
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
st.plotly_chart(pcoa_genus, use_container_width=True)

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
    file_name=f'sample_info_allbiomes.csv',
    mime='text/csv',
)

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')