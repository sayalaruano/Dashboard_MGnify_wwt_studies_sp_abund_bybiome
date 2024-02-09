# Web app
import streamlit as st
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
import plotly.express as px
import plotly.graph_objects as go
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
        taxonomic_levels = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Species']
        for level in taxonomic_levels:
            if level in abund_table.columns:
                abund_table = abund_table.drop(columns=level)
        
        # Set the genus column as index
        abund_table = abund_table.set_index('Genus')

    return abund_table

# Function to download the data
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

# Add a title and info about the app
st.title('Summary and EDA of waste water treatment studies from Mgnify combined by biomes')
st.header('Plots to summarize all studies')

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
st.subheader("Sampling countries")

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
        legend_title=dict(text='Biome', font=dict(size=24)),
        legend=dict(font=dict(size=20))
    )

st.plotly_chart(pie_plot_biomes, use_container_width=True)

# Get the unique biomes removeing nan values
biomes = st.session_state.studies_data["biomes"].unique()
biomes = [b for b in biomes if str(b) != 'nan']

# Selectbox to choose a biome
st.subheader("Select a biome to create the plots")
biome = st.selectbox("Biome:", biomes) 

# Replace the ":" and " " characters by "_" in the biome name
biome = biome.replace(":", "_").replace(" ", "_")

# Create PCoA plots for all studies at genus level
# Load the merged abundance table
merged_df_genus = pd.read_csv(f"Abundance_tables/{biome}_merged_all_abund_tables_genus.csv", index_col=0)

# Extract the study_id column and remove it from the DataFrame
study_id = merged_df_genus['study_id']
merged_df_genus = merged_df_genus.drop(columns=['study_id'])

# Transpose the DataFrame
merged_df_genus_transp = merged_df_genus.T

# Extract analyses names as numpy arrays
analyses_names = list(merged_df_genus.index.values)

# Convert abundance table to numpy array
abund_table_mat_genus = merged_df_genus.to_numpy()

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
bc_pcoa_genus_data.index = merged_df_genus_transp.columns

# Add study_id column to the PCoA df
bc_pcoa_genus_data['study_id'] = study_id

# Add biome column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'biomes']], on='study_id')

# Add country column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'sampling_country']], on='study_id')

# Add type of data column to the PCoA df
bc_pcoa_genus_data = bc_pcoa_genus_data.merge(st.session_state.studies_data[['study_id', 'experiment_type']], on='study_id')

# Add biome in the study id column
bc_pcoa_genus_data['study_id'] = bc_pcoa_genus_data['study_id'].str.cat(bc_pcoa_genus_data['biomes'], sep=' - ')

# Get explained variance ratio
explained_var_ratio = bc_pcoa_genus.proportion_explained

# Create a boxplot to show the distribution of the distances within and between the studies
st.subheader(f'Distribution of the Bray-Curtis distances at Genus level')
# Initialize an empty DataFrame to store distances and corresponding study IDs
distances_df = pd.DataFrame()

# Loop over each unique study ID
for study in np.unique(study_id):
    # Get the indices of samples belonging to this study
    indices = np.where(study_id == study)[0]

    # Extract the submatrix of distances for these samples
    submatrix = bc_mat_genus[np.ix_(indices, indices)]

    # Convert the submatrix to a 1D array of distances
    dist_within = squareform(submatrix)

    # Filter out zero distances (self-comparisons)
    dist_within = dist_within[dist_within != 0]

    # Add these distances to the DataFrame
    temp_df = pd.DataFrame({'Distance': dist_within, 'Study': study})
    distances_df = pd.concat([distances_df, temp_df], ignore_index=True)

# Create the violin plot for each study
violin_plot_genus = px.violin(distances_df, y='Distance', x='Study', box=True,
                color='Study', color_discrete_sequence=px.colors.qualitative.Dark24)

# Rotate x-axis labels to display them vertically
violin_plot_genus.update_layout(
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

st.plotly_chart(violin_plot_genus, use_container_width=True)
st.warning('Some studies have just a few samples, so they are not shown in this plot', icon="⚠️")

# Create a boxplot to show the distribution of the distances within and between the studies
st.subheader(f'Top 5 genera per study')
# Add study_id back to merged_df_genus
merged_df_genus['study_id'] = study_id

# Initialize an empty DataFrame for the top 5 genera data
top_genera_df = pd.DataFrame()

# Loop over each study to find the top 5 genera
for study in merged_df_genus['study_id'].unique():
    study_data = merged_df_genus[merged_df_genus['study_id'] == study]
    top_genera = study_data.drop(columns='study_id').sum().nlargest(5)
    total_abundance = top_genera.sum()
    
    # Normalize the abundance for relative comparison
    top_genera_relative = (top_genera / total_abundance) * 100

    temp_df = pd.DataFrame({
        'Study': study,
        'Genus': top_genera_relative.index,
        'Relative Abundance': top_genera_relative.values
    })
    top_genera_df = pd.concat([top_genera_df, temp_df])

# Create a stacked bar chart
top_genera_plot = px.scatter(top_genera_df, x='Study', y='Relative Abundance', text='Genus', color = 'Study',
             category_orders={"Genus": top_genera_df['Genus'].unique()}, opacity=0.8,
             color_discrete_sequence=px.colors.qualitative.Dark24)

# Adjust text label size and position
for trace in top_genera_plot.data:
    trace.textposition = 'middle center'
    trace.marker.size = 7
    trace.textfont = dict(size=16)  # Increase the text size

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
    showlegend=False
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

# Generate a palette with 35 unique colors
# palette = cc.glasbey_bw_minc_20_maxl_70
# Convert the colorcet palette to HEX format
palette_hex = ['#' + ''.join([f'{int(c*255):02x}' for c in rgb]) for rgb in cc.glasbey_bw_minc_20]

# Plot PCoA colored by biome
st.subheader(f"PCoA plot (Bray Curtis distance) of the analyses from all studies at Genus level")

# Dropdown menu for selecting the color variable
color_option = st.selectbox("Select a variable to color by:", 
                            ('Biomes', 'Study ID and Biome', 'Sampling country', 'Data type'))

# Create a function to update the figure
def update_figure(selected_variable):
    if selected_variable == 'Biomes':
        return 'biomes'
    elif selected_variable == 'Study ID and Biome':
        return 'study_id'
    elif selected_variable == 'Sampling country':
        return 'sampling_country'
    elif selected_variable == 'Data type':
        return 'experiment_type'

# Select colors based on the unique values of the selected variable
color_var = update_figure(color_option)
unique_values = bc_pcoa_genus_data[color_var].nunique()
selected_palette = palette_hex[:unique_values]

# Make the plot
pcoa_genus = px.scatter(bc_pcoa_genus_data, x='PC1', y='PC2', opacity=0.8, color=color_var,
                            hover_data=['study_id'], color_discrete_sequence=selected_palette)

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

# Add info on the sidebar
st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')