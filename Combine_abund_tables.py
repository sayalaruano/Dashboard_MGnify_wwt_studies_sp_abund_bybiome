#%%
import pandas as pd
import glob
import os

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

# Function to load data for a specific study
def load_study_data(all_data, selected_study):
    # Filter data for the selected study
    study_info = all_data[all_data['study_id'] == selected_study]
    
    # Load sample information for the study
    sample_info = pd.read_csv(f"Samples_metadata/{selected_study}/{selected_study}_samples_metadata.csv")
    
    return sample_info

#%%
# Separate the studies by their biomes
# Load the summary table of all studies
studies_data = pd.read_csv('wwt_studies_more_median_taxa_bybiome.csv')

# Remove rows with nan values in all columns
studies_data = studies_data.dropna(how='all')

# Extract the study ids from the csv file
wwt_studies_allsamples_list = studies_data["study_id"].tolist()

# Remove nan values from the list
# wwt_studies_more10samples_list = [x for x in wwt_studies_more10samples_list if str(x) != 'nan']
wwt_studies_allsamples_list = [x for x in wwt_studies_allsamples_list if str(x) != 'nan']

# Count the number of studies per biome
biome_count = studies_data["biomes"].value_counts()

# Rename the entries in the biomes column
studies_data['biomes'] = studies_data['biomes'].replace({
    "root:Engineered:Wastewater:Nutrient removal:Biological phosphorus removal:Activated sludge": "root:Engineered:Wastewater:Activated Sludge",
    "root:Engineered:Wastewater:Nutrient removal:Dissolved organics (anaerobic)": "root:Engineered:Wastewater",
    "root:Engineered:Wastewater:Nutrient removal:Nitrogen removal": "root:Engineered:Wastewater",
    "root:Engineered:Wastewater:Industrial wastewater:Petrochemical": "root:Engineered:Wastewater:Industrial wastewater",
    "root:Engineered:Wastewater:Industrial wastewater:Agricultural wastewater": "root:Engineered:Wastewater:Industrial wastewater",
    "root:Engineered:Wastewater:Activated Sludge, root:Engineered:Wastewater:Industrial wastewater": "root:Engineered:Wastewater:Activated Sludge"
})

# Count the number of studies per biome
biome_count_new = studies_data["biomes"].value_counts()

# Get the lists of study ids for each biome
biome_list = studies_data["biomes"].unique()

# Create a dictionary to store the study ids for each biome
biome_dict = {}

# Iterate through the biomes and store the study ids in the dictionary
for biome in biome_list:
    biome_dict[biome] = studies_data[studies_data["biomes"] == biome]["study_id"].tolist()

# Remove the "root:Engineered:" prefix from the biome names in the dictionary
biome_dict = {biome.replace("root:Engineered:", ""): study_ids for biome, study_ids in biome_dict.items()}

# Replace ":" and speaces with "_" in the biome names
biome_dict = {biome.replace(":", "_").replace(" ", "_"): study_ids for biome, study_ids in biome_dict.items()}

#%%
# Iterate through the biomes in the dictionary
for biome, study_ids in biome_dict.items():
    # Initialize empty DataFrames for phylum and genus merged data
    merged_df_phylum = pd.DataFrame()
    merged_df_genus = pd.DataFrame()
    merged_df_sample_metadata = pd.DataFrame()

    # Process phylum and genus tables for each biome
    for selected_study in study_ids:
        # Process phylum table
        abund_table_phylum = load_abund_table(selected_study, phylum=True)
        if abund_table_phylum is not None:
            abund_table_phylum = preprocess_abund_table(abund_table_phylum, phylum=True)

            # Add a row with the study ID for all samples 
            study_id_row_phylum = pd.Series(selected_study, index=abund_table_phylum.columns)
            abund_table_phylum = pd.concat([pd.DataFrame([study_id_row_phylum]), abund_table_phylum])

            # Add the data from the current study to the merged DataFrame
            merged_df_phylum = pd.concat([merged_df_phylum, abund_table_phylum], axis=1)

        # Process genus table
        abund_table_genus = load_abund_table(selected_study, phylum=False)
        if abund_table_genus is not None:
            abund_table_genus = preprocess_abund_table(abund_table_genus, phylum=False)

            # Add a row with the study ID for all samples
            study_id_row_genus = pd.Series(selected_study, index=abund_table_genus.columns)
            abund_table_genus = pd.concat([pd.DataFrame([study_id_row_genus]), abund_table_genus])
            
            # Add the data from the current study to the merged DataFrame
            merged_df_genus = pd.concat([merged_df_genus, abund_table_genus], axis=1)
        
        # Load metadata for the study
        sample_metadata_df = load_study_data(studies_data, selected_study)
        if sample_metadata_df is not None:
            # Add a column with the study ID for all samples 
            sample_metadata_df['study_id'] = selected_study

            # Add the data from the current study to the merged DataFrame
            merged_df_sample_metadata = pd.concat([merged_df_sample_metadata, sample_metadata_df], axis=0)

    # Fill NaN values with 0, transpose the DataFrames and expor them
    if not merged_df_phylum.empty:
        merged_df_phylum = merged_df_phylum.fillna(0).T.rename(columns={0: 'study_id'})
        merged_df_phylum.index.name = 'assembly_run_ids'
        merged_df_phylum.to_csv(f'Abundance_tables/{biome}_merged_all_abund_tables_phylum.csv')

    if not merged_df_genus.empty:
        merged_df_genus = merged_df_genus.fillna(0).T.rename(columns={0: 'study_id'})
        merged_df_genus.index.name = 'assembly_run_ids'
        merged_df_genus.to_csv(f'Abundance_tables/{biome}_merged_all_abund_tables_genus.csv')
    
    if not merged_df_sample_metadata.empty:
        merged_df_sample_metadata.to_csv(f'Samples_metadata/{biome}_merged_samples_metadata.csv', index=False)

#%%
# Reset index of merged_df_genus
merged_df_genus_biplot = merged_df_genus.reset_index()
merged_df_genus_biplot = merged_df_genus_biplot.drop(columns=['index'])
bc_pcoa_genus.samples = bc_pcoa_genus.samples.reset_index()
bc_pcoa_genus.samples = bc_pcoa_genus.samples.drop(columns=['index'])

# Compute the projection of descriptors into a PCoA matrix
bc_pcoa_biplot = pcoa_biplot(bc_pcoa_genus, merged_df_genus_biplot)

# Extract biplot scores for features
biplot_scores = bc_pcoa_biplot.features

st.write(biplot_scores)

# Calculate the magnitude of loadings for each species
biplot_scores['Magnitude'] = np.sqrt(biplot_scores['PC1']**2 + biplot_scores['PC2']**2)

# Sort the species based on the magnitude and select the top 5
top5_species = biplot_scores.nlargest(5, 'Magnitude')

# Create a color map for your categories
category_colors = px.colors.qualitative.Dark24
color_map = {category: category_colors[i % len(category_colors)] 
             for i, category in enumerate(bc_pcoa_genus_data[update_figure(color_option)].unique())}

# Map your categories to colors
mapped_colors = bc_pcoa_genus_data[update_figure(color_option)].map(color_map)

# Create a new figure and add scatter points for the PCoA data
fig = go.Figure()

# Add the original scatter points from the PCoA data
fig.add_trace(go.Scatter(
    x=bc_pcoa_genus_data['PC1'],
    y=bc_pcoa_genus_data['PC2'],
    mode='markers',
    marker=dict(color=mapped_colors, size=6),
    text=bc_pcoa_genus_data['study_id'],  # hover text
    hoverinfo='text'
))

# Add arrows for top 5 species
for species in top5_species.index:
    fig.add_trace(go.Scatter(
        x=[0, top5_species.loc[species, 'PC1']], 
        y=[0, top5_species.loc[species, 'PC2']],
        mode='lines+markers+text',
        text=[None, species],
        textposition="top center",
        line=dict(color='green', width=2),
        showlegend=False
    ))

# Update layout with titles and axis labels
fig.update_layout(
    title=f'PCoA plot (Bray Curtis distance) of the analyses from all studies at Genus level',
    xaxis=dict(title=f'PCo1 ({explained_var_ratio[0]:.2%})'),
    yaxis=dict(title=f'PCo2 ({explained_var_ratio[1]:.2%})'),
    legend_title=dict(text=color_option)
)

# Show the plot in Streamlit
st.plotly_chart(fig, use_container_width=True)