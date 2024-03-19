#%%
import pandas as pd
import numpy as np
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
        # List of expected taxonomic levels
        expected_tax_levels = ['kingdom', 'phylum']
    
    else:
        # List of expected taxonomic levels
        expected_tax_levels = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Family_Genus']

        # Add missing taxonomic levels with empty values
        for i, tax_level in enumerate(expected_tax_levels):
            if tax_level not in abund_table.columns:
                abund_table.insert(i, tax_level, "")

    # Add the taxonomic columns to a different df and drop them from the original abundance table
    taxonomic_df = pd.DataFrame()
    for tax_level in expected_tax_levels:
        if tax_level in abund_table.columns:
            taxonomic_df[tax_level] = abund_table[tax_level]
            abund_table = abund_table.drop(columns=[tax_level])
            
    # Define the column names for this format
    taxonomic_df.columns  = expected_tax_levels

    return abund_table, taxonomic_df

# Function to load data for a specific study
def load_sample_data(selected_study):    
    # Load sample information for the study
    sample_info = pd.read_csv(f"Samples_metadata/{selected_study}/{selected_study}_samples_metadata.csv")
    
    return sample_info

def create_consensus_column(df, column_name):
    # Filter columns that have the same name as 'column_name'
    columns = df.filter(like=column_name).columns

    # Use the first non-NaN value found across the columns for each row
    consensus_column = df[columns].apply(lambda x: next((item for item in x if pd.notna(item)), np.nan), axis=1)

    return consensus_column

def update_taxname_for_conflicts(abund_table, tax_rank):
    """
    Update the tax_rank column in the abundance table to include higher_tax_rank values for entries with 
    the same tax_rank but different higher_tax_rank.
    The updated name is created by appending the higher taxonomic rank to the tax_rank name for conflicting entries.
    Input: abund_table (DataFrame) - DataFrame with the abundance table for the study
           tax_rank (str) - taxonomic rank to preprocess
    Output: abund_table (DataFrame) - DataFrame with the abundance table for the study after updating the tax_rank column
    """
    # Find the index of the tax_rank column
    tax_rank_index = abund_table.columns.get_loc(tax_rank)

    # Get the higher taxonomic rank's name
    higher_tax_rank = abund_table.columns[tax_rank_index - 1]

    # Create a temporary column to check for conflicts without affecting the original tax_rank column
    abund_table['temp_tax_rank'] = abund_table[tax_rank]
    
    # Group by the tax_rank and check for duplicates in the higher_tax_rank within each group
    for tax_name, group in abund_table.groupby(tax_rank):
        if len(group[higher_tax_rank].unique()) > 1:  # More than one unique higher_tax_rank for this tax_rank
            # Update the temporary tax_rank column for these entries
            abund_table.loc[group.index, 'temp_tax_rank'] = group[tax_rank] + '_' + group[higher_tax_rank].astype(str)

    # Remove the duplicate rows based on higher_tax_rank and temp_tax_rank
    abund_table.drop_duplicates(subset=[higher_tax_rank, 'temp_tax_rank'], keep='last', inplace=True)

    # Replace the tax_rank column with the updated values from temp_tax_rank
    abund_table[tax_rank] = abund_table['temp_tax_rank']
    abund_table.drop(columns=['temp_tax_rank'], inplace=True)

    # Replace "_nan" with empty string
    abund_table[tax_rank] = abund_table[tax_rank].str.replace('_nan', '')

    return abund_table

# Separate the studies by their biomes
# Load the summary table of all studies
studies_data = pd.read_csv('wwt_studies_more_median_taxa_bybiome.csv')

# Remove rows with nan values in all columns
studies_data = studies_data.dropna(how='all')

# Extract the study ids from the csv file
wwt_studies_allsamples_list = studies_data["study_id"].tolist()

# Remove nan values from the list
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

# Initialize global lists
global_abund_dfs_phylum = []
global_abund_dfs_genus = []
global_taxa_dfs_phylum = []
global_taxa_dfs_genus = []
global_metadata_dfs = []
global_studies_per_sample = []

# Process and merge data for each biome
for biome, study_ids in biome_dict.items():
    abund_dfs_phylum = []
    abund_dfs_genus = []
    taxa_dfs_phylum = []
    taxa_dfs_genus = []
    metadata_dfs = []
    studies_per_sample = []

    for selected_study in study_ids:
        # Process phylum-level data
        abund_table_phylum = load_abund_table(selected_study, phylum=True)
        if abund_table_phylum is not None:
            # Preprocess the abundance table
            processed_abund_phylum, processed_taxa_phylum = preprocess_abund_table(abund_table_phylum, phylum=True)
            
            # Set the index to the phylum column for both dataframes
            processed_taxa_phylum.set_index('phylum', inplace=True)
            processed_abund_phylum.set_index(processed_taxa_phylum.index, inplace=True)
            
            # Append the processed dataframes to the lists
            abund_dfs_phylum.append(processed_abund_phylum)
            taxa_dfs_phylum.append(processed_taxa_phylum)

        # Process genus-level data
        abund_table_genus = load_abund_table(selected_study, phylum=False)
        if abund_table_genus is not None:
            # Preprocess the abundance table
            processed_abund_genus, processed_taxa_genus = preprocess_abund_table(abund_table_genus, phylum=False)
            
            # Set the index to the Family_Genus column for both dataframes
            processed_taxa_genus.set_index('Family_Genus', inplace=True)
            processed_abund_genus.set_index(processed_taxa_genus.index, inplace=True)

            # Add a row with the study ID for all samples 
            study_id_row_genus = pd.Series(selected_study, index=processed_abund_genus.columns)
            studies_per_sample.append(study_id_row_genus)

            # Append the processed dataframes to the lists
            abund_dfs_genus.append(processed_abund_genus)
            taxa_dfs_genus.append(processed_taxa_genus)

        # Load and process sample metadata
        sample_metadata_df = load_sample_data(selected_study)
        if sample_metadata_df is not None:
            sample_metadata_df['study_id'] = selected_study
            metadata_dfs.append(sample_metadata_df)

    # Concatenate and aggregate phylum-level abundance data
    if abund_dfs_phylum:
        # Concatenate data from all studies, fill NaNs with 0, and convert to integers
        merged_abund_df_phylum = pd.concat(abund_dfs_phylum, axis = 1)
        merged_abund_df_phylum = merged_abund_df_phylum.fillna(0).astype(int)

        # Move index to a column
        merged_abund_df_phylum.reset_index(inplace=True)

        # Create a new index with OTU IDs
        otu_index = ['OTU' + str(i+1) for i in range(len(merged_abund_df_phylum))]
        merged_abund_df_phylum.index = otu_index
        merged_abund_df_phylum.index.name = 'OTU'

        # Save the DataFrame to a csv file
        merged_abund_df_phylum.to_csv(f'Abundance_tables/Merged_tables/{biome}/{biome}_merged_abund_tables_phylum.csv')

    # Concatenate and aggregate genus-level abundance data
    if abund_dfs_genus:
        # Concatenate data from all studies, fill NaNs with 0, and convert to integers
        merged_abund_df_genus = pd.concat(abund_dfs_genus, axis = 1)
        merged_abund_df_genus = merged_abund_df_genus.fillna(0).astype(int)

        # Create a new index with OTU IDs
        otu_index = ['OTU' + str(i+1) for i in range(len(merged_abund_df_genus))]
        merged_abund_df_genus.index = otu_index
        merged_abund_df_genus.index.name = 'OTU'

        # Save the DataFrame to a csv file
        merged_abund_df_genus.to_csv(f'Abundance_tables/Merged_tables/{biome}/{biome}_merged_abund_tables_genus.csv')

        # Concatenate the study IDs for each sample
        studies_per_sample_genus_df = pd.concat(studies_per_sample, axis=0)

        # Reset the index and rename the columns
        studies_per_sample_genus_df = studies_per_sample_genus_df.reset_index()
        studies_per_sample_genus_df.columns = ['assembly_run_ids', 'study_id']

        # Save the DataFrame to a csv file
        studies_per_sample_genus_df.to_csv(f'Abundance_tables/Merged_tables/{biome}/{biome}_studies_per_sample.csv', index=False)

    # Concatenate and deduplicate phylum-level taxonomic data
    if taxa_dfs_phylum:
        merged_taxa_df_phylum = pd.concat(taxa_dfs_phylum, axis=1)

        # Identify duplicated column names
        duplicated_columns = merged_taxa_df_phylum.columns[merged_taxa_df_phylum.columns.duplicated()]

        # Create a dictionary to hold the consensus columns
        consensus_columns_phylum = {}
                
        # Create a consensus column for each duplicated column and store it in the dictionary
        for column in duplicated_columns.unique():
            consensus_columns_phylum[column + '_consensus'] = create_consensus_column(merged_taxa_df_phylum, column)

        # Convert the dictionary to a DataFrame and concatenate it with the original DataFrame
        consensus_df = pd.DataFrame(consensus_columns_phylum)
        merged_taxa_df_phylum = pd.concat([merged_taxa_df_phylum, consensus_df], axis=1)

        # Drop the original duplicated columns
        merged_taxa_df_phylum = merged_taxa_df_phylum.drop(columns=duplicated_columns)

        # Remove '_consensus' from the column names
        merged_taxa_df_phylum.columns = [col.replace('_consensus', '') for col in merged_taxa_df_phylum.columns]

        # Move index to a column
        merged_taxa_df_phylum.reset_index(inplace=True)
        
        # Create a new index with OTU IDs
        otu_index = ['OTU' + str(i+1) for i in range(len(merged_taxa_df_phylum))]
        merged_taxa_df_phylum.index = otu_index
        merged_taxa_df_phylum.index.name = 'OTU'

        # Save the DataFrame to a csv file
        merged_taxa_df_phylum.to_csv(f'Abundance_tables/Merged_tables/{biome}/{biome}_merged_taxa_tables_phylum.csv')

    # Concatenate and deduplicate genus-level taxonomic data
    if taxa_dfs_genus:
        merged_taxa_df_genus = pd.concat(taxa_dfs_genus, axis=1)

        # Identify duplicated column names
        duplicated_columns = merged_taxa_df_genus.columns[merged_taxa_df_genus.columns.duplicated()]

        # Create a dictionary to hold the consensus columns
        consensus_columns_genus = {}
                
        # Create a consensus column for each duplicated column and store it in the dictionary
        for column in duplicated_columns.unique():
            consensus_columns_genus[column + '_consensus'] = create_consensus_column(merged_taxa_df_genus, column)

        # Convert the dictionary to a DataFrame and concatenate it with the original DataFrame
        consensus_df = pd.DataFrame(consensus_columns_genus)
        merged_taxa_df_genus = pd.concat([merged_taxa_df_genus, consensus_df], axis=1)

        # Drop the original duplicated columns
        merged_taxa_df_genus = merged_taxa_df_genus.drop(columns=duplicated_columns)

        # Remove '_consensus' from the column names
        merged_taxa_df_genus.columns = [col.replace('_consensus', '') for col in merged_taxa_df_genus.columns]

        # Apply function to update taxonomic names for conflicts
        merged_taxa_df_genus = update_taxname_for_conflicts(merged_taxa_df_genus, 'Genus')

        # Create a new index with OTU IDs
        otu_index = ['OTU' + str(i+1) for i in range(len(merged_taxa_df_genus))]
        merged_taxa_df_genus.index = otu_index
        merged_taxa_df_genus.index.name = 'OTU'

        # Save the DataFrame to a csv file
        merged_taxa_df_genus.to_csv(f'Abundance_tables/Merged_tables/{biome}/{biome}_merged_taxa_tables_genus.csv')

    # Concatenate sample metadata
    if metadata_dfs:
        merged_metadata_df = pd.concat(metadata_dfs)
        merged_metadata_df.to_csv(f'Samples_metadata/Merged_tables/{biome}_merged_samples_metadata.csv', index=False)

    # Append the biome-specific data to the global lists
    global_abund_dfs_phylum.extend(abund_dfs_phylum)
    global_abund_dfs_genus.extend(abund_dfs_genus)
    global_taxa_dfs_phylum.extend(taxa_dfs_phylum)
    global_taxa_dfs_genus.extend(taxa_dfs_genus)
    global_metadata_dfs.extend(metadata_dfs)
    global_studies_per_sample.extend(studies_per_sample)


# Concatenate and aggregate phylum-level abundance data for all biomes
if global_abund_dfs_phylum:
    merged_global_abund_df_phylum = pd.concat(global_abund_dfs_phylum, axis=1)
    merged_global_abund_df_phylum = merged_global_abund_df_phylum.fillna(0).astype(int)
    otu_index = ['OTU' + str(i+1) for i in range(len(merged_global_abund_df_phylum))]
    merged_global_abund_df_phylum.index = otu_index
    merged_global_abund_df_phylum.index.name = 'OTU'
    merged_global_abund_df_phylum.reset_index(inplace=True)
    merged_global_abund_df_phylum.to_csv('Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_abund_tables_phylum.csv', index=False)

# Concatenate and aggregate genus-level abundance data for all biomes
if global_abund_dfs_genus:
    merged_global_abund_df_genus = pd.concat(global_abund_dfs_genus, axis=1)
    merged_global_abund_df_genus = merged_global_abund_df_genus.fillna(0).astype(int)
    otu_index = ['OTU' + str(i+1) for i in range(len(merged_global_abund_df_genus))]
    merged_global_abund_df_genus.index = otu_index
    merged_global_abund_df_genus.index.name = 'OTU'
    merged_global_abund_df_genus.reset_index(inplace=True)
    merged_global_abund_df_genus.to_csv('Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_abund_tables_genus.csv', index=False)

# Concatenate and deduplicate phylum-level taxonomic data for all biomes
if global_taxa_dfs_phylum:
    merged_global_taxa_df_phylum = pd.concat(global_taxa_dfs_phylum).drop_duplicates()
    merged_global_taxa_df_phylum.reset_index(inplace=True, drop=True)
    merged_global_taxa_df_phylum.to_csv('Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_taxa_tables_phylum.csv', index=False)

# Concatenate and deduplicate genus-level taxonomic data for all biomes
if global_taxa_dfs_genus:
    merged_global_taxa_df_genus = pd.concat(global_taxa_dfs_genus).drop_duplicates()
    merged_global_taxa_df_genus.reset_index(inplace=True, drop=True)
    merged_global_taxa_df_genus.to_csv('Abundance_tables/Merged_tables/All_biomes/all_biomes_merged_taxa_tables_genus.csv', index=False)

# Concatenate sample metadata for all biomes
if global_metadata_dfs:
    merged_global_metadata_df = pd.concat(global_metadata_dfs).drop_duplicates()
    merged_global_metadata_df.reset_index(inplace=True, drop=True)
    merged_global_metadata_df.to_csv('Samples_metadata/Merged_tables/all_biomes_merged_samples_metadata.csv', index=False)

# Concatenate the studies per sample data for all biomes
if global_studies_per_sample:
    # Concatenate all study_id_row_genus DataFrames
    merged_global_studies_per_sample_df = pd.concat(global_studies_per_sample, axis=0)

    # Reset the index and rename the columns
    merged_global_studies_per_sample_df = merged_global_studies_per_sample_df.reset_index()
    merged_global_studies_per_sample_df.columns = ['assembly_run_ids', 'study_id']
    merged_global_studies_per_sample_df.to_csv('Abundance_tables/Merged_tables/All_biomes/all_biomes_studies_per_sample.csv', index=False)
#%%