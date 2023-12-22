#%%
import os
import glob
import pandas as pd
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.tree import nj
from skbio import TreeNode
from skbio import DistanceMatrix
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "vscode"

# Function to load abundance table for a specific study
def load_abund_table(selected_study):
    # Set the file name 
    folder_path = f"{selected_study}/"
    filename = os.path.join(folder_path, f"{selected_study}*phylum_taxonomy*.tsv")

    # Use glob to find files that match the pattern
    file = glob.glob(filename)

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


# Load and preprocess abundance table for a specific study
abund_table = load_abund_table('MGYS00005846')
abund_table = preprocess_abund_table(abund_table)

# Transpose abundance table
abund_table = abund_table.set_index('phylum').T

# Extract OTU names and sample names as numpy arrays
otu_names = list(abund_table.columns.values)
sample_names = list(abund_table.index.values)

# Convert abundance table to numpy array
abund_table = abund_table.to_numpy()
#%%
# Obtain bray-curtis distance matrix
bc_mat = beta_diversity("braycurtis", abund_table, sample_names)

# Run PCoA
bc_pcoa = pcoa(bc_mat)

# Extract the data to plot the PcoA
bc_pcoa_data = bc_pcoa.samples[['PC1', 'PC2']]

# Create a 2d plot with plotly
fig = px.scatter(bc_pcoa_data, x='PC1', y='PC2', hover_name=bc_pcoa_data.index)
fig.show()

#%%
# Create phylogenetic tree
dm = DistanceMatrix(bc_mat, sample_names)
tree = nj(dm)
rooted_tree = tree.root_at_midpoint()

# Obtain the weighted unifrac distance matrix
weighunifrac_mat = beta_diversity("weighted_unifrac", abund_table, tree = rooted_tree,
                                    ids = sample_names, otu_ids = otu_names)

#%%
import pandas as pd
# Load data 
study_info_all = pd.read_csv('Mgnify_studies_wwt_shot_metag_assembly_more10samples.csv')

# Select the study of interest
# study_info_filt = study_info[study_info['study_id'] == 'MGYS00005846']

# Create list with 2 studies of interest
studies = ['MGYS00005846', 'MGYS00005847']

# Create empty df
study_info_filt = pd.DataFrame()

# Create a df to summarize the study information
for study in studies:
    print(study)
    
    study_info = study_info_all[study_info_all['study_id'] == study]

    study_info = study_info.T

    study_info.columns = [study]

    study_info = study_info.drop(index='study_id')

    #study_info_filt[study] = study_info[study]

    print(type(study_info))

    # study_info_filt = study_info_filt.append(study_info)

# Create a df to summarize the study information
# summ_df = pd.DataFrame(study_info_filt).T

# Remove study_id row
# summ_df = summ_df.drop(index='study_id')

# Replace indices and column names
#summ_df.index = ['Description', 'Bioproject ID', 'Centre name', 'Number of samples', 
 #                'Biome','Experiment type', 'Pipeline version', 'Sampling country']

# summ_df.columns = ['Summary metrics']

# %%
