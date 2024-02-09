# Web app
import streamlit as st
import pandas as pd

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

# Function to load all data
@st.cache_data
def load_studies_data(filename):
    all_data = pd.read_csv(filename)

    # Rename the entries in the biomes column
    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Nutrient removal:Dissolved organics (anaerobic)", 
                                                                                "root:Engineered:Wastewater:Nutrient removal")
    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Nutrient removal:Biological phosphorus removal:Activated sludge", 
                                                                                "root:Engineered:Wastewater:Nutrient removal")
    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Nutrient removal:Nitrogen removal", 
                                                                                "root:Engineered:Wastewater:Nutrient removal")

    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Industrial wastewater:Petrochemical", 
                                                                                "root:Engineered:Wastewater:Industrial wastewater")
    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Industrial wastewater:Agricultural wastewater", 
                                                                                "root:Engineered:Wastewater:Industrial wastewater")

    all_data["biomes"] = all_data["biomes"].replace("root:Engineered:Wastewater:Activated Sludge, root:Engineered:Wastewater:Industrial wastewater", 
                                                                                "root:Engineered:Wastewater:Activated Sludge")
    
    # Remove the "root:Engineered:" part of the biomes column entries
    all_data["biomes"] = all_data["biomes"].str.replace("root:Engineered:", "")

    # Drop Unnamed columns
    all_data = all_data.drop(all_data.columns[all_data.columns.str.contains('unnamed',case = False)],axis = 1)

    return all_data

# Loading the Mgnify waste water treatment studies data
# This will be done only once and all the pages will have access to the data
if 'studies_data' not in st.session_state:
    st.session_state.studies_data = load_studies_data("wwt_studies_more_median_taxa_bybiome.csv")

# Add a title and info about the app
st.title('Summary and EDA of waste water treatment studies from Mgnify combined by biomes')

with st.expander('About this app'):
    st.write('''
    This app shows a summary and teh exploratory data analysis of the waste water treatment studies from Mgnify with more than 10 samples. It includes informations about the studies and its samples.
    
    **Credits**
    - Developed by [Sebastián Ayala Ruano](https://sayalaruano.github.io/).
    - The data was obtained using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI).
      ''')

st.subheader('Welcome!')
st.info('Look at the summary of all studies, the comparison between studies, or a detailed information for every study', icon='👈')

st.sidebar.header('Data')
st.sidebar.write('The data used in this app was obtained from [Mgnify](https://www.ebi.ac.uk/metagenomics/api/latest/) using the scripts available in this [GitHub repository](https://github.com/Multiomics-Analytics-Group/Retrieve_info_MGnifyAPI)')

st.sidebar.header('Code availability')
st.sidebar.write('The code for this project is available under the [MIT License](https://mit-license.org/) in this [GitHub repo](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies). If you use or modify the source code of this project, please provide the proper attributions for this work.')

st.sidebar.header('Contact')
st.sidebar.write('If you have any comments or suggestions about this work, please [create an issue](https://github.com/sayalaruano/Dashboard_MGnify_wwt_studies/issues/new) in the GitHub repository of this project.')
    


