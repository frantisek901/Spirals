#!/usr/bin/env python
# coding: utf-8

# # Stepwise Data Fitting
# 
# This is an attempt to incrementally analyze the model as we complicate it slowly by adding new dynamics (except for identity which has been left out here). 
# 
# 
# ## Models and Simulations
# Every simulation run has 1000 agents that represent humans, and the richest model has 10 media agents as well as a variable media influence factor. For the first three steps, every parameter combination is run with 30 random initializations. Steps 4 and 5 have been run with 5 random initializations.
# 
# The following models were tried:
# 
# 1. Simple Hegselmann-Krause with full network.
# 2. Simple Hegselmann-Krause with scale-free network.
# 3. Hegselmann-Krause with normally distributed open-mindedness and scale-free network.
# 3. Hegselmann-Krause with normally distributed open-mindedness, normally distributed initial opinion, and scale-free network.
# 5. Hegselmann-Krause with normally distributed open-mindedness, normally distributed initial opinion, media agents ($MedN = 10$), variable media distribution parameters, variable media influence factor, and scale-free network.
# 6. Hegselmann-Krause with normally distributed open-mindedness, normally distributed initial opinion, media agents ($MedN = 10$), variable media distribution parameters, variable media influence factor, silence mechanisms with memory, and scale-free network.
# 
# **Notes:** 
# 1. The media influence factor  $MedInf \epsilon [0, 1]$ and is a weighting factor for media-sourced influence. 
# So $MedInf = 0$ implies no influence from media sources, while $MedInf = 1$ implies that media-sourced influence is as strong as other people's opinions.
# 2. The scale-free network parameter is set to 1.
# 3. For Step 4, initial opinion distribution parameters were fixed to values of $\mu_{InitialPopOp} = 0$ and $\sigma_{InitialPopOp} = 0.8$ after pilot simulations showed significantly good fits for these values.
# 4. Media opinion positions are drawn from a deterministic set of positions approximating the normal pdf, and have the parameters $\mu_{MediaOp}$ and $\sigma_{MediaOp}$.
# 
# | Model    | Parameters | Random Initializations | Total Simulations |
# | -------- | ------- | ------- | ------- |
# | Step 1  | $\varepsilon$ | 30 | | 
# | Step 2  | $\varepsilon$ | 30 | | 
# | Step 3  | $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$| 30 | | 
# | Step 4  | $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$, $\sigma_{InitialOp}$| 5 | | 
# | Step 5  | $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$,  $\mu_{MediaOp}$, $\sigma_{MediaOp}$, $MedInf$| 5 | | 
# | Step 6  | $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$,  $\mu_{MediaOp}$, $\sigma_{MediaOp}$, $MedInf$, $Silence_{\tau}$,  $Silence_{\delta_{0}}$,  $Silence_{\alpha}$| 2 | | 
# 
# ## Survey Data
# 
# <>
# 
# ### Parlemeter Data - Initial fitting patterns
# 
# <>
# 
# ### Dynamic Data - 
# 
# <>
# 
# ## Fitting Strategy
# 
# <>
# 
# 
# 

# # Preprocessing and Saving Simulation Data

# ## Globals (SET EACH TIME)

# In[1]:


current_runs_title = "15.09.25"

# Set below to True if rerunning after saving appropriately in //data//preprocessed//<current_runs_title>
load_preprocessed_data = True

# Set below to True if there is a need to preprocess data (such as while running for the first time for new data)
preprocess_data = True

# Set below to True if running for the first time or need to preprocess and save for some other reason
save_preprocessed_data = True

save_fits = True

# set a smaller processing unit if it benefits you, else keep it to range(1, 5)
steps_to_process =  range(1, 7)

print("Running stepwiseDataFitting_Cluster.py")
print(f"The following steps will be processed: \n")
for i in (steps_to_process):
    print(f"Step {i} \n")

# For the folder names for each step
step_input_folder_name = []
step_output_folder_name = []

# Specific sub-folder names (folder name within data/cluster/endsim)
step_input_folder_name.append("Step1") 
step_input_folder_name.append("Step2") 
step_input_folder_name.append("Step3") 
step_input_folder_name.append("Step4_NEWEST") 
step_input_folder_name.append("Step5_New") 
step_input_folder_name.append("Step6_New") 

step_output_folder_name = step_input_folder_name

step_titles = ["Simple HK with Full Network", "Simple HK with Scale-Free Network", "Heterogenous Openness", "Heterogenous Openness with normal opinion distributions", "Media and Influence", "Media with Silence"]


# In[2]:


import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import seaborn as sns
from scipy.spatial.distance import jensenshannon


# ## Paths

# In[3]:


# Where the preprocessed data will be saved
full_preprocessed_data_folder_path = os.path.join("data", "preprocessed", "cluster", current_runs_title)


# Wheer the output plots will be saved
plots_folder_path = os.path.join("analysis", "plots", "cluster", current_runs_title)


# ## Issue - Step 4 has a LOT OF redundancies and no media
# 
# I accidentally ran this step with 0 media agents, so we have a great number of repeat simulations and no media in this step.
# We will rectify this by adding a step 5 with media and rerunning the simulations with the new setting.
# 
# However, for the time being we need to remove the extra simulations.

# In[4]:


# # Define the regex pattern to extract parameters
# pattern = r"EndSim_epsM(-?[\d\.]+)_epsSD(-?[\d\.]+)___OpD([a-zA-Z-]+)_OpM(-?[\d\.]+)_OpSD(-?[\d\.]+)___Net([a-zA-Z-]+)___NAgents(\d+)___RS(\d+)__MedInfF([\d\.]+)___MedD([a-zA-Z-]+)_MedN(\d+)_MedM(-?[\d\.]+)_MedSD(-?[\d\.]+)"

# def extract_eps_params(filename):
#     """Extract epsM and epsSD parameters from filename"""
#     base_name = os.path.splitext(filename)[0]
#     match = re.search(pattern, base_name)
#     if match:
#         epsM = float(match.group(1))
#         epsSD = float(match.group(2))
#         RS = int(match.group(8))
#         return (epsM, epsSD, RS)
#     return None

# # Alternative: Move files instead of deleting (safer approach)
# def move_redundant_files(directory, backup_dir="redundant_backup"):
#     """Move redundant files to a backup directory instead of deleting"""
#     # Create backup directory if it doesn't exist
#     backup_path = os.path.join(directory, backup_dir)
#     os.makedirs(backup_path, exist_ok=True)

#     # Get all CSV files
#     csv_files = [f for f in os.listdir(directory) if f.endswith('.csv')]

#     unique_combinations = {}
#     files_to_keep = []
#     files_to_move = []

#     # Process each file
#     for filename in csv_files:
#         params = extract_eps_params(filename)
#         if params:
#             epsM, epsSD, RS = params
#             combination = (epsM, epsSD, RS)

#             if combination not in unique_combinations:
#                 unique_combinations[combination] = filename
#                 files_to_keep.append(filename)
#             else:
#                 files_to_move.append(filename)
#         else:
#             print(f"Warning: Could not parse parameters from {filename}")

#     # Move redundant files
#     moved_count = 0
#     for filename in files_to_move:
#         src_path = os.path.join(directory, filename)
#         dest_path = os.path.join(backup_path, filename)
#         try:
#             os.rename(src_path, dest_path)
#             moved_count += 1
#             # print(f"Moved to backup: {filename}")
#         except Exception as e:
#             print(f"Error moving {filename}: {e}")

#     # Print summary
#     print(f"\nSummary:")
#     print(f"Total files processed: {len(csv_files)}")
#     print(f"Unique (epsM, epsSD, RS) combinations: {len(unique_combinations)}")
#     print(f"Files kept: {len(files_to_keep)}")
#     print(f"Files moved to backup: {moved_count}")

#     return files_to_keep, files_to_move


# move_redundant_files("data/cluster/endsim/step4_5reps")


# ##  Reading and Preprocessing Simulation Data

# In[5]:


# Regex pattern for all files
pattern = (
    r"EndSim_epsM(-?[\d\.]+)_epsSD(-?[\d\.]+)"
    r"___OpD([a-zA-Z-]+)_OpM(-?[\d\.]+)_OpSD(-?[\d\.]+)"
    r"___Net([a-zA-Z-]+)___NAgents(\d+)___RS(\d+)"
    r"__MedInfF([\d\.]+)___MedD([a-zA-Z-]+)_MedN(\d+)_MedM(-?[\d\.]+)_MedSD(-?[\d\.]+)"
    r"(?:_Silence_Alpha(-?[\d\.]+)_Silence_Tau(-?[\d\.]+)_Silence_Delta0(-?[\d\.]+)_SilenceByBoundary\?([a-zA-Z]+))?"
)

def parse_filename(filename):
    """Extract parameters from filename using regex with optional silence parameters"""
    # Remove .csv extension before parsing
    base_name = os.path.splitext(filename)[0]
    match = re.search(pattern, base_name)

    print("Now parsing file: \n")
    print(filename)
    print("\n")

    if match:
        result = {
            'epsM': float(match.group(1)),
            'epsSD': float(match.group(2)),
            'OpD': match.group(3),
            'OpM': float(match.group(4)),
            'OpSD': float(match.group(5)),
            'NetworkType': match.group(6),
            'NAgents': int(match.group(7)),
            'RandomSeed': int(match.group(8)),
            'MedInfF': float(match.group(9)),
            'MedD': match.group(10),
            'MedN': int(match.group(11)),
            'MedM': float(match.group(12)),
            'MedSD': float(match.group(13))
        }
        print(f"MedSD should be {match.group(13)}")

        # Check if silence parameters are present (group 14 is the start of silence block)
        if match.group(14) is not None:  # The entire silence block exists
            # Group indices after the silence block:
            # group 14: Silence_Alpha value (e.g., "0.79")
            # group 15: Silence_Tau value (e.g., "5")  
            # group 16: Silence_Delta0 value (e.g., "0")
            # group 17: SilenceByBoundary value (e.g., "false")
            
            print(f"Silence_Alpha value string: {match.group(14)}")
            print(f"Silence_Tau value string: {match.group(15)}")
            print(f"Silence_Delta0 value string: {match.group(16)}")
            print(f"SilenceByBoundary value string: {match.group(17)}")
            
            result.update({
                'Silence_Alpha': float(match.group(14)),  # Was group 15, now 14
                'Silence_Tau': float(match.group(15)),   # Was group 16, now 15
                'Silence_Delta0': float(match.group(16)), # Was group 17, now 16
                'SilenceByBoundary': match.group(17).lower() == 'true'  # Was group 18, now 17
            })
        else:
            # Set default values for silence parameters if not present
            result.update({
                'Silence_Alpha': None,
                'Silence_Tau': None,
                'Silence_Delta0': None,
                'SilenceByBoundary': None
            })

        return result
    return None

    
def extract_opinion_sections(content):
    """Extract media and individual opinion sections from file content"""
    # Split content into sections using separator lines
    sections = re.split(r'-{10,}', content)

    # We expect at least 4 sections:
    # 0: Metadata
    # 1: Network connections (skip)
    # 2: Media opinions
    # 3: Individual opinions
    if len(sections) < 4:
        return None, None

    # Extract media opinions section
    media_section = sections[2].strip()
    # Extract individual opinions section
    individual_section = sections[3].strip()

    return media_section, individual_section

def parse_opinion_section(section):
    """Parse a single opinion section into a DataFrame"""
    if not section:
        return pd.DataFrame()

    # Clean section by removing any empty lines and brackets
    cleaned_lines = []
    for line in section.split('\n'):
        if line.strip() and not line.startswith('[') and not line.startswith('--'):
            cleaned_lines.append(line)

    if not cleaned_lines:
        return pd.DataFrame()

    # Recreate CSV content
    csv_content = "\n".join(cleaned_lines)

    # Parse as CSV
    try:
        return pd.read_csv(StringIO(csv_content))
    except:
        return pd.DataFrame()

def process_simulation_file(filepath):
    """Process a single simulation file with section handling"""
    # Extract parameters from filename
    filename = os.path.basename(filepath)
    params = parse_filename(filename)
    if not params:
        print(f"Skipping {filename} - pattern not matched")
        return None

    try:
        # Read entire file content
        with open(filepath, 'r') as f:
            content = f.read()
        # Extract opinion sections
        media_section, individual_section = extract_opinion_sections(content)

        # Parse media opinions
        media_df = parse_opinion_section(media_section)

        media_opinions = []
        # print(media_df.columns)

        if not media_df.empty and ' opinion)' in media_df:
            media_opinions = media_df[' opinion)'].apply(
                lambda x: float(str(x).strip('[]'))).values

        # Parse individual opinions
        individual_df = parse_opinion_section(individual_section)
        individual_opinions = []
        if not individual_df.empty and 'opinion' in individual_df:
            individual_opinions = individual_df['opinion'].apply(
                lambda x: float(str(x).strip('[]'))).values


        individual_initial_opinions = []
        if not individual_df.empty and 'initialOpinion' in individual_df:
            individual_initial_opinions = individual_df['initialOpinion'].apply(
                lambda x: float(str(x).strip('[]'))).values

        # Add to params
        params['individual_opinions'] = np.array(individual_opinions)
        params['individual_initial_opinions'] = np.array(individual_initial_opinions)
        params['media_opinions'] = np.array(media_opinions)

        # Add static parameters
        params['NAgents'] = 1000
        params['MedN'] = 20
        params['NetType'] = "scale-free"
        params['ScaleFreeDegree'] = 1

        return params

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")
        return None

def load_all_simulations(directory):
    """Process all simulation files in a directory"""
    all_data = []
    processed_count = 0
    skipped_count = 0
    error_count = 0

    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            filepath = os.path.join(directory, filename)
            result = process_simulation_file(filepath)
            if result:
                all_data.append(result)
                processed_count += 1
            else:
                skipped_count += 1
            if(processed_count%1000 == 0):
                print(f"Processed {processed_count} sims for next step.")
    print(f"\n\nSummary: Processed={processed_count}, Skipped={skipped_count}")
    return pd.DataFrame(all_data) if all_data else pd.DataFrame()

def get_simulations_for_step(stepNo):
    simulations_df = load_all_simulations(os.path.join("data", "cluster", "endsim", step_input_folder_name[stepNo - 1]))

    if not simulations_df.empty:
        print(f"\nSuccessfully loaded {len(simulations_df)} simulations for Step {stepNo}")
        print("Parameter ranges:")
        for param in ['epsM', 'epsSD', 'OpM', 'OpSD', 'MedM', 'MedSD', 'MedInfF']:
            values = simulations_df[param]
            print(f"{param}: min={values.min()}, max={values.max()}")

        # Show example opinion data
        sample = simulations_df.iloc[0]
        print(f"\nSample simulation opinions:")
        print(f"Individual opinions: {len(sample['individual_opinions'])} values")
        print(f"Media opinions: {len(sample['media_opinions'])} values")
    else:
        print(f"No simulations loaded for Step {stepNo}. Please check file patterns and data quality.")

    return(simulations_df)


# ## Loading up Simulation Data into DataFrames

# In[6]:


# To Load simulation data
simulation_df_step = []    # The main list of Dataframes in which we will be loading the data



if(preprocess_data):
    print("-------------------- PREPROCESSING DATA --------------------------------")

    # Running all sims and storing them as DF's
    for stepNo in steps_to_process:
        simulation_df_step.append(get_simulations_for_step(stepNo))
        print(f"-------------------- PREPROCESSED STEP #{stepNo} DATA --------------------------------")

    
print(f"-------------------- PREPROCESSED ALL STEPS DATA --------------------------------")


# In[7]:


simulation_df_step


# ## Saving Preprocessed Simulation Dataframes as CSV

# In[8]:


if(save_preprocessed_data and preprocess_data):
    print("-------------------- SAVING PREPROCESSED DATA --------------------------------")

    if not os.path.exists(full_preprocessed_data_folder_path):
        os.makedirs(full_preprocessed_data_folder_path)
    for stepNo in steps_to_process:
        file_name = f"simulated_data_Step{stepNo}.csv"
        file_path = os.path.join(full_preprocessed_data_folder_path, file_name)
        simulation_df_step[stepNo - 1].to_csv(file_path, index=False)
    print("-------------------- PREPROCESSED DATA SAVED --------------------------------")



# ## Loading Preprocessed Simulation Dataframes from CSV

# In[9]:


if(load_preprocessed_data):
    print("-------------------- LOADING PREPROCESSED DATA --------------------------------")

    del simulation_df_step
    simulation_df_step = [] # Don't do this unless we are loading subsequently, otherwise we will have an emptied list and much sadness to go with.
    for stepNo in steps_to_process:
        file_name = f"simulated_data_Step{stepNo}.csv"
        file_path = os.path.join(full_preprocessed_data_folder_path, file_name)
        df = pd.read_csv(file_path)

        # Need to convert the arrays to numpy arrays
        individual_opinions = []
        if not df.empty and 'individual_opinions' in df:
            df['individual_opinions'] = df['individual_opinions'].apply(
                lambda x: np.fromstring(x.strip('[]').replace('\n', ' '), sep=' ', dtype=float))

        individual_initial_opinions = []
        if not df.empty and 'individual_initial_opinions' in df:
            df['individual_initial_opinions'] = df['individual_initial_opinions'].apply(
                lambda x: np.fromstring(x.strip('[]').replace('\n', ' '), sep=' ', dtype=float))

        simulation_df_step.append(df) # append to the stepwise dataframe list
        print(f"---- Preprocessed data loaded for Step #{stepNo} -----")


# In[10]:


# simulation_df_step[5]['Silence_Tau']


# ## For Binning Simulation Data

# In[11]:
print("-------------------- BINNING SIMULATION DATA --------------------------------")


# Define bin edges for 10-point scale mapping
bin_edges = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# Create bin interval labels for display
bin_intervals = [
    "[-1.0, -0.8)",
    "[-0.8, -0.6)",
    "[-0.6, -0.4)",
    "[-0.4, -0.2)",
    "[-0.2, 0.0)",
    "[0.0, 0.2)",
    "[0.2, 0.4)",
    "[0.4, 0.6)",
    "[0.6, 0.8)",
    "[0.8, 1.0]"
]

def bin_opinions(opinions):
    """Bin opinions into 10-point scale categories"""
    # Digitize opinions (returns bin indices 1-10)
    bin_indices = np.digitize(opinions, bin_edges[1:-1], right=False)
    return bin_indices

def bin_distribution(opinions):
    """Compute normalized distribution across bins"""
    counts, _ = np.histogram(opinions, bins=bin_edges)
    return counts / counts.sum()


# In[12]:


# Adding binned columns to DF
# for  s in steps_to_process:
#     print(s)
for stepNo in steps_to_process:
    simulations_df = simulation_df_step[stepNo - 1]

    # Add binned opinions and distributions to DataFrame
    simulations_df['binned_opinions'] = simulations_df['individual_opinions'].apply(
        lambda x: bin_opinions(np.array(x))
    )

    simulations_df['binned_distribution'] = simulations_df['individual_opinions'].apply(
        lambda x: bin_distribution(np.array(x))
    )

    # Add bin interval labels as a constant column
    simulations_df['bin_intervals'] = [bin_intervals] * len(simulations_df)

    # Club the two central bins to account for anchoring in the real data
    simulations_df['binned_distributions_5_6clubbed'] = simulations_df['binned_distribution'].apply(
        lambda dist: np.concatenate([dist[:4], [dist[4] + dist[5]], dist[6:]])
    )

    # Repeating the same for the initial opinion distribution

    # Add binned opinions and distributions to DataFrame
    simulations_df['binned_initial_opinions'] = simulations_df['individual_initial_opinions'].apply(
        lambda x: bin_opinions(np.array(x))
    )

    simulations_df['binned_initial_distribution'] = simulations_df['individual_initial_opinions'].apply(
        lambda x: bin_distribution(np.array(x))
    )

    # Club the two central bins to account for anchoring in the real data
    simulations_df['binned_initial_distributions_5_6clubbed'] = simulations_df['binned_initial_distribution'].apply(
        lambda dist: np.concatenate([dist[:4], [dist[4] + dist[5]], dist[6:]])
    )



    # Display the updated DataFrame
    # print("\nDataFrame with binned opinions:")
    # print(simulations_df[['epsM', 'epsSD', 'OpM', 'OpSD', 'MedM', 'MedSD', 
    #                       'binned_opinions', 'binned_distribution', 'bin_intervals']].head())


# 

# In[13]:


# Adding a Simulation Index column to the dataframes


# Initialize cumulative index counter
cumulative_index = 0


for stepNo in steps_to_process:
    # Get the current DataFrame
    simulations_df = simulation_df_step[stepNo - 1]

    # Add model simulation index (starts at 1 for each DataFrame)
    simulations_df['model_simulation_index'] = range(1, len(simulations_df) + 1)

    # Add cumulative simulation index (continues from previous DataFrame)
    simulations_df['cumulative_simulation_index'] = range(cumulative_index + 1, 
                                                         cumulative_index + len(simulations_df) + 1)

    # Update cumulative index for next DataFrame
    cumulative_index += len(simulations_df)

    # Update the DataFrame in the list
    simulation_df_step[stepNo - 1] = simulations_df    



# In[14]:


for stepNo in steps_to_process:
    simulations_df = simulation_df_step[stepNo - 1]
    print(f"Step {stepNo}:")
    print(f"  Model indices: {simulations_df['model_simulation_index'].min()} to {simulations_df['model_simulation_index'].max()}")
    print(f"  Cumulative indices: {simulations_df['cumulative_simulation_index'].min()} to {simulations_df['cumulative_simulation_index'].max()}")
    print(f"  Number of simulations: {len(simulations_df)}")
    print()


# In[ ]:





# ## Pooled Histograms for each step

# In[15]:

print("-------------------- NOT GENERATING POOLED HISTOGRAMS --------------------------------")


# for stepNo in steps_to_process:
#     simulations_df = simulation_df_step[stepNo - 1]
#     # Create pooled histogram of all opinions
#     all_opinions = np.concatenate(simulations_df['individual_opinions'].values)
#     pooled_counts, _ = np.histogram(all_opinions, bins=bin_edges)
#     pooled_distribution = pooled_counts / pooled_counts.sum()

#     # Plot the pooled histogram
#     plt.figure(figsize=(12, 6))
#     plt.bar(range(1, 11), pooled_distribution, color='skyblue', edgecolor='black')
#     plt.xticks(range(1, 11), bin_intervals, rotation=45, ha='right')
#     plt.xlabel('Opinion Bins')
#     plt.ylabel('Proportion of Opinions')
#     plt.title('Pooled Opinion Distribution Across All Simulations in ' +  step_titles[stepNo - 1] + " (Step " +  str(stepNo)+ ")")
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.tight_layout()
#     plt.savefig('pooled_opinion_distribution.png', dpi=300)
    # plt.show()


# ## Generating Individual histograms for specific simulations

# In[16]:


print("-------------------- NOT GENERATING INDIVIDUAL HISTOGRAMS --------------------------------")


# import os
# import matplotlib.pyplot as plt
# import numpy as np
# from pathlib import Path

# # Define bin edges and labels
# bin_edges = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# bin_intervals = [
#     "[-1.0, -0.8)", "[-0.8, -0.6)", "[-0.6, -0.4)", "[-0.4, -0.2)", 
#     "[-0.2, 0.0)", "[0.0, 0.2)", "[0.2, 0.4)", "[0.4, 0.6)", 
#     "[0.6, 0.8)", "[0.8, 1.0]"
# ]

# def save_simulation_histogram(row, index, stepNum, histograms_save_path):
#     # histograms_save_path = os.join(histograms_save_path, step_output_folder_name[stepNum - 1])
#     """Create and save histogram for a single simulation with both initial and final opinions"""
#     # Extract parameters for title and filename based on step number
#     if stepNum == 1 or stepNum == 2:
#         # Only epsM for steps 1 and 2
#         filename = f"sim_{index:04d}_epsM{row['epsM']:.2f}.png".replace(".", "dot")
#         title = f"Opinion Distribution (Step {stepNum}): εμ={row['epsM']:.2f}"
#     elif stepNum == 3:
#         # epsM and epsSD for step 3
#         filename = f"sim_{index:04d}_epsM{row['epsM']:.2f}_epsSD{row['epsSD']:.2f}.png".replace(".", "dot")
#         title = f"Opinion Distribution (Step {stepNum}): εμ={row['epsM']:.2f}, εσ={row['epsSD']:.2f}"
#     else:
#         # All parameters for step 4 and beyond
#         filename = (
#             f"sim_{index:04d}_"
#             f"epsM{row['epsM']:.2f}_epsSD{row['epsSD']:.2f}_"
#             f"OpM{row['OpM']:.2f}_OpSD{row['OpSD']:.2f}_"
#             f"MedM{row['MedM']:.2f}_MedSD{row['MedSD']:.2f}__"
#             f"MedInf{row['MedInfF']:.2f}.png"
#         ).replace(".", "dot")
#         title = (
#             f"Opinion Distribution (Step {stepNum}): "
#             f"εμ={row['epsM']:.2f}, εσ={row['epsSD']:.2f}, "
#             f"Medμ={row['MedM']:.2f}, Medσ={row['MedSD']:.2f}, "
#             f"MedInfF={row['MedInfF']:.2f}"
#         )

#     # Create figure
#     plt.figure(figsize=(12, 7))

#     # Plot histograms for both initial and final opinions
#     n_initial, _, _ = plt.hist(
#         row['individual_initial_opinions'],
#         bins=bin_edges,
#         density=False,
#         alpha=0.5,
#         color='steelblue',
#         edgecolor='black',
#         label='Initial Opinions'
#     )

#     n_final, bins, patches = plt.hist(
#         row['individual_opinions'], 
#         bins=bin_edges, 
#         density=False, 
#         alpha=0.7, 
#         color='green',
#         edgecolor='black',
#         label='Final Opinions'
#     )



#     # Add vertical lines at media positions for reference
#     if(stepNo == 5):
#         for medop in row['media_opinions']:
#                 plt.axvline(medop, color='red', linestyle='--', alpha=0.5, label='Media Opinion' if medop == row['media_opinions'][0] else "")

#     # Format plot
#     plt.title(title, fontsize=14, pad=20)
#     plt.xlabel('Opinion Value', fontsize=12)
#     plt.ylabel('Frequency', fontsize=12)
#     plt.xticks(bin_edges, rotation=45)
#     plt.xlim(-1.05, 1.05)
#     # plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.legend()

#     # Add bin labels
#     bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#     # for i in range(len(bin_intervals)):
#     #     plt.text(
#     #         bin_centers[i], max(n_final[i], n_initial[i]) + 0.01, 
#     #         f"F:{n_final[i]}\nI:{n_initial[i]}\n{bin_intervals[i]}", 
#     #         ha='center', va='bottom', fontsize=8
#     #     )

#     # Add statistics box
#     stats_text = (
#         f"N: {len(row['individual_opinions'])}\n"
#         f"Initial Mean: {np.mean(row['individual_initial_opinions']):.3f}\n"
#         f"Final Mean: {np.mean(row['individual_opinions']):.3f}\n"
#         f"Change: {np.mean(row['individual_opinions']) - np.mean(row['individual_initial_opinions']):.3f}"
#     )

#     plt.gcf().text(0.50, 0.85, stats_text, fontsize=10, 
#                    bbox=dict(facecolor='white', alpha=0.5))

#     # Save and close
#     plt.tight_layout()
#     plt.savefig(os.path.join(histograms_save_path, filename), dpi=150)
#     plt.close()
#     return filename

# # Function to generate histogram by cumulative index
# def generate_histogram_by_cumulative_index(cumulative_index, histograms_save_path):
#     """Generate histogram for a simulation based on its cumulative index"""
#     for stepNo in steps_to_process:
#         simulations_df = simulation_df_step[stepNo - 1]
#         if cumulative_index in simulations_df['cumulative_simulation_index'].values:
#             row = simulations_df[simulations_df['cumulative_simulation_index'] == cumulative_index].iloc[0]
#             model_index = row['model_simulation_index']
#             return save_simulation_histogram(row, cumulative_index, stepNo, histograms_save_path)
#     print(f"Error: Cumulative index {cumulative_index} not found.")
#     return None

# # Function to generate histogram by step and model index
# def generate_histogram_by_step_model_index(stepNo, model_index, histograms_save_path):
#     """Generate histogram for a simulation based on step number and model index"""
#     if stepNo > len(simulation_df_step):
#         print(f"Error: Step number {stepNo} is out of range.")
#         return None

#     simulations_df = simulation_df_step[stepNo - 1]
#     if model_index not in simulations_df['model_simulation_index'].values:
#         print(f"Error: Model index {model_index} not found in step {stepNo}.")
#         return None

#     row = simulations_df[simulations_df['model_simulation_index'] == model_index].iloc[0]
#     cumulative_index = row['cumulative_simulation_index']
#     return save_simulation_histogram(row, cumulative_index, stepNo, histograms_save_path)

# # Setup directories for each step
# for stepNo in steps_to_process:
#     # Setup Path
#     histograms_save_path = os.path.join(plots_folder_path, "selected_sims_histograms")
#     if not os.path.exists(histograms_save_path):
#         os.makedirs(histograms_save_path)

#     # You can now call either:
#     # generate_histogram_by_cumulative_index(1, histograms_save_path)
#     # or
#     generate_histogram_by_step_model_index(4, 1, histograms_save_path)


# # Load Survey Data.

# In[17]:
print("-------------------- LOADING SURVEY DATA --------------------------------")


def parse_year(year_val):
    """Convert year string like '2022S', '2022A', '2022W' to numeric year."""
    if isinstance(year_val, str):
        if year_val.endswith("S"):
            return int(year_val[:-1]) - 0.5
        elif year_val.endswith(("A", "W")):
            return int(year_val[:-1])
        else:
            return int(year_val)  # just in case it's string but plain year
    return int(year_val)  # if it's already an int


# Set up paths
folder_path = "EU_dataset"
file_title = "Cleaned_Parlemeter_Data-LeftRight"
file_name = file_title + ".xlsx"
survey_file = os.path.join(folder_path, file_name)

# sheets = pd.ExcelFile(survey_file).sheet_names[1:]  # Skip 'RowHeaders'

## As xlrd doesnt support xlsx
sheets = pd.ExcelFile(survey_file, engine='openpyxl').sheet_names[1:]  # Skip 'RowHeaders'


survey_data = []

for sheet in sheets:
    country = sheet.split("-")[0]  # Extract country code (e.g., "BE")
    df_sheet = pd.read_excel(survey_file, sheet_name=sheet, header=None, engine='openpyxl')

    # Years are in the first row
    years = df_sheet.iloc[0, :].tolist()

    # Extract counts for positions 1-10 (skip % rows and DK/Refusal)
    counts_matrix = []
    for i in range(10):  # Positions 1-10
        row_index = 1 + 2 * i  # Count rows are at indices 1,3,5,...,19
        counts_matrix.append(df_sheet.iloc[row_index, :].tolist())

    # Process each year
    for j, year in enumerate(years):
        year_num = parse_year(year)  # use parsed numeric year
        counts = [float(counts_matrix[i][j]) for i in range(10)]
        total_valid = sum(counts)
        distribution = [c / total_valid for c in counts]
        pooled_counts = counts[:4] + [counts[4] + counts[5]] + counts[6:]
        pooled_distribution = [c / total_valid for c in pooled_counts]

        survey_data.append({
            "country": country,
            "year": year_num,
            "binned_distribution": distribution,
            "binned_distribution_5_6clubbed": pooled_distribution,
            "total_respondents": total_valid
        })
survey_df = pd.DataFrame(survey_data)
print("-------------------- SURVEY DATA LOADED --------------------------------")


# In[18]:


survey_df


# ## Fitting Data for all Steps

# ### Old Code (fits and saves data)

# In[19]:

print("-------------------- PREPARING TO FIT --Fitting Data for all Steps --------------------------------")
# Ensure vectors are numpy arrays of floats
def to_array(x):
    return np.array(x, dtype=float)

# List of dataframes for the fits (now combined initial and final)
jsfit_dfs_stepwise = []
best_jsfits_stepwise = []

# Define all possible parameters including silence variables
base_params = [
    'cumulative_simulation_index', 'epsM', 'epsSD', 'OpM', 'OpSD', 
    'MedM', 'RandomSeed', 'MedSD', 'MedInfF'
]

silence_vars = [
    "Silence_Alpha",
    "Silence_Tau", 
    "Silence_Delta0",
    "SilenceByBoundary"
]

all_possible_params = base_params + silence_vars

for stepNo in steps_to_process:
    simulations_df = simulation_df_step[stepNo - 1]
    this_fit_data = []

    # Identify which parameters actually exist in this step's data
    existing_params = []
    missing_params = []

    for param in all_possible_params:
        if param in simulations_df.columns:
            existing_params.append(param)
        else:
            missing_params.append(param)

    print(f"\nStep {stepNo}: Parameter Analysis")
    print(f"  Existing parameters: {len(existing_params)}")
    print(f"  Missing parameters: {len(missing_params)}")

    if missing_params:
        print(f"  Missing: {missing_params}")

    # Check specifically for silence variables
    existing_silence = [var for var in silence_vars if var in existing_params]
    print(f"  Silence variables found: {existing_silence}")

    # Precompute survey distributions once
    survey_dists = []
    for _, surv_row in survey_df.iterrows():
        survey_dists.append({
            'country': surv_row['country'],
            'year': surv_row['year'],
            'distribution': to_array(surv_row['binned_distribution_5_6clubbed'])
        })

    for sim_idx, sim_row in simulations_df.iterrows():
        if(sim_idx%10000 == 0): 
            print(f"  Step {stepNo} has reached index {sim_idx}")
        sim_dist_final = to_array(sim_row['binned_distributions_5_6clubbed'])
        sim_dist_initial = to_array(sim_row['binned_initial_distributions_5_6clubbed'])

        # Extract simulation parameters - only those that exist
        sim_params = {}
        for param in existing_params:
            # Handle column name variations
            if param == 'cumulative_simulation_index' and param not in sim_row:
                # Try alternative name
                if 'cum_sim_index' in sim_row:
                    sim_params['cum_sim_index'] = sim_row['cum_sim_index']
                else:
                    sim_params['cum_sim_index'] = sim_idx
            else:
                sim_params[param] = sim_row[param]

        for surv_data in survey_dists:
            surv_dist = surv_data['distribution']

            # Jensen–Shannon divergence for both initial and final
            jsd_final = jensenshannon(sim_dist_final, surv_dist)
            jsd_initial = jensenshannon(sim_dist_initial, surv_dist)

            # Create a single entry with both distances
            fit_entry = {
                **sim_params,  # Unpack simulation parameters
                'country': surv_data['country'],
                'year': surv_data['year'],
                'distance_final': jsd_final,
                'distance_initial': jsd_initial,
                'distance_delta': jsd_final - jsd_initial
            }

            this_fit_data.append(fit_entry)

    # Create a single DataFrame for this step
    thisfit_df = pd.DataFrame(this_fit_data)
    jsfit_dfs_stepwise.append(thisfit_df)

    # Show sample of what was included
    print(f"Step {stepNo} final DataFrame shape: {thisfit_df.shape}")
    print(f"Step {stepNo} columns: {list(thisfit_df.columns)}")

    if existing_silence:
        sample_silence = thisfit_df[existing_silence].iloc[0].to_dict()
        print(f"Step {stepNo} sample silence values: {sample_silence}")

    # Extract best fits for both initial and final
    best_fit_final = thisfit_df.sort_values('distance_final').groupby(['country', 'year']).head(5)
    best_fit_initial = thisfit_df.sort_values('distance_initial').groupby(['country', 'year']).head(5)

    # Combine best fits into a single DataFrame with a type indicator
    best_fit_final['distance_type'] = 'final'
    best_fit_initial['distance_type'] = 'initial'

    best_fit_combined = pd.concat([best_fit_final, best_fit_initial], ignore_index=True)
    best_jsfits_stepwise.append(best_fit_combined)

print("\n" + "="*50)
print("PROCESSING SUMMARY")
print("="*50)
for i, df in enumerate(jsfit_dfs_stepwise):
    step_no = i + 1
    silence_included = [var for var in silence_vars if var in df.columns]
    print(f"Step {step_no}: {len(df):,} rows | Silence vars: {len(silence_included)}")
    if silence_included:
        unique_vals = {var: df[var].unique() for var in silence_included}
        print(f"  Unique values: {unique_vals}")


# ## Alternative processing pipeline for memory efficiency
# 
# ### Block 1

# In[ ]:


# import numpy as np
# import pandas as pd
# from scipy.spatial.distance import jensenshannon

# # Ensure vectors are numpy arrays of floats
# def to_array(x):
#     return np.array(x, dtype=float)

# # List of dataframes for the fits (now combined initial and final)
# jsfit_dfs_stepwise = []
# best_jsfits_stepwise = []

# # Define all possible parameters including silence variables
# base_params = [
#     'cumulative_simulation_index', 'epsM', 'epsSD', 'OpM', 'OpSD', 
#     'MedM', 'RandomSeed', 'MedSD', 'MedInfF'
# ]

# silence_vars = [
#     "Silence_Alpha",
#     "Silence_Tau", 
#     "Silence_Delta0",
#     "SilenceByBoundary"
# ]

# all_possible_params = base_params + silence_vars

# # Precompute survey distributions once for all steps
# print("Precomputing survey distributions...")
# survey_dists = []
# for _, surv_row in survey_df.iterrows():
#     survey_dists.append({
#         'country': surv_row['country'],
#         'year': surv_row['year'],
#         'distribution': to_array(surv_row['binned_distribution_5_6clubbed'])
#     })
# num_surveys = len(survey_dists)
# print(f"Precomputed {num_surveys} survey distributions")

# for stepNo in steps_to_process:
#     print(f"\n{'='*60}")
#     print(f"PROCESSING STEP {stepNo}")
#     print('='*60)

#     simulations_df = simulation_df_step[stepNo - 1]

#     # Identify which parameters actually exist in this step's data
#     existing_params = []
#     missing_params = []

#     for param in all_possible_params:
#         if param in simulations_df.columns:
#             existing_params.append(param)
#         else:
#             missing_params.append(param)

#     print(f"Parameter Analysis:")
#     print(f"  Existing parameters: {len(existing_params)}")
#     print(f"  Missing parameters: {len(missing_params)}")

#     if missing_params:
#         print(f"  Missing: {missing_params}")

#     # Check specifically for silence variables
#     existing_silence = [var for var in silence_vars if var in existing_params]
#     print(f"  Silence variables found: {existing_silence}")

#     # Process data in chunks to avoid memory issues
#     CHUNK_SIZE = 250  # Adjust based on your memory constraints
#     all_chunks_data = []
#     total_simulations = len(simulations_df)

#     print(f"\nProcessing {total_simulations:,} simulations in chunks of {CHUNK_SIZE}...")
#     print(f"Each simulation will be matched with {num_surveys} surveys")
#     print(f"Expected total rows: {total_simulations * num_surveys:,}")

#     # Process simulations in chunks
#     for chunk_start in range(0, total_simulations, CHUNK_SIZE):
#         chunk_end = min(chunk_start + CHUNK_SIZE, total_simulations)
#         chunk_df = simulations_df.iloc[chunk_start:chunk_end]

#         chunk_data = []

#         for sim_idx, sim_row in chunk_df.iterrows():
#             # Progress reporting
#             if (sim_idx - chunk_start) % 1000 == 0 and (sim_idx - chunk_start) > 0:
#                 print(f"  Step {stepNo}: Processed {sim_idx:,} simulations")

#             sim_dist_final = to_array(sim_row['binned_distributions_5_6clubbed'])
#             sim_dist_initial = to_array(sim_row['binned_initial_distributions_5_6clubbed'])

#             # Extract simulation parameters - only those that exist
#             sim_params = {}
#             for param in existing_params:
#                 # Handle column name variations
#                 if param == 'cumulative_simulation_index' and param not in sim_row:
#                     # Try alternative name
#                     if 'cum_sim_index' in sim_row:
#                         sim_params['cum_sim_index'] = sim_row['cum_sim_index']
#                     else:
#                         sim_params['cum_sim_index'] = sim_idx
#                 else:
#                     sim_params[param] = sim_row[param]

#             # Process all surveys for this simulation
#             for surv_data in survey_dists:
#                 surv_dist = surv_data['distribution']

#                 # Jensen–Shannon divergence for both initial and final
#                 jsd_final = jensenshannon(sim_dist_final, surv_dist)
#                 jsd_initial = jensenshannon(sim_dist_initial, surv_dist)

#                 # Create entry without building massive dictionary list
#                 entry = {
#                     **sim_params,  # Unpack simulation parameters
#                     'country': surv_data['country'],
#                     'year': surv_data['year'],
#                     'distance_final': jsd_final,
#                     'distance_initial': jsd_initial,
#                     'distance_delta': jsd_final - jsd_initial
#                 }
#                 chunk_data.append(entry)

#         # Store chunk data
#         all_chunks_data.append(chunk_data)
#         print(f"  Completed chunk {chunk_start:,}-{chunk_end:,} ({len(chunk_data):,} rows)")

#     # Store chunked data for this step (we'll create DataFrame in Block 2)
#     jsfit_dfs_stepwise.append({
#         'stepNo': stepNo,
#         'chunks': all_chunks_data,
#         'existing_params': existing_params,
#         'existing_silence': existing_silence,
#         'num_simulations': total_simulations,
#         'num_surveys': num_surveys
#     })

#     print(f"\nStep {stepNo} completed processing.")
#     print(f"Total chunks: {len(all_chunks_data)}")


# In[21]:


# print(all_chunks_data)


# In[22]:


# import gc
# import pandas as pd
# import os
# import json

# print("\n" + "="*60)
# print("DIRECT PARQUET WRITING WITH PROGRESS SAVING")
# print("="*60)

# # Create output directory
# output_dir = "step_data"
# os.makedirs(output_dir, exist_ok=True)

# # Try to load existing progress
# progress_file = f"{output_dir}/progress.json"
# if os.path.exists(progress_file):
#     with open(progress_file, 'r') as f:
#         step_data_info = json.load(f)
#     print(f"Loaded existing progress for {len(step_data_info)} steps")
# else:
#     step_data_info = []

# # Process only steps that haven't been saved yet
# processed_steps = {info['stepNo'] for info in step_data_info}

# for step_idx, step_data in enumerate(jsfit_dfs_stepwise):
#     stepNo = step_data['stepNo']

#     # Skip if already processed
#     if stepNo in processed_steps:
#         print(f"\nStep {stepNo}: Already processed, skipping...")
#         continue

#     all_chunks_data = step_data['chunks']
#     existing_silence = step_data['existing_silence']

#     print(f"\nStep {stepNo}: Processing {len(all_chunks_data)} chunks...")

#     # Process chunks in batches to avoid memory issues
#     batch_size = 10  # Process 10 chunks at a time
#     all_batches = []

#     for batch_start in range(0, len(all_chunks_data), batch_size):
#         batch_end = min(batch_start + batch_size, len(all_chunks_data))
#         batch_chunks = all_chunks_data[batch_start:batch_end]

#         print(f"  Processing batch {batch_start//batch_size + 1}/{(len(all_chunks_data) + batch_size - 1)//batch_size}...")

#         # Combine chunks in this batch
#         batch_rows = []
#         for chunk in batch_chunks:
#             batch_rows.extend(chunk)

#         # Convert to DataFrame
#         batch_df = pd.DataFrame(batch_rows)
#         all_batches.append(batch_df)

#         # Clear memory
#         del batch_rows, batch_df
#         gc.collect()

#     # Combine all batches
#     print(f"  Combining {len(all_batches)} batches...")
#     if len(all_batches) == 1:
#         thisfit_df = all_batches[0]
#     else:
#         thisfit_df = pd.concat(all_batches, ignore_index=True)

#     print(f"  Final DataFrame: {thisfit_df.shape}")

#     # Save to parquet
#     parquet_path = f"{output_dir}/step_{stepNo}.parquet"
#     print(f"  Saving to {parquet_path}...")
#     thisfit_df.to_parquet(parquet_path, index=False, compression='snappy')

#     file_size_mb = os.path.getsize(parquet_path) / (1024 * 1024)

#     print(f"  Saved! Size: {file_size_mb:.2f} MB")

#     # Store metadata
#     step_data_info.append({
#         'stepNo': stepNo,
#         'parquet_path': parquet_path,
#         'total_rows': len(thisfit_df),
#         'existing_silence': existing_silence,
#         'file_size_mb': file_size_mb
#     })

#     # Save progress after each step
#     with open(progress_file, 'w') as f:
#         json.dump(step_data_info, f, indent=2)

#     print(f"  Progress saved")

#     # Clear memory
#     del all_batches, thisfit_df
#     jsfit_dfs_stepwise[step_idx]['chunks'] = None
#     gc.collect()

# print(f"\n{'='*60}")
# print("ALL STEPS PROCESSED")
# print('='*60)
# print(f"Data saved to: {output_dir}/")
# print(f"Progress saved to: {progress_file}")


# In[ ]:


# import pandas as pd
# import pyarrow.dataset as ds
# import pyarrow.compute as pc
# import gc

# print("\n" + "="*60)
# print("PROCESSING BEST FITS FROM PARQUET")
# print("="*60)

# new_jsfit_dfs_stepwise = []
# new_best_jsfits_stepwise = []

# # Get unique country-year pairs
# unique_country_years = list(survey_df[['country', 'year']].drop_duplicates().itertools(index=False, name=None))
# print(f"Found {len(unique_country_years)} unique country-year pairs")

# for info in step_data_info:
#     stepNo = info['stepNo']
#     parquet_path = info['parquet_path']
#     total_rows = info['total_rows']

#     print(f"\nStep {stepNo}: Processing {total_rows:,} rows from parquet...")

#     # OPTION 1: Use pyarrow dataset for efficient filtering (most memory efficient)
#     print("  Using pyarrow dataset for efficient processing...")

#     dataset = ds.dataset(parquet_path, format='parquet')

#     # Process each country-year separately to avoid memory issues
#     best_fits_rows = []

#     for i, (country, year) in enumerate(unique_country_years):
#         if i % 10 == 0:
#             print(f"    Processing country-year {i+1}/{len(unique_country_years)}: {country}, {year}")

#         # Filter for this specific country-year using pyarrow
#         filter_expr = (ds.field('country') == country) & (ds.field('year') == year)

#         # Convert filtered data to pandas
#         filtered_table = dataset.to_table(filter=filter_expr)

#         if filtered_table.num_rows > 0:
#             filtered_df = filtered_table.to_pandas()

#             # Get top 5 for final distance
#             top_final = filtered_df.nsmallest(5, 'distance_final')
#             top_final['distance_type'] = 'final'

#             # Get top 5 for initial distance
#             top_initial = filtered_df.nsmallest(5, 'distance_initial')
#             top_initial['distance_type'] = 'initial'

#             best_fits_rows.append(top_final)
#             best_fits_rows.append(top_initial)

#         # Clear memory
#         del filtered_table, filtered_df
#         gc.collect()

#     # Combine all best fits
#     if best_fits_rows:
#         best_fit_combined = pd.concat(best_fits_rows, ignore_index=True)
#     else:
#         best_fit_combined = pd.DataFrame()

#     new_best_jsfits_stepwise.append(best_fit_combined)
#     print(f"  Best fits extracted: {len(best_fit_combined):,} rows")

#     # OPTION 2: Load full DataFrame for this step (only if you REALLY need it)
#     # Use memory mapping with chunks
#     print(f"  Loading full DataFrame with memory mapping...")

#     # Read in chunks to avoid memory issues
#     chunk_size = 100000  # Adjust based on available memory
#     df_chunks = []

#     for chunk in pd.read_parquet(parquet_path, chunksize=chunk_size):
#         df_chunks.append(chunk)
#         print(f"    Loaded chunk {len(df_chunks)}: {len(chunk):,} rows")

#         # If we have too many chunks, start combining incrementally
#         if len(df_chunks) >= 10:
#             combined = pd.concat(df_chunks, ignore_index=True)
#             df_chunks = [combined]  # Replace with combined chunk
#             print(f"    Combined to single chunk: {len(combined):,} rows")
#             gc.collect()

#     # Final combination
#     if df_chunks:
#         thisfit_df = pd.concat(df_chunks, ignore_index=True)
#     else:
#         thisfit_df = pd.DataFrame()

#     print(f"  Final DataFrame shape: {thisfit_df.shape}")
#     new_jsfit_dfs_stepwise.append(thisfit_df)

#     # Clear memory
#     del df_chunks, thisfit_df
#     gc.collect()

# # Replace the original lists
# jsfit_dfs_stepwise = new_jsfit_dfs_stepwise
# best_jsfits_stepwise = new_best_jsfits_stepwise

# print("\n" + "="*60)
# print("PROCESSING COMPLETE")
# print("="*60)
# for i, df in enumerate(jsfit_dfs_stepwise):
#     step_no = i + 1
#     print(f"Step {step_no}: {len(df):,} rows | Best fits: {len(best_jsfits_stepwise[i]):,} rows")


# In[ ]:


# type(jsfit_dfs_stepwise[0])


# In[ ]:


# import csv
# import os

# def save_large_list_to_csv(data_list, filename):
#     """Save large list of dictionaries directly to CSV without creating DataFrame"""
#     if not data_list:
#         print("No data to save")
#         return

#     # Get column names from first dictionary
#     fieldnames = list(data_list[0].keys())

#     with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
#         writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#         writer.writeheader()

#         # Write in chunks to manage memory
#         chunk_size = 10000
#         for i in range(0, len(data_list), chunk_size):
#             chunk = data_list[i:i + chunk_size]
#             writer.writerows(chunk)
#             print(f"Saved chunk {i//chunk_size + 1}/{(len(data_list)-1)//chunk_size + 1}")

#     print(f"Successfully saved {len(data_list)} rows to {filename}")

# # Usage - save your this_fit_data directly
# save_large_list_to_csv(this_fit_data, "step6_fit_data.csv")

# # Then append to your list (just the filename for reference)
# jsfit_dfs_stepwise.append("step6_fit_data.csv")  # Store filename instead of DataFrame


# ## Adding Silence param to step 6 data

# In[ ]:


# import pandas as pd

# # pick the relevant dfs
# sim_df = simulation_df_step[5]
# jsfit_df = jsfit_dfs_stepwise[5]

# # select only the columns we need from sim_df
# cols_to_copy = [
#     "cumulative_simulation_index",
#     "Silence_Alpha",
#     "Silence_Tau",
#     "Silence_Delta0",
#     "SilenceByBoundary",
# ]

# sim_subset = sim_df[cols_to_copy]

# # merge: match cum_sim_index in jsfit_df to cumulative_simulation_index in sim_subset
# merged = jsfit_df.merge(
#     sim_subset,
#     how="left",  # keep all rows in jsfit_df
#     left_on="cum_sim_index",
#     right_on="cumulative_simulation_index",
# )

# # drop duplicate key column if you don’t want it twice
# merged = merged.drop(columns=["cumulative_simulation_index"])

# # replace back in list if you want
# jsfit_dfs_stepwise[5] = merged


# In[ ]:


# print(type(jsfit_dfs_stepwise[0]))


# ## Saving Fits

# In[ ]:


save_fits = True
if(save_fits):
    if not os.path.exists(full_preprocessed_data_folder_path):
        os.makedirs(full_preprocessed_data_folder_path)
    for stepNo in steps_to_process:
        file_name = f"JS_fits_Step{stepNo}.csv"
        file_path = os.path.join(full_preprocessed_data_folder_path, file_name)
        jsfit_dfs_stepwise[stepNo - 1].to_csv(file_path, index=False)


# In[ ]:





# ## Plotting JSD v/s Epsilon
# 
# ### Steps 1 and 2

# In[ ]:


best_fit_threshold = 0.05 # Set this to filter the best fits

# Set up plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("colorblind")

# Create output directory for plots
plot_output_dir = os.path.join(plots_folder_path, "Fit_Plots")
os.makedirs(plot_output_dir, exist_ok=True)

# Process steps 1 and 2 (index 0 and 1 in our list)
for stepNo in steps_to_process:
    step_idx = stepNo - 1  # Convert to 1-based indexing
    fit_df = jsfit_dfs_stepwise[step_idx]

    print(f"Processing Step {stepNo} with {len(fit_df)} data points")

    # Plot 1: Final distance against epsM
    plt.figure(figsize=(10, 6))

    # Group by epsM and calculate mean and standard error
    grouped = fit_df.groupby('epsM')['distance_final'].agg(['mean', 'std', 'count'])
    grouped['se'] = grouped['std'] / np.sqrt(grouped['count'])


    grouped_initial = fit_df.groupby('epsM')['distance_initial'].agg(['mean', 'std', 'count'])
    grouped_initial['se']  =  grouped_initial['std'] / np.sqrt(grouped_initial['count'])


    plt.errorbar(grouped.index, grouped['mean'], yerr=grouped['se'], 
                 fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8, label="End of Sim divergence")

    #Initial for comparison
    plt.errorbar(grouped_initial.index, grouped_initial['mean'], yerr=grouped_initial['se'], 
                 fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8, label="Beginning of Sim divergence")

    plt.xlabel('ε (Openness)', fontsize=14)
    plt.ylabel('End of Sim Jensen-Shannon Divergence', fontsize=14)
    plt.title(f'Final Fit Quality vs Openness (ε)\n{step_titles[stepNo - 1]}', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(os.path.join(plot_output_dir, f'step{stepNo}_final_vs_epsM.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Final distance against epsM (only improved fits)
    plt.figure(figsize=(10, 6))

    # Filter for improved fits (final distance < initial distance)
    improved_fits = fit_df[fit_df['distance_delta'] < 0]

    if len(improved_fits) > 0:
        # Group by epsM and calculate statistics
        grouped_improved = improved_fits.groupby('epsM')['distance_final'].agg(['mean', 'std', 'count'])
        grouped_improved['se'] = grouped_improved['std'] / np.sqrt(grouped_improved['count'])

        plt.errorbar(grouped_improved.index, grouped_improved['mean'], yerr=grouped_improved['se'], 
                     fmt='s-', capsize=5, capthick=2, linewidth=2, markersize=8, color='green', label = "Improved fits only")

        # Add original for comparison (light gray)
        plt.errorbar(grouped.index, grouped['mean'], yerr=grouped['se'], 
                     fmt='o--', capsize=3, capthick=1, linewidth=1, markersize=5, 
                     color='gray', alpha=0.5, label='All simulations')

        plt.xlabel('ε (Openness)', fontsize=14)
        plt.ylabel('End of Sim Jensen-Shannon Divergence', fontsize=14)
        plt.title(f'Improved Fits Only (Final < Initial)\n{step_titles[stepNo - 1]}', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plot_output_dir, f'step{stepNo}_improved_final_vs_epsM.png'), dpi=300, bbox_inches='tight')
        plt.close()
    else:
        print(f"No improved fits found for Step {stepNo}")
        plt.close()

    # Plot 3: Delta distance against epsM
    plt.figure(figsize=(10, 6))

    # Group by epsM and calculate mean and standard error for delta
    grouped_delta = fit_df.groupby('epsM')['distance_delta'].agg(['mean', 'std', 'count'])
    grouped_delta['se'] = grouped_delta['std'] / np.sqrt(grouped_delta['count'])

    # Create the plot
    plt.errorbar(grouped_delta.index, grouped_delta['mean'], yerr=grouped_delta['se'], 
                 fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8, color='purple')

    # Add horizontal line at zero for reference
    plt.axhline(y=0, color='red', linestyle='--', alpha=0.7, linewidth=1)

    plt.xlabel('ε (Openness)', fontsize=14)
    plt.ylabel('Δ Divergence (Final - Initial)', fontsize=14)
    plt.title(f'Fit Improvement vs ε (Negative values = Improvement)\n{step_titles[stepNo - 1]}', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_output_dir, f'step{stepNo}_delta_vs_epsM.png'), dpi=300, bbox_inches='tight')
    plt.close()


    # Plot 4: Only the best of the best fits.
    best_fits = fit_df[fit_df['distance_final'] < best_fit_threshold]
    print("Best fits dataframe: \n")
    print(best_fits)
    if len(best_fits) > 0:
        # Group by epsM and calculate statistics
        grouped_best = best_fits.groupby('epsM')['distance_final'].agg(['mean', 'std', 'count'])
        grouped_best['se'] = grouped_best['std'] / np.sqrt(grouped_best['count'])


        grouped_best_initial = best_fits.groupby('epsM')['distance_initial'].agg(['mean', 'std', 'count'])
        grouped_best_initial['se'] = grouped_best_initial['std'] / np.sqrt(grouped_best_initial['count'])

        plt.errorbar(grouped_best.index, grouped_best['mean'], yerr=grouped_best['se'], 
                     fmt='s-', capsize=5, capthick=2, linewidth=2, markersize=8, color='green', label = "End of Sim Divergence")

        # Initial for comparison
        plt.errorbar(grouped_best.index, grouped_best['mean'], yerr=grouped_best['se'], 
                 fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8, label="Beginning of Sim divergence")

        plt.xlabel('ε (Openness)', fontsize=14)
        plt.ylabel('End of Sim Jensen-Shannon Divergence', fontsize=14)
        plt.title(f'Best Fits Only (JS < {best_fit_threshold})\n{step_titles[stepNo - 1]}', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plot_output_dir, f'step{stepNo}_best_final_vs_epsM.png'), dpi=300, bbox_inches='tight')
        plt.close()
    else:
        print(f"No good fits found for Step {stepNo}")
        plt.close()



    # Additional analysis: Print summary statistics
    print(f"\nStep {stepNo} Summary Statistics:")
    print(f"Total simulations: {fit_df['cum_sim_index'].nunique()}")
    print(f"Total country-year combinations: {fit_df[['country', 'year']].drop_duplicates().shape[0]}")
    print(f"Mean final distance: {fit_df['distance_final'].mean():.4f}")
    print(f"Mean initial distance: {fit_df['distance_initial'].mean():.4f}")
    print(f"Mean delta: {fit_df['distance_delta'].mean():.4f}")
    print(f"Percentage improved: {(len(improved_fits) / len(fit_df) * 100):.1f}%")

    # Best epsM values
    best_epsM_final = grouped.sort_values('mean').index[0]
    best_epsM_delta = grouped_delta.sort_values('mean').index[0]
    print(f"Best epsM for final distance: {best_epsM_final:.3f}")
    print(f"Best epsM for improvement: {best_epsM_delta:.3f}")

print(f"\nAll plots saved to: {plot_output_dir}")


# In[ ]:





# ## Color Plots

# In[ ]:


from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Create output directory for plots
plot_output_dir = os.path.join(plots_folder_path, "Colorplots")
best_fit_threshold = 0.05

# Greek symbols for parameter names
param_labels = {
    'epsM': r'$\mu_{\varepsilon}$',
    'epsSD': r'$\sigma_{\varepsilon}$', 
    'OpSD': r'$\sigma_{InitialOp}$',
    'MedM': r'$\mu_{MediaOp}$',
    'MedSD': r'$\sigma_{MediaOp}$',
    'MedInfF': 'Media Influence'
}

# Global style settings
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 8,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 18
})

# Global flags
SHOW_GRID = True  # Set to True to enable grid

# Common color normalization for all plots across all steps
def get_global_color_norm(variant, all_steps_data):
    """Get common color normalization across all steps for a variant"""
    all_values = []

    for step_data in all_steps_data:
        if step_data is not None and len(step_data) > 0:
            if variant == 'mean_jsd':
                values = step_data.groupby(['epsM', 'epsSD'])['distance_final'].mean()
            elif variant == 'mean_jsd_improved':
                improved = step_data[step_data['distance_delta'] < 0]
                values = improved.groupby(['epsM', 'epsSD'])['distance_final'].mean()
            elif variant == 'mean_jsd_best_fit':
                best_fit = step_data[step_data['distance_final'] < best_fit_threshold]
                values = best_fit.groupby(['epsM', 'epsSD'])['distance_final'].mean()
            elif variant == 'best_jsd':
                values = step_data.groupby(['epsM', 'epsSD'])['distance_final'].min()

            if len(values) > 0:
                all_values.extend(values.values)

    if not all_values:
        return colors.Normalize(vmin=0, vmax=1)

    return colors.Normalize(vmin=min(all_values), vmax=max(all_values))

# Function to create categorical heatmap
def create_categorical_heatmap(pivot_table, xlabel, ylabel, title, norm, cmap='viridis_r', ax=None):
    """Create a categorical heatmap with proper discrete parameter handling"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Get the unique parameter values
    x_categories = pivot_table.columns.tolist()
    y_categories = pivot_table.index.tolist()

    # Create the heatmap
    im = ax.imshow(pivot_table.values, origin='lower', aspect='auto',
                   norm=norm, cmap=cmap)

    # Set ticks and labels at cell centers
    ax.set_xticks(np.arange(len(x_categories)))
    ax.set_yticks(np.arange(len(y_categories)))
    ax.set_xticklabels([f"{x:.2f}" for x in x_categories])
    ax.set_yticklabels([f"{y:.2f}" for y in y_categories])

    # Add proper grid lines if enabled
    if SHOW_GRID:
        # Add grid lines at cell boundaries
        ax.set_xticks(np.arange(len(x_categories) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(y_categories) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linestyle='-', linewidth=1.5, alpha=0.8)
        ax.tick_params(which="minor", length=0)

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)

    return im

# Function to process data for a variant
def process_variant_data(df, variant):
    """Process data for different variants"""
    if variant == 'mean_jsd':
        return df.groupby(['epsM', 'epsSD'])['distance_final'].mean().reset_index()
    elif variant == 'mean_jsd_improved':
        improved = df[df['distance_delta'] < 0]
        return improved.groupby(['epsM', 'epsSD'])['distance_final'].mean().reset_index()
    elif variant == 'mean_jsd_best_fit':
        best_fit = df[df['distance_final'] < best_fit_threshold]
        return best_fit.groupby(['epsM', 'epsSD'])['distance_final'].mean().reset_index()
    elif variant == 'best_jsd':
        return df.groupby(['epsM', 'epsSD'])['distance_final'].min().reset_index()

# Collect all data for global normalization
all_steps_data = jsfit_dfs_stepwise[2:5]  # Steps 3, 4, 5

# Precompute global norms for each variant
global_norms = {}
variants = ['mean_jsd', 'mean_jsd_improved', 'mean_jsd_best_fit', 'best_jsd']

for variant in variants:
    global_norms[variant] = get_global_color_norm(variant, all_steps_data)
    print(f"{variant}: vmin={global_norms[variant].vmin:.4f}, vmax={global_norms[variant].vmax:.4f}")

# Step 3: epsM and epsSD only
print("Processing Step 3...")
step3_dir = os.path.join(plot_output_dir, "step3")
os.makedirs(step3_dir, exist_ok=True)

step3_data = jsfit_dfs_stepwise[2] if len(jsfit_dfs_stepwise) > 2 else None

if step3_data is not None:
    # Get unique parameter values for the entire step
    unique_epsM = sorted(step3_data['epsM'].unique())
    unique_epsSD = sorted(step3_data['epsSD'].unique())

    for variant in variants:
        variant_data = process_variant_data(step3_data, variant)

        if len(variant_data) > 0:
            # Create pivot table and ensure all parameter combinations are represented
            pivot = variant_data.pivot_table(index='epsSD', columns='epsM', values='distance_final')
            pivot = pivot.reindex(index=unique_epsSD, columns=unique_epsM)

            fig, ax = plt.subplots(figsize=(12, 9))
            norm = global_norms[variant]

            im = create_categorical_heatmap(pivot, param_labels['epsM'], param_labels['epsSD'],
                                          f'Step 3: {variant.replace("_", " ").title()}', norm, ax=ax)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Jensen-Shannon Distance', fontsize=14)

            # Save best_jsd directly in step folder, others in subfolders
            if variant == 'best_jsd':
                plt.savefig(os.path.join(step3_dir, f'step3_{variant}.png'), 
                           dpi=300, bbox_inches='tight')
            else:
                variant_dir = os.path.join(step3_dir, variant)
                os.makedirs(variant_dir, exist_ok=True)
                plt.savefig(os.path.join(variant_dir, f'step3_{variant}.png'), 
                           dpi=300, bbox_inches='tight')

            plt.close()

# Step 4: epsM, epsSD, and OpSD
print("Processing Step 4...")
step4_dir = os.path.join(plot_output_dir, "step4")
os.makedirs(step4_dir, exist_ok=True)

step4_data = jsfit_dfs_stepwise[3] if len(jsfit_dfs_stepwise) > 3 else None

if step4_data is not None:
    # Get unique parameter values for the entire step
    unique_epsM = sorted(step4_data['epsM'].unique())
    unique_epsSD = sorted(step4_data['epsSD'].unique())
    unique_OpSD = sorted(step4_data['OpSD'].unique())

    for variant in variants:
        for OpSD_val in unique_OpSD:
            OpSD_data = step4_data[step4_data['OpSD'] == OpSD_val]
            variant_data = process_variant_data(OpSD_data, variant)

            if len(variant_data) > 0:
                # Create pivot table and ensure all parameter combinations are represented
                pivot = variant_data.pivot_table(index='epsSD', columns='epsM', values='distance_final')
                pivot = pivot.reindex(index=unique_epsSD, columns=unique_epsM)

                fig, ax = plt.subplots(figsize=(12, 9))
                norm = global_norms[variant]

                im = create_categorical_heatmap(pivot, param_labels['epsM'], param_labels['epsSD'],
                                              f'Step 4: {variant.replace("_", " ").title()}\n{param_labels["OpSD"]} = {OpSD_val:.2f}',
                                              norm, ax=ax)

                # Add colorbar
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Jensen-Shannon Distance', fontsize=14)

                # Save best_jsd directly in step folder, others in subfolders
                if variant == 'best_jsd':
                    plt.savefig(os.path.join(step4_dir, f'step4_OpSD_{OpSD_val:.2f}_{variant}.png'), 
                               dpi=300, bbox_inches='tight')
                else:
                    variant_dir = os.path.join(step4_dir, variant)
                    os.makedirs(variant_dir, exist_ok=True)
                    plt.savefig(os.path.join(variant_dir, f'step4_OpSD_{OpSD_val:.2f}_{variant}.png'), 
                               dpi=300, bbox_inches='tight')

                plt.close()

# Step 5: epsM, epsSD, and OpSD with media
print("Processing Step 5...")
step5_dir = os.path.join(plot_output_dir, "step5")
os.makedirs(step5_dir, exist_ok=True)

step5_data = jsfit_dfs_stepwise[4] if len(jsfit_dfs_stepwise) > 4 else None

if step5_data is not None:
    # Get unique parameter values for the entire step
    unique_epsM = sorted(step5_data['epsM'].unique())
    unique_epsSD = sorted(step5_data['epsSD'].unique())
    unique_MedM = sorted(step5_data['MedM'].unique())
    unique_MedSD = sorted(step5_data['MedSD'].unique())
    unique_MedInfF = sorted(step5_data['MedInfF'].unique())

    for variant in variants:
        for MedInfF_val in unique_MedInfF:
            MedInfF_data = step5_data[step5_data['MedInfF'] == MedInfF_val]

            # Create subplot grid
            fig, axes = plt.subplots(len(unique_MedSD), len(unique_MedM), 
                                    figsize=(6*len(unique_MedM), 5*len(unique_MedSD)),
                                    sharex=True, sharey=True)

            if len(unique_MedSD) == 1 and len(unique_MedM) == 1:
                axes = np.array([[axes]])
            elif len(unique_MedSD) == 1:
                axes = axes[np.newaxis, :]
            elif len(unique_MedM) == 1:
                axes = axes[:, np.newaxis]

            # Get global norm for this variant
            norm = global_norms[variant]

            # Find the first valid subplot to create colorbar
            first_valid_im = None

            for i, MedSD_val in enumerate(unique_MedSD):
                for j, MedM_val in enumerate(unique_MedM):
                    ax = axes[i, j]

                    # Filter data for this parameter combination
                    param_data = MedInfF_data[
                        (MedInfF_data['MedSD'] == MedSD_val) & 
                        (MedInfF_data['MedM'] == MedM_val)
                    ]

                    if len(param_data) > 0:
                        variant_data = process_variant_data(param_data, variant)
                        if len(variant_data) > 0:
                            # Create pivot table and ensure all parameter combinations are represented
                            pivot = variant_data.pivot_table(index='epsSD', columns='epsM', values='distance_final')
                            pivot = pivot.reindex(index=unique_epsSD, columns=unique_epsM)

                            im = create_categorical_heatmap(pivot, '', '', '', norm, ax=ax)
                            if first_valid_im is None:
                                first_valid_im = im

                            # Add parameter labels to subplots
                            if i == len(unique_MedSD) - 1:
                                ax.set_xlabel(param_labels['epsM'], fontsize=12)
                            if j == 0:
                                ax.set_ylabel(param_labels['epsSD'], fontsize=12)

                            # Add MedM and MedSD values to subplot title
                            ax.set_title(f'{param_labels["MedM"]}={MedM_val:.2f}\n{param_labels["MedSD"]}={MedSD_val:.2f}', 
                                       fontsize=11)
                    else:
                        ax.axis('off')

            # Add overall title
            fig.suptitle(f'Step 5: {variant.replace("_", " ").title()}\n'
                        f'{param_labels["MedInfF"]} = {MedInfF_val:.2f}', 
                        fontsize=18, y=0.95)

            # Add single colorbar for the entire figure
            if first_valid_im is not None:
                cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
                cbar = fig.colorbar(first_valid_im, cax=cbar_ax)
                cbar.set_label('Jensen-Shannon Distance', fontsize=14)

            plt.tight_layout()
            plt.subplots_adjust(top=0.9, bottom=0.1, right=0.9)

            # Save best_jsd directly in step folder, others in subfolders
            if variant == 'best_jsd':
                plt.savefig(os.path.join(step5_dir, f'step5_MedInfF_{MedInfF_val:.2f}_{variant}.png'), 
                           dpi=300, bbox_inches='tight')
            else:
                variant_dir = os.path.join(step5_dir, variant)
                os.makedirs(variant_dir, exist_ok=True)
                plt.savefig(os.path.join(variant_dir, f'step5_MedInfF_{MedInfF_val:.2f}_{variant}.png'), 
                           dpi=300, bbox_inches='tight')

            plt.close()

print("All color plots generated successfully!")
print(f"Plots saved to: {plot_output_dir}")


# In[ ]:


# Test block for Step 5 with categorical treatment of parameters
def create_categorical_heatmap(pivot_table, xlabel, ylabel, title, norm, cmap='viridis_r', ax=None):
    """Create a categorical heatmap with proper discrete parameter handling"""

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Get the unique parameter values
    x_categories = pivot_table.columns.tolist()
    y_categories = pivot_table.index.tolist()

    # Add grid lines between cells
    ax.set_xticks(np.arange(len(x_categories)+1)-0.5, minor=True)
    ax.set_yticks(np.arange(len(y_categories)+1)-0.5, minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", size=0)


    # Create the heatmap
    im = ax.imshow(pivot_table.values, origin='lower', aspect='auto',
                   norm=norm, cmap=cmap)

    # Set ticks and labels
    ax.set_xticks(np.arange(len(x_categories)))
    ax.set_yticks(np.arange(len(y_categories)))
    ax.set_xticklabels([f"{x:.2f}" for x in x_categories])
    ax.set_yticklabels([f"{y:.2f}" for y in y_categories])

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)

    return im

# Test with a specific variant and MedInfF value from Step 5
test_variant = 'best_jsd'
test_MedInfF = step5_data['MedInfF'].iloc[0]  # Use the first value for testing

# Get unique parameter values for the entire step
unique_epsM = sorted(step5_data['epsM'].unique())
unique_epsSD = sorted(step5_data['epsSD'].unique())
unique_MedM = sorted(step5_data['MedM'].unique())
unique_MedSD = sorted(step5_data['MedSD'].unique())

# Filter data for the test MedInfF value
test_data = step5_data[step5_data['MedInfF'] == test_MedInfF]

# Create subplot grid
fig, axes = plt.subplots(len(unique_MedSD), len(unique_MedM), 
                        figsize=(6*len(unique_MedM), 5*len(unique_MedSD)),
                        sharex=True, sharey=True)

if len(unique_MedSD) == 1 and len(unique_MedM) == 1:
    axes = np.array([[axes]])
elif len(unique_MedSD) == 1:
    axes = axes[np.newaxis, :]
elif len(unique_MedM) == 1:
    axes = axes[:, np.newaxis]

# Get global norm for this variant
norm = global_norms[test_variant]

# Find the first valid subplot to create colorbar
first_valid_im = None

for i, MedSD_val in enumerate(unique_MedSD):
    for j, MedM_val in enumerate(unique_MedM):
        ax = axes[i, j]

        # Filter data for this parameter combination
        param_data = test_data[
            (test_data['MedSD'] == MedSD_val) & 
            (test_data['MedM'] == MedM_val)
        ]

        if len(param_data) > 0:
            variant_data = process_variant_data(param_data, test_variant)
            if len(variant_data) > 0:
                # Create pivot table and ensure all parameter combinations are represented
                pivot = variant_data.pivot_table(index='epsSD', columns='epsM', values='distance_final')

                # Reindex to include all possible parameter values
                pivot = pivot.reindex(index=unique_epsSD, columns=unique_epsM)

                im = create_categorical_heatmap(pivot, '', '', '', norm, ax=ax)
                if first_valid_im is None:
                    first_valid_im = im

                # Add parameter labels to subplots
                if i == len(unique_MedSD) - 1:
                    ax.set_xlabel(param_labels['epsM'], fontsize=12)
                if j == 0:
                    ax.set_ylabel(param_labels['epsSD'], fontsize=12)

                # Add MedM and MedSD values to subplot title
                ax.set_title(f'{param_labels["MedM"]}={MedM_val:.2f}\n{param_labels["MedSD"]}={MedSD_val:.2f}', 
                           fontsize=11)
        else:
            ax.axis('off')

# Add overall title
fig.suptitle(f'Step 5: {test_variant.replace("_", " ").title()}\n'
            f'{param_labels["MedInfF"]} = {test_MedInfF:.2f}', 
            fontsize=18, y=0.95)

# Add single colorbar for the entire figure
if first_valid_im is not None:
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(first_valid_im, cax=cbar_ax)
    cbar.set_label('Jensen-Shannon Distance', fontsize=14)

plt.tight_layout()
plt.subplots_adjust(top=0.9, bottom=0.1, right=0.9)

# Save the test plot
test_dir = os.path.join(plot_output_dir, "test")
os.makedirs(test_dir, exist_ok=True)
plt.savefig(os.path.join(test_dir, f'step5_test_{test_variant}_MedInfF_{test_MedInfF:.2f}.png'), 
           dpi=300, bbox_inches='tight')
plt.close()

print("Test plot created successfully!")


# # 

# # Step 6 - Media with Silence

# In[ ]:


# jsfit_dfs_stepwise[5]


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib import colors

# Focus only on Step 6
step6_data = jsfit_dfs_stepwise[5] if len(jsfit_dfs_stepwise) > 5 else None

if step6_data is not None:
    # Set up output directory
    plot_output_dir = os.path.join(plots_folder_path, "Colorplots_Step6_Percentage")
    os.makedirs(plot_output_dir, exist_ok=True)

    # Get unique parameter values for Step 6
    unique_epsM = sorted(step6_data['epsM'].unique())
    unique_epsSD = sorted(step6_data['epsSD'].unique())
    unique_MedM = sorted(step6_data['MedM'].unique())
    unique_MedSD = sorted(step6_data['MedSD'].unique())
    unique_Silence_Alpha = sorted(step6_data['Silence_Alpha'].dropna().unique())
    unique_Silence_Tau = sorted(step6_data['Silence_Tau'].dropna().unique())
    unique_Silence_Delta0 = sorted(step6_data['Silence_Delta0'].dropna().unique())
    unique_SilenceByBoundary = sorted(step6_data['SilenceByBoundary'].dropna().unique())

    print(f"Step 6 Parameter Ranges:")
    print(f"  epsM: {unique_epsM}")
    print(f"  epsSD: {unique_epsSD}")
    print(f"  MedM: {unique_MedM}")
    print(f"  MedSD: {unique_MedSD}")
    print(f"  Silence_Alpha: {unique_Silence_Alpha}")
    print(f"  Silence_Tau: {unique_Silence_Tau}")
    print(f"  Silence_Delta0: {unique_Silence_Delta0}")
    print(f"  SilenceByBoundary: {unique_SilenceByBoundary}")

    # Parameter labels
    param_labels = {
        'epsM': r'$\mu_{\varepsilon}$',
        'epsSD': r'$\sigma_{\varepsilon}$', 
        'MedM': r'$\mu_{MediaOp}$',
        'MedSD': r'$\sigma_{MediaOp}$',
        'Silence_Alpha': r'$\alpha_{Silence}$',
        'Silence_Tau': r'$\tau_{Silence}$',
        'Silence_Delta0': r'$\delta_0_{Silence}$',
        'SilenceByBoundary': 'SilenceByBoundary'
    }

    def calculate_percentage_below_threshold(df):
        """Calculate percentage of simulations below threshold for each epsM/epsSD combination"""
        if len(df) == 0:
            return pd.DataFrame()

        result = df.groupby(['epsM', 'epsSD']).apply(
            lambda x: (x['distance_final'] <= best_fit_threshold).mean() * 100
        ).reset_index()
        result.columns = ['epsM', 'epsSD', 'pct_below_threshold']

        return result

    def create_categorical_heatmap(pivot_table, xlabel, ylabel, title, ax):
        """Create a categorical heatmap for percentage below threshold"""
        # Get the unique parameter values
        x_categories = pivot_table.columns.tolist()
        y_categories = pivot_table.index.tolist()

        # Create the heatmap with percentage colormap
        im = ax.imshow(pivot_table.values, origin='lower', aspect='auto',
                      cmap='RdYlBu', vmin=0, vmax=100)

        # Set ticks and labels at cell centers
        ax.set_xticks(np.arange(len(x_categories)))
        ax.set_yticks(np.arange(len(y_categories)))
        ax.set_xticklabels([f"{x:.2f}" for x in x_categories], fontsize=8)
        ax.set_yticklabels([f"{y:.2f}" for y in y_categories], fontsize=8)

        # Add value annotations
        for i in range(len(y_categories)):
            for j in range(len(x_categories)):
                value = pivot_table.iloc[i, j]
                if not np.isnan(value):
                    ax.text(j, i, f'{value:.0f}%', 
                           ha='center', va='center', 
                           fontsize=7, fontweight='bold',
                           color='white' if value > 50 else 'black')

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11)

        return im

    # Generate plots for all silence parameter combinations
    total_combinations = (len(unique_Silence_Alpha) * len(unique_Silence_Tau) * 
                         len(unique_Silence_Delta0) * len(unique_SilenceByBoundary))
    current_combination = 0

    for silence_alpha in unique_Silence_Alpha:
        for silence_tau in unique_Silence_Tau:
            for silence_delta0 in unique_Silence_Delta0:
                for silence_boundary in unique_SilenceByBoundary:
                    current_combination += 1
                    print(f"Processing combination {current_combination}/{total_combinations}: "
                          f"SA={silence_alpha}, ST={silence_tau}, SDO={silence_delta0}, SBB={silence_boundary}")

                    # Filter data for this silence parameter combination
                    silence_data = step6_data[
                        (step6_data['Silence_Alpha'] == silence_alpha) &
                        (step6_data['Silence_Tau'] == silence_tau) &
                        (step6_data['Silence_Delta0'] == silence_delta0) &
                        (step6_data['SilenceByBoundary'] == silence_boundary)
                    ]

                    if len(silence_data) == 0:
                        print(f"  No data for this combination, skipping...")
                        continue

                    # Create subplot grid for Media parameters
                    fig, axes = plt.subplots(len(unique_MedSD), len(unique_MedM), 
                                            figsize=(4*len(unique_MedM), 3.5*len(unique_MedSD)))

                    if len(unique_MedSD) == 1 and len(unique_MedM) == 1:
                        axes = np.array([[axes]])
                    elif len(unique_MedSD) == 1:
                        axes = axes[np.newaxis, :]
                    elif len(unique_MedM) == 1:
                        axes = axes[:, np.newaxis]

                    # Find the first valid subplot to create colorbar
                    first_valid_im = None

                    for i, med_sd in enumerate(unique_MedSD):
                        for j, med_m in enumerate(unique_MedM):
                            ax = axes[i, j]

                            # Filter data for this media parameter combination
                            media_data = silence_data[
                                (silence_data['MedSD'] == med_sd) & 
                                (silence_data['MedM'] == med_m)
                            ]

                            if len(media_data) > 0:
                                # Calculate percentage below threshold
                                percentage_data = calculate_percentage_below_threshold(media_data)

                                if len(percentage_data) > 0:
                                    # Create pivot table
                                    pivot = percentage_data.pivot_table(index='epsSD', columns='epsM', 
                                                                       values='pct_below_threshold')
                                    pivot = pivot.reindex(index=unique_epsSD, columns=unique_epsM)

                                    im = create_categorical_heatmap(pivot, '', '', '', ax)
                                    if first_valid_im is None:
                                        first_valid_im = im

                                    # Add Media parameter labels to subplots
                                    if i == len(unique_MedSD) - 1:
                                        ax.set_xlabel(param_labels['epsM'], fontsize=9)
                                    if j == 0:
                                        ax.set_ylabel(param_labels['epsSD'], fontsize=9)

                                    # Add MedM and MedSD values to subplot title
                                    ax.set_title(f'{param_labels["MedM"]}={med_m:.2f}\n{param_labels["MedSD"]}={med_sd:.2f}', 
                                               fontsize=9)
                            else:
                                ax.axis('off')
                                ax.text(0.5, 0.5, 'No Data', 
                                       ha='center', va='center', transform=ax.transAxes, fontsize=10)

                    # Create comprehensive title with all parameters
                    title_params = [
                        f"{param_labels['Silence_Alpha']}={silence_alpha:.2f}",
                        f"{param_labels['Silence_Tau']}={silence_tau:.2f}",
                        f"{param_labels['Silence_Delta0']}={silence_delta0:.2f}",
                        f"{param_labels['SilenceByBoundary']}={silence_boundary}"
                    ]

                    fig.suptitle(f'Step 6: Percentage of Simulations with JSD ≤ {best_fit_threshold}\n' +
                                ' | '.join(title_params), 
                                fontsize=14, y=0.98)

                    # Add single colorbar for the entire figure
                    if first_valid_im is not None:
                        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
                        cbar = fig.colorbar(first_valid_im, cax=cbar_ax)
                        cbar.set_label('Percentage Below Threshold (%)', fontsize=10)

                    # Create filename with all parameters
                    filename_parts = [
                        f"step6",
                        f"SA{silence_alpha:.2f}",
                        f"ST{silence_tau:.2f}",
                        f"SDO{silence_delta0:.2f}",
                        f"SBB{silence_boundary}"
                    ]
                    filename = "_".join(filename_parts) + ".png"

                    plt.tight_layout()
                    plt.subplots_adjust(top=0.92, bottom=0.1, right=0.9)
                    plt.savefig(os.path.join(plot_output_dir, filename), 
                               dpi=300, bbox_inches='tight')
                    plt.close()

                    print(f"  Saved: {filename}")

    print(f"\nStep 6 colorplots completed! Saved to: {plot_output_dir}")

    # Generate a summary statistics file
    summary_stats = []
    for silence_alpha in unique_Silence_Alpha:
        for silence_tau in unique_Silence_Tau:
            for silence_delta0 in unique_Silence_Delta0:
                for silence_boundary in unique_SilenceByBoundary:
                    silence_data = step6_data[
                        (step6_data['Silence_Alpha'] == silence_alpha) &
                        (step6_data['Silence_Tau'] == silence_tau) &
                        (step6_data['Silence_Delta0'] == silence_delta0) &
                        (step6_data['SilenceByBoundary'] == silence_boundary)
                    ]

                    if len(silence_data) > 0:
                        overall_percentage = (silence_data['distance_final'] <= best_fit_threshold).mean() * 100
                        summary_stats.append({
                            'Silence_Alpha': silence_alpha,
                            'Silence_Tau': silence_tau,
                            'Silence_Delta0': silence_delta0,
                            'SilenceByBoundary': silence_boundary,
                            'Overall_Percentage': overall_percentage,
                            'Total_Simulations': len(silence_data)
                        })

    if summary_stats:
        summary_df = pd.DataFrame(summary_stats)
        summary_df.to_csv(os.path.join(plot_output_dir, "step6_summary_statistics.csv"), index=False)
        print(f"Summary statistics saved to: step6_summary_statistics.csv")

else:
    print("Step 6 data not found in jsfit_dfs_stepwise")


# In[ ]:


step6_data


# # FIT SUMMARY ACROSS STEPS BROKEN DOWN BY PARAMERTERS

# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns

# Set up the output directory structure
plot_output_dir = os.path.join(plots_folder_path, "Percent_Fit_By_Param")
threshold_folder = f"threshold_{best_fit_threshold:.2f}".replace('.', 'dot')
threshold_dir = os.path.join(plot_output_dir, threshold_folder)
os.makedirs(threshold_dir, exist_ok=True)

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("tab10")

# Define parameters for each step
step_parameters = {
    1: ['epsM'],
    2: ['epsM'],
    3: ['epsM', 'epsSD'],
    4: ['epsM', 'epsSD', 'OpSD'],
    5: ['epsM', 'epsSD', 'MedM', 'MedSD', 'MedInfF'],
    6: ['epsM', 'epsSD', 'MedM', 'MedSD', 'Silence_Alpha', 'Silence_Tau', 'Silence_Delta0', 'SilenceByBoundary']
}

# Parameter labels for plotting
param_labels = {
    'epsM': r'$\mu_{\varepsilon}$',
    'epsSD': r'$\sigma_{\varepsilon}$', 
    'OpSD': r'$\sigma_{InitialOp}$',
    'MedM': r'$\mu_{MediaOp}$',
    'MedSD': r'$\sigma_{MediaOp}$',
    'MedInfF': 'Media Influence',
    'Silence_Alpha': r'$\alpha_{Silence}$',
    'Silence_Tau': r'$\tau_{Silence}$',
    'Silence_Delta0': r'$\delta0_{Silence}$',
    'SilenceByBoundary': 'Silence by Boundary'
}

def calculate_fit_percentage_by_param(step_data, param_name):
    """Calculate percentage of simulations that fit for each value of a parameter"""
    if param_name not in step_data.columns:
        return None

    # Group by parameter value and calculate fit percentage
    result = step_data.groupby(param_name).apply(
        lambda x: (x['distance_final'] <= best_fit_threshold).mean() * 100
    ).reset_index()
    result.columns = [param_name, 'fit_percentage']

    # Also calculate number of simulations per parameter value
    counts = step_data.groupby(param_name).size().reset_index(name='count')
    result = result.merge(counts, on=param_name)

    return result

def create_parameter_plot(step_data, param_name, step_no, ax):
    """Create a plot showing fit percentage for a specific parameter"""
    data = calculate_fit_percentage_by_param(step_data, param_name)

    if data is None or len(data) == 0:
        ax.text(0.5, 0.5, f'No data for {param_name}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title(f"{param_labels.get(param_name, param_name)} - No Data", fontsize=14)
        return

    # Sort by parameter value for better visualization
    data = data.sort_values(param_name)

    # Create the plot
    if data[param_name].dtype in [np.float64, np.int64]:
        # Numerical parameter - use line plot
        ax.plot(data[param_name], data['fit_percentage'], 
                marker='o', linewidth=2, markersize=6, markerfacecolor='white')
        ax.set_xlabel(param_labels.get(param_name, param_name), fontsize=12)
    else:
        # Categorical parameter - use bar plot
        x_pos = np.arange(len(data))
        bars = ax.bar(x_pos, data['fit_percentage'], alpha=0.7)
        ax.set_xlabel(param_labels.get(param_name, param_name), fontsize=12)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'{x:.2f}' if isinstance(x, (int, float)) else str(x) 
                           for x in data[param_name]], rotation=45)

        # Add value labels on bars
        for bar, percentage, count in zip(bars, data['fit_percentage'], data['count']):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{percentage:.1f}%\n(n={count})', ha='center', va='bottom', fontsize=8)

    ax.set_ylabel('Fit Percentage (%)', fontsize=12)
    # ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3, linestyle='--')

    # Add threshold line and better title
    # ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    ax.set_title(f"{param_labels.get(param_name, param_name)}\nStep {step_no}", fontsize=14)

    # Add legend for numerical plots
    if data[param_name].dtype in [np.float64, np.int64]:
        ax.legend()

# Process each step
for step_no in steps_to_process:
    step_idx = step_no + 1
    step_data = jsfit_dfs_stepwise[step_idx]

    if step_data is None or len(step_data) == 0:
        print(f"Skipping Step {step_no} - no data")
        continue

    print(f"Processing Step {step_no}...")

    # Create step directory
    step_dir = os.path.join(threshold_dir, f"step{step_no}")
    os.makedirs(step_dir, exist_ok=True)

    # Get parameters for this step
    if step_no in step_parameters:
        params = step_parameters[step_no]
    else:
        # Fallback: use all available parameters except distance columns and indices
        exclude_cols = ['distance_final', 'distance_initial', 'distance_delta', 
                       'sim_index', 'country', 'year', 'survey_index']
        params = [col for col in step_data.columns if col not in exclude_cols 
                 and step_data[col].notna().any()]

    print(f"  Parameters to analyze: {params}")

    # Calculate overall fit percentage for this step
    overall_fit_percentage = (step_data['distance_final'] <= best_fit_threshold).mean() * 100
    total_simulations = len(step_data)

    print(f"  Overall fit percentage: {overall_fit_percentage:.1f}%")
    print(f"  Total simulations: {total_simulations}")

    # Create individual plots for each parameter
    for param in params:
        if param not in step_data.columns:
            print(f"    Skipping {param} - not in data")
            continue

        # Create individual plot for this parameter
        fig, ax = plt.subplots(figsize=(10, 6))
        create_parameter_plot(step_data, param, step_no, ax)

        # Add overall statistics to the plot
        fig.text(0.02, 0.02, 
                f'Overall fit: {overall_fit_percentage:.1f}% | Total simulations: {total_simulations} | Threshold: {best_fit_threshold}',
                fontsize=10, style='italic')

        plt.tight_layout()
        plt.savefig(os.path.join(step_dir, f"step{step_no}_{param}_fit_percentage.png"), 
                   dpi=300, bbox_inches='tight')
        plt.close()

    # Create a summary figure with all parameters for this step
    if len(params) > 0:
        n_cols = min(3, len(params))
        n_rows = (len(params) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))

        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)

        # Flatten axes array for easy iteration
        axes_flat = axes.flatten()

        for i, param in enumerate(params):
            if i < len(axes_flat):
                create_parameter_plot(step_data, param, step_no, axes_flat[i])
            else:
                break

        # Hide unused subplots
        for i in range(len(params), len(axes_flat)):
            axes_flat[i].axis('off')

        # Add overall title
        fig.suptitle(f'Step {step_no}: Fit Percentage by Parameter\n(Threshold = {best_fit_threshold})', 
                    fontsize=16, y=0.98)

        # Add overall statistics
        fig.text(0.02, 0.02, 
                f'Overall fit percentage: {overall_fit_percentage:.1f}% | Total simulations: {total_simulations}',
                fontsize=10, style='italic')

        plt.tight_layout()
        plt.subplots_adjust(top=0.92, bottom=0.08)
        plt.savefig(os.path.join(step_dir, f"step{step_no}_summary_fit_percentage.png"), 
                   dpi=300, bbox_inches='tight')
        plt.close()

print(f"\nAnalysis complete! Plots saved to: {threshold_dir}")

# Create a comprehensive summary across all steps
print("\n=== COMPREHENSIVE SUMMARY ACROSS ALL STEPS ===")
summary_data = []

for step_idx in range(len(jsfit_dfs_stepwise)):
    step_no = step_idx + 1
    step_data = jsfit_dfs_stepwise[step_idx]

    if step_data is None or len(step_data) == 0:
        continue

    fit_percentage = (step_data['distance_final'] <= best_fit_threshold).mean() * 100
    total_sims = len(step_data)

    summary_data.append({
        'Step': step_no,
        'Fit_Percentage': fit_percentage,
        'Total_Simulations': total_sims,
        'Parameters': ', '.join(step_parameters.get(step_no, ['Unknown']))
    })

if summary_data:
    summary_df = pd.DataFrame(summary_data)
    print("\nStep Summary:")
    print(summary_df.to_string(index=False))

    # Create summary plot across steps
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Fit percentage by step
    ax1.bar(summary_df['Step'], summary_df['Fit_Percentage'], color='skyblue', alpha=0.7)
    ax1.set_xlabel('Step', fontsize=12)
    ax1.set_ylabel('Fit Percentage (%)', fontsize=12)
    ax1.set_title('Fit Percentage by Step', fontsize=14)
    ax1.grid(True, alpha=0.3)

    # Add value labels on bars
    for i, (step, percentage) in enumerate(zip(summary_df['Step'], summary_df['Fit_Percentage'])):
        ax1.text(step, percentage + 1, f'{percentage:.1f}%', 
                ha='center', va='bottom', fontsize=10)

    # Simulation count by step
    ax2.bar(summary_df['Step'], summary_df['Total_Simulations'], color='lightcoral', alpha=0.7)
    ax2.set_xlabel('Step', fontsize=12)
    ax2.set_ylabel('Number of Simulations', fontsize=12)
    ax2.set_title('Simulation Count by Step', fontsize=14)
    ax2.grid(True, alpha=0.3)

    # Add value labels on bars
    for i, (step, count) in enumerate(zip(summary_df['Step'], summary_df['Total_Simulations'])):
        ax2.text(step, count + max(summary_df['Total_Simulations']) * 0.01, 
                f'{count:,}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(threshold_dir, "cross_step_summary.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()


# In[ ]:


jsfit_dfs_stepwise[5]

