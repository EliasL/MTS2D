import sys
from pathlib import Path
from settings import settings
import pandas as pd
import matplotlib.pyplot as plt

# Add Management to sys.path (used to import files)
sys.path.append(str(Path(__file__).resolve().parent.parent / 'Management'))

# Now we can import from Management
from configGenerator import ConfigGenerator

# We want to plot several runs ontop of each other
# In order to find the folders, we use the name generator in ConfigGenerator

seeds = range(0,11)
configs = ConfigGenerator.generateOverSeeds(seeds, nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.00001, maxLoad=1)

outPath = "/media/elias/T7 Sheild/output/"
macroPath = f"/{settings['MACRODATAFILE']}.csv"
filePaths = [outPath + config.generate_name(False) + macroPath for config in configs]

# Initialize a list to store DataFrames for each run
dfs = []

# Set up the plot for all runs and the average
plt.figure(figsize=(6, 5))

# Loop through each file path to read CSV data
for (i,filePath) in enumerate(filePaths):
    # Read the CSV data into a DataFrame
    df = pd.read_csv(filePath, usecols=['Load', 'Avg. energy'], dtype={'Load': 'float64', 'Avg. energy': 'float64'})
    dfs.append(df)
    
    # Plot individual run
    plt.plot(df['Load'], df['Avg. energy'], label=f'seed: {i}', alpha=0.5) 

# Concatenate all DataFrames for average calculation
combined_df = pd.concat(dfs)

# Group by Load and calculate mean of Avg. energy for the average plot
average_df = combined_df.groupby('Load')['Avg. energy'].mean().reset_index()

# Plot combined average with a distinct color and label
plt.plot(average_df['Load'], average_df['Avg. energy'], label='Average', color='black', linewidth=2)

# Customize the plot
plt.xlabel('Load')
plt.ylabel('Average Energy')
plt.title('Energy vs. Load + Average')
plt.legend()
plt.grid(True)

# Save the plot
visuals_folder = str(Path(__file__).resolve().parent.parent / 'Visuals')
plt.savefig(visuals_folder+f"/RunCollections/seeds0_{len(seeds)}_{configs[0].generate_name(False)}.pdf")