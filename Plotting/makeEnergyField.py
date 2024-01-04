from matplotlib import pyplot as plt
import numpy as np
from icecream import ic
import os
import csv


from settings import settings

def makeEnergyField(path, csv_file):
    ic("Plotting energy field...")

    dataPath = path + settings["DATAFOLDERPATH"]
    filePath = dataPath + csv_file

    if not os.path.exists(filePath):
        print(f"No file found at: {filePath}")
        return

    # Reading data from CSV
    data = np.genfromtxt(filePath, delimiter=',')
    
    # Replace NaN and -NaN with infinity
    data[np.isnan(data)] = np.inf

    # Extracting x, y, and energy values
    x_vals = data[:, 0]
    y_vals = data[:, 1]
    energies = data[:, 2]

    # Assuming equal spacing and regular grid
    grid_size_x = len(np.unique(x_vals))
    grid_size_y = len(np.unique(y_vals))
    energy_grid = energies.reshape((grid_size_y, grid_size_x)).transpose()


    # Create the plot
    plt.figure(figsize=(8, 6))

    # Set the minimum and maximum values for the color bar
    min_energy = 0  # Replace with your desired minimum value
    max_energy = 30  # Replace with your desired maximum value

    plt.imshow(energy_grid, cmap='viridis', origin='lower', vmin=min_energy, vmax=max_energy)


    # Adjusting ticks
    plt.xticks(np.linspace(0, grid_size_x - 1, 5), np.linspace(x_vals.min(), x_vals.max(), 5).round(2))
    plt.yticks(np.linspace(0, grid_size_y - 1, 5), np.linspace(y_vals.min(), y_vals.max(), 5).round(2))

    plt.colorbar(label='Energy')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Energy Distribution Heatmap')

    output_pdf_path = path + "energy_field.pdf"
    plt.savefig(output_pdf_path, format='pdf')


if __name__ == "__main__":
    # Replace 'your_pvd_file.pvd' with the path to your .pvd file
    makeEnergyField('build/output/testing/','energy_grid.csv')