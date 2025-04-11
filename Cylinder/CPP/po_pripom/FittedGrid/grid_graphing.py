# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 19:57:49 2025

@author: vithr
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

case = 2
bound = "Custom"
N = 400
iterations = 30000
Re = 10

# Initialize a list to store rows
data = []
file_path = 'MeshgridMaxDens25.000000MinDens2.000000'+'.txt'
# Open and read the file
with open(file_path, 'r') as file:
    # # Read the first line to get the double value
    # first_line = file.readline().strip()
    # time = float(first_line)
    for line in file:
        # Clean the line to remove extra spaces and newlines
        line = line.strip()
        # Split the line by ']' and process each element in brackets
        elements = line.split('] [')
        row = []
        for element in elements:
            # Remove the surrounding '[' and ']' characters
            element = element.strip('[]')
            # Convert the string elements into float and store them as numpy array
            row.append(np.array([float(x) for x in element.split()]))
        # Append the row (which is a list of numpy arrays) to data
        data.append(row)

# Convert the list of rows into a 2D numpy array (of numpy arrays, each of 5 values)
matrix = np.empty((len(data), len(data[0])), dtype=object)

# Assign the numpy arrays into the final 2D array
for i in range(len(data)):
    for j in range(len(data[i])):
        matrix[i, j] = data[i][j]


# Create three 2D matrices to store the first, second, and third subelements
x_coords = np.empty_like(matrix, dtype=np.float64)
y_coords = np.empty_like(matrix, dtype=np.float64)

# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(matrix.shape[0]):  # Iterate over rows
    for j in range(matrix.shape[1]):  # Iterate over columns
        x_coords[i, j] = matrix[i, j][0]
        y_coords[i, j] = matrix[i, j][1]

# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
x_vals = np.unique(x_coords)
y_vals = np.unique(y_coords)

# Create a 2D meshgrid based on x_vals and y_vals
X, Y = np.meshgrid(x_vals, y_vals)

# Plot the grid
plt.figure(figsize=(8, 6))

# Plot the grid lines for both x and y directions
for y in Y:  # Iterate over the rows (y-coordinates)
    plt.plot(x_vals, y, color='black', linestyle='-', linewidth=0.3)  # Plot horizontal lines

for x in X.T:  # Iterate over the columns (x-coordinates), transpose X for this
    plt.plot(x, y_vals, color='black', linestyle='-', linewidth=0.3)  # Plot vertical lines

plt.title("Meshgrid of Coordinates")
plt.xlabel("X Coordinates")
plt.ylabel("Y Coordinates")
plt.grid(False)
# Set the x and y limits
plt.xlim(min(x_vals), max(x_vals))  # Set the x-axis limit based on x_vals
plt.ylim(min(y_vals), max(y_vals))  # Set the y-axis limit based on y_vals
plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal aspect ratio
plt.show()

