import numpy as np
import matplotlib.pyplot as plt

case = 3
N = 160
iterations = 120000
bound = "DNN"
Re = 2

# Initialize a list to store rows
data = []
file_path = 'flowdata_cylinder'+str(case)+'_nondimensional_CPP_N'+str(N)+'_iter'+str(iterations)+'_CFL05_beta1.000000.txt'
# Open and read the file
with open(file_path, 'r') as file:
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
past = np.empty((len(data), len(data[0])), dtype=object)

# Assign the numpy arrays into the final 2D array
for i in range(len(data)):
    for j in range(len(data[i])):
        past[i, j] = data[i][j]


# Create three 2D matrices to store the first, second, and third subelements
pressures = np.empty_like(past, dtype=np.float64)
u_speeds = np.empty_like(past, dtype=np.float64)
v_speeds = np.empty_like(past, dtype=np.float64)
x_coords = np.empty_like(past, dtype=np.float64)
y_coords = np.empty_like(past, dtype=np.float64)

# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(past.shape[0]):  # Iterate over rows
    for j in range(past.shape[1]):  # Iterate over columns
        pressures[i, j] = past[i, j][0]  # First subelement
        u_speeds[i, j] = past[i, j][1]  # Second subelement
        v_speeds[i, j] = past[i, j][2]  # Third subelement
        x_coords[i, j] = past[i, j][3]
        y_coords[i, j] = past[i, j][4]

# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
x_vals = np.unique(x_coords)
y_vals = np.unique(y_coords)

# Generate a 2D grid based on unique x and y coordinates
X, Y = np.meshgrid(x_vals, y_vals)

# Initialize an empty grid for u_speeds
u_speeds_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    u_speeds_grid[y_index, x_index] = u_speeds[i]
    
# Initialize an empty grid for u_speeds
v_speeds_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    v_speeds_grid[y_index, x_index] = v_speeds[i]

# Plot streamlines
plt.figure(figsize=(10, 8))
plt.streamplot(X, Y, u_speeds_grid, v_speeds_grid, color='black', linewidth=0.5, density=4, integration_direction='both')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Proudnice při ' + bound + ', Re = ' + str(Re) + ', N = ' + str(N))
plt.figtext(0.5, 0.04, 'Dochází k přenulování tlaku uvnitř válce', ha='center', va='center', fontsize=10)
plt.show()

