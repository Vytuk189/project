import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

case = 2
bound = "Custom"
N = 400
iterations = 50000
Re = 10

# Initialize a list to store rows
data = []
file_path = 'flowdata_channel_nondimensional_N100_iter220000_CFL05_beta1.000000.txt'
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



# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])

# Plot the grid
plt.figure(figsize=(10, 5))

# Plot the grid lines for both x and y directions
for y in Y:  # Iterate over the rows (y-coordinates)
    plt.plot(x_vals, y, color='black', linestyle='-', linewidth=0.7)  # Plot horizontal lines

for x in X.T:  # Iterate over the columns (x-coordinates), transpose X for this
    plt.plot(x, y_vals, color='black', linestyle='-', linewidth=0.7)  # Plot vertical lines

N = len(x_vals)
M = len(y_vals)
DN = 10
plt.title("Síť " + str(N-3) + "x" + str(M-1), fontsize=18, pad=12) 
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(False)
# Set the x and y limits
plt.xlim(min(x_vals), max(x_vals))  # Set the x-axis limit based on x_vals
plt.ylim(min(y_vals), max(y_vals))  # Set the y-axis limit based on y_vals
plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal aspect ratio






# Initialize an empty grid for u_speeds
u_speeds_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    u_speeds_grid[y_index, x_index] = u_speeds[i]

# Plot using pcolormesh
plt.figure(figsize=(10, 5))
plt.pcolormesh(X, Y, u_speeds_grid, cmap='jet', shading='auto')
plt.colorbar(label='U')
plt.title('L.-F.: Horizontální rychlost U při síti ' + str(N-3) + "x" + str(M-1), fontsize=16, pad=12)
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional

# Plot the second mesh, only where binary_data is 0, set color to white
cmap_binary = ListedColormap(["none", "white"])
#plt.pcolormesh(X, Y, mask, shading='auto', cmap=cmap_binary)

plt.show()


# Initialize an empty grid for u_speeds
v_speeds_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    v_speeds_grid[y_index, x_index] = v_speeds[i]

# Plot using pcolormesh
plt.figure(figsize=(10, 5))
plt.pcolormesh(X, Y, v_speeds_grid, cmap='jet', shading='auto')
plt.colorbar(label='V')
plt.title('L.-F.: Vertikální rychlost V při síti ' + str(N-3) + "x" + str(M-1), fontsize=16, pad=12)
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
# Plot the second mesh, only where binary_data is 0, set color to white
cmap_binary = ListedColormap(["none", "white"])
#plt.pcolormesh(X, Y, mask, shading='auto', cmap=cmap_binary)

plt.show()

# Initialize an empty grid for u_speeds
pressures_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    pressures_grid[y_index, x_index] = pressures[i]

# Plot using pcolormesh
plt.figure(figsize=(10, 5))
plt.pcolormesh(X, Y, pressures_grid, cmap='jet', shading='auto')
plt.colorbar(label='P')
plt.title('L.-F.: Tlak P při síti ' + str(N-3) + "x" + str(M-1), fontsize=16, pad=12)
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
# Plot the second mesh, only where binary_data is 0, set color to white
cmap_binary = ListedColormap(["none", "white"])
#plt.pcolormesh(X, Y, mask, shading='auto', cmap=cmap_binary)
plt.show()


# Load the data
data_time = np.loadtxt('times_channel_nondimensional_CPP_N100_iter220000_CFL05_beta1.000000.txt')
data_res = np.loadtxt('residues_channel_nondimensional_CPP_N100_iter220000_CFL05_beta1.000000.txt')

# Assuming the file has three columns:
var1 = data_res[:, 0]
var2 = data_res[:, 1]
var3 = data_res[:, 2]

# Create an iteration number vector (starting from 1)
iterations = np.arange(1, data_res.shape[0] + 1)

# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot the residues
ax1.plot(iterations, var1, label='P', color="red")
ax1.plot(iterations, var2, label='U', color="blue")
ax1.plot(iterations, var3, label='V', color="green")

# Set the y-axis to logarithmic scale
ax1.set_yscale('log')
ax1.set_xlim(0, 125000)

# Labeling the plot
ax1.set_xlabel('Iterace', fontsize = 14)
ax1.set_ylabel('Reziduum', fontsize = 14)
plt.title(r'L.-F.: Rezidua pro $\,\xi$=0.01', fontsize=16, pad=12)
ax1.legend()
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

# Create a second x-axis for simulation time
ax2 = ax1.twiny()  # Create a second axis that shares the same y-axis
ax2.set_xlim(ax1.get_xlim())  # Ensure both axes have the same range for iterations
ax2.set_xlabel('Simulační čas [s]', fontsize = 12)
# Set x-ticks at intervals (optional: adjust the step size as needed)
ax2.set_xlim(0, 125000)

# Round the simulation time to one decimal place for tick labels
rounded_time = np.round(data_time[::20000], 1)  # Round the time values

# Set the tick labels for simulation time with rounded values
ax2.set_xticklabels(rounded_time)  # Set the rounded simulation time values

# Display the plot
plt.show()




# Define the function u = f(x)
def f(x):
    deltaP = -1.6
    rho = 10
    Re = 10
    L = 2
    H = 1
    return deltaP/(2*rho*L)*Re*H*(x**2-x)  # Example function


# Generate y values (ranging from 0 to 1)
y_values = np.linspace(0, 1, 100)

# Calculate corresponding u values using the function
u_values = f(y_values)

# Create the plot
plt.plot(u_values, y_values, linestyle='--', color='blue', label="teor." )
#

u_speeds_plot = u_speeds[-1, :]
plt.scatter(u_speeds_plot, y_vals, color='red', marker='x', label="numer.")

plt.xlim(-0.005, 0.16)

# Label the axes
plt.xlabel('Horizontální rychlost U')
plt.ylabel('y')

# Show the plot
plt.title('Porovnání teoretického a numerického řešení')
plt.legend()
plt.grid(True)
plt.show()

print(max(u_speeds_plot))









