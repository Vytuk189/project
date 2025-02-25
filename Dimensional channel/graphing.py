import numpy as np
import matplotlib.pyplot as plt

# Load the residues list from the text file
residues = []
with open('residues.txt', 'r') as file:
    for line in file:
        # Assuming each line in the file is a string representation of a NumPy array
        residues.append(np.fromstring(line.strip(), sep=' '))

# Load the past array from the binary file
past = np.load('past_array.npy', allow_pickle=True)

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

# Plot using pcolormesh
plt.figure(figsize=(10, 5))
plt.pcolormesh(X, Y, u_speeds_grid, cmap='jet', shading='auto')
plt.colorbar(label='U')
plt.title('Horizontální rychlost U')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
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
plt.title('Vertikální rychlost V')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
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
plt.title('Tlak P')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
plt.show()


# Replace 'data.txt' with your filename
data = np.loadtxt('residues.txt')


# Assuming the file has three columns:
var1 = data[:, 0]
var2 = data[:, 1]
var3 = data[:, 2]

# Create an iteration number vector (starting from 1)
iterations = np.arange(1, data.shape[0] + 1)

# Plot the variables
plt.figure(figsize=(10, 6))
plt.plot(iterations, var1, label='P Reziduum', color="red")
plt.plot(iterations, var2, label='U Reziduum', color="blue")
plt.plot(iterations, var3, label='V Reziduum', color="green")

# Set y-axis to logarithmic scale
plt.yscale('log')
plt.xlim(0, 300000)

# Labeling the plot
plt.xlabel('Iterace')
plt.ylabel('Reziduum')
plt.title('Hodnota reziduí')
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Display the plot
plt.show()





# Define the function u = f(x)
def f(x):
    deltaP = 1.6
    mu = 0.63
    L = 2
    H = 1
    return -deltaP/(2*mu*L)*(x**2-H*x)  # Example function


# Generate y values (ranging from 0 to 1)
y_values = np.linspace(0, 1, 100)

# Calculate corresponding u values using the function
u_values = f(y_values)

# Create the plot
plt.plot(u_values, y_values, linestyle='--', color='blue', label="teor." )
#
plt.scatter(u_speeds[-1, :], y_vals, color='red', marker='x', label="numer.")

plt.xlim(-0.005, 0.16)

# Label the axes
plt.xlabel('Horizontální rychlost U')
plt.ylabel('y')

# Show the plot
plt.title('Porovnání teoretického a numerického řešení')
plt.legend()
plt.grid(True)
plt.show()







