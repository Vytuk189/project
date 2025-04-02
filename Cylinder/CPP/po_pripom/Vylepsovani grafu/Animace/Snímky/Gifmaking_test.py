

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.animation as animation

# Parameters
N = 50
Re = 80
start_iter = 0
max_iter = 100000  # maximum iteration value
step = 100         # step between iterations

# Create a list of file names for every 100 iterations
file_names = [
    'flowdata_cylinder_nondimensional_N'+str(N)+'_Re'+str(Re)+'.000000_iter'+str(iter_val)+'.txt'
    for iter_val in range(0, max_iter+1, step)
]

def load_data(file_path):
    """
    Reads the text file and returns:
        time: simulation time (float)
        X, Y: mesh grid coordinates (2D arrays)
        u_speeds_grid, v_speeds_grid: velocity fields (2D arrays)
        mask: binary mask for the cylinder region (2D boolean array)
    """
    data = []
    with open(file_path, 'r') as file:
        # First line contains the time value
        time = float(file.readline().strip())
        for line in file:
            line = line.strip()
            # Split the line by ']' and process each bracketed element
            elements = line.split('] [')
            row = []
            for element in elements:
                element = element.strip('[]')
                row.append(np.array([float(x) for x in element.split()]))
            data.append(row)

    # Convert list to 2D numpy array of objects
    past = np.empty((len(data), len(data[0])), dtype=object)
    for i in range(len(data)):
        for j in range(len(data[i])):
            past[i, j] = data[i][j]

    # Prepare arrays for each subelement
    pressures = np.empty_like(past, dtype=np.float64)
    u_speeds = np.empty_like(past, dtype=np.float64)
    v_speeds = np.empty_like(past, dtype=np.float64)
    x_coords = np.empty_like(past, dtype=np.float64)
    y_coords = np.empty_like(past, dtype=np.float64)
    structure = np.empty_like(past, dtype=np.float64)

    for i in range(past.shape[0]):
        for j in range(past.shape[1]):
            pressures[i, j] = past[i, j][0]
            u_speeds[i, j] = past[i, j][1]
            v_speeds[i, j] = past[i, j][2]
            x_coords[i, j] = past[i, j][3]
            y_coords[i, j] = past[i, j][4]
            structure[i, j] = past[i, j][5]

    # Create grid from unique x and y coordinates
    x_vals = np.unique(x_coords)
    y_vals = np.unique(y_coords)
    X, Y = np.meshgrid(x_vals, y_vals)

    # Map structure and velocity data to a grid
    structure_grid = np.zeros((len(y_vals), len(x_vals)))
    u_speeds_grid = np.zeros((len(y_vals), len(x_vals)))
    v_speeds_grid = np.zeros((len(y_vals), len(x_vals)))
    
    # Flatten the coordinate arrays for indexing
    for i in range(len(x_coords)):
        x_index = np.searchsorted(x_vals, x_coords[i])
        y_index = np.searchsorted(y_vals, y_coords[i])
        structure_grid[y_index, x_index] = structure[i]
        u_speeds_grid[y_index, x_index] = u_speeds[i]
        v_speeds_grid[y_index, x_index] = v_speeds[i]

    # Create a mask where structure is zero (assumed to be the fluid region)
    mask = structure_grid == 0
    

    return time, X, Y, u_speeds_grid, v_speeds_grid, mask

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 5))
cmap_binary = ListedColormap(["none", "black"])

def update(frame):
    file_path = file_names[frame]
    try:
        time, X, Y, u_speeds_grid, v_speeds_grid, mask = load_data(file_path)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return

    ax.clear()  # Clear the previous frame

    # Plot streamlines
    ax.streamplot(
        X, Y, u_speeds_grid, v_speeds_grid,
        color='black', linewidth=0.7, density=1.5,
        integration_direction='both', arrowstyle='-', minlength=0.2
    )
    
    # Overlay the cylinder (mask)
    ax.pcolormesh(X, Y, mask, shading='auto', cmap=cmap_binary)
    
    ax.set_title(f'Time = {time:.3f}')
    ax.axis('off')
    print("On file: " + file_path)

# Create the animation: frames = number of files
ani = animation.FuncAnimation(fig, update, frames=len(file_names), repeat=False)

# Save the animation as a GIF
ani.save('flow_animation.gif', writer='pillow', dpi=200, fps=20)

plt.show()
