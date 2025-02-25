import numpy as np
import matplotlib.pyplot as plt

def inCircle(x, y, center_x, center_y, R):
    
    fluid = 1
    
    distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    if distance <= R:
        fluid = 0
                
    return fluid

#Podmínky dle domluvy
u_max = 1 #m/s
L = 2 #Length of channel
H = 1 #Height of channel
Re = 10

nu = u_max*H/Re

mu = 0.63 #oil

rho = mu/nu


#beta = u_max*(rho)**0.5
beta = 1

D = np.array([[0, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

Dbeta = np.array([[1/(rho*beta**2), 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

invDbeta = np.linalg.inv(Dbeta)


N = 400 #No. of points in x direction

h = L/N #Space step
M = int(H/h) #No. of points in y direction

#Defining obstacle - circle

R = 0.2
center_x = L/3
center_y = H/2


u_in = u_max
p_in = 1.6
p_out = 0

# Create a matrix with (N+3) rows and (M+1) columns, each element is an array of 3 values
initial_matrix = np.empty((N+3, M+1), dtype=object)

# Assign each element as an array of 5 values
# p, u, v, x, y
        
#Rozliti jen tlaku
for i in range(N+3):
    for j in range(M+1):
        x = -h+i*h
        y = j*h
        fluid = inCircle(x, y, center_x, center_y, R)
        initial_matrix[i, j] = np.array([x, y, fluid])
        
        
# Create three 2D matrices to store the first, second, and third subelements
inside = np.empty_like(initial_matrix, dtype=np.float64)
x_coords = np.empty_like(initial_matrix, dtype=np.float64)
y_coords = np.empty_like(initial_matrix, dtype=np.float64)

# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(initial_matrix.shape[0]):  # Iterate over rows
    for j in range(initial_matrix.shape[1]):  # Iterate over columns
        inside[i, j] = initial_matrix[i, j][2]  # Third subelement
        x_coords[i, j] = initial_matrix[i, j][0]
        y_coords[i, j] = initial_matrix[i, j][1]

# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
x_vals = np.unique(x_coords)
y_vals = np.unique(y_coords)

# Generate a 2D grid based on unique x and y coordinates
X, Y = np.meshgrid(x_vals, y_vals)

# Initialize an empty grid for u_speeds
inside_grid = np.zeros((len(y_vals), len(x_vals)))

# Map the u_speeds values to their correct coordinates
for i in range(len(x_coords)):
    x_index = np.searchsorted(x_vals, x_coords[i])
    y_index = np.searchsorted(y_vals, y_coords[i])
    inside_grid[y_index, x_index] = inside[i]

# Plot using pcolormesh
plt.figure(figsize=(10, 5))
plt.pcolormesh(X, Y, inside_grid, cmap='jet', shading='auto')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure axes are proportional
plt.show()


