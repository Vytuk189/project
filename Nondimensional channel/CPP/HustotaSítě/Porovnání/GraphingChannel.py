import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

# Initialize a list to store rows
data = []
file_path = 'N20.txt'
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


u_speeds20 = np.empty_like(past, dtype=np.float64)
y_coords = np.empty_like(past, dtype=np.float64)


# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(past.shape[0]):  # Iterate over rows
    for j in range(past.shape[1]):  # Iterate over columns
        u_speeds20[i, j] = past[i, j][1]  # Second subelement
        y_coords[i, j] = past[i, j][4]
        
# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
y_vals20 = np.unique(y_coords)








# Initialize a list to store rows
data = []
file_path = 'N50.txt'
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


u_speeds50 = np.empty_like(past, dtype=np.float64)
y_coords = np.empty_like(past, dtype=np.float64)


# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(past.shape[0]):  # Iterate over rows
    for j in range(past.shape[1]):  # Iterate over columns
        u_speeds50[i, j] = past[i, j][1]  # Second subelement
        y_coords[i, j] = past[i, j][4]
        
# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
y_vals50 = np.unique(y_coords)



# Initialize a list to store rows
data = []
file_path = 'N100.txt'
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


u_speeds100 = np.empty_like(past, dtype=np.float64)
y_coords = np.empty_like(past, dtype=np.float64)


# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(past.shape[0]):  # Iterate over rows
    for j in range(past.shape[1]):  # Iterate over columns
        u_speeds100[i, j] = past[i, j][1]  # Second subelement
        y_coords[i, j] = past[i, j][4]
        
# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
y_vals100 = np.unique(y_coords)




# Initialize a list to store rows
data = []
file_path = 'N200.txt'
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


u_speeds200 = np.empty_like(past, dtype=np.float64)
y_coords = np.empty_like(past, dtype=np.float64)


# Loop over each element in the past matrix to extract the corresponding subelements
for i in range(past.shape[0]):  # Iterate over rows
    for j in range(past.shape[1]):  # Iterate over columns
        u_speeds200[i, j] = past[i, j][1]  # Second subelement
        y_coords[i, j] = past[i, j][4]
        
# Assuming 'u_speeds', 'x_coords', and 'y_coords' arrays are already available

# Create the grid using x_coords and y_coords
y_vals200 = np.unique(y_coords)








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
plt.plot(u_values, y_values, linestyle='--', color='black', label="Exakt." )
#

u_speeds_plot20 = u_speeds20[-1, :]
plt.scatter(u_speeds_plot20, y_vals20, color='orange', marker='^', label="20x10", s=8)


u_speeds_plot50 = u_speeds50[-1, :]
plt.scatter(u_speeds_plot50, y_vals50, color='red', marker='o', label="50x25", s=5)


u_speeds_plot100 = u_speeds100[-1, :]
plt.scatter(u_speeds_plot100, y_vals100, color='purple', marker='v', label="100x50", s=5)

u_speeds_plot200 = u_speeds200[-1, :]
plt.scatter(u_speeds_plot200, y_vals200, color='blue', marker='x', label="200x0", s=5)




plt.xlim(0, 0.10)
plt.ylim(0, 1)

# Label the axes
plt.xlabel('Horizontální rychlost U')
plt.ylabel('y')

# Show the plot
plt.title('Exaktní a numerické řešení')
plt.legend(fontsize=8)
plt.grid(True)

plt.show()


# Create the plot
plt.plot(u_values, y_values, linestyle='--', color='black', label="Exakt." )
#


u_speeds_plot20 = u_speeds20[-1, :]
plt.scatter(u_speeds_plot20, y_vals20, color='orange', marker='^', label="20x10", s=8)


u_speeds_plot50 = u_speeds50[-1, :]
plt.scatter(u_speeds_plot50, y_vals50, color='red', marker='o', label="50x25", s=5)


u_speeds_plot100 = u_speeds100[-1, :]
plt.scatter(u_speeds_plot100, y_vals100, color='purple', marker='v', label="100x50", s=5)

u_speeds_plot200 = u_speeds200[-1, :]
plt.scatter(u_speeds_plot200, y_vals200, color='blue', marker='x', label="200x100", s=5)




plt.xlim(0.08, 0.10)
plt.ylim(0.2, 0.8)

# Label the axes
plt.xlabel('Horizontální rychlost U')
plt.ylabel('y')

# Show the plot
plt.title('Exaktní a numerické řešení (detail)')
plt.legend(fontsize=8)
plt.grid(True)
# Format x-axis to show three decimal places
plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.3f}'))

plt.show()



print(max(u_speeds_plot20))
print(max(u_speeds_plot50))
print(max(u_speeds_plot100))
print(max(u_speeds_plot200))




