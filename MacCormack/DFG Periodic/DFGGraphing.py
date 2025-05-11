import numpy as np
import matplotlib.pyplot as plt


# # Initialize the list to store differences and times
# differences_and_times = []

# # Loop through the file names from Counter0.txt to Counter9990.txt with a step of 10
# for counter in range(0, 10000, 10):
#     print(counter)
#     # Construct the file path dynamically based on the counter
#     file_path = f'DFGPeriodic_DataPart_Counter{counter}.txt'
    
#     # Initialize a list to store rows
#     data = []
    
#     # Open and read the file
#     with open(file_path, 'r') as file:
#         # Read the first line to get the double value
#         first_line = file.readline().strip()
#         time = float(first_line)  # Extract the time
        
#         # Process the rest of the file
#         for line in file:
#             line = line.strip()  # Clean the line to remove extra spaces and newlines
#             elements = line.split('] [')  # Split the line by ']' and process each element in brackets
#             row = []
#             for element in elements:
#                 element = element.strip('[]')  # Remove the surrounding '[' and ']' characters
#                 # Convert the string elements into float and store them as a numpy array
#                 row.append(np.array([float(x) for x in element.split()]))
#             data.append(row)
    
#     # Convert the list of rows into a 2D numpy array (of numpy arrays, each of 5 values)
#     past = np.empty((len(data), len(data[0])), dtype=object)
    
#     # Assign the numpy arrays into the final 2D array
#     for i in range(len(data)):
#         for j in range(len(data[i])):
#             past[i, j] = data[i][j]
    
#     # Create a matrix to store the pressures
#     pressures = np.empty_like(past, dtype=np.float64)
    
#     # Extract the pressures from the past matrix
#     for i in range(past.shape[0]):  # Iterate over rows
#         for j in range(past.shape[1]):  # Iterate over columns
#             pressures[i, j] = past[i, j][0]  # First subelement is the pressure
    
#     # Get the specific pressure values at indices (41, 83) and (89, 83)
#     p1 = pressures[41][83]
#     p2 = pressures[89][83]
    
#     # Calculate the difference
#     diff = p1 - p2
    
#     # Store the result as a tuple (difference, time)
#     differences_and_times.append((diff, time))
    
# # Extract the times and differences
# times = [entry[1] for entry in differences_and_times]
# differences = [entry[0] for entry in differences_and_times]

# # Subtract the first time from all the time values to make time start at 0
# time_offset = times[0]
# times = [t - time_offset for t in times]

# # Plot the difference over time
# plt.plot(times, differences, marker='o', linestyle='-', color='b')
# plt.xlabel('Time (s)')
# plt.ylabel('Pressure Difference')
# plt.title('Pressure Difference Over Time')
# plt.grid(True)
# plt.show()

# # Save the differences and times into a text file
# with open('pressure_differences.txt', 'w') as output_file:
#     output_file.write("Time (s)\tPressure Difference\n")
#     for time, diff in zip(times, differences):
#         output_file.write(f"{time:.6f}\t{diff:.6f}\n")

# print("Data has been saved to 'pressure_differences.txt'")


loaded_times = []
loaded_differences = []

# Open the file to read the data
with open('pressure_differences.txt', 'r') as input_file:
    # Skip the header line
    header = input_file.readline()
    
    # Read each subsequent line and extract the time and pressure difference values
    for line in input_file:
        # Split the line by tab character to separate time and difference
        time_str, diff_str = line.split('\t')
        
        # Convert the values from string to float
        time = float(time_str)
        diff = float(diff_str)
        
        # Append the values to the lists
        loaded_times.append(time)
        loaded_differences.append(diff)

# Print out the loaded data
print("Loaded times:", loaded_times)
print("Loaded pressure differences:", loaded_differences)

times = loaded_times
differences = loaded_differences

# Subtract the first time from all the time values to make time start at 0
time_offset = times[0]
times = [t - time_offset for t in times]

# Plot the difference over time
plt.plot(times, differences, marker='o', linestyle='-', color='r', markersize=1)
plt.xlabel('Čas (s)')
plt.ylabel(r'$\Delta p$', fontsize=14)
plt.title('Rozdíl tlaků')
plt.grid(True)
# Set x-axis limits
plt.xlim(0, 0.4)  # Adjust these limits as per your data

plt.show()



