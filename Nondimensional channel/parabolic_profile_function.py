import numpy as np
import matplotlib.pyplot as plt

# Function for parabolic profile
def parabolic_function(y, delta_p, L, rho, Re, H):
    u = -delta_p*Re*H*(y**2-y)/(2*L*rho)
    return u 

# Horizontal speed vector of parabolic profile, creates the vector of y_values
def parabolic_profile1(delta_p, L, rho, Re, H, y_count):
    
    
    # Create numpy array of M elements ranging from 0 to 1
    y_values = np.linspace(0, 1, y_count)
    
    # Calculate corresponding u values using the parabolic function
    u_values = parabolic_function(y_values, delta_p, L, rho, Re, H)
    
    return u_values

# Horizontal speed vector of parabolic profile, takes in a vector of coordinates
def parabolic_profile2(delta_p, L, rho, Re, H, y_values):
    
    # Calculate corresponding u values using the parabolic function
    u_values = parabolic_function(y_values, delta_p, L, rho, Re, H)
    
    return u_values

y_count = 50

# Create numpy array of M elements ranging from 0 to 1
y_values = np.linspace(0, 1, y_count)

u_values = parabolic_profile2(1.6, 2, 10, 10, 1, y_values)
# Display the results
print("y values:", y_values)
print("u values:", u_values)
