#This is a MacCormack solver for 2D Poissel flow of incompressible Newtonian fluid
#between parallel plates

#The flow is solved for non-dimensional pressure and speeds

import numpy as np
import copy

def parabolic_function(y, delta_p, L, rho, Re, H):
    u = -delta_p*Re*H*(y**2-y)/(2*L*rho)
    return u 

# A function that returns 0 if the point is in circle, 1 if the point is in fluid
def inCircle(x, y, center_x, center_y, R):
    
    fluid = 1
    
    distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    if distance <= R:
        fluid = 0
                
    return fluid


#Flow and channel conditions
rho = 10 #Density
#Square control area
L = 2 #Length of area
H = L/2 #Height of area
Re = 1

# 1/Re takes the same place in the formula with nondimensional values,
# as nu in the formula with dimensional values.
# The solver was written for dimensional values, this line changes it to nondimensional
nu = 1/Re


beta = 1

D = np.array([[0, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

Dbeta = np.array([[1/(beta**2), 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

invDbeta = np.linalg.inv(Dbeta)

iterations = 100000
N = 200 #No. of points in x direction
M = int(N/2) #No. of points in y direction
h = L/N #Space step


#Defining obstacle - circle
R = L/20
center_x = L/3
center_y = H/2

#Nondimensional values
u_in = 1
p_in = 0

# Create a matrix with (N+3) rows and (M+1) columns, each element is an array of 3 values
initial_matrix = np.empty((N+3, M+1), dtype=object)

# Assign each element as an array of 6 values
# p, u, v, x, y, 
        
#Initial conditions
for i in range(N+3):
    for j in range(M+1):
        x = -h+i*h
        y = j*h
        
        #   Initializes the area uniformly
        p = p_in
        u = u_in
        v = 0
        
        # Finds whether point is in fluid (1) or in object (0)
        fluid = inCircle(x, y, center_x, center_y, R)
                   

        initial_matrix[i, j] = np.array([p, u, v, x, y, fluid])
        

#Object
out_of_object = np.zeros_like(initial_matrix)
for i in range(N+3):
    for j in range(M+1):
        out_of_object[i, j] = initial_matrix[i,j][5]


     
residues = []
counter = 0

#Next step:
past = copy.deepcopy(initial_matrix)
predict = copy.deepcopy(past)
deltasPredict= np.zeros_like(initial_matrix)
deltasCorrect= np.zeros_like(initial_matrix)

# while residuum < :

while counter < iterations:
    u_speeds = np.array([[sub_array[1] for sub_array in row] for row in past])
    v_speeds = np.array([[sub_array[2] for sub_array in row] for row in past])
    
    
    umax_now = np.max(u_speeds)
    vmax_now = np.max(v_speeds)
    
    pressures = np.array([[sub_array[0] for sub_array in row] for row in past])
    p_max = np.max(pressures)
    p_min = np.min(pressures)
    
    
    
    tau = 0.5 * 1 / ( (umax_now + (umax_now**2 + beta**2)**0.5)/h + (vmax_now + (vmax_now**2 + beta**2)**0.5)/h + 2*nu*(1/(h**2) + 1/(h**2))  )
    
    suma = np.array([0,0,0])
    
    # #Predictor section
    
    # for i in range(1,N+2):
    #     for j in range(1,M):
    #         a_this = out_of_object[i,j]
    #         p_this = a_this*past[i, j][0]
    #         u_this = a_this*past[i, j][1]
    #         v_this = a_this*past[i, j][2]
            
    #         a_left = out_of_object[i-1, j]
    #         p_left = a_left*past[i-1, j][0]
    #         u_left = a_left*past[i-1, j][1]
    #         v_left = a_left*past[i-1, j][2]
            
    #         a_right = out_of_object[i+1, j]
    #         p_right = a_right*past[i+1, j][0]
    #         u_right = a_right*past[i+1, j][1]
    #         v_right = a_right*past[i+1, j][2]
            
    #         a_up = out_of_object[i, j+1]
    #         p_up = a_up*past[i, j+1][0]
    #         u_up = a_up*past[i, j+1][1]
    #         v_up = a_up*past[i, j+1][2]
            
    #         a_down = out_of_object[i, j-1]
    #         p_down = a_down*past[i, j-1][0]
    #         u_down = a_down*past[i, j-1][1]
    #         v_down = a_down*past[i, j-1][2]
            
    #         W_this = np.array([p_this, u_this, v_this])
    #         W_left = np.array([p_left, u_left, v_left])
    #         W_right = np.array([p_right, u_right, v_right])
    #         W_up = np.array([p_up, u_up, v_up])
    #         W_down = np.array([p_down, u_down, v_down])
            
            
            
    #         F_this = np.array([u_this, u_this**2 + p_this/rho, u_this*v_this])
    #         F_left = np.array([u_left, u_left**2 + p_left/rho, u_left*v_left])
            
    #         G_this = np.array([v_this, v_this*u_this , v_this**2 + p_this/rho])
    #         G_down = np.array([v_down, v_down*u_down, v_down**2 + p_down/rho])
            
    #         temp = nu*np.dot(D, (W_right - 2*W_this + W_left)/(h**2) + (W_up - 2*W_this + W_down)/(h**2) )
    #         deltaWpredict = np.dot(invDbeta, -(F_this-F_left)/h - (G_this-G_down)/h  +  temp )
            
    #         deltasPredict[i, j] = deltaWpredict         
            
    #         W_predict = W_this + tau*deltaWpredict
            
    #         #Matrix of predictor values
    #         predict[i,j][0] = W_predict[0]
    #         predict[i,j][1] = W_predict[1]
    #         predict[i,j][2] = W_predict[2] 
            
    #         #Pressure at walls
    #         for k in range(0,N+3):
    #             predict[k, 0][0] = predict[k, 1][0]
    #             predict[k,-1][0] = predict[k,-2][0]
                
    #         #Left and right boundary
    #         for m in range(0, M+1):
    #             predict[0, m][1] = predict[1, m][1]
    #             predict[-1, m][1] = predict[-2, m][1]
    #             predict[0, m][2] = predict[1, m][2]
    #             predict[-1, m][2] = predict[-2, m][2]
            
            
                    
            
    # #Corrector section

    # for i in range(1,N+2):
    #     for j in range(1,M):       
    #         a_this = out_of_object[i, j]
    #         p_this = a_this*predict[i, j][0]
    #         u_this = a_this*predict[i, j][1]
    #         v_this = a_this*predict[i, j][2]
            
    #         a_left = out_of_object[i-1, j]
    #         p_left = a_left*predict[i-1, j][0]
    #         u_left = a_left*predict[i-1, j][1]
    #         v_left = a_left*predict[i-1, j][2]
            
    #         a_right = out_of_object[i+1, j]
    #         p_right = a_right*predict[i+1, j][0]
    #         u_right = a_right*predict[i+1, j][1]
    #         v_right = a_right*predict[i+1, j][2]
            
    #         a_up = out_of_object[i, j+1]
    #         p_up = a_up*predict[i, j+1][0]
    #         u_up = a_up*predict[i, j+1][1]
    #         v_up = a_up*predict[i, j+1][2]
            
    #         a_down = out_of_object[i, j-1]
    #         p_down = a_down*predict[i, j-1][0]
    #         u_down = a_down*predict[i, j-1][1]
    #         v_down = a_down*predict[i, j-1][2]
            
    #         W_this = np.array([p_this, u_this, v_this])
    #         W_left = np.array([p_left, u_left, v_left])
    #         W_right = np.array([p_right, u_right, v_right])
    #         W_up = np.array([p_up, u_up, v_up])
    #         W_down = np.array([p_down, u_down, v_down])
            
             
    #         F_this = np.array([u_this, u_this**2 + p_this/rho, u_this*v_this])
    #         F_right = np.array([u_right, u_right**2 + p_right/rho, u_right*v_right])
            
    #         G_this = np.array([v_this, v_this*u_this, v_this**2 + p_this/rho])
    #         G_up = np.array([v_up, v_up*u_up, v_up**2 + p_up/rho])
            
    #         temp = nu*np.dot(D, (W_right - 2*W_this + W_left)/(h**2) + (W_up - 2*W_this + W_down)/(h**2)  )
    #         deltaWcorrect = np.dot(invDbeta, -(F_right-F_this)/h - (G_up-G_this)/h  + temp   )
    #         deltasCorrect[i,j] = deltaWcorrect
    
    #Predictor section
    
    for i in range(1,N+2):
        for j in range(1,M):
            p_this = past[i, j][0]
            u_this = past[i, j][1]
            v_this = past[i, j][2]
            
            p_left = past[i-1, j][0]
            u_left = past[i-1, j][1]
            v_left = past[i-1, j][2]
            
            p_right = past[i+1, j][0]
            u_right = past[i+1, j][1]
            v_right = past[i+1, j][2]
            
            p_up = past[i, j+1][0]
            u_up = past[i, j+1][1]
            v_up = past[i, j+1][2]
            
            p_down = past[i, j-1][0]
            u_down = past[i, j-1][1]
            v_down = past[i, j-1][2]
            
            W_this = np.array([p_this, u_this, v_this])
            W_left = np.array([p_left, u_left, v_left])
            W_right = np.array([p_right, u_right, v_right])
            W_up = np.array([p_up, u_up, v_up])
            W_down = np.array([p_down, u_down, v_down])
            
            
            
            F_this = np.array([u_this, u_this**2 + p_this/rho, u_this*v_this])
            F_left = np.array([u_left, u_left**2 + p_left/rho, u_left*v_left])
            
            G_this = np.array([v_this, v_this*u_this , v_this**2 + p_this/rho])
            G_down = np.array([v_down, v_down*u_down, v_down**2 + p_down/rho])
            
            temp = nu*np.dot(D, (W_right - 2*W_this + W_left)/(h**2) + (W_up - 2*W_this + W_down)/(h**2) )
            deltaWpredict = np.dot(invDbeta, -(F_this-F_left)/h - (G_this-G_down)/h  +  temp )
            
            deltasPredict[i, j] = deltaWpredict         
            
            W_predict = W_this + tau*deltaWpredict
            
            #Matrix of predictor values
            predict[i,j][0] = W_predict[0]
            predict[i,j][1] = W_predict[1]
            predict[i,j][2] = W_predict[2] 
            
    #Pressure at walls
    for k in range(0,N+3):
        predict[k, 0][0] = predict[k, 1][0]
        predict[k,-1][0] = predict[k,-2][0]
        
    #Left and right boundary
    for m in range(0, M+1):
        predict[0, m][1] = predict[1, m][1]
        predict[-1, m][1] = predict[-2, m][1]
        predict[0, m][2] = predict[1, m][2]
        predict[-1, m][2] = predict[-2, m][2]
        
                    
            
    #Corrector section

    for i in range(1,N+2):
        for j in range(1,M):       
            p_this = predict[i, j][0]
            u_this = predict[i, j][1]
            v_this = predict[i, j][2]
            
            p_left = predict[i-1, j][0]
            u_left = predict[i-1, j][1]
            v_left = predict[i-1, j][2]
            
            p_right = predict[i+1, j][0]
            u_right = predict[i+1, j][1]
            v_right = predict[i+1, j][2]
            
            p_up = predict[i, j+1][0]
            u_up = predict[i, j+1][1]
            v_up = predict[i, j+1][2]
            
            p_down = predict[i, j-1][0]
            u_down = predict[i, j-1][1]
            v_down = predict[i, j-1][2]
            
            W_this = np.array([p_this, u_this, v_this])
            W_left = np.array([p_left, u_left, v_left])
            W_right = np.array([p_right, u_right, v_right])
            W_up = np.array([p_up, u_up, v_up])
            W_down = np.array([p_down, u_down, v_down])
            
             
            F_this = np.array([u_this, u_this**2 + p_this/rho, u_this*v_this])
            F_right = np.array([u_right, u_right**2 + p_right/rho, u_right*v_right])
            
            G_this = np.array([v_this, v_this*u_this, v_this**2 + p_this/rho])
            G_up = np.array([v_up, v_up*u_up, v_up**2 + p_up/rho])
            
            temp = nu*np.dot(D, (W_right - 2*W_this + W_left)/(h**2) + (W_up - 2*W_this + W_down)/(h**2)  )
            deltaWcorrect = np.dot(invDbeta, -(F_right-F_this)/h - (G_up-G_this)/h  + temp   )
            deltasCorrect[i,j] = deltaWcorrect
            
            
            deltaW = 0.5*(deltasPredict[i,j]+deltasCorrect[i,j])
            suma = suma + deltaW*deltaW
            a = out_of_object[i,j]
            past[i,j][0] = a*(past[i,j][0] + tau*deltaW[0])
            past[i,j][1] = a*(past[i,j][1] + tau*deltaW[1])
            past[i,j][2] = a*(past[i,j][2] + tau*deltaW[2])
            
    
            
    # Top and bottom border
    # Enforcing Neuman conditions 
    for i in range(0,N+3):
        #Pressure/U speed/V speed dont change between the bottom (or top) two layers
        past[i, 0][0] = past[i, 1][0]
        past[i,-1][0] = past[i,-2][0]
        past[i, 0][1] = past[i, 1][1]
        past[i,-1][1] = past[i,-2][1]
        past[i, 0][2] = past[i, 1][2]
        past[i,-1][2] = past[i,-2][2]
        
        
        
    # Left and right boundary
    # Dirichlet conditions are enforced from initial setup - 
    # the scheme loop doesnt affect borders
    
    # # Enforcing Neuman conditions
    # for j in range(0, M+1):
    #     past[0, j][0] = past[1, j][0]
    #     past[-1, j][0] = past[-2, j][0]
        
    #     past[0, j][1] = past[1, j][1]
    #     past[-1, j][1] = past[-2, j][1]
        
    #     past[0, j][2] = past[1, j][2]
    #     past[-1, j][2] = past[-2, j][2]
        
    
    
    residue = np.sqrt(suma/tau)/(N*M)
    residues.append(residue)
    print("***************************************************************")
    print("Iteration: " + str(counter))
    print("U max = " + str(umax_now))
    print("V max = " + str(vmax_now))
    print("P max = " + str(p_max) + " P min = " + str(p_min))
    print("Tau = " + str(tau))
    print("Residuum p: " + str(residue[0]))
    print("Residuum u: " + str(residue[1]))
    print("Residuum v: " + str(residue[2]))
    counter = counter + 1
    
# Save the residues list to a text file
np.savetxt('residues_cylinder3'+str(N)+'_iter'+str(iterations)+'.txt', residues)

# Save the past array to a binary file
np.save('dataflow_cylinder3_N'+str(N)+'_iter'+str(iterations)+'.npy', past)
        


