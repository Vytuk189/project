import numpy as np
import copy



#Podmínky dle domluvy
u_max = 1 #m/s
L = 2 #Length of channel
H = 1 #Height of channel
Re = 10
#Predelat na bezrozmerove

nu = u_max*H/Re

mu = 0.63 #oil

rho = mu/nu


#beta = u_max*(rho)**0.5
#beta =
beta = 1

D = np.array([[0, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

Dbeta = np.array([[1/(rho*beta**2), 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

invDbeta = np.linalg.inv(Dbeta)


N = 40 #No. of points in x direction

h = L/N #Space step
M = int(H/h) #No. of points in y direction

u_in = 1
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
        p = 0
        u = u_in
        v = 0
        if i == 0:
            p = p_in
        if i == N+2:  # accesses last element
            p = p_out
        if j == 0 or j == M:  # no-slip condition
            u = 0
            v = 0

        initial_matrix[i, j] = np.array([p, u, v, x, y])



     
residues = []
counter = 0

#Next step:
past = copy.deepcopy(initial_matrix)
predict = copy.deepcopy(past)
deltasPredict= np.zeros_like(initial_matrix)
deltasCorrect= np.zeros_like(initial_matrix)

# while residuum < :

while counter < 500000:
    u_speeds = np.array([[sub_array[1] for sub_array in row] for row in past])
    v_speeds = np.array([[sub_array[2] for sub_array in row] for row in past])
    
    
    umax_now = np.max(u_speeds)
    vmax_now = np.max(v_speeds)
    
    pressures = np.array([[sub_array[0] for sub_array in row] for row in past])
    p_max = np.max(pressures)
    p_min = np.min(pressures)
    
    
    
    tau = 0.1 * 1 / ( (umax_now + (umax_now**2 + beta**2)**0.5)/h + (vmax_now + (vmax_now**2 + beta**2)**0.5)/h + 2*nu*(1/(h**2) + 1/(h**2))  )
    
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
            
            
            # #Pressure at walls
            # predict[:, 0][0] = predict[:, 1][0]
            # predict[:,-1][0] = predict[:,-2][0]
                
            
            # # Left and right boundaries:
            # predict[0, :][1] = predict[1, :][1]
            # predict[-1, :][1] = predict[-2, :][1]
            # predict[0, :][2] = predict[1, :][2]
            # predict[-1, :][2] = predict[-2, :][2]
                    
            
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
            
    
    suma = np.array([0,0,0])
    for i in range(1,N+2):
        for j in range(1,M):  
            deltaW = 0.5*(deltasPredict[i,j]+deltasCorrect[i,j])
            suma = suma + deltaW*deltaW
            past[i,j][0] = past[i,j][0] + tau*deltaW[0]
            past[i,j][1] = past[i,j][1] + tau*deltaW[1]
            past[i,j][2] = past[i,j][2] + tau*deltaW[2]
            
    #Pressure at walls
    for i in range(0,N+3):
        past[i, 0][0] = past[i, 1][0]
        past[i,-1][0] = past[i,-2][0]
        
    #Left and right boundary
    for j in range(0, M+1):
        past[0, j][1] = past[1, j][1]
        past[-1, j][1] = past[-2, j][1]
        
        past[0, j][2] = past[1, j][2]
        past[-1, j][2] = past[-2, j][2]
        
    
    
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
np.savetxt('residues.txt', residues)

# Save the past array to a binary file
np.save('past_array.npy', past)
        
