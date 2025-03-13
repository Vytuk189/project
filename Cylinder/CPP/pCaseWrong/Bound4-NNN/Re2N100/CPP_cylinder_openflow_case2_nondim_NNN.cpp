/*
This is a MacCormack solver for 2D Poissel flow of incompressible Newtonian fluid around a cylinder in an openflow

The flow is solved for non-dimensional pressure and speeds

The boundary conditions in this specific case are Neumann (Case NNN)

The flow is viscous - cylinder is placed in the middle and the flow distortion should be more or less symetrical (wont be completely because of the assymetry in the boundary conditions)
*/


#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>


std::vector<double> scaleVector(double scale, std::vector<double> vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scale*vec[i];
    }
    return result;
}

std::vector<double> addVectors(std::vector<double> vec1, std::vector<double> vec2) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

std::vector<double> subtractVectors(std::vector<double> vec1, std::vector<double> vec2) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

std::vector<double> multiplyMatrixWithVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
   
   // Initialize the result vector with the correct size (number of rows in the matrix)
    std::vector<double> result(matrix.size(), 0.0);

    // Perform matrix-vector multiplication
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}


// Flow and channel conditions
const double rho = 1.0; // Density
const double L = 2.0; // Length of channel
const double H = 1.5; // Height of channel
const double Re = 2;

/* 1/Re takes the same place in the formula with nondimensional values,
     as nu in the formula with dimensional values.
    The solver was written for dimensional values, this line changes it to nondimensional
*/
const double nu = 1.0 / Re;

const double beta = 1.0;

std::vector<std::vector<double>> D ={{1, 0, 0},
                                    {0, 1, 0},
                                    {0, 0, 1}};

std::vector<std::vector<double>> Dbeta = {{1/(beta*beta), 0, 0},
                                          {     0       , 1, 0},
                                          {     0       , 0, 1}};

std::vector<std::vector<double>> invDbeta = {{ beta*beta,  0, 0},
                                             {  0       ,  1, 0},
                                             {  0       ,  0, 1}};                                       


const double cfl = 0.5;
const int iter = 800000; // Number of iterations
const int N = 100; // Number of points in x direction
const int M = static_cast<int>(H / (L / N)); // Number of points in y direction
const double h = L / N; // Space step

const double u_in = 1.0;
const double p_in = 0.0;
const double p_out = 0.0;

const double center_x = L/2;
const double center_y = H/2;
const double R = L/20;

using Matrix = std::vector<std::vector<std::vector<double>>>;
// Create a matrix with (N+3) rows and (M+1) columns, each element is an array of 6 values
// # Assign each element as an array of 6 values - p, u, v, x, y, in/out of circle
Matrix createInitialMatrix() {
    Matrix initial_matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(6, 0.0)));

    double distance = 0.0;
    double fluid = 1.0;
    
    for (int i = 0; i < N + 3; ++i) {
        for (int j = 0; j < M + 1; ++j) {
            double x = -h + i * h;
            double y = j * h;

            // Initializes the area uniformly
            double p = p_in;
            double u = u_in;
            double v = 0.;

            fluid = 1.0;

            // Calculate the distance between the point (x, y) and the center (center_x, center_y)
            distance = std::sqrt((x - center_x) * (x - center_x) + (y - center_y) * (y - center_y));

            // Check if the point is within the circle
            if (distance <= R) {
                fluid = 0.0;  // The point is inside the circle
            }
        
            
            initial_matrix[i][j][0] = p;
            initial_matrix[i][j][1] = u;
            initial_matrix[i][j][2] = v;
            initial_matrix[i][j][3] = x;
            initial_matrix[i][j][4] = y;
            initial_matrix[i][j][5] = fluid;
        }
    }
    
    return initial_matrix;
}

void saveMatrix(const Matrix& past, int iter_number) {
    std::string filename = "flowdata_cylinder2_nondimensional_CPP_N" + std::to_string(N) + "_iter" + std::to_string(iter_number) + "_CFL05_beta"+std::to_string(beta)+".txt";
    std::ofstream file(filename);
    for (const auto& row : past) {
        for (const auto& col : row) {
            file << "[ ";
            for (double val : col) {
                file << val << " ";
            }
            file << "] ";
        }
       file << std::endl; 
    }

     // Close the file
    file.close();
    std::cout << "Data saved to " << filename << std::endl;

}

void saveResidues(const std::vector<std::vector<double>>& residues, int iter_number) {
    std::string filename = "residues_cylinder2_nondimensional_CPP_N" + std::to_string(N) + "_iter" + std::to_string(iter_number) + "_CFL05_beta"+std::to_string(beta)+".txt";
    std::ofstream file(filename);

    // Iterate through the matrix and write each element on a new line
    for (const auto& row : residues) {
        for (const auto& element : row) {
            file << element << " ";  // Write the element to the file, each on a new line
        } file << std::endl;
    }

    // Close the file
    file.close();
    std::cout << "Residues saved to " << filename << std::endl;
}


int main(){

    Matrix initial_matrix = createInitialMatrix();

    Matrix past = initial_matrix;
    Matrix predict = past;
    Matrix deltasPredict = Matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(3, 0.0)));
    Matrix deltasCorrect = Matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(3, 0.0)));

    std::vector<std::vector<double>> residues;


    int counter = 0;
    double a = 1; //Will store whether a point is in fluid or in circle
    while (counter < iter) {
    
        double umax_now = 0.0;
        double vmax_now = 0.0;
        double p_max = -1e10;
        double p_min = 1e10;
        
        // Extract u_speeds, v_speeds, pressures to find max values
        for (int i = 0; i < N + 3; ++i) {
            for (int j = 0; j < M + 1; ++j) {
                double p = past[i][j][0];
                double u = past[i][j][1];
                double v = past[i][j][2];
                umax_now = std::max(umax_now, std::abs(u));
                vmax_now = std::max(vmax_now, std::abs(v));
                p_max = std::max(p_max, p);
                p_min = std::min(p_min, p);
            }
        }
        
        double tau = cfl * 1.0 / (((umax_now + std::sqrt(umax_now * umax_now + beta * beta)) / h) +
                                    ((vmax_now + std::sqrt(vmax_now * vmax_now + beta * beta)) / h) + 
                                    2.0 * nu * (1.0 / (h * h) + 1.0 / (h * h)));
    
        std::vector<double> sums = {0.0, 0.0, 0.0};
        
        // Predictor section
        for (int i = 1; i < N + 2; ++i) {
            for (int j = 1; j < M; ++j) {
                // Extract data for the current grid point and its neighbors
                double p_this = past[i][j][0];
                double u_this = past[i][j][1];
                double v_this = past[i][j][2];
                
                double p_left = past[i - 1][j][0];
                double u_left = past[i - 1][j][1];
                double v_left = past[i - 1][j][2];
                
                double p_right = past[i + 1][j][0];
                double u_right = past[i + 1][j][1];
                double v_right = past[i + 1][j][2];
                
                double p_up = past[i][j + 1][0];
                double u_up = past[i][j + 1][1];
                double v_up = past[i][j + 1][2];
                
                double p_down = past[i][j - 1][0];
                double u_down = past[i][j - 1][1];
                double v_down = past[i][j - 1][2];


                std::vector<double> W_this = {p_this, u_this, v_this};
                std::vector<double> W_left = {p_left, u_left, v_left};
                std::vector<double> W_right = {p_right, u_right, v_right};
                std::vector<double> W_up = {p_up, u_up, v_up};
                std::vector<double> W_down = {p_down, u_down, v_down};
                
                std::vector<double> F_this = {u_this, u_this * u_this + p_this / rho, u_this * v_this};
                std::vector<double> F_left = {u_left, u_left * u_left + p_left / rho, u_left * v_left};
                
                std::vector<double> G_this = {v_this, v_this * u_this, v_this * v_this + p_this / rho};
                std::vector<double> G_down = {v_down, v_down * u_down, v_down * v_down + p_down / rho};
                

                //////////////////
                // Calculating the delta values
                std::vector<double> deltaWpredict(3, 0.0);
                // Do the matrix multiplication and other computations here for deltaWpredict
                double inv_h = 1/h;
                double inv_h_square = 1/(h*h);
                std::vector<double> brack_hor = scaleVector(inv_h_square, subtractVectors( addVectors (W_right, W_left), scaleVector(2.0, W_this)  ) );
                std::vector<double> brack_ver = scaleVector(inv_h_square, subtractVectors( addVectors (W_up, W_down), scaleVector(2.0, W_this)  ) );

                std::vector<double> temp = scaleVector(nu, multiplyMatrixWithVector(D, addVectors(brack_hor, brack_ver)));

                std::vector<double> brack = scaleVector(-inv_h, addVectors(subtractVectors(F_this, F_left), subtractVectors(G_this, G_down)));
                deltaWpredict = multiplyMatrixWithVector(invDbeta, addVectors(temp, brack));
                deltasPredict[i][j] = deltaWpredict;

                //std::vector<double> W_predict = addVectors(W_predict, scaleVector(tau, deltaWpredict));

                a = initial_matrix[i][j][5];
                predict[i][j][0] = past[i][j][0] + tau*deltaWpredict[0];
                predict[i][j][1] = a*(past[i][j][1] + tau*deltaWpredict[1]);
                predict[i][j][2] = a*(past[i][j][2] + tau*deltaWpredict[2]);
            }
        }
        
         // Top, bottom and right border
        // Enforcing Neuman Conditions
        for (int k = 0; k < N + 3; ++k) {
            // Pressure/U speed/V speed dont change between the bottom (or top) two layers
            predict[k][0][0] = predict[k][1][0];          // Bottom boundary condition
            predict[k][M][0] = predict[k][M-1][0];      // Top  boundary condition
            predict[k][0][1] = predict[k][1][1];          // Bottom boundary condition
            predict[k][M][1] = predict[k][M-1][1];      // Top  boundary condition
            predict[k][0][2] = predict[k][1][2];          // Bottom boundary condition
            predict[k][M][2] = predict[k][M-1][2];      // Top  boundary condition
        }

        for (int j = 0; j < M + 1; ++j) {
            // Pressure/U speed/V speed dont change between the last layer and the layer before it
            predict[N+2][j][0] = predict[N+1][j][0];          
            predict[N+2][j][1] = predict[N+1][j][1];          
            predict[N+2][j][2] = predict[N+1][j][2]; 
            
            // Neumann condition on the left border
            // Pressure doesnt change between the first and second layer
            predict[0][j][0] = predict[1][j][0];
            predict[0][j][1] = predict[1][j][1];
            predict[0][j][2] = predict[1][j][2];
        }





        // Corrector section
        for (int i = 1; i < N + 2; ++i) {
            for (int j = 1; j < M; ++j) {
                // Corrector section similar to the predictor
                // Calculate deltasCorrect here, similar to deltasPredict
                // Extract data for the current grid point and its neighbors
                double p_this = predict[i][j][0];
                double u_this = predict[i][j][1];
                double v_this = predict[i][j][2];
                
                double p_left = predict[i - 1][j][0];
                double u_left = predict[i - 1][j][1];
                double v_left = predict[i - 1][j][2];
                
                double p_right = predict[i + 1][j][0];
                double u_right = predict[i + 1][j][1];
                double v_right = predict[i + 1][j][2];
                
                double p_up = predict[i][j + 1][0];
                double u_up = predict[i][j + 1][1];
                double v_up = predict[i][j + 1][2];
                
                double p_down = predict[i][j - 1][0];
                double u_down = predict[i][j - 1][1];
                double v_down = predict[i][j - 1][2];


                std::vector<double> W_this = {p_this, u_this, v_this};
                std::vector<double> W_left = {p_left, u_left, v_left};
                std::vector<double> W_right = {p_right, u_right, v_right};
                std::vector<double> W_up = {p_up, u_up, v_up};
                std::vector<double> W_down = {p_down, u_down, v_down};
                
                std::vector<double> F_this = {u_this, u_this * u_this + p_this / rho, u_this * v_this};
                std::vector<double> F_right = {u_right, u_right * u_right + p_right / rho, u_right * v_right};
                
                std::vector<double> G_this = {v_this, v_this * u_this, v_this * v_this + p_this / rho};
                std::vector<double> G_up = {v_up, v_up * u_up, v_up * v_up + p_up / rho};

                double inv_h = 1/h;
                double inv_h_square = 1/(h*h);
                std::vector<double> brack_hor = scaleVector(inv_h_square, subtractVectors( addVectors (W_right, W_left), scaleVector(2.0, W_this)  ) );
                std::vector<double> brack_ver = scaleVector(inv_h_square, subtractVectors( addVectors (W_up, W_down), scaleVector(2.0, W_this)  ) );

                std::vector<double> temp = scaleVector(nu, multiplyMatrixWithVector(D, addVectors(brack_hor, brack_ver)));

                std::vector<double> brack = scaleVector(-inv_h, addVectors(subtractVectors(F_right, F_this), subtractVectors(G_up, G_this)));
                std::vector<double> deltaWCorrect = multiplyMatrixWithVector(invDbeta, addVectors(temp, brack));

                deltasCorrect[i][j] = deltaWCorrect;

                // Compute residue and save values
                
                std::vector<double> deltaW = scaleVector(0.5, addVectors(deltasPredict[i][j], deltasCorrect[i][j])); // Replace with actual deltaW calculation
                sums[0] += deltaW[0] * deltaW[0];
                sums[1] += deltaW[1] * deltaW[1];
                sums[2] += deltaW[2] * deltaW[2];
                a = initial_matrix[i][j][5];
                past[i][j][0] = past[i][j][0] + tau*deltaW[0];
                past[i][j][1] = a*(past[i][j][1] + tau*deltaW[1]);
                past[i][j][2] = a*(past[i][j][2] + tau*deltaW[2]);


            }
        }

        
         // Top, bottom and left and right border
        // Enforcing Neuman Conditions
        for (int k = 0; k < N + 3; ++k) {
            // Pressure/U speed/V speed dont change between the bottom (or top) two layers
            past[k][0][0] = past[k][1][0];          // Bottom boundary condition
            past[k][M][0] = past[k][M-1][0];      // Top  boundary condition
            past[k][0][1] = past[k][1][1];          // Bottom boundary condition
            past[k][M][1] = past[k][M-1][1];      // Top  boundary condition
            past[k][0][2] = past[k][1][2];          // Bottom boundary condition
            past[k][M][2] = past[k][M-1][2];      // Top  boundary condition
        }

        for (int j = 0; j < M + 1; ++j) {
            // Pressure/U speed/V speed dont change between the last layer and the layer before it
            past[N+2][j][0] = past[N+1][j][0];          
            past[N+2][j][1] = past[N+1][j][1];          
            past[N+2][j][2] = past[N+1][j][2];     

            // Neumann condition on the left border
            // Pressure/U speed/V doesnt change between the first and second layer
            past[0][j][0] = past[1][j][0];
            past[0][j][1] = past[1][j][1];
            past[0][j][2] = past[1][j][2];

        }



        


        
        double residual_p = std::sqrt(sums[0] / tau)/(N*M);
        double residual_u = std::sqrt(sums[1] / tau)/(N*M);
        double residual_v = std::sqrt(sums[2] / tau)/(N*M);

        residues.push_back({residual_p, residual_u, residual_v});
        
        if(counter % 10 == 0){
        std::cout << "***************************************************************\n";
        std::cout << "Case DNN: N=" << N << ", M=" << M << ", L=" << L << ", H=" << H << ", Re=" << Re << ", beta=" << beta <<  "\n";
        std::cout << "Iteration: " << counter << "/" << iter << "\n";
        std::cout << "U max: " << umax_now << "\n";
        std::cout << "V max: " << vmax_now << "\n";
        std::cout << "P max: " << p_max << "          P min:" << p_min << "\n";
        std::cout << "Tau: " << tau << "\n";
        std::cout << "Residuum u: " << residual_u << "\n";
        std::cout << "Residuum v: " << residual_v << "\n";
        std::cout << "Residuum p: " << residual_p << "\n";
        }

        if(counter % 1000 == 0){
            saveMatrix(past, counter);
            saveResidues(residues, counter);
        }


        counter++;
    }
    

    saveMatrix(past, iter);
    saveResidues(residues, iter);

    return 0;

}