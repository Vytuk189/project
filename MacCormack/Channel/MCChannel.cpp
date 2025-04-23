/*
This is a MacCormack solver for 2D Poissel flow of incompressible Newtonian fluid
between parallel plates

The flow is solved for non-dimensional pressure and speeds
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
const double rho = 10.0; // Density
const double L = 2.0; // Length of channel
const double H = 1.0; // Height of channel
const double Re = 200;

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
const double p_in = 1.6;
const double p_out = 0.0;

using Matrix = std::vector<std::vector<std::vector<double>>>;

Matrix createInitialMatrix() {
    Matrix initial_matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(5, 0.0)));
    
    for (int i = 0; i < N + 3; ++i) {
        for (int j = 0; j < M + 1; ++j) {
            double x = -h + i * h;
            double y = j * h;
            double p = 0.0;
            double u = u_in;
            double v = 0.;
            
            if (i == 0) p = p_in;
            if (i == N + 2) p = p_out;
            if (j == 0 || j == M) { // no-slip condition
                u = 0.0;
                v = 0.0;
            }
            
            initial_matrix[i][j][0] = p;
            initial_matrix[i][j][1] = u;
            initial_matrix[i][j][2] = v;
            initial_matrix[i][j][3] = x;
            initial_matrix[i][j][4] = y;
        }
    }
    
    return initial_matrix;
}

void saveMatrix(const Matrix& past, int counter) {
    std::string filename = "flowdata_channel_nondimensional_N" + std::to_string(N) + "_iter" + std::to_string(counter) + "_CFL05_beta"+std::to_string(beta)+".txt";
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

void saveResidues(const std::vector<std::vector<double>>& residues, int counter) {
    std::string filename = "residues_channel_nondimensional_CPP_N" + std::to_string(N) + "_iter" + std::to_string(counter) + "_CFL05_beta"+std::to_string(beta)+".txt";
    std::ofstream file(filename);

    // Iterate through the matrix and write each element on a new line
    for (const auto& row : residues) {
        for (const auto& element : row) {
            file << element << " ";  // Write the element to the file, each on a new line
        } file << std::endl;
    }

    // Close the file
    file.close();
    std::cout << "Residues saved to " << "residues_CPP.txt" << std::endl;
}

void saveTimes(const std::vector<double>& times, int counter) {
    std::string filename = "times_channel_nondimensional_CPP_N" + std::to_string(N) + "_iter" + std::to_string(counter) + "_CFL05_beta"+std::to_string(beta)+".txt";
    std::ofstream file(filename);

    // Iterate through the matrix and write each element on a new line
    for (const auto& time : times) {
        file << time << std::endl;  // Write the element to the file, each on a new line
    }

    // Close the file
    file.close();
    std::cout << "Times saved to " << filename << std::endl;
}


int main(){

    Matrix initial_matrix = createInitialMatrix();

    Matrix past = initial_matrix;
    Matrix predict = past;

    std::vector<std::vector<double>> residues;
    std::vector<double> times;


    int counter = 0;
    double time = 0;
    

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
                double delta_x_left = past[i][j][3] - past[i - 1][j][3];
                
                double p_right = past[i + 1][j][0];
                double u_right = past[i + 1][j][1];
                double v_right = past[i + 1][j][2];
                double delta_x_right = past[i+1][j][3] - past[i][j][3];
                
                double p_up = past[i][j + 1][0];
                double u_up = past[i][j + 1][1];
                double v_up = past[i][j + 1][2];
                double delta_y_up = past[i][j + 1][4] - past[i][j][4];
                
                double p_down = past[i][j - 1][0];
                double u_down = past[i][j - 1][1];
                double v_down = past[i][j - 1][2];
                double delta_y_down = past[i][j][4] - past[i][j-1][4];


                
                double central_first_u_x = 0.5*((u_right-u_this)/delta_x_right+(u_this-u_left)/delta_x_left);
                double central_first_v_x = 0.5*((v_right-v_this)/delta_x_right+(v_this-v_left)/delta_x_left);
                double central_first_p_x = 0.5*((p_right-p_this)/delta_x_right+(p_this-p_left)/delta_x_left);
                double backward_u_x = (u_this-u_left)/delta_x_left;
                double backward_v_x = (v_this-v_left)/delta_x_left;
                double backward_p_x = (p_this-p_left)/delta_x_left;
                double central_second_u_x = (2*delta_x_left/(delta_x_left+delta_x_right)*u_right - 2*u_this + 2*delta_x_right/(delta_x_left+delta_x_right)*u_left)/(delta_x_right*delta_x_left);
                double central_second_v_x = (2*delta_x_left/(delta_x_left+delta_x_right)*v_right - 2*v_this + 2*delta_x_right/(delta_x_left+delta_x_right)*v_left)/(delta_x_right*delta_x_left);

    
                
                double central_first_u_y = 0.5*((u_up-u_this)/delta_y_up+(u_this-u_down)/delta_y_down);
                double central_first_v_y = 0.5*((v_up-v_this)/delta_y_up+(v_this-v_down)/delta_y_down);
                double central_first_p_y = 0.5*((p_up-p_this)/delta_y_up+(p_this-p_down)/delta_y_down);
                double backward_u_y = (u_this-u_down)/delta_y_down;
                double backward_v_y = (v_this-v_down)/delta_y_down;
                double backward_p_y = (p_this-p_down)/delta_y_down;
                double central_second_u_y = (2*delta_y_down/(delta_y_down+delta_y_up)*u_up - 2*u_this + 2*delta_y_up/(delta_y_down+delta_y_up)*u_down)/(delta_y_down*delta_y_up);
                double central_second_v_y = (2*delta_y_down/(delta_y_down+delta_y_up)*v_up - 2*v_this + 2*delta_y_up/(delta_y_down+delta_y_up)*v_down)/(delta_y_down*delta_y_up);
               
       

                predict[i][j][0] = p_this - beta*beta*rho*tau*(backward_u_x+backward_v_y);
                predict[i][j][1] = u_this + tau*(-u_this*backward_u_x - v_this*backward_u_y - 1/rho*backward_p_x + nu*(central_second_u_x + central_second_u_y));
                predict[i][j][2] = v_this + tau*(-u_this*backward_v_x - v_this*backward_v_y - 1/rho*backward_p_y + nu*(central_second_v_x + central_second_v_y));


            }
        }

        //Boundary conditions for predict matrix
        // Pressure at walls
        for (int k = 0; k < N + 3; ++k) {
            predict[k][0][0] = predict[k][1][0];          // Bottom wall boundary condition
            predict[k][M][0] = predict[k][M-1][0];      // Top wall boundary condition
        }

        // Left and right boundary conditions
        for (int m = 0; m <= M; ++m) {
            predict[0][m][1] = predict[1][m][1];          // Left boundary for second component
            predict[N + 2][m][1] = predict[N + 1][m][1];  // Right boundary for second component

            predict[0][m][2] = predict[1][m][2];          // Left boundary for third component
            predict[N + 2][m][2] = predict[N + 1][m][2];  // Right boundary for third component
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
                double delta_x_left = predict[i][j][3] - predict[i - 1][j][3];
                
                double p_right = predict[i + 1][j][0];
                double u_right = predict[i + 1][j][1];
                double v_right = predict[i + 1][j][2];
                double delta_x_right = predict[i+1][j][3] - predict[i][j][3];
                
                double p_up = predict[i][j + 1][0];
                double u_up = predict[i][j + 1][1];
                double v_up = predict[i][j + 1][2];
                double delta_y_up = predict[i][j + 1][4] - predict[i][j][4];
                
                double p_down = predict[i][j - 1][0];
                double u_down = predict[i][j - 1][1];
                double v_down = predict[i][j - 1][2];
                double delta_y_down = predict[i][j][4] - predict[i][j-1][4];

                double p_past = past[i][j][0];
                double u_past = past[i][j][1];
                double v_past = past[i][j][2];
                
                double central_first_u_x = 0.5*((u_right-u_this)/delta_x_right+(u_this-u_left)/delta_x_left);
                double central_first_v_x = 0.5*((v_right-v_this)/delta_x_right+(v_this-v_left)/delta_x_left);
                double central_first_p_x = 0.5*((p_right-p_this)/delta_x_right+(p_this-p_left)/delta_x_left);
                double forward_u_x = (u_right-u_this)/delta_x_right;
                double forward_v_x = (v_right-v_this)/delta_x_right;
                double forward_p_x = (p_right-p_this)/delta_x_right;
                double central_second_u_x = (2*delta_x_left/(delta_x_left+delta_x_right)*u_right - 2*u_this + 2*delta_x_right/(delta_x_left+delta_x_right)*u_left)/(delta_x_right*delta_x_left);
                double central_second_v_x = (2*delta_x_left/(delta_x_left+delta_x_right)*v_right - 2*v_this + 2*delta_x_right/(delta_x_left+delta_x_right)*v_left)/(delta_x_right*delta_x_left);

                double central_first_u_y = 0.5*((u_up-u_this)/delta_y_up+(u_this-u_down)/delta_y_down);
                double central_first_v_y = 0.5*((v_up-v_this)/delta_y_up+(v_this-v_down)/delta_y_down);
                double central_first_p_y = 0.5*((p_up-p_this)/delta_y_up+(p_this-p_down)/delta_y_down);
                double forward_u_y = (u_up-u_this)/delta_y_up;
                double forward_v_y = (v_up-v_this)/delta_y_up;
                double forward_p_y = (p_up-p_this)/delta_y_up;
                double central_second_u_y = (2*delta_y_down/(delta_y_down+delta_y_up)*u_up - 2*u_this + 2*delta_y_up/(delta_y_down+delta_y_up)*u_down)/(delta_y_down*delta_y_up);
                double central_second_v_y = (2*delta_y_down/(delta_y_down+delta_y_up)*v_up - 2*v_this + 2*delta_y_up/(delta_y_down+delta_y_up)*v_down)/(delta_y_down*delta_y_up);
               
       

                double new_p = 0.5*( p_this + p_past ) - beta*beta*rho*tau/2*(forward_u_x+forward_v_y);
                double new_u = 0.5*(u_this + u_past) + tau/2*(-u_this*forward_u_x - v_this*forward_u_y - 1/rho*forward_p_x + nu*(central_second_u_x + central_second_u_y));
                double new_v = 0.5*(v_this + v_past) + tau/2*(-u_this*forward_v_x - v_this*forward_v_y - 1/rho*forward_p_y + nu*(central_second_v_x + central_second_v_y));


                sums[0] += (new_p - past[i][j][0])*(new_p - past[i][j][0]);
                sums[1] += (new_u - past[i][j][1])*(new_u - past[i][j][1]);
                sums[2] += (new_v - past[i][j][2])*(new_v - past[i][j][2]);

                past[i][j][0] = new_p;
                past[i][j][1] = new_u;
                past[i][j][2] = new_v;

            }
        }

        //Boundary conditions for past matrix
        // Pressure at walls
        for (int k = 0; k < N + 3; ++k) {
            past[k][0][0] = past[k][1][0];          // Bottom wall boundary condition
            past[k][M][0] = past[k][M-1][0];      // Top wall boundary condition
        }

        // Left and right boundary conditions
        for (int m = 0; m <= M; ++m) {
            past[0][m][1] = past[1][m][1];          // Left boundary for second component
            past[N + 2][m][1] = past[N + 1][m][1];  // Right boundary for second component

            past[0][m][2] = past[1][m][2];          // Left boundary for third component
            past[N + 2][m][2] = past[N + 1][m][2];  // Right boundary for third component
        }
        


        
        double residual_p = std::sqrt(sums[0] / tau)/(N*M);
        double residual_u = std::sqrt(sums[1] / tau)/(N*M);
        double residual_v = std::sqrt(sums[2] / tau)/(N*M);


        residues.push_back({residual_p, residual_u, residual_v});
        times.push_back(time);
        
        if(counter%100 == 0) {
            std::cout << "***************************************************************\n";
            std::cout << "Case: N=" << N << ", beta=" << beta << "\n";
            std::cout << "Iteration: " << counter << "/" << iter << "\n";
            std::cout << "U max: " << umax_now << "\n";
            std::cout << "V max: " << vmax_now << "\n";
            std::cout << "P max: " << p_max << "          P min:" << p_min << "\n";
            std::cout << "Time: " << time << "\n";
            std::cout << "Tau: " << tau << "\n";
            std::cout << "Residuum u: " << residual_u << "\n";
            std::cout << "Residuum v: " << residual_v << "\n";
            std::cout << "Residuum p: " << residual_p << "\n";
        }

        time = time + tau;

        if(counter%10000 == 0) {
            saveMatrix(past, counter);
            saveResidues(residues, counter);
            saveTimes(times, counter);
        }

        
        counter++;
    }
    

    return 0;

}