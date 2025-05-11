/*
This is a MacCormack solver for 2D Poissel flow of incompressible Newtonian fluid above a cavity

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
const double L = 2; // Length of channel
const double H = 1.5; // Height of channel
const double Re = 10;

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


const double cfl = 0.1;
const int iter = 800000; // Number of iterations
const int N = 100; // Number of points in x direction
const int M = static_cast<int>(H / (L / N)); // Number of points in y direction
const double h = L / N; // Space step

const double u_in = 0.1;
const double p_in = 0;
const double p_out = 0.0;

const double side_cavity = 0.5; // Length of the side of the square cavity
const double x_topleft_corner = (L-side_cavity)/2; // Top left corner of the cavity
const double y_topleft_corner = side_cavity;
const double x_bottomright_corner = x_topleft_corner+side_cavity;
const double y_bottomright_corner = 0;

using Matrix = std::vector<std::vector<std::vector<double>>>;

// Create a matrix with (N+3) rows and (M+1) columns, each element is an array of 6 values
// # Assign each element as an array of 6 values - p, u, v, x, y, in/out of circle
Matrix createInitialMatrix() {
    Matrix initial_matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(6, 0.0)));

    double distance = 0.0;
    double fluid = 1.0;
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = i * h;
            double y = j * h;

            // Initializes the area uniformly
            double p = p_in;
            double u = 0;
            double v = 0.;

            fluid = 1.0;

            // Determine whether the point is in fluid or in wall
            if(x <= x_topleft_corner && y<=y_topleft_corner){
                fluid = 0;
            }

            if(x >= x_bottomright_corner && y<=y_topleft_corner){
                fluid = 0;
            }

            if(y==0){
                fluid = 0;
            }
        
            
            initial_matrix[i][j][0] = p;
            initial_matrix[i][j][1] = u;
            initial_matrix[i][j][2] = v;
            initial_matrix[i][j][3] = x;
            initial_matrix[i][j][4] = y;
            initial_matrix[i][j][5] = fluid;
        }
    }

    double y_wall_actual;
    int j_wall_actual;
    double y_wall_top = initial_matrix[0][M-1][4];

    for (int j = 0; j < M; j++){
        if(initial_matrix[0][j][5]==1){
            y_wall_actual = initial_matrix[0][j][4];
            j_wall_actual = j;
            break;
        }
    }

    double H_actual = initial_matrix[0][M-1][4] - initial_matrix[0][j_wall_actual][4];

    for(int i = 0; i<N; i++ ){
        for(int j = j_wall_actual; j<M; j++){
            double y_temp = initial_matrix[i][j][4] - y_wall_actual;
            double u = -4*u_in*(y_temp*y_temp-y_temp*H_actual)/(H_actual*H_actual);
            initial_matrix[i][j][1] = u;
        }
    }
 
    return initial_matrix;
}

void saveMatrix(const Matrix& past, int iter_number) {
    std::string filename = "flowdata_cavity_N" + std::to_string(N) + "_Re"+std::to_string(Re)+"_iter" + std::to_string(iter_number) + ".txt";
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

void saveMatrixWithTime(const Matrix& past, int iter_number, double time) {
    std::string filename = "flowdata_cavity_N" + std::to_string(N) + "_Re"+std::to_string(Re)+"_iter" + std::to_string(iter_number) + ".txt";
    std::ofstream file(filename);
    file << time << std::endl;
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
    std::string filename = "residues_cavity_N" + std::to_string(N) + "_Re"+std::to_string(Re)+"_iter" + std::to_string(iter_number) + ".txt";
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

    std::vector<std::vector<double>> residues;
    std::vector<double> times;
    saveMatrix(initial_matrix, 0);
    

    int counter = 0;
    double time = 0;
    double a = 1; //Will store whether a point is in fluid or in circle
    double cfl = 0.1;

    double umax_now = 0.0;
    double vmax_now = 0.0;
    double p_max = -1e10;
    double p_min = 1e10;
        
    while (counter < iter) {

        umax_now = 0.0;
        vmax_now = 0.0;
        p_max = -1e10;
        p_min = 1e10;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
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
                                    2.0 * nu * (1.0 / (h* h) + 1.0 / (h*h)));

       //tau = 0.0001;
    
        std::vector<double> sums = {0.0, 0.0, 0.0};
        
        // Predictor section
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < M - 1; ++j) {

                //checks to see if point is in fluid or in wall
                a = initial_matrix[i][j][5];
                

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
                double central_second_p_x = (2*delta_x_left/(delta_x_left+delta_x_right)*p_right - 2*p_this + 2*delta_x_right/(delta_x_left+delta_x_right)*p_left)/(delta_x_right*delta_x_left);
    
                
                double central_first_u_y = 0.5*((u_up-u_this)/delta_y_up+(u_this-u_down)/delta_y_down);
                double central_first_v_y = 0.5*((v_up-v_this)/delta_y_up+(v_this-v_down)/delta_y_down);
                double central_first_p_y = 0.5*((p_up-p_this)/delta_y_up+(p_this-p_down)/delta_y_down);
                double backward_u_y = (u_this-u_down)/delta_y_down;
                double backward_v_y = (v_this-v_down)/delta_y_down;
                double backward_p_y = (p_this-p_down)/delta_y_down;
                double central_second_u_y = (2*delta_y_down/(delta_y_down+delta_y_up)*u_up - 2*u_this + 2*delta_y_up/(delta_y_down+delta_y_up)*u_down)/(delta_y_down*delta_y_up);
                double central_second_v_y = (2*delta_y_down/(delta_y_down+delta_y_up)*v_up - 2*v_this + 2*delta_y_up/(delta_y_down+delta_y_up)*v_down)/(delta_y_down*delta_y_up);
                double central_second_p_y = (2*delta_y_down/(delta_y_down+delta_y_up)*p_up - 2*p_this + 2*delta_y_up/(delta_y_down+delta_y_up)*p_down)/(delta_y_down*delta_y_up);
       

                predict[i][j][0] = p_this + beta*beta*tau*(central_second_p_x + central_second_p_y - rho*(backward_u_x+backward_v_y));
                predict[i][j][1] = a*(u_this + tau*(-u_this*backward_u_x - v_this*backward_u_y - 1/rho*backward_p_x + nu*(central_second_u_x + central_second_u_y)));
                predict[i][j][2] = a*(v_this + tau*(-u_this*backward_v_x - v_this*backward_v_y - 1/rho*backward_p_y + nu*(central_second_v_x + central_second_v_y)));



            }
        }


        //Boundary conditions for predict matrix
        // Pressure at walls
        for (int k = 0; k < N; ++k) {
            predict[k][0][0] = predict[k][1][0];          // Bottom wall boundary condition
            predict[k][M-1][0] = predict[k][M-2][0];      // Top wall boundary condition
            //Dirichlet conditions are enforced from initial conditions
        }

        // Left and right boundary conditions
        for (int m = 0; m < M; ++m) {
            predict[0][m][0] = predict[1][m][0];          // Left Neuman for pressure
            predict[N-1][m][1] = predict[N - 2][m][1];  // Right Neuman for speeds
            predict[N-1][m][2] = predict[N - 2][m][2];  
        }


        // Corrector section
        for (int i = 1; i < N -1; ++i) {
            for (int j = 1; j < M-1; ++j) {

                //checks to see if point is in fluid or in wall
                a = initial_matrix[i][j][5];
                
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
                double central_second_p_x = (2*delta_x_left/(delta_x_left+delta_x_right)*p_right - 2*p_this + 2*delta_x_right/(delta_x_left+delta_x_right)*p_left)/(delta_x_right*delta_x_left);


                double central_first_u_y = 0.5*((u_up-u_this)/delta_y_up+(u_this-u_down)/delta_y_down);
                double central_first_v_y = 0.5*((v_up-v_this)/delta_y_up+(v_this-v_down)/delta_y_down);
                double central_first_p_y = 0.5*((p_up-p_this)/delta_y_up+(p_this-p_down)/delta_y_down);
                double forward_u_y = (u_up-u_this)/delta_y_up;
                double forward_v_y = (v_up-v_this)/delta_y_up;
                double forward_p_y = (p_up-p_this)/delta_y_up;
                double central_second_u_y = (2*delta_y_down/(delta_y_down+delta_y_up)*u_up - 2*u_this + 2*delta_y_up/(delta_y_down+delta_y_up)*u_down)/(delta_y_down*delta_y_up);
                double central_second_v_y = (2*delta_y_down/(delta_y_down+delta_y_up)*v_up - 2*v_this + 2*delta_y_up/(delta_y_down+delta_y_up)*v_down)/(delta_y_down*delta_y_up);
                double central_second_p_y = (2*delta_y_down/(delta_y_down+delta_y_up)*p_up - 2*p_this + 2*delta_y_up/(delta_y_down+delta_y_up)*p_down)/(delta_y_down*delta_y_up);
       

                double new_p = p_this + beta*beta*tau*(central_second_p_x + central_second_p_y - rho*(forward_u_x+forward_v_y));
                double new_u = a*(0.5*(u_this + u_past) + tau/2*(-u_this*forward_u_x - v_this*forward_u_y - 1/rho*forward_p_x + nu*(central_second_u_x + central_second_u_y)));
                double new_v = a*(0.5*(v_this + v_past) + tau/2*(-u_this*forward_v_x - v_this*forward_v_y - 1/rho*forward_p_y + nu*(central_second_v_x + central_second_v_y)));


                sums[0] += (new_p - past[i][j][0])*(new_p - past[i][j][0]);
                sums[1] += (new_u - past[i][j][1])*(new_u - past[i][j][1]);
                sums[2] += (new_v - past[i][j][2])*(new_v - past[i][j][2]);

                past[i][j][0] = new_p;
                past[i][j][1] = new_u;
                past[i][j][2] = new_v;
                
                


            }
        }

        //Boundary conditions for predict matrix
        // Pressure at walls
        for (int k = 0; k < N; ++k) {
            past[k][0][0] = past[k][1][0];          // Bottom wall boundary condition
            past[k][M-1][0] = past[k][M-2][0];      // Top wall boundary condition
            //Dirichlet conditions are enforced from initial conditions
        }

        // Left and right boundary conditions
        for (int m = 0; m < M; ++m) {
            past[0][m][0] = past[1][m][0];          // Left Neuman for pressure
            past[N-1][m][1] = past[N - 2][m][1];  // Right Neuman for speeds
            past[N-1][m][2] = past[N - 2][m][2];  
        }


        
        double residual_p = std::sqrt(sums[0])/(N*M);
        double residual_u = std::sqrt(sums[1])/(N*M);
        double residual_v = std::sqrt(sums[2])/(N*M);

        residues.push_back({residual_p, residual_u, residual_v});
        
        time = time + tau;
        times.push_back(time);
        if(counter % 100 == 0){
        std::cout << "***************************************************************\n";
        std::cout << "Case: N=" << N << ", M=" << M << ", L=" << L << ", H=" << H << ", Re=" << Re << ", nu=" << nu <<  "\n";
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

 


        if(counter % 10000 == 0){
            //saveMatrix(past, counter);
            saveMatrixWithTime(past, counter, time);
        }

        if(counter % 50000 == 0){
            saveResidues(residues, counter);
        }
       



        counter++;
    }

    return 0;

}