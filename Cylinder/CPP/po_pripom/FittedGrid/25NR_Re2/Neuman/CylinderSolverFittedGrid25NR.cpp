/*
This is a MacCormack solver for 2D Poissel flow of incompressible Newtonian fluid around a cylinder in an openflow

The flow is solved for non-dimensional pressure and speeds

The boundary conditions in this specific case are give the supervisors instruction - Zero pressure gradient at left, top, bottom boundary, dirichlet on the right boundary
                                                                                   - Dirichlet velocity on the intake, zero velocity gradient everywhere else

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
const double L = 2; // Length of channel
const double H = 1; // Height of channel
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



const int iter = 800000; // Number of iterations

const double u_in = 1;
const double p_in = 0.0;
const double p_out = 0.0;

const double center_x = 0.45*L;
const double center_y = 0.5*H;
const double R = 0.05*H;


using Matrix = std::vector<std::vector<std::vector<double>>>;

const double max_node_density = 25; //Pocet uzlu na R
const double min_node_density = 0.5; //outlet
const double min_node_density_inlet = 0.5; //inlet
const double min_node_density_walls = 1;
const double delta_x_min = R/max_node_density;
const double delta_x_max = R/min_node_density;
const double delta_x_max_inlet = R/min_node_density_inlet;
const double delta_y_min = R/max_node_density;
const double delta_y_max = R/min_node_density_walls;
const int left_zone_count = 10;
const int right_zone_count = 10;
const int vert_zone_count = 10;


const double x_dense_left = center_x - 1.25*R;
const double x_dense_right = center_x + 1.75*R;
const double y_dense_bottom = center_y - 1.25*R;
const double y_dense_top = center_y + 1.25*R;

Matrix CreateMeshFreeFlow() {



    // Creation of vector of x coordinates
    std::vector<double> x_coordinates;
    double x_temp = 0;
    int x_zone_number_left = 0;
    std::cout << "Navrzeny nejmensi krok x: " << std::to_string(delta_x_min) << std::endl;
    std::cout << "Navrzeny nejmensi krok y: " << std::to_string(delta_x_min) << std::endl;
    

    double delta_x = delta_x_max_inlet;
    while(x_temp < x_dense_left) {
        x_coordinates.push_back(x_temp);
        if(x_temp < (x_zone_number_left+1)*x_dense_left/left_zone_count){
            x_temp = x_temp + delta_x;
        }
        else {
            x_zone_number_left++;
            delta_x = delta_x_max_inlet - x_zone_number_left*(delta_x_max_inlet-delta_x_min)/left_zone_count;
            x_temp = x_temp + delta_x;
        }
    }

    while(x_temp < x_dense_right) {
        x_coordinates.push_back(x_temp);
        x_temp = x_temp + delta_x_min;
    }

    
    int x_zone_number_right = 0;
    delta_x = delta_x_min + (x_zone_number_right+1)*(delta_x_max-delta_x_min)/right_zone_count;

    while(x_temp < L-0.5*delta_x_max) {
        x_coordinates.push_back(x_temp);
        if(x_temp < x_dense_right + (x_zone_number_right+1)*(L-x_dense_right)/right_zone_count){
            x_temp = x_temp + delta_x;
        }
        else {
            x_zone_number_right++;
            delta_x = delta_x_min + (x_zone_number_right+1)*(delta_x_max-delta_x_min)/right_zone_count;
            x_temp = x_temp + delta_x;
        }
    }

    x_coordinates.push_back(L);

    std::cout << "Počet bodů v x směru je: " << x_coordinates.size() << std::endl;

    // Creation of vector of y coordinates
    std::vector<double> y_coordinates;
    double y_temp = 0;
    int y_zone_number_bottom = 0;

    double delta_y = delta_y_max;
    while(y_temp < y_dense_bottom) {
        y_coordinates.push_back(y_temp);
        if(y_temp < (y_zone_number_bottom+1)*y_dense_bottom/vert_zone_count){
            y_temp = y_temp + delta_y;
        }
        else {
            y_zone_number_bottom++;
            delta_y = delta_y_max - (y_zone_number_bottom)*(delta_y_max-delta_y_min)/vert_zone_count;
            y_temp = y_temp + delta_y;
        }
    }

    while(y_temp < y_dense_top) {
        y_coordinates.push_back(y_temp);
        y_temp = y_temp + delta_y_min;
    }

    int y_zone_number_top = 0;
    delta_y = delta_y_min + (y_zone_number_top+1)*(delta_y_max-delta_y_min)/vert_zone_count;

    while(y_temp < H-0.5*delta_y_max) {
        y_coordinates.push_back(y_temp);
        if(y_temp < y_dense_top + (y_zone_number_top+1)*(H-y_dense_top)/vert_zone_count){
            y_temp = y_temp + delta_y;
        }
        else {
            y_zone_number_top++;
            delta_y = delta_y_min + (y_zone_number_top)*(delta_y_max-delta_y_min)/vert_zone_count;
            y_temp = y_temp + delta_y;
        }
    }

    y_coordinates.push_back(H);

    std::cout << "Počet bodů v y směru je: " << y_coordinates.size() << std::endl;

    int N = x_coordinates.size();
    int M = y_coordinates.size();

    double min_delta_x = 1000000;
    double min_delta_y = 1000000;

    Matrix coordinates_matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(2, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = x_coordinates[i];
            double y = y_coordinates[j];
            if(x < min_delta_x){
                min_delta_x = x;
            }
            if(y < min_delta_y){
                min_delta_y = y;
            }

            coordinates_matrix[i][j][0] = x;
            coordinates_matrix[i][j][1] = y;
        }
    }

    

    return coordinates_matrix;
}
void saveMeshGrid(const Matrix& past) {
    std::string filename = "MeshgridMaxDens"+std::to_string(max_node_density)+"MinDens"+std::to_string(min_node_density)+".txt";
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

// Create a matrix with (N) rows and (M) columns, each element is an array of 6 values
// # Assign each element as an array of 6 values - p, u, v, x, y, in/out of circle
Matrix createInitialMatrix(const Matrix& coordinates) {


    int N = coordinates.size();
    int M = coordinates[0].size();

    
    Matrix initial_matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(6, 0.0)));

    double distance = 0.0;
    double fluid = 1.0;
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = coordinates[i][j][0];
            double y = coordinates[i][j][1];

            // Initializes the area uniformly
            double p = p_in;
            //Creates a parabolic profile at the inlet
            double u = 1;
            double v = 0.;

            fluid = 1.0;

            // Calculate the distance between the point (x, y) and the center (center_x, center_y)
            distance = std::sqrt((x - center_x) * (x - center_x) + (y - center_y) * (y - center_y));

            // Check if the point is within the circle
            if (distance <= R) {
                fluid = 0.0;  // The point is inside the circle
                u = 0.0; //
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
    std::string filename = "flowdata_cylinder_NR"+std::to_string(max_node_density)+"_Re"+std::to_string(Re)+"_Iter"+ std::to_string(iter_number) + ".txt";
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
    std::string filename = "flowdata_cylinder_NR"+std::to_string(max_node_density)+"_Re"+std::to_string(Re)+"_Iter"+ std::to_string(iter_number) + ".txt";
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
    std::string filename = "residues_cylinder_NR"+std::to_string(max_node_density)+"_Re"+std::to_string(Re)+"_Iter"+ std::to_string(iter_number) + ".txt";
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

    Matrix coordinates = CreateMeshFreeFlow();

    Matrix initial_matrix = createInitialMatrix(coordinates);
    int N = coordinates.size();
    int M = coordinates[0].size();

    double min_delta_x = 100000000;
    double min_delta_y = 100000000;
    for (int i = 1; i < N; ++i) {
        for (int j = 1; j < M; ++j) {
           

            double x_right = coordinates[i][j][0];
            double y_up = coordinates[i][j][1];
            double x_this = coordinates[i-1][j][0];
            double y_this = coordinates[i][j-1][1];
            if(x_right - x_this < min_delta_x){
                min_delta_x = x_right - x_this;
            }
            if(y_up - y_this < min_delta_y){
                min_delta_y = y_up - y_this;
                std::cout << std::endl;
                std::cout << y_up << std::endl;
            }
        }
    }


    Matrix past = initial_matrix;
    Matrix predict = past;
    Matrix deltasPredict = Matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(3, 0.0)));
    Matrix deltasCorrect = Matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(3, 0.0)));

    std::vector<std::vector<double>> residues;
    saveMatrix(initial_matrix, 0);
    

    int counter = 0;
    double time = 0;
    double a = 1; //Will store whether a point is in fluid or in circle
    double cfl = 0.9;

    double umax_now = 0.0;
    double vmax_now = 0.0;
    double p_max = -1e10;
    double p_min = 1e10;

    std::cout << "Nejmensi x krok je: " << min_delta_x << std::endl;
    std::cout << "Nejmensi y krok je: " << min_delta_y << std::endl;
        
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
        
        
        double tau = cfl * 1.0 / (((umax_now + std::sqrt(umax_now * umax_now + beta * beta)) / min_delta_x) +
                                    ((vmax_now + std::sqrt(vmax_now * vmax_now + beta * beta)) / min_delta_y) + 
                                    2.0 * nu * (1.0 / (min_delta_y* min_delta_y) + 1.0 / (min_delta_y *min_delta_y)));

       //tau = 0.0001;
    
        std::vector<double> sums = {0.0, 0.0, 0.0};
        
        // Predictor section
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < M - 1; ++j) {
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
                double inv_delta_x_left = 1/delta_x_left;
                double inv_delta_x_right = 1/delta_x_right;
                double inv_delta_y_down = 1/delta_y_down;
                double inv_delta_y_up = 1/delta_y_up;

                double inv_delta_x_dot = inv_delta_x_left*inv_delta_x_right;
                double inv_delta_y_dot = inv_delta_y_up*inv_delta_y_down;

                std::vector<double> brack_W_hor_r = scaleVector(2*delta_x_left/(delta_x_left+delta_x_right), W_right);
                std::vector<double> brack_W_hor_l = scaleVector(2*delta_x_right/(delta_x_left+delta_x_right), W_left);
                std::vector<double> brack_hor = scaleVector(inv_delta_x_dot, subtractVectors( addVectors (brack_W_hor_r, brack_W_hor_l), scaleVector(2.0, W_this)  ) );


                std::vector<double> brack_W_ver_u = scaleVector(2*delta_y_down/(delta_y_down+delta_y_up), W_up);
                std::vector<double> brack_W_ver_d = scaleVector(2*delta_y_up/(delta_y_down+delta_y_up), W_down);
                std::vector<double> brack_ver = scaleVector(inv_delta_y_dot, subtractVectors( addVectors (brack_W_ver_u, brack_W_ver_d), scaleVector(2.0, W_this)  ) );

                std::vector<double> temp = scaleVector(nu, multiplyMatrixWithVector(D, addVectors(brack_hor, brack_ver)));

                std::vector<double> brack_F = scaleVector(-inv_delta_x_left, subtractVectors(F_this, F_left));
                std::vector<double> brack_G = scaleVector(-inv_delta_y_down, subtractVectors(G_this, G_down));
                std::vector<double> brack = addVectors(brack_F, brack_G);
                deltaWpredict = multiplyMatrixWithVector(invDbeta, addVectors(temp, brack));
                deltasPredict[i][j] = deltaWpredict;

                //std::vector<double> W_predict = addVectors(W_predict, scaleVector(tau, deltaWpredict));

                a = initial_matrix[i][j][5];
                predict[i][j][0] = past[i][j][0] + tau*deltaWpredict[0];
                predict[i][j][1] = a*(past[i][j][1] + tau*deltaWpredict[1]);
                predict[i][j][2] = a*(past[i][j][2] + tau*deltaWpredict[2]);
            }
        }
        
         // Top and bottom wall - no slip zero velocity is passed on from the initial conditons (the scheme loop doesnt affect borders)
        for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the bottom (or top) two layers
            predict[k][0][0] = predict[k][1][0];          // Bottom boundary condition
            predict[k][M-1][0] = predict[k][M-2][0];      // Top  boundary condition
            // Velocity doesnt change between the bottom (or top) two layers
            predict[k][0][1] = predict[k][1][1];          // Bottom boundary condition
            predict[k][M-1][1] = predict[k][M-2][1];      // Top  boundary condition
            predict[k][0][2] = predict[k][1][2];          // Bottom boundary condition
            predict[k][M-1][2] = predict[k][M-2][2];      // Top  boundary condition
        }

        for (int j = 0; j < M; ++j) {
            // Neuman condition for speed at the outlet
            //U speed/V speed dont change between the last layer and the layer before it   
            predict[N-1][j][1] = predict[N-2][j][1];          
            predict[N-1][j][2] = predict[N-2][j][2]; 
            
            // Neumann condition for pressure on the inlet
            // Pressure doesnt change between the first and second layer
            predict[0][j][0] = predict[1][j][0];
        }







        // Corrector section
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < M - 1; ++j) {
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


                std::vector<double> W_this = {p_this, u_this, v_this};
                std::vector<double> W_left = {p_left, u_left, v_left};
                std::vector<double> W_right = {p_right, u_right, v_right};
                std::vector<double> W_up = {p_up, u_up, v_up};
                std::vector<double> W_down = {p_down, u_down, v_down};
                
                std::vector<double> F_this = {u_this, u_this * u_this + p_this / rho, u_this * v_this};
                std::vector<double> F_right = {u_right, u_right * u_right + p_right / rho, u_right * v_right};
                
                std::vector<double> G_this = {v_this, v_this * u_this, v_this * v_this + p_this / rho};
                std::vector<double> G_up = {v_up, v_up * u_up, v_up * v_up + p_up / rho};

                double inv_delta_x_left = 1/delta_x_left;
                double inv_delta_x_right = 1/delta_x_right;
                double inv_delta_y_down = 1/delta_y_down;
                double inv_delta_y_up = 1/delta_y_up;

                double inv_delta_x_dot = inv_delta_x_left*inv_delta_x_right;
                double inv_delta_y_dot = inv_delta_y_up*inv_delta_y_down;

                std::vector<double> brack_W_hor_r = scaleVector(2*delta_x_left/(delta_x_left+delta_x_right), W_right);
                std::vector<double> brack_W_hor_l = scaleVector(2*delta_x_right/(delta_x_left+delta_x_right), W_left);
                std::vector<double> brack_hor = scaleVector(inv_delta_x_dot, subtractVectors( addVectors (brack_W_hor_r, brack_W_hor_l), scaleVector(2.0, W_this)  ) );


                std::vector<double> brack_W_ver_u = scaleVector(2*delta_y_down/(delta_y_down+delta_y_up), W_up);
                std::vector<double> brack_W_ver_d = scaleVector(2*delta_y_up/(delta_y_down+delta_y_up), W_down);
                std::vector<double> brack_ver = scaleVector(inv_delta_y_dot, subtractVectors( addVectors (brack_W_ver_u, brack_W_ver_d), scaleVector(2.0, W_this)  ) );

                std::vector<double> temp = scaleVector(nu, multiplyMatrixWithVector(D, addVectors(brack_hor, brack_ver)));

                std::vector<double> brack_F = scaleVector(-inv_delta_x_right, subtractVectors(F_right, F_this));
                std::vector<double> brack_G = scaleVector(-inv_delta_y_up, subtractVectors(G_up, G_this));
                std::vector<double> brack = addVectors(brack_F, brack_G);
                std::vector<double> deltaWCorrect = multiplyMatrixWithVector(invDbeta, addVectors(temp, brack));

                deltasCorrect[i][j] = deltaWCorrect;

                // Compute residue and save values
                

                //Check if here isnt the reason for why residues are not going down??? the difference comes from the edge of the circle
                std::vector<double> deltaW = scaleVector(0.5, addVectors(deltasPredict[i][j], deltasCorrect[i][j]));

                a = initial_matrix[i][j][5];

                double old_p = past[i][j][0];
                double old_u = past[i][j][1];
                double old_v = past[i][j][2];

                double new_p = past[i][j][0] + tau*deltaW[0];
                double new_u = a*(past[i][j][1] + tau*deltaW[1]);
                double new_v = a*(past[i][j][2] + tau*deltaW[2]);

                sums[0] += (new_p - old_p)*(new_p - old_p);
                sums[1] += (new_u - old_u)*(new_u - old_u);
                sums[2] += (new_v - old_v)*(new_v - old_v);


                past[i][j][0] = new_p;
                past[i][j][1] = new_u;
                past[i][j][2] = new_v;


            }
        }

        
         // Top and bottom wall - no slip zero velocity is passed on from the initial conditons (the scheme loop doesnt affect borders)
         for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the bottom (or top) two layers
            past[k][0][0] = past[k][1][0];          // Bottom boundary condition
            past[k][M-1][0] = past[k][M-2][0];      // Top  boundary condition
            // Velocity doesnt change between the bottom (or top) two layers
            past[k][0][1] = past[k][1][1];          // Bottom boundary condition
            past[k][M-1][1] = past[k][M-2][1];      // Top  boundary condition
            past[k][0][2] = past[k][1][2];          // Bottom boundary condition
            past[k][M-1][2] = past[k][M-2][2];      // Top  boundary condition
        }

        for (int j = 0; j < M; ++j) {
            // Neuman condition for speed at the outlet
            //U speed/V speed dont change between the last layer and the layer before it   
            past[N-1][j][1] = past[N-2][j][1];          
            past[N-1][j][2] = past[N-2][j][2]; 
            
            // Neumann condition for pressure on the inlet
            // Pressure doesnt change between the first and second layer
            past[0][j][0] = past[1][j][0];
        }


        


        
        double residual_p = std::sqrt(sums[0])/(N*M);
        double residual_u = std::sqrt(sums[1])/(N*M);
        double residual_v = std::sqrt(sums[2])/(N*M);

        residues.push_back({residual_p, residual_u, residual_v});
        
        time = time + tau;
        if(counter % 10 == 0){
        std::cout << "***************************************************************\n";
        std::cout << "Case Custom BC: N=" << N << ", M=" << M << ", L=" << L << ", H=" << H << ", Re=" << Re << ", beta=" << beta <<  "\n";
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

 


            if(counter % 1000 == 0){
                //saveMatrix(past, counter);
                saveMatrixWithTime(past, counter, time);
            }

            if(counter % 10000 == 0){
                saveResidues(residues, counter);
            }
       



        counter++;
    }

    return 0;

}