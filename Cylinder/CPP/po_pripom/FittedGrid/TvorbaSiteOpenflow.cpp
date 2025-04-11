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

using Matrix = std::vector<std::vector<std::vector<double>>>;


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

void printVector(const std::vector<double>& vec) {
    for(int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}


// Flow and channel conditions
const double rho = 1.0; // Density
const double L = 2.0; // Length of channel
const double H = 1; // Height of channel
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

const double u_in = 1.0;
const double p_in = 0.0;
const double p_out = 0.0;


const double cfl = 0.5;
const int iter = 800000; // Number of iterations
const int N = 400; // Number of points in x direction
const int M = static_cast<int>(H / (L / N)); // Number of points in y direction
const double h = L / N; // Space step


const double max_node_density = 25; //Pocet uzlu na R
const double min_node_density = 2;
const int left_zone_count = 5;
const int right_zone_count = 10;
const int vert_zone_count = 5;

const double center_x = 0.25*L;
const double center_y = 0.5*H;
const double R = 0.125*H;

const double x_dense_left = center_x - 1.5*R;
const double x_dense_right = center_x + 2*R;
const double y_dense_bottom = center_y - 1.5*R;
const double y_dense_top = center_y + 1.5*R;

Matrix TvorbaSiteFreeFlow() {

    // Creation of vector of x coordinates
    std::vector<double> x_coordinates;
    double x_temp = 0;
    int x_zone_number_left = 0;
    double delta_x_min = R/max_node_density;
    double delta_x_max = R/min_node_density;

    double delta_x = delta_x_max;
    while(x_temp < x_dense_left) {
        x_coordinates.push_back(x_temp);
        if(x_temp < (x_zone_number_left+1)*x_dense_left/left_zone_count){
            x_temp = x_temp + delta_x;
        }
        else {
            x_zone_number_left++;
            delta_x = delta_x_max - x_zone_number_left*(delta_x_max-delta_x_min)/left_zone_count;
            x_temp = x_temp + delta_x;
        }
    }

    x_temp = x_dense_left;
    while(x_temp <= x_dense_right-delta_x_min) {
        x_coordinates.push_back(x_temp);
        x_temp = x_temp + delta_x_min;
    }

    x_temp = x_dense_right;
    int x_zone_number_right = 0;
    delta_x = delta_x_min + (x_zone_number_right+1)*(delta_x_max-delta_x_min)/right_zone_count;

    while(x_temp < L) {
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
    double delta_y_min = R/max_node_density;
    double delta_y_max = R/min_node_density;

    double delta_y = delta_y_max;
    while(y_temp < y_dense_bottom) {
        y_coordinates.push_back(y_temp);
        if(y_temp < (y_zone_number_bottom+1)*y_dense_bottom/vert_zone_count){
            y_temp = y_temp + delta_y;
        }
        else {
            y_zone_number_bottom++;
            delta_y = delta_y_max - y_zone_number_bottom*(delta_y_max-delta_y_min)/vert_zone_count;
            y_temp = y_temp + delta_y;
        }
    }

    y_temp = y_dense_bottom;
    while(y_temp < y_dense_top) {
        y_coordinates.push_back(y_temp);
        y_temp = y_temp + delta_y_min;
    }

    y_temp = y_dense_top;
    int y_zone_number_top = 0;
    delta_y = delta_y_min + (y_zone_number_top+1)*(delta_y_max-delta_y_min)/vert_zone_count;

    while(y_temp < H) {
        y_coordinates.push_back(y_temp);
        if(y_temp < y_dense_top + (y_zone_number_top+1)*(H-y_dense_top)/vert_zone_count){
            y_temp = y_temp + delta_y;
        }
        else {
            y_zone_number_top++;
            delta_y = delta_y_min + (y_zone_number_top+1)*(delta_y_max-delta_y_min)/vert_zone_count;
            y_temp = y_temp + delta_y;
        }
    }

    y_coordinates.push_back(H);

    std::cout << "Počet bodů v y směru je: " << y_coordinates.size() << std::endl;

    int N = x_coordinates.size();
    int M = y_coordinates.size();

    Matrix coordinates_matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(2, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = x_coordinates[i];
            double y = y_coordinates[j];

            coordinates_matrix[i][j][0] = x;
            coordinates_matrix[i][j][1] = y;
        }
    }

    std::cout << coordinates_matrix.size() << std::endl;
    std::cout << coordinates_matrix[0].size() << std::endl;

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



int main() {

    Matrix coordinates = TvorbaSiteFreeFlow();
    saveMeshGrid(coordinates);

    return 0;
}

