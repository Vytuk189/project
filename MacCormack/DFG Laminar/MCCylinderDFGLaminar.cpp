/*
This is a modified Lax-Friedrichs solver for 2D Poissel flow of incompressible Newtonian fluid around a cylinder in an openflow

It is using the non-conservative NS equation

The flow is solved for non-dimensional pressure and speeds

The flow is viscous - cylinder is placed in the middle and the flow distortion should be more or less symetrical (wont be completely because of the assymetry in the boundary conditions
*/


#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

// Flow and channel conditions
const double rho = 1.0; // Density
const double L = 2.2; // Length of channel
const double H = 0.41; // Height of channel
const double Re = 20;

/* 1/Re takes the same place in the formula with nondimensional values,
     as nu in the formula with dimensional values.
    The solver was written for dimensional values, this line changes it to nondimensional
*/





const double beta = 1.0;
const double zeta = 0.001;

                      



const int iter = 2000000; // Number of iterations

const double u_in = 0.3;
const double p_in = 0.0;
const double p_out = 0.0;

const double center_x = 0.2;
const double center_y = 0.2;
const double R = 0.05;


const double nu = 0.001;


using Matrix = std::vector<std::vector<std::vector<double>>>;

const double max_node_density = 25; //Pocet uzlu na R
const double delta_x_min = R/max_node_density;
const double delta_y_min = R/max_node_density;

const double rate_of_step_increase = 1.05;
const double rate_of_step_increase_vert = 1.01;


const double x_dense_left = center_x - 1.5*R;
const double x_dense_right = center_x + 1.75*R;
const double y_dense_bottom = center_y - 1.25*R;
const double y_dense_top = center_y + 1.25*R;

Matrix CreateMeshFreeFlow() {
    std::cout << "Nejmenší krok by měl být " << std::to_string(delta_x_min) << std::endl;

    // Creation of vector of x coordinates
    std::vector<double> x_coordinates;

    double x_temp = x_dense_left;

    double delta_x = delta_x_min;
    while(x_temp < x_dense_right) {
        x_coordinates.push_back(x_temp);
        x_temp = x_temp + delta_x;
    }

    while(x_temp < L-delta_x/2) {
        x_coordinates.push_back(x_temp);
        delta_x = rate_of_step_increase*delta_x;
        x_temp = x_temp + delta_x;
    }

    x_coordinates.push_back(L);

    delta_x = rate_of_step_increase*delta_x_min;

    x_temp = x_dense_left - delta_x;
    
    while(x_temp > 0+delta_x/2){
        x_coordinates.insert(x_coordinates.begin(), x_temp);
        delta_x = rate_of_step_increase*delta_x;
        x_temp = x_temp - delta_x;
    }

    x_coordinates.insert(x_coordinates.begin(), 0);

    std::cout << "Počet bodů v x směru je: " << x_coordinates.size() << std::endl;

    // Creation of vector of y coordinates
    std::vector<double> y_coordinates;

    double y_temp = y_dense_bottom;

    double delta_y = delta_y_min;

    while(y_temp < y_dense_top) {
        y_coordinates.push_back(y_temp);
        y_temp = y_temp + delta_y;
    }

    while(y_temp < H-delta_y/2) {
        y_coordinates.push_back(y_temp);
        delta_y = rate_of_step_increase_vert*delta_y;
        y_temp = y_temp + delta_y;
    }

    y_coordinates.push_back(H);

    delta_y = rate_of_step_increase_vert*delta_y_min;

    y_temp = y_dense_bottom - delta_y;
    
    while(y_temp > 0+delta_y/2){
        y_coordinates.insert(y_coordinates.begin(), y_temp);
        delta_y =rate_of_step_increase_vert*delta_y;
        y_temp = y_temp - delta_y;
    }

    y_coordinates.insert(y_coordinates.begin(), 0);



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
    std::string filename = "MeshgridMaxDens"+std::to_string(max_node_density)+".txt";
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


    double closest_to_center_x;
    double closest_to_center_y;
    double diff_x = 100000;
    double diff_y = 100000;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = coordinates[i][j][0];
            double y = coordinates[i][j][1];

            double diff_x_now = std::abs(x-center_x);
            if(diff_x_now < diff_x){
                closest_to_center_x = x;
                diff_x = diff_x_now;
            }
           
            double diff_y_now = std::abs(y-center_y);
            if(diff_y_now < diff_y){
                closest_to_center_y = y;
                diff_y = diff_y_now;
            }
        }
    }


    double distance = 0.0;
    double fluid = 1.0;
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = coordinates[i][j][0];
            double y = coordinates[i][j][1];

            // Initializes the area uniformly
            double p = p_in;
            //Creates a parabolic profile at the inlet
            double u = 4*u_in*y*(0.41-y)/(0.41*0.41);
            double v = 0.;

            if(j == 0){
                u = 0;
            }
            if(j == M-1){
                u = 0;
            }


            fluid = 1.0;

            // Calculate the distance between the point (x, y) and the center (center_x, center_y)
            distance = std::sqrt((x - closest_to_center_x) * (x - closest_to_center_x) + (y - closest_to_center_y) * (y - closest_to_center_y));

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

    //Setting up the cylinder within the matrix
    std::vector<std::vector<double>> structure_matrix(N, std::vector<double>(M, 0.0));
    double fluid_this;
    double fluid_left;
    double fluid_right;
    double fluid_up;
    double fluid_down;

    for (int i = 1; i < N-1; ++i) {
        for (int j = 1; j < M-1; ++j) {
        
            fluid_this = initial_matrix[i][j][5];
            fluid_left = initial_matrix[i-1][j][5];
            fluid_right = initial_matrix[i+1][j][5];
            fluid_up = initial_matrix[i][j+1][5];
            fluid_down = initial_matrix[i][j-1][5];

            if(fluid_this+fluid_left+fluid_right+fluid_down+fluid_up < 1){
                structure_matrix[i][j] = std::nan("");
            }
            else{
                structure_matrix[i][j] = fluid_this;
            }

        }
    }

    //Borders are all fluid
    for (int i = 0; i < N; ++i) { 
        structure_matrix[i][0] = 1;
        structure_matrix[i][M-1] = 1;  
    }
    for (int j = 1; j < M-1; ++j) { 
        structure_matrix[0][j] = 1;
        structure_matrix[N-1][j] = 1;  
    }

    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            fluid_this = structure_matrix[i][j];
            if(std::isnan(fluid_this)){
                initial_matrix[i][j][0] = std::nan("");
                initial_matrix[i][j][1] = std::nan("");
                initial_matrix[i][j][2] = std::nan("");
            }
            if(fluid_this == 0){
                initial_matrix[i][j][1] = 0;
                initial_matrix[i][j][2] = 0;
            }
            initial_matrix[i][j][5] = structure_matrix[i][j];
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

void saveTimes(const std::vector<double>& times, int counter) {
    std::string filename = "times_channel_nondimensional_CPP_DR" + std::to_string(max_node_density) + "_iter" + std::to_string(counter) + "_CFL05_beta"+std::to_string(beta)+".txt";
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

    std::vector<std::vector<double>> residues;
    std::vector<double> times;
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

                //checks to see if point is in fluid or in wall
                a = initial_matrix[i][j][5];
                
                // Extract data for the current grid point and its neighbors
                //Skips to next point if this point is inside the wall
                if(std::isnan(a)){
                    continue;
                }
                

                if(a == 1) {
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

                
                continue;
                }
                
                // Enforces the wall boundary conditions when point is on the wall
                if(a == 0) {
                    predict[i][j][1] = 0;
                    predict[i][j][2] = 0;

                    double a_this = a;
                    double a_left = initial_matrix[i-1][j][5];
                    double a_right = initial_matrix[i+1][j][5];
                    double a_up = initial_matrix[i][j+1][5];
                    double a_down = initial_matrix[i][j-1][5];

                    double x_this = initial_matrix[i][j][3];
                    double y_this = initial_matrix[i][j][4];

                    double nx = (x_this - center_x)/std::sqrt((x_this - center_x)*(x_this - center_x) +  (y_this - center_y)*(y_this - center_y)); //cos theta
                    double ny = (y_this - center_y)/std::sqrt((x_this - center_x)*(x_this - center_x) +  (y_this - center_y)*(y_this - center_y)); //sin theta

                    double new_p = 0;

                    // 1st Quadrant
                    if (nx >= 0 && ny >= 0 && ny <= 1){
                        double p_right = predict[i+1][j][0];
                        if(a_right == 0){
                            p_right = predict[i+1][j+1][0];
                        }

                        double p_up = predict[i][j+1][0];
                        if(a_up == 0){
                            p_up = predict[i+1][j+1][0];
                        }

                        double delta_x_right = predict[i+1][j][3] - predict[i][j][3];
                        double delta_y_up = predict[i][j+1][4] - predict[i][j][4];

                        //Add calculation
                        new_p = 1/(nx/delta_x_right+ny/delta_y_up)*(nx/delta_x_right*p_right + ny/delta_y_up*p_up);
                    }

                    // 2nd Quadrant
                    if (-1 <= nx && nx <= 0 && 0 <= ny){
                        double p_left = predict[i-1][j][0];
                        if(a_left == 0){
                            p_left = predict[i-1][j+1][0];
                        }

                        double p_up = predict[i][j+1][0];
                        if(a_up == 0){
                            p_up = predict[i-1][j+1][0];
                        }
                       
                        //Add calculation
                        double delta_x_left = predict[i][j][3] - predict[i - 1][j][3];
                        double delta_y_up = predict[i][j+1][4] - predict[i][j][4];

                        new_p = 1/(nx/delta_x_left-ny/delta_y_up)*(nx/delta_x_left*p_left - ny/delta_y_up*p_up);
                
                    }

                    // 3rd Quadrant
                    if (nx <= 0 && -1 <= ny && ny <= 0){
                        double p_left = predict[i-1][j][0];
                        if(a_left == 0){
                            p_left = predict[i-1][j-1][0];
                        }
                   
                        double p_down = predict[i][j-1][0];
                        if(a_down == 0){
                            p_down = predict[i-1][j-1][0];
                        }


                        //Add calculation
                        double delta_x_left = predict[i][j][3] - predict[i-1][j][3];
                        double delta_y_down = predict[i][j][4] - predict[i][j-1][4];

                        //Add calculation
                        new_p = 1/(nx/delta_x_left+ny/delta_y_down)*(nx/delta_x_left*p_left + ny/delta_y_down*p_down);
                    }

                    // 4th Quadrant
                    if (0 <= nx && nx < 1 && ny < 0){
                        double p_right = predict[i+1][j][0];
                        if(a_right == 0){
                            p_right = predict[i+1][j-1][0];
                        }
                        
                        double p_down = predict[i][j-1][0];
                        if(a_down == 0){
                            p_down = predict[i+1][j-1][0];
                        }
                        

                        //Add calculation
                        double delta_x_right = predict[i+1][j][3] - predict[i][j][3];
                        double delta_y_down = predict[i][j][4] - predict[i][j-1][4];

                        //Add calculation
                        new_p = 1/(-nx/delta_x_right+ny/delta_y_down)*(-nx/delta_x_right*p_right + ny/delta_y_down*p_down);
                    }



                    predict[i][j][0] = new_p;
                    continue;
                }


            }
        }


        
         // Top and bottom wall - no slip zero velocity is passed on from the initial conditons (the scheme loop doesnt affect borders)
        for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the bottom (or top) two layers
            predict[k][0][0] = predict[k][1][0];          // Bottom boundary condition
            predict[k][M-1][0] = predict[k][M-2][0];      // Top  boundary condition
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

                //checks to see if point is in fluid or in wall
                a = initial_matrix[i][j][5];
                
                // Extract data for the current grid point and its neighbors
                //Skips to next point if this point is inside the wall
                if(std::isnan(a)){
                    continue;
                }
                

                if(a == 1) {
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
                
                continue;
                }


                // Enforces the wall boundary conditions when point is on the wall
                if(a == 0) {
                past[i][j][1] = 0;
                past[i][j][2] = 0;

                double a_this = a;
                double a_left = initial_matrix[i-1][j][5];
                double a_right = initial_matrix[i+1][j][5];
                double a_up = initial_matrix[i][j+1][5];
                double a_down = initial_matrix[i][j-1][5];

                double x_this = initial_matrix[i][j][3];
                double y_this = initial_matrix[i][j][4];

                double nx = (x_this - center_x)/std::sqrt((x_this - center_x)*(x_this - center_x) +  (y_this - center_y)*(y_this - center_y)); //cos theta
                double ny = (y_this - center_y)/std::sqrt((x_this - center_x)*(x_this - center_x) +  (y_this - center_y)*(y_this - center_y)); //sin theta

                double new_p = 0;



                // 1st Quadrant
                if (nx >= 0 && ny >= 0 && ny <= 1){
                    double p_right = past[i+1][j][0];
                    if(a_right == 0){
                        p_right = past[i+1][j+1][0];
                    }

                    double p_up = past[i][j+1][0];
                    if(a_up == 0){
                        p_up = past[i+1][j+1][0];
                    }

                    double delta_x_right = past[i+1][j][3] - past[i][j][3];
                    double delta_y_up = past[i][j+1][4] - past[i][j][4];
                    //std::cout << "1st Quadrant" << std::endl;

                    //Add calculation
                    new_p = 1/(nx/delta_x_right+ny/delta_y_up)*(nx/delta_x_right*p_right + ny/delta_y_up*p_up);
                }

                // 2nd Quadrant
                if (-1 <= nx && nx <= 0 && 0 <= ny){
                    //std::cout << "2nd Quadrant" << std::endl;
                    double p_left = past[i-1][j][0];
                    if(a_left == 0){
                        p_left = past[i-1][j+1][0];
                    }

                    double p_up = past[i][j+1][0];
                    if(a_up == 0){
                        p_up = past[i-1][j+1][0];
                    }
                    
                    //Add calculation
                    double delta_x_left = past[i][j][3] - past[i - 1][j][3];
                    double delta_y_up = past[i][j+1][4] - past[i][j][4];

                    new_p = 1/(nx/delta_x_left-ny/delta_y_up)*(nx/delta_x_left*p_left - ny/delta_y_up*p_up);
            
                }

                // 3rd Quadrant
                if (nx <= 0 && -1 <= ny && ny <= 0){
                    //std::cout << "3rd Quadrant" << std::endl;
                    double p_left = past[i-1][j][0];
                    if(a_left == 0){
                        p_left = past[i-1][j-1][0];
                    }
                
                    double p_down = past[i][j-1][0];
                    if(a_down == 0){
                        p_down = past[i-1][j-1][0];
                    }


                    //Add calculation
                    double delta_x_left = past[i][j][3] - past[i-1][j][3];
                    double delta_y_down = past[i][j][4] - past[i][j-1][4];

                    //Add calculation
                    new_p = 1/(nx/delta_x_left+ny/delta_y_down)*(nx/delta_x_left*p_left + ny/delta_y_down*p_down);
                }

                // 4th Quadrant
                if (0 <= nx && nx < 1 && ny < 0){
                    //std::cout << "4th Quadrant" << std::endl;
                    double p_right = past[i+1][j][0];
                    if(a_right == 0){
                        p_right = past[i+1][j-1][0];
                    }
                    
                    double p_down = past[i][j-1][0];
                    if(a_down == 0){
                        p_down = past[i+1][j-1][0];
                    }
                    

                    //Add calculation
                    double delta_x_right = past[i+1][j][3] - past[i][j][3];
                    double delta_y_down = past[i][j][4] - past[i][j-1][4];

                    //Add calculation
                    new_p = 1/(-nx/delta_x_right+ny/delta_y_down)*(-nx/delta_x_right*p_right + ny/delta_y_down*p_down);
                }

                double old_p_here = past[i][j][0];
                past[i][j][0] = new_p;

                sums[0] += (past[i][j][0] - old_p_here)*(past[i][j][0] - old_p_here);
                continue;
                }
                


            }
        }

            // Top and bottom wall - no slip zero velocity is passed on from the initial conditons (the scheme loop doesnt affect borders)
        for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the bottom (or top) two layers
            past[k][0][0] = past[k][1][0];          // Bottom boundary condition
            past[k][M-1][0] = past[k][M-2][0];      // Top  boundary condition
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
        times.push_back(time);
        if(counter % 10 == 0){
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

 


        if(counter % 5000 == 0){
            //saveMatrix(past, counter);
            saveMatrixWithTime(past, counter, time);
        }

        if(counter % 100000 == 0){
            saveResidues(residues, counter);
            saveTimes(times, counter);
        }
       



        counter++;
    }

    return 0;

}