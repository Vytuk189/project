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


const double cfl = 0.8;
const int iter = 800000; // Number of iterations
const int N = 100; // Number of points in x direction
const int M = static_cast<int>(H / (L / N)); // Number of points in y direction
const double h = L / N; // Space step

const double u_in = 1;
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





    //Setting up walls within the matrix
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

    
    //Left and right boundaries (without corners)
    for (int j = 1; j<M-1; ++j){

        int i = 0;

        fluid_this = initial_matrix[i][j][5];
        fluid_right = initial_matrix[i+1][j][5];
        fluid_up = initial_matrix[i][j+1][5];
        fluid_down = initial_matrix[i][j-1][5];

        if(fluid_this+fluid_right+fluid_down+fluid_up < 1){
            structure_matrix[i][j] = std::nan("");
        }
        else{
            structure_matrix[i][j] = fluid_this;
        }

        i = N-1;

        fluid_this = initial_matrix[i][j][5];
        fluid_left = initial_matrix[i-1][j][5];
        fluid_up = initial_matrix[i][j+1][5];
        fluid_down = initial_matrix[i][j-1][5];

        if(fluid_this+fluid_left+fluid_down+fluid_up < 1){
            structure_matrix[i][j] = std::nan("");
        }
        else{
            structure_matrix[i][j] = fluid_this;
        }
    }


    // Top and bottom boundaries (without corners)
    for (int i = 1; i<N-1; ++i){

        int j = 0;

        fluid_this = initial_matrix[i][j][5];
        fluid_left = initial_matrix[i-1][j][5];
        fluid_right = initial_matrix[i+1][j][5];
        fluid_up = initial_matrix[i][j+1][5];

        if(fluid_this+fluid_right+fluid_left+fluid_up < 1){
            structure_matrix[i][j] = std::nan("");
        }
        else{
            structure_matrix[i][j] = fluid_this;
        }

        j = M-1;

        fluid_this = initial_matrix[i][j][5];
        fluid_left = initial_matrix[i-1][j][5];
        fluid_right = initial_matrix[i+1][j][5];
        fluid_down = initial_matrix[i][j-1][5];

        if(fluid_this+fluid_left+fluid_down+fluid_right < 1){
            structure_matrix[i][j] = std::nan("");
        }
        else{
            structure_matrix[i][j] = fluid_this;
        }
    }



    //Bottom left corner
    fluid_this = initial_matrix[0][0][5];
    fluid_right = initial_matrix[1][0][5];
    fluid_up = initial_matrix[0][1][5];

    if(fluid_this+fluid_right+fluid_up < 1){
        structure_matrix[0][0] = std::nan("");
    }
    else{
        structure_matrix[0][0] = fluid_this;
    }

    //Bottom right corner
    fluid_this = initial_matrix[N-1][0][5];
    fluid_left = initial_matrix[N-2][0][5];
    fluid_up = initial_matrix[N-1][1][5];

    if(fluid_this+fluid_left+fluid_up < 1){
        structure_matrix[N-1][0] = std::nan("");
    }
    else{
        structure_matrix[N-1][0] = fluid_this;
    }

    //Top right corner
    fluid_this = initial_matrix[N-1][M-1][5];
    fluid_left = initial_matrix[N-2][M-1][5];
    fluid_down = initial_matrix[N-1][M-2][5];

    if(fluid_this+fluid_left+fluid_down < 1){
        structure_matrix[N-1][M-1] = std::nan("");
    }
    else{
        structure_matrix[N-1][M-1] = fluid_this;
    }

    //Top left corner
    fluid_this = initial_matrix[0][M-1][5];
    fluid_right = initial_matrix[1][M-1][5];
    fluid_down = initial_matrix[0][M-2][5];

    if(fluid_this+fluid_right+fluid_down < 1){
        structure_matrix[0][M-1] = std::nan("");
    }
    else{
        structure_matrix[0][M-1] = fluid_this;
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
    saveMatrixWithTime(initial_matrix, 0, 0);
    
    Matrix past = initial_matrix;
    Matrix predict = past;
    Matrix deltasPredict = Matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(3, 0.0)));
    Matrix deltasCorrect = Matrix(N, std::vector<std::vector<double>>(M, std::vector<double>(3, 0.0)));

    std::vector<std::vector<double>> residues;


    int counter = 0;
    double time = 0;
    double a = 1; //Will store whether a point is in fluid or in structure
    double cfl = 0.5;

    double umax_now = 0.0;
    double vmax_now = 0.0;
    double p_max = -1e10;
    double p_min = 1e10;
    while (counter < 800000) {

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
                                    2.0 * nu * (1.0 / (h*h) + 1.0 / (h *h)));

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
                // Enforces the wall boundary conditions when point is on the wall
                if(a == 0) {
                    predict[i][j][1] = 0;
                    predict[i][j][2] = 0;

                    double a_this = a;
                    double a_left = initial_matrix[i-1][j][5];
                    double a_right = initial_matrix[i+1][j][5];
                    double a_up = initial_matrix[i][j+1][5];
                    double a_down = initial_matrix[i][j-1][5];

                    double sum_of_pressures = 0;
                    int count_of_p = 0;
                    if(a_left == 1){
                        sum_of_pressures = sum_of_pressures + past[i-1][j][0];
                        count_of_p++;
                    }
                    if(a_right == 1){
                        sum_of_pressures = sum_of_pressures + past[i+1][j][0];
                        count_of_p++;
                    }
                    if(a_up == 1){
                        sum_of_pressures = sum_of_pressures + past[i][j+1][0];
                        count_of_p++;
                    }
                    if(a_down == 1){
                        sum_of_pressures = sum_of_pressures + past[i][j-1][0];
                        count_of_p++;
                    }

                    predict[i][j][0] = sum_of_pressures/count_of_p;
                    continue;
                }


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

                
                predict[i][j][0] = past[i][j][0] + tau*deltaWpredict[0];
                predict[i][j][1] = past[i][j][1] + tau*deltaWpredict[1];
                predict[i][j][2] = past[i][j][2] + tau*deltaWpredict[2];
            }
        }
        
         // Top and bottom wall - Neumann Zero Conditions for pressure
        for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the top two layers
            predict[k][0][0] = predict[k][1][0];      // Bottom  boundary condition
            predict[k][M-1][0] = predict[k][M-2][0];      // Top  boundary condition
            // Dirichlet zero for velocity from initial conditions
        }

        // Inlet and outlet
        for (int j = 0; j < M; ++j) {

            //Dirichlet for speed on the inlet is enforced from initial conditions

            // Neuman condition for speed at the outlet
            //U speed/V speed dont change between the last layer and the layer before it   
            predict[N-1][j][1] = predict[N-2][j][1];          
            predict[N-1][j][2] = predict[N-2][j][2]; 
            
            // Neumann condition for pressure on the inlet
            // Pressure doesnt change between the first and second layer
            predict[0][j][0] = predict[1][j][0];

            //Dirichlet zero condition for pressure at the outlet is enforced from initial conditions

        }







        // Corrector section
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < M - 1; ++j) {
                // Corrector section similar to the predictor
                // Calculate deltasCorrect here, similar to deltasPredict
                // Extract data for the current grid point and its neighbors

                //checks to see if point is in fluid or in wall
                a = initial_matrix[i][j][5];
                
                // Extract data for the current grid point and its neighbors
                //Skips to next point if this point is inside the wall
                if(std::isnan(a)){
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

                    double sum_of_pressures = 0;
                    int count_of_p = 0;
                    if(a_left == 1){
                        sum_of_pressures = sum_of_pressures + predict[i-1][j][0];
                        count_of_p++;
                    }
                    if(a_right == 1){
                        sum_of_pressures = sum_of_pressures + predict[i+1][j][0];
                        count_of_p++;
                    }
                    if(a_up == 1){
                        sum_of_pressures = sum_of_pressures + predict[i][j+1][0];
                        count_of_p++;
                    }
                    if(a_down == 1){
                        sum_of_pressures = sum_of_pressures + predict[i][j-1][0];
                        count_of_p++;
                    }

                    double old_p_here = past[i][j][0];
                    past[i][j][0] = sum_of_pressures/count_of_p;

                    sums[0] += (past[i][j][0] - old_p_here)*(past[i][j][0] - old_p_here);

                    continue;
                }



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

                past[i][j][0] = new_p;
                past[i][j][1] = new_u;
                past[i][j][2] = new_v;


                sums[0] += (new_p - old_p)*(new_p - old_p);
                sums[1] += (new_u - old_u)*(new_u - old_u);
                sums[2] += (new_v - old_v)*(new_v - old_v);


            }
        }


        // Top wall - Neumann Zero Conditions for pressure
        for (int k = 0; k < N; ++k) {
            // Pressure doesnt change between the top two layers
            past[k][0][0] = past[k][1][0];      // Bottom  boundary condition
            past[k][M-1][0] = past[k][M-2][0];      // Top  boundary condition
            // Dirichlet zero for velocity from initial conditions
        }

        // Inlet and outlet
        for (int j = 0; j < M; ++j) {

            //Dirichlet for speed on the inlet is enforced from initial conditions

            // Neuman condition for speed at the outlet
            //U speed/V speed dont change between the last layer and the layer before it   
            past[N-1][j][1] = past[N-2][j][1];          
            past[N-1][j][2] = past[N-2][j][2]; 
            
            // Neumann condition for pressure on the inlet
            // Pressure doesnt change between the first and second layer
            past[0][j][0] = past[1][j][0];

            //Dirichlet zero condition for pressure at the outlet is enforced from initial conditions

        }

        

        
        double residual_p = std::sqrt(sums[0])/(N*M);
        double residual_u = std::sqrt(sums[1])/(N*M);
        double residual_v = std::sqrt(sums[2])/(N*M);

        residues.push_back({residual_p, residual_u, residual_v});
        
        time = time + tau;
        if(counter % 10 == 0){
        std::cout << "***************************************************************\n";
        std::cout << "Case Custom BC: N=" << N << ", M=" << M << ", L=" << L << ", H=" << H << ", Re=" << Re << ", nu=" << nu <<  "\n";
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