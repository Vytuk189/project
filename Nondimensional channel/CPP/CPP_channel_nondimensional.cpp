#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Flow and channel conditions
const double rho = 10.0; // Density
const double L = 2.0; // Length of channel
const double H = 1.0; // Height of channel
const double Re = 10.0;

const double nu = 1.0 / Re;

const double beta = 1.0;

const int N = 200; // Number of points in x direction
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
            double v = 0.0;
            
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

void saveResidues(const std::vector<double>& residues) {
    std::ofstream file("residues.txt");
    for (const auto& residue : residues) {
        file << residue << std::endl;
    }
}

void savePast(const Matrix& past) {
    std::ofstream file("past_array.txt");
    for (const auto& row : past) {
        for (const auto& col : row) {
            for (double val : col) {
                file << val << " ";
            }
            file << std::endl;
        }
    }
}

int main() {
    Matrix initial_matrix = createInitialMatrix();
    Matrix past = initial_matrix;
    Matrix predict = past;
    Matrix deltasPredict = Matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(3, 0.0)));
    Matrix deltasCorrect = Matrix(N + 3, std::vector<std::vector<double>>(M + 1, std::vector<double>(3, 0.0)));
    
    std::vector<double> residues;
    int counter = 0;
    
    while (counter < 400000) {
        if (counter % 100 == 0) {
            double umax_now = 0.0;
            double vmax_now = 0.0;
            double p_max = -1e10;
            double p_min = 1e10;
            
            // Extract u_speeds, v_speeds, pressures to find max values
            for (int i = 0; i < N + 3; ++i) {
                for (int j = 0; j < M + 1; ++j) {
                    double u = past[i][j][1];
                    double v = past[i][j][2];
                    double p = past[i][j][0];
                    umax_now = std::max(umax_now, u);
                    vmax_now = std::max(vmax_now, v);
                    p_max = std::max(p_max, p);
                    p_min = std::min(p_min, p);
                }
            }
            
            double tau = 0.1 * 1.0 / (((umax_now + std::sqrt(umax_now * umax_now + beta * beta)) / h) +
                                     ((vmax_now + std::sqrt(vmax_now * vmax_now + beta * beta)) / h) + 
                                     2.0 * nu * (1.0 / (h * h) + 1.0 / (h * h)));
        }
        
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
                
                // Velocity and pressure differences
                std::vector<double> F_this = {u_this, u_this * u_this + p_this / rho, u_this * v_this};
                std::vector<double> F_left = {u_left, u_left * u_left + p_left / rho, u_left * v_left};
                
                std::vector<double> G_this = {v_this, v_this * u_this, v_this * v_this + p_this / rho};
                std::vector<double> G_down = {v_down, v_down * u_down, v_down * v_down + p_down / rho};
                
                // Calculating the delta values
                std::vector<double> deltaWpredict(3, 0.0);
                // Do the matrix multiplication and other computations here for deltaWpredict
                
                // Update the predictor matrix
                for (int k = 0; k < 3; ++k) {
                    predict[i][j][k] = past[i][j][k] + tau * deltaWpredict[k];
                }
            }
        }
        
        // Corrector section
        for (int i = 1; i < N + 2; ++i) {
            for (int j = 1; j < M; ++j) {
                // Corrector section similar to the predictor
                // Calculate deltasCorrect here, similar to deltasPredict
            }
        }
        
        // Compute residue and save values
        double residue = 0.0;
        for (int i = 1; i < N + 2; ++i) {
            for (int j = 1; j < M; ++j) {
                std::vector<double> deltaW = {0.0, 0.0, 0.0}; // Replace with actual deltaW calculation
                residue += deltaW[0] * deltaW[0] + deltaW[1] * deltaW[1] + deltaW[2] * deltaW[2];
            }
        }
        
        residue = std::sqrt(residue / (tau * N * M));
        residues.push_back(residue);
        
        std::cout << "***************************************************************\n";
        std::cout << "Iteration: " << counter << "\n";
        std::cout << "Residue: " << residue << "\n";
        counter++;
    }
    
    // Save the results
    saveResidues(residues);
    savePast(past);
    
    return 0;
}
