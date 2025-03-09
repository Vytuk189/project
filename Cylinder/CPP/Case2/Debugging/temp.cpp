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

int main() {
    // Initialize a sample vector for testing
    std::vector<double> vec1 = {1.0, 2.0, 3.0};
    std::vector<double> vec2 = {4.0, 5.0, 6.0};
    
    // Scale the vector by a factor of 2
    double scale = 2.0;
    std::vector<double> scaledVec = scaleVector(scale, vec1);
    
    std::cout << "Scaled vector (scale by " << scale << "): ";
    for (double value : scaledVec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Add the vectors
    std::vector<double> addedVec = addVectors(vec1, vec2);
    
    std::cout << "Added vector: ";
    for (double value : addedVec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Subtract the vectors
    std::vector<double> subtractedVec = subtractVectors(vec1, vec2);
    
    std::cout << "Subtracted vector: ";
    for (double value : subtractedVec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Initialize a sample matrix for testing
    std::vector<std::vector<double>> matrix = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };

    // Multiply the matrix with the vector
    std::vector<double> multipliedVec = multiplyMatrixWithVector(matrix, vec1);
    
    std::cout << "Matrix-vector multiplication result: ";
    for (double value : multipliedVec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
