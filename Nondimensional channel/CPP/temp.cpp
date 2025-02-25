#include <iostream>
#include <fstream>
#include <vector>

// Function to save a 2D matrix to a text file
void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    // Open the file in write mode
    std::ofstream file(filename);

    // Check if the file is open
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    // Iterate through the matrix and write each element on a new line
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            file << element << " ";  // Write the element to the file, each on a new line
        } file << std::endl;
    }

    // Close the file
    file.close();
    std::cout << "Matrix saved to " << filename << std::endl;
}

int main() {
    // Example matrix (2D vector)
    std::vector<std::vector<double>> residues = {
        {1.1, 2.2, 3.3},
        {4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9}
    };

    // Save the matrix to a file named "matrix.txt"
    saveMatrixToFile(residues, "matrix.txt");

    return 0;
}
