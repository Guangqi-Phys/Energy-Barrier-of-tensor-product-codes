#include "../include/energy_barrier.hpp"
#include "../include/generate_codeword.hpp"
#include "../include/tensor_product.hpp"
#include "../include/energy_barrier_exhaust.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <climits>

// Helper function to convert string codeword to vector<int>
vector<int> stringToVector(const string& s) {
    vector<int> v;
    v.reserve(s.length());
    for(char c : s) {
        v.push_back(c - '0');
    }
    return v;
}

// Helper function to print matrix
void printMatrix(const vector<vector<int>>& matrix, const string& name) {
    cout << name << " (" << matrix.size() << " x " 
         << (matrix.empty() ? 0 : matrix[0].size()) << "):" << endl;
    for(const auto& row : matrix) {
        for(int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int main() {
    // Example parity-check matrices
    vector<vector<int>> H1 = {
        {1,1,0},
        {0,1,1},
        {1,0,1}
    };

    vector<vector<int>> H2 = {
        {1,1,0},
        {0,1,1},
        {1,0,1}
    };

    // 1. Print input matrices and their properties
    cout << "=== Analysis of H1 ===" << endl;
    printMatrix(H1, "H1");
    int d1 = computeMinimumDistance(H1);
    cout << "Minimum distance of H1: " << d1 << endl;
    
    // Compute energy barrier of H1
    vector<string> codewords1 = computeAllCodewordsGF2(H1);
    int minBarrier1 = INT_MAX;
    string minBarrierCodeword1;
    
    for(const auto& cw : codewords1) {
        if(cw.find('1') == string::npos) continue; // Skip zero codeword
        vector<int> codewordVec = stringToVector(cw);
        int barrier = computeEnergyBarrier(H1, codewordVec);
        if(barrier >= 0 && barrier < minBarrier1) {
            minBarrier1 = barrier;
            minBarrierCodeword1 = cw;
        }
    }
    cout << "Energy barrier of H1: " << minBarrier1 << endl;
    cout << "Achieved by codeword: " << minBarrierCodeword1 << "\n\n";

    cout << "=== Analysis of H2 ===" << endl;
    printMatrix(H2, "H2");
    int d2 = computeMinimumDistance(H2);
    cout << "Minimum distance of H2: " << d2 << endl;
    
    // Compute energy barrier of H2
    vector<string> codewords2 = computeAllCodewordsGF2(H2);
    int minBarrier2 = INT_MAX;
    string minBarrierCodeword2;
    
    for(const auto& cw : codewords2) {
        if(cw.find('1') == string::npos) continue; // Skip zero codeword
        vector<int> codewordVec = stringToVector(cw);
        int barrier = computeEnergyBarrier(H2, codewordVec);
        if(barrier >= 0 && barrier < minBarrier2) {
            minBarrier2 = barrier;
            minBarrierCodeword2 = cw;
        }
    }
    cout << "Energy barrier of H2: " << minBarrier2 << endl;
    cout << "Achieved by codeword: " << minBarrierCodeword2 << "\n\n";

    // Build and analyze tensor product code H3
    cout << "=== Analysis of Tensor Product Code H3 ===" << endl;
    vector<vector<int>> H3 = buildTensorProductParityCheck(H1, H2);
    printMatrix(H3, "H3");
    int d3 = computeMinimumDistance(H3);
    cout << "Minimum distance of H3: " << d3 << endl;
    
    // Compute energy barrier of H3
    vector<string> codewords3 = computeAllCodewordsGF2(H3);
    int minBarrier3 = INT_MAX;
    string minBarrierCodeword3;
    
    // First, get the total number of codewords to process
    int total_codewords = 0;
    for(const auto& cw : codewords3) {
        if(cw.find('1') != string::npos) { // Count non-zero codewords
            total_codewords++;
        }
    }

    // Initialize counter for processed codewords
    int processed = 0;

    // Process codewords with progress bar
    for(const auto& cw : codewords3) {
        if(cw.find('1') == string::npos) continue; // Skip zero codeword
        
        // Update and display progress
        processed++;
        float percentage = (100.0 * processed) / total_codewords;
        cout << "\rProgress: [";
        int pos = 50 * processed / total_codewords;
        for (int i = 0; i < 50; ++i) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(percentage) << "%" << flush;

        // Original processing code
        vector<int> codewordVec = stringToVector(cw);
        int barrier = computeEnergyBarrier(H3, codewordVec);
        if(barrier >= 0 && barrier < minBarrier3) {
            minBarrier3 = barrier;
            minBarrierCodeword3 = cw;
        }
    }
    cout << endl; // New line after progress bar is complete
    cout << "Energy barrier of H3: " << minBarrier3 << endl;
    cout << "Achieved by codeword: " << minBarrierCodeword3 << "\n\n";

    // Print summary
    cout << "=== Summary ===" << endl;
    cout << "H1: distance = " << d1 << ", energy barrier = " << minBarrier1 << endl;
    cout << "H2: distance = " << d2 << ", energy barrier = " << minBarrier2 << endl;
    cout << "H3 (tensor product): distance = " << d3 << ", energy barrier = " << minBarrier3 << endl;

    return 0;
}

