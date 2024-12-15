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

    // 1. Print input matrices
    printMatrix(H1, "H1");
    printMatrix(H2, "H2");

    // 2. Build tensor product parity-check matrix H3
    vector<vector<int>> H3 = buildTensorProductParityCheck(H1, H2);
    printMatrix(H3, "H3 (Tensor Product)");

    // 3. Get all codewords of H3
    vector<string> codewords = computeAllCodewordsGF2(H3);
    
    cout << "Found " << codewords.size() << " codewords for H3:\n";
    for(const auto& cw : codewords) {
        cout << cw << "\n";
    }
    cout << endl;

    // 4. Compute energy barrier for each non-zero codeword
    int minBarrier = INT_MAX;
    string minBarrierCodeword;

    for(const auto& cw : codewords) {
        // Skip the zero codeword (all zeros)
        if(cw.find('1') == string::npos) continue;
        
        // Convert string codeword to vector<int>
        vector<int> codewordVec = stringToVector(cw);
        
        // Compute energy barrier for this codeword
        int barrier = computeEnergyBarrier(H3, codewordVec);
        
        cout << "Energy barrier for codeword " << cw << ": " << barrier << endl;
        
        if(barrier >= 0 && barrier < minBarrier) {
            minBarrier = barrier;
            minBarrierCodeword = cw;
        }
    }

    // 5. Output the minimum energy barrier
    if(minBarrier == INT_MAX) {
        cout << "No valid energy barriers found." << endl;
    } else {
        cout << "\nMinimum energy barrier of the tensor product code: " << minBarrier << endl;
        cout << "Achieved by codeword: " << minBarrierCodeword << endl;
    }

    return 0;
}
