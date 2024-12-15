#include "../include/energy_barrier.hpp"
#include "../include/generate_codeword.hpp"
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

int main() {
    // Example parity-check matrix H
vector<vector<int>> H = {
    {1,0,0,1,0,0,0,0,0},
    {0,1,0,0,1,0,0,0,0},
    {0,0,1,0,0,1,0,0,0},
    {0,0,0,1,0,0,1,0,0},
    {0,0,0,0,1,0,0,1,0},
    {0,0,0,0,0,1,0,0,1},
    {1,0,0,0,0,0,1,0,0},
    {0,1,0,0,0,0,0,1,0},
    {0,0,1,0,0,0,0,0,1},
    {1,1,0,0,0,0,0,0,0},
    {0,1,1,0,0,0,0,0,0},
    {1,0,1,0,0,0,0,0,0},
    {0,0,0,1,1,0,0,0,0},
    {0,0,0,0,1,1,0,0,0},
    {0,0,0,1,0,1,0,0,0},
    {0,0,0,0,0,0,1,1,0},
    {0,0,0,0,0,0,0,1,1},
    {0,0,0,0,0,0,1,0,1}
};

    // 1. Get all codewords using generate_codeword functions
    vector<string> codewords = computeAllCodewordsGF2(H);
    
    cout << "Found " << codewords.size() << " codewords:\n";
    for(const auto& cw : codewords) {
        cout << cw << "\n";
    }
    cout << endl;

    // 2. Compute energy barrier for each non-zero codeword
    int minBarrier = INT_MAX;
    string minBarrierCodeword;

    for(const auto& cw : codewords) {
        // Skip the zero codeword (all zeros)
        if(cw.find('1') == string::npos) continue;
        
        // Convert string codeword to vector<int>
        vector<int> codewordVec = stringToVector(cw);
        
        // Compute energy barrier for this codeword
        int barrier = computeEnergyBarrierExhaustive(H, codewordVec);
        
        cout << "Energy barrier for codeword " << cw << ": " << barrier << endl;
        
        if(barrier >= 0 && barrier < minBarrier) {
            minBarrier = barrier;
            minBarrierCodeword = cw;
        }
    }

    // 3. Output the minimum energy barrier
    if(minBarrier == INT_MAX) {
        cout << "No valid energy barriers found." << endl;
    } else {
        cout << "\nMinimum energy barrier of the code: " << minBarrier << endl;
        cout << "Achieved by codeword: " << minBarrierCodeword << endl;
    }

    return 0;
}
