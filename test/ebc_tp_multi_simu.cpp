#include "../include/energy_barrier.hpp"
#include "../include/generate_codeword.hpp"
#include "../include/tensor_product.hpp"
#include "../include/energy_barrier_exhaust.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>


using namespace std;

// Function to display progress bar
void showProgress(int current, int total) {
    const int barWidth = 50;
    float progress = (float)current / total;
    int pos = barWidth * progress;

    cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")" << flush;
}

// Function to run one simulation
bool runSingleSimulation(int m1, int n1, int m2, int n2, int w,
                        vector<vector<int>>& H1, vector<vector<int>>& H2,
                        int& d1, int& E1, int& d2, int& E2, int& E3) {
    
    // Generate random parity check matrices
    H1 = generateRandomParityCheckMatrix(m1, n1, w);
    H2 = generateRandomParityCheckMatrix(m2, n2, w);

    // Verify matrix constraints
    if (!verifyMatrixConstraints(H1, w) || !verifyMatrixConstraints(H2, w)) {
        return false;
    }

    // Compute distances
    d1 = computeMinimumDistance(H1);
    d2 = computeMinimumDistance(H2);
    if (d1 <= 0 || d2 <= 0) return false;  // Invalid codes

    // Get a non-zero codeword for each code
    vector<int> codewords1 = findSingleCodeword(H1);
    vector<int> codewords2 = findSingleCodeword(H2);
    
    // Convert first non-zero codeword to vector<int>
    // vector<int> c1, c2;
    // for (const string& cw : codewords1) {
    //     if (hammingWeight(cw) > 0) {
    //         c1.resize(cw.length());
    //         for (size_t i = 0; i < cw.length(); ++i) c1[i] = cw[i] - '0';
    //         break;
    //     }
    // }
    // for (const string& cw : codewords2) {
    //     if (hammingWeight(cw) > 0) {
    //         c2.resize(cw.length());
    //         for (size_t i = 0; i < cw.length(); ++i) c2[i] = cw[i] - '0';
    //         break;
    //     }
    // }
    // if (c1.empty() || c2.empty()) return false;

    // Compute energy barriers
    E1 = computeEnergyBarrier(H1, codewords1);
    E2 = computeEnergyBarrier(H2, codewords2);
    if (E1 < 0 || E2 < 0) return false;  // Invalid barriers

    // Compute tensor product
    vector<vector<int>> H3 = buildTensorProductParityCheck(H1, H2);
    
    // Get a non-zero codeword for H3
    vector<int> codewords3 = buildTensorProductCodeword(codewords1, codewords2);
    // vector<int> c3;
    // for (const string& cw : codewords3) {
    //     if (hammingWeight(cw) > 0) {
    //         c3.resize(cw.length());
    //         for (size_t i = 0; i < cw.length(); ++i) c3[i] = cw[i] - '0';
    //         break;
    //     }
    // }
    if (codewords3.empty()) return false;

    // Compute energy barrier for tensor product code
    E3 = computeEnergyBarrier(H3, codewords3);
    return (E3 >= 0);
}

int main() {
    const int nn = 100;  // Number of simulations
    bool foundCounterexample = false;

    cout << "Starting simulation with " << nn << " iterations...\n";

    for (int iter = 0; iter < nn; ++iter) {
        showProgress(iter, nn);

        // Random dimensions in range [6,9]
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dim_dist(5, 6);
        
        int m1 = dim_dist(gen);
        int n1 = dim_dist(gen);
        int m2 = dim_dist(gen);
        int n2 = dim_dist(gen);
        const int w = 3;

        vector<vector<int>> H1, H2;
        int d1, E1, d2, E2, E3;

        if (runSingleSimulation(m1, n1, m2, n2, w, H1, H2, d1, E1, d2, E2, E3)) {
            int min_bound = min(d1 * E2, E1 * d2);
            
            if (E3 < min_bound) {
                foundCounterexample = true;
                cout << "\nFound counterexample in iteration " << iter + 1 << ":\n";
                cout << "H1: " << m1 << "x" << n1 << " matrix, d1=" << d1 << ", E1=" << E1 << "\n";
                cout << "H2: " << m2 << "x" << n2 << " matrix, d2=" << d2 << ", E2=" << E2 << "\n";
                cout << "H3 (tensor product): E3=" << E3 << "\n";
                cout << "min(d1*E2, E1*d2)=" << min_bound << "\n";
                break;
            }
        }
    }

    showProgress(nn, nn);  // Complete the progress bar
    cout << "\n\n";

    if (!foundCounterexample) {
        cout << "For all simulations, E3 >= min(d1*E2, E1*d2) was satisfied.\n";
    }

    return 0;
}

