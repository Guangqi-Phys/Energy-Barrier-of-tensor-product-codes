#include "../include/energy_barrier.hpp"
#include "../include/generate_codeword.hpp"
#include "../include/tensor_product.hpp"
#include "../include/energy_barrier_exhaust.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "/opt/homebrew/opt/libomp/include/omp.h"


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
                        vector<vector<int>>& H3,
                        int& d1, int& E1, int& d2, int& E2, int& E3,
                        vector<int>& codewords1,
                        vector<int>& codewords2,
                        vector<int>& codewords3) {
    try {
        cout << "\nDebug: Starting new simulation with dimensions: " 
             << "m1=" << m1 << ", n1=" << n1 
             << ", m2=" << m2 << ", n2=" << n2 << ", w=" << w << endl;

        // Generate random parity check matrices
        try {
            cout << "Debug: Generating random matrices..." << endl;
            H1 = generateRandomParityCheckMatrix(m1, n1, w);
            H2 = generateRandomParityCheckMatrix(m2, n2, w);
        } catch (const exception& e) {
            cout << "Error in generating matrices: " << e.what() << endl;
            return false;
        }

        // Verify matrix constraints
        try {
            cout << "Debug: Verifying matrix constraints..." << endl;
            if (!verifyMatrixConstraints(H1, w) || !verifyMatrixConstraints(H2, w)) {
                cout << "Matrix constraints verification failed" << endl;
                return false;
            }
        } catch (const exception& e) {
            cout << "Error in verifying matrix constraints: " << e.what() << endl;
            return false;
        }

        // Find codewords
        try {
            cout << "Debug: Finding codewords..." << endl;
            codewords1 = findSingleCodeword(H1);
            codewords2 = findSingleCodeword(H2);
            
            if (codewords1.empty() || codewords2.empty()) {
                cout << "Failed to find valid codewords" << endl;
                return false;
            }
        } catch (const exception& e) {
            cout << "Error in finding codewords: " << e.what() << endl;
            return false;
        }

        // Compute Hamming distances
        try {
            cout << "Debug: Computing Hamming distances of codewords..." << endl;
            d1 = 0;
            for (int bit : codewords1) {
                if (bit == 1) d1++;
            }
            d2 = 0;
            for (int bit : codewords2) {
                if (bit == 1) d2++;
            }
            
            if (d1 <= 0 || d2 <= 0) {
                cout << "Invalid distances found" << endl;
                return false;
            }
        } catch (const exception& e) {
            cout << "Error in computing distances: " << e.what() << endl;
            return false;
        }

        // Compute energy barriers
        try {
            cout << "Debug: Computing energy barriers..." << endl;
            E1 = computeEnergyBarrier(H1, codewords1);
            E2 = computeEnergyBarrier(H2, codewords2);
            
            if (E1 < 0 || E2 < 0) {
                cout << "Invalid energy barriers found" << endl;
                return false;
            }
        } catch (const exception& e) {
            cout << "Error in computing energy barriers: " << e.what() << endl;
            return false;
        }

        // Build tensor product
        try {
            cout << "Debug: Building tensor product..." << endl;
            H3 = buildTensorProductParityCheck(H1, H2);
            codewords3 = buildTensorProductCodeword(codewords1, codewords2);
            
            if (codewords3.empty()) {
                cout << "Failed to build tensor product codeword" << endl;
                return false;
            }
        } catch (const exception& e) {
            cout << "Error in building tensor product: " << e.what() << endl;
            return false;
        }

        // Compute tensor product energy barrier
        try {
            cout << "Debug: Computing tensor product energy barrier..." << endl;
            E3 = computeEnergyBarrier(H3, codewords3);
        } catch (const exception& e) {
            cout << "Error in computing tensor product energy barrier: " << e.what() << endl;
            return false;
        }

        return (E3 >= 0);
    } catch (const exception& e) {
        cout << "Unexpected error in simulation: " << e.what() << endl;
        return false;
    }
}

int main() {
    const int nn = 100;  // Number of simulations
    bool foundCounterexample = false;
    int successCount = 0;  // Add this counter

    cout << "Starting simulation with " << nn << " iterations...\n";

    // Use thread-local random generators
    random_device rd;
    
    // Create a mutex for synchronized console output
    #pragma omp declare reduction(|| : bool : omp_out = omp_out || omp_in)

    #pragma omp parallel reduction(+:successCount)  // Add reduction for successCount
    {
        mt19937 local_gen(rd() + omp_get_thread_num()); // Different seed for each thread
        uniform_int_distribution<> dim_dist(4, 7);

        #pragma omp for reduction(||:foundCounterexample)
        for (int iter = 0; iter < nn; ++iter) {
            // Synchronize progress bar updates
            #pragma omp critical
            {
                showProgress(iter, nn);
            }

            int m1 = dim_dist(local_gen);
            int n1 = dim_dist(local_gen);
            int m2 = dim_dist(local_gen);
            int n2 = dim_dist(local_gen);
            const int w = 3;

            vector<vector<int>> H1, H2, H3;
            vector<int> codewords1, codewords2, codewords3;
            int d1, E1, d2, E2, E3;

            auto start = chrono::steady_clock::now();
            bool success = false;
            
            try {
                success = runSingleSimulation(m1, n1, m2, n2, w, H1, H2, H3, 
                                            d1, E1, d2, E2, E3,
                                            codewords1, codewords2, codewords3);
                
                auto current = chrono::steady_clock::now();
                if (chrono::duration_cast<chrono::seconds>(current - start).count() > 5) {
                    #pragma omp critical
                    {
                        cout << "\nIteration " << iter << " timed out, skipping...\n";
                    }
                    continue;
                }
            } catch (const exception& e) {
                #pragma omp critical
                {
                    cout << "\nError in iteration " << iter << ": " << e.what() << ", skipping...\n";
                }
                continue;
            }

            if (success) {
                successCount++;  // Increment on successful simulation
                int min_bound = min(d1 * E2, E1 * d2);
                
                if (E3 < min_bound - 2) {
                    foundCounterexample = true;
                    #pragma omp critical
                    {
                        cout << "\nFound counterexample in iteration " << iter + 1 << ":\n";
                        cout << "H1: " << m1 << "x" << n1 << " matrix, d1=" << d1 << ", E1=" << E1 << "\n";
                        // Print H1
                        cout << "H1 matrix:\n";
                        for (const auto& row : H1) {
                            for (int val : row) {
                                cout << val << " ";
                            }
                            cout << "\n";
                        }
                        cout << "Codeword 1: ";
                        for (int val : codewords1) {
                            cout << val << " ";
                        }
                        cout << "\n\n";

                        cout << "H2: " << m2 << "x" << n2 << " matrix, d2=" << d2 << ", E2=" << E2 << "\n";
                        // Print H2
                        cout << "H2 matrix:\n";
                        for (const auto& row : H2) {
                            for (int val : row) {
                                cout << val << " ";
                            }
                            cout << "\n";
                        }
                        cout << "Codeword 2: ";
                        for (int val : codewords2) {
                            cout << val << " ";
                        }
                        cout << "\n\n";

                        cout << "H3 (tensor product): E3=" << E3 << "\n";
                        // Print H3
                        cout << "H3 matrix:\n";
                        for (const auto& row : H3) {
                            for (int val : row) {
                                cout << val << " ";
                            }
                            cout << "\n";
                        }
                        cout << "Codeword 3: ";
                        for (int val : codewords3) {
                            cout << val << " ";
                        }
                        cout << "\n\n";

                        cout << "min(d1*E2, E1*d2)=" << min_bound << "\n";
                    }
                    #pragma omp cancel for
                }
            }
        }
    }

    showProgress(nn, nn);
    cout << "\n\nSuccessful simulations: " << successCount << "/" << nn << "\n";

    if (!foundCounterexample) {
        cout << "For all simulations, E3 >= min(d1*E2, E1*d2) was satisfied.\n";
    }

    return 0;
}

