#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string>
using namespace std;

// Helper function to count 1-bits in an integer (mod 2)
inline int xorBit(int a, int b) {
    return (a + b) & 1;
}

// Convert a binary vector to a string like "0101"
string vectorToString(const vector<int>& vec) {
    string s;
    s.reserve(vec.size());
    for (int bit : vec) {
        s.push_back(bit ? '1' : '0');
    }
    return s;
}

// Gaussian Elimination (over GF(2)) to find the RREF of H
// Returns a tuple: (H in RREF, pivotCols, rank)
tuple<vector<vector<int>>, vector<int>, int> 
gaussianEliminationGF2(const vector<vector<int>>& H_in) {
    // Make a local copy so we don't destroy the input
    vector<vector<int>> H(H_in);
    int rows = (int)H.size();
    int cols = (rows == 0) ? 0 : (int)H[0].size();

    int pivotRow = 0;
    vector<int> pivotCols;  // store pivot column indices

    for (int col = 0; col < cols && pivotRow < rows; col++) {
        // Find pivot (row with a '1' in the current column)
        int pivotCandidate = -1;
        for (int r = pivotRow; r < rows; r++) {
            if (H[r][col] == 1) {
                pivotCandidate = r;
                break;
            }
        }
        if (pivotCandidate == -1) {
            continue; // no pivot in this column
        }

        // Swap pivotCandidate row into pivotRow
        if (pivotCandidate != pivotRow) {
            swap(H[pivotRow], H[pivotCandidate]);
        }
        pivotCols.push_back(col);

        // Eliminate below pivotRow
        for (int r = pivotRow + 1; r < rows; r++) {
            if (H[r][col] == 1) {
                // XOR row pivotRow into row r
                for (int c = col; c < cols; c++) {
                    H[r][c] ^= H[pivotRow][c];
                }
            }
        }
        pivotRow++;
    }
    int rank = pivotRow;

    // Convert to Reduced Row Echelon Form (RREF) by eliminating above pivots
    for (int i = rank - 1; i >= 0; i--) {
        int pivotCol = pivotCols[i];
        // pivot is at row i, col pivotCol
        // Eliminate above
        for (int r = i - 1; r >= 0; r--) {
            if (H[r][pivotCol] == 1) {
                for (int c = pivotCol; c < cols; c++) {
                    H[r][c] ^= H[i][c];
                }
            }
        }
    }

    return make_tuple(H, pivotCols, rank);
}

/*
 * ComputeAllCodewordsGF2(H):
 *   Input: H is an ℓ x n parity-check matrix over GF(2).
 *   Output: A vector of binary strings, each representing one codeword in the null space of H.
 */
vector<string> computeAllCodewordsGF2(const vector<vector<int>>& H) {
    if (H.empty()) {
        // Edge case: no parity checks -> all 2^n strings are valid
        // But we don't know n from H. Return empty or handle differently.
        return {"0"}; // minimal placeholder
    }

//    int rows = (int)H.size();
    int cols = (int)H[0].size();

    // 1) Compute RREF of H
    auto [RREF, pivotCols, rank] = gaussianEliminationGF2(H);

    // pivotCols: columns where we have pivots
    // freeCols: columns that are not in pivotCols
    vector<int> freeCols;
    {
        vector<int> isPivot(cols, 0);
        for (int pc : pivotCols) {
            isPivot[pc] = 1;
        }
        for (int c = 0; c < cols; c++) {
            if (!isPivot[c]) freeCols.push_back(c);
        }
    }

    // The dimension of the code = # of free columns = k
    int k = (int)freeCols.size();

    // If the code is the zero code (rank = n), then only codeword is the zero vector
    if (k == 0) {
        // Check if the zero vector indeed satisfies Hx=0 (it does).
        return {string(cols, '0')};
    }

    // 2) For each free column, create a basis vector that sets that free col=1, others=0,
    //    then solve for pivot columns using RREF relationships.

    vector<vector<int>> basis; // each basis vector is length 'cols'
    basis.reserve(k);

    for (int i = 0; i < k; i++) {
        vector<int> v(cols, 0);
        // set the i-th free column to 1
        int fc = freeCols[i];
        v[fc] = 1;

        // Solve pivot variables from RREF
        // If row i has a pivot in pivotCols[i], that row implies pivot variable = sum of free variables
        for (int pivot_i = 0; pivot_i < (int)pivotCols.size(); pivot_i++) {
            int pcol = pivotCols[pivot_i];
            // The pivot is in row pivot_i. The row is RREF[pivot_i].
            int sumFree = 0;
            for (int c : freeCols) {
                if (RREF[pivot_i][c] == 1) {
                    sumFree ^= v[c];
                }
            }
            // pivot variable = sumFree to satisfy Hx=0
            v[pcol] = sumFree;
        }

        basis.push_back(v);
    }

    // 3) Enumerate all linear combinations of the k basis vectors
    vector<string> allCodewords;
    allCodewords.reserve((size_t)1 << k);

    // We'll generate every 0/1 combination of these k basis vectors
    for (int mask = 0; mask < (1 << k); mask++) {
        vector<int> codeword(cols, 0);
        // Combine basis vectors according to bits in 'mask'
        for (int b = 0; b < k; b++) {
            if ((mask >> b) & 1) {
                // XOR basis[b] into codeword
                for (int c = 0; c < cols; c++) {
                    codeword[c] ^= basis[b][c];
                }
            }
        }
        allCodewords.push_back(vectorToString(codeword));
    }

    // Optional: sort allCodewords lexicographically
    sort(allCodewords.begin(), allCodewords.end());
    
    return allCodewords;
}

// Add this function to src/generate_codeword.cpp

/*
 * Compute the Hamming weight (number of 1s) of a binary string
 */
int hammingWeight(const string& s) {
    int weight = 0;
    for (char c : s) {
        if (c == '1') weight++;
    }
    return weight;
}

/*
 * Compute the minimum Hamming distance of the code
 * The minimum distance is the minimum weight of any non-zero codeword
 * Returns:
 * - The minimum Hamming distance
 * - Returns -1 if the code contains only the zero codeword
 */
int computeMinimumDistance(const vector<vector<int>>& H) {
    // Get all codewords
    vector<string> codewords = computeAllCodewordsGF2(H);
    
    // Find minimum weight of non-zero codewords
    int minWeight = -1;  // -1 indicates no non-zero codeword found yet
    
    for (const string& codeword : codewords) {
        int weight = hammingWeight(codeword);
        if (weight > 0) {  // Skip the zero codeword
            if (minWeight == -1 || weight < minWeight) {
                minWeight = weight;
            }
        }
    }
    
    return minWeight;  // Will be -1 if only zero codeword exists
}

// Example main
// int main(){
//     // Example usage: H is a 2x4 parity-check matrix
//     // H = [ 1 1 0 1
//     //       0 1 1 1 ]
//     // Over GF(2)
//     vector<vector<int>> H = {
//         {1,1,0,0},
//         {0,1,1,0},
//         {0,0,1,1}
//     };

//     vector<string> codewords = computeAllCodewordsGF2(H);

//     cout << "All codewords of the code defined by H:\n";
//     for(const auto &cw : codewords){
//         cout << cw << "\n";
//     }
//     return 0;
// }
