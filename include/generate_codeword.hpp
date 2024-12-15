#ifndef GENERATE_CODEWORD_HPP
#define GENERATE_CODEWORD_HPP

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string>

using namespace std;

// Helper function to count 1-bits in an integer (mod 2)
inline int xorBit(int a, int b);

// Convert a binary vector to a string like "0101"
string vectorToString(const vector<int>& vec);

// Gaussian Elimination (over GF(2)) to find the RREF of H
// Returns a tuple: (H in RREF, pivotCols, rank)
tuple<vector<vector<int>>, vector<int>, int> 
gaussianEliminationGF2(const vector<vector<int>>& H_in);

/*
 * ComputeAllCodewordsGF2(H):
 *   Input: H is an ℓ x n parity-check matrix over GF(2).
 *   Output: A vector of binary strings, each representing one codeword in the null space of H.
 */
vector<string> computeAllCodewordsGF2(const vector<vector<int>>& H);

/*
 * Compute the Hamming weight (number of 1s) of a binary string
 */
int hammingWeight(const string& s);

/*
 * Compute the minimum Hamming distance of the code
 * The minimum distance is the minimum weight of any non-zero codeword
 * Returns:
 * - The minimum Hamming distance
 * - Returns -1 if the code contains only the zero codeword
 */
int computeMinimumDistance(const vector<vector<int>>& H);

#endif // GENERATE_CODEWORD_HPP