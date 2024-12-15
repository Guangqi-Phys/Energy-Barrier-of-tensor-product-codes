#ifndef ENERGY_BARRIER_HPP
#define ENERGY_BARRIER_HPP

#include <vector>
#include <unordered_map>
#include <queue>
#include <string>

// Function to compute the syndrome H*x^T over GF(2) and return its Hamming weight.
int energyOfState(const std::vector<std::vector<int>>& H, const std::vector<int>& x);

// Function to compute the minimal energy barrier from the zero codeword to c_target by single-bit flips.
int computeEnergyBarrier(const std::vector<std::vector<int>>& H, const std::vector<int>& c_target);



#endif // ENERGY_BARRIER_HPP