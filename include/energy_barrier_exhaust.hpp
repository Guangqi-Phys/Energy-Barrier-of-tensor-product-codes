#ifndef ENERGY_BARRIER_EXHAUST_HPP
#define ENERGY_BARRIER_EXHAUST_HPP

#include <vector>

/*
 * Compute the number of violated parity checks for state x
 * given parity-check matrix H (ℓ x n). 
 * E(x) = Hamming weight of H*x^T over GF(2).
 */
int energyOfStateex(const std::vector<std::vector<int>>& H, const std::vector<int>& x);

/*
 * Exhaustively try ALL single-bit-flip paths from 0^n to c_target.
 * Track the minimal possible peak energy (barrier).
 *
 * This is exponential in the worst case (potentially exploring 
 * many paths if we allow revisits). We prune whenever we revisit 
 * a state with a worse or equal barrier than already found.
 *
 * Parameters:
 * H - Parity-check matrix (ℓ x n)
 * c_target - target codeword in {0,1}^n
 *
 * Returns:
 * The minimal energy barrier, or -1 if c_target is not reachable
 */
int computeEnergyBarrierExhaustive(
    const std::vector<std::vector<int>>& H,
    const std::vector<int>& c_target
);

/*
 * Helper function for recursive path exploration in brute force approach
 */
void exploreAllPaths(const std::vector<std::vector<int>>& H, 
                     const std::vector<int>& c_target,
                     std::vector<int>& current_state,
                     std::vector<bool>& visited,
                     int current_max_energy,
                     int& global_min_barrier);

/*
 * Compute the minimal energy barrier by exploring all possible paths
 * Warning: This is exponential time complexity - only suitable for small inputs
 */
int computeEnergyBarrierBruteForce(const std::vector<std::vector<int>>& H, 
                                  const std::vector<int>& c_target);

#endif // ENERGY_BARRIER_EXHAUST_HPP