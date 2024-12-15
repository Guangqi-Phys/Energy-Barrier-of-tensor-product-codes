#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <functional>
#include <climits>
using namespace std;

/*
 * Compute the number of violated parity checks for state x
 * given parity-check matrix H (ℓ x n). 
 * E(x) = Hamming weight of H*x^T over GF(2).
 */
int energyOfStateex(const vector<vector<int>>& H, const vector<int>& x) {
    int rows = (int)H.size();
    int cols = rows ? (int)H[0].size() : 0;
    int countViolated = 0;
    for(int r = 0; r < rows; r++){
        int dot = 0;
        for(int c = 0; c < cols; c++){
            dot ^= (H[r][c] & x[c]);
        }
        if(dot == 1) countViolated++;
    }
    return countViolated;
}

/*
 * Exhaustively try ALL single-bit-flip paths from 0^n to c_target.
 * Track the minimal possible peak energy (barrier).
 *
 * This is exponential in the worst case (potentially exploring 
 * many paths if we allow revisits). We prune whenever we revisit 
 * a state with a worse or equal barrier than already found.
 */
int computeEnergyBarrierExhaustive(
    const vector<vector<int>>& H,       // Parity-check matrix (ℓ x n)
    const vector<int>& c_target        // target codeword in {0,1}^n
){
    int n = (int)c_target.size();
    
    // Quick check if c_target is the all-zero codeword
    bool allZero = true;
    for(int b : c_target) {
        if(b == 1){ allZero = false; break; }
    }
    if(allZero) {
        return 0; // trivial barrier
    }

    // We'll store the minimum barrier found for each visited state to prune paths
    // key: string representation of the n-bit state, value: best (lowest) barrier so far
    unordered_map<string,int> bestBarrierForState;

    // Convert a state vector<int> to a string "0101..."
    auto vecToString = [&](const vector<int>& v){
        string s; s.reserve(v.size());
        for(int bit : v) s.push_back(bit ? '1':'0');
        return s;
    };

    // A global variable (or captured reference) to store the best barrier found
    // for a path that reaches c_target.
    int globalMinBarrier = INT_MAX;

    // Depth-first search (DFS) recursion
    // currentState: current bit configuration
    // currentBarrier: the highest energy encountered so far along the path
    function<void(vector<int>,int)> dfs = [&](vector<int> state, int currentBarrier){
        // If we've reached c_target, update globalMinBarrier
        if(state == c_target) {
            globalMinBarrier = min(globalMinBarrier, currentBarrier);
            return;
        }

        // Prune if currentBarrier is already worse than globalMinBarrier
        if(currentBarrier >= globalMinBarrier) {
            return;
        }

        // Explore single-bit flips from 'state'
        for(int i = 0; i < n; i++){
            vector<int> nextState = state;
            nextState[i] ^= 1; // flip bit i
            int eNext = energyOfStateex(H, nextState);
            int nextBarrier = max(currentBarrier, eNext);

            // If we haven't visited nextState or found a better barrier now:
            string key = vecToString(nextState);
            if(!bestBarrierForState.count(key) || bestBarrierForState[key] > nextBarrier){
                bestBarrierForState[key] = nextBarrier;
                dfs(nextState, nextBarrier);
            }
        }
    };

    // Start from zero state
    vector<int> zeroState(n, 0);
    int e0 = energyOfStateex(H, zeroState); // usually 0 if zeroState is a valid codeword
    bestBarrierForState[vecToString(zeroState)] = e0;

    // Launch DFS
    dfs(zeroState, e0);

    // If globalMinBarrier is still INT_MAX, it means c_target wasn't reached
    if(globalMinBarrier == INT_MAX) {
        // Possibly c_target is not a valid codeword or unreachable via single-bit flips
        cerr << "Exhaustive search: c_target not reached. Possibly invalid codeword." << endl;
        return -1;
    }

    return globalMinBarrier;
}

/*
 * Helper function for recursive path exploration
 */
void exploreAllPaths(const vector<vector<int>>& H, 
                     const vector<int>& c_target,
                     vector<int>& current_state,
                     vector<bool>& visited,
                     int current_max_energy,
                     int& global_min_barrier) {
    
    // If we reached the target state, update the minimum barrier if needed
    if(current_state == c_target) {
        global_min_barrier = min(global_min_barrier, current_max_energy);
        return;
    }
    
    int n = (int)current_state.size();
    
    // Try flipping each bit
    for(int i = 0; i < n; i++) {
        // Create next state by flipping bit i
        current_state[i] ^= 1;
        
        // Convert state to string for visited checking
        string state_key;
        for(int b : current_state) state_key.push_back(b ? '1' : '0');
        
        // If we haven't visited this state
        if(!visited[stoi(state_key, nullptr, 2)]) {
            visited[stoi(state_key, nullptr, 2)] = true;
            
            // Calculate energy of new state
            int new_energy = energyOfStateex(H, current_state);
            int new_max_energy = max(current_max_energy, new_energy);
            
            // Only continue if this path could potentially give us a better barrier
            if(new_max_energy < global_min_barrier) {
                exploreAllPaths(H, c_target, current_state, visited, 
                              new_max_energy, global_min_barrier);
            }
            
            visited[stoi(state_key, nullptr, 2)] = false;
        }
        
        // Restore the bit
        current_state[i] ^= 1;
    }
}

/*
 * Compute the minimal energy barrier by exploring all possible paths
 * Warning: This is exponential time complexity - only suitable for small inputs
 */
int computeEnergyBarrierBruteForce(const vector<vector<int>>& H, 
                                  const vector<int>& c_target) {
    int n = (int)c_target.size();
    
    // Check trivial case
    bool isAllZero = true;
    for(int bit : c_target) if(bit == 1) { isAllZero = false; break; }
    if(isAllZero) return 0;
    
    // Start from zero state
    vector<int> current_state(n, 0);
    
    // Track visited states to avoid cycles
    vector<bool> visited(1 << n, false); // 2^n possible states
    string init_state(n, '0');
    visited[0] = true; // mark zero state as visited
    
    // Initialize with worst possible barrier
    int global_min_barrier = n + 1;
    
    // Start the recursive exploration
    exploreAllPaths(H, c_target, current_state, visited, 
                   energyOfStateex(H, current_state), global_min_barrier);
    
    return (global_min_barrier == n + 1) ? -1 : global_min_barrier;
}

// ------------------- Example usage -------------------
// int main(){
//     // Example parity-check matrix H (3 checks x 4 bits, for instance)
//     vector<vector<int>> H = {
//         {1,1,0,1},
//         {0,1,1,0},
//         {1,0,1,0}
//     };
//     // Some example target codeword c_target in {0,1}^4
//     vector<int> c_target = {1,1,0,0}; // Make sure it's actually in the code if we want a successful path.

//     int barrier = computeEnergyBarrierExhaustive(H, c_target);
//     if(barrier >= 0){
//         cout << "Energy barrier (exhaustive) from 0^n to c_target: " << barrier << endl;
//     } else {
//         cout << "No valid path found or c_target not in the code." << endl;
//     }

//     return 0;
// }
