#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <string>
#include <algorithm>
using namespace std;

/*
 * Compute the syndrome H*x^T over GF(2) and return its Hamming weight.
 * H is an ℓ x n matrix, x is length n. The syndrome is length ℓ.
 * E(x) = number of parity checks violated.
 */
int energyOfState(const vector<vector<int>>& H, const vector<int>& x) {
    int rows = (int)H.size();
    int cols = (rows > 0) ? (int)H[0].size() : 0;
    int countViolated = 0;

    for(int r = 0; r < rows; r++) {
        // Compute row r dot product with x over GF(2)
        int dot = 0;
        for(int c = 0; c < cols; c++){
            dot ^= (H[r][c] & x[c]); 
        }
        // If dot == 1, parity check r is violated
        if(dot == 1) countViolated++;
    }
    return countViolated;
}

/*
 * Compute the minimal energy barrier from the zero codeword (all 0's) 
 * to c_target by single-bit flips. 
 * The energy barrier is the minimal possible max(E(x_t)) along any path.
 * 
 * We'll use a best-first search with a priority queue 
 * keyed by the current path's max energy so far.
 *
 * H: parity-check matrix (ℓ x n)
 * c_target: length n codeword (which must satisfy H*c_target^T = 0 for it to be in the code).
 *
 * Return: minimal energy barrier as an integer.
 */
int computeEnergyBarrier(const vector<vector<int>>& H, const vector<int>& c_target) {
    int n = (int)c_target.size();
    // Check trivial case
    bool isAllZero = true;
    for(int bit : c_target) if(bit == 1) { isAllZero = false; break; }
    if(isAllZero) {
        // If c_target is the zero vector, barrier is obviously 0
        return 0;
    }

    // We'll represent states as bitstrings in an integer form if n <= ~32, 
    // but let's do a vector<int> approach for clarity.
    // The search space can be huge, so be mindful for large n.

    // A structure to hold (peakSoFar, stateVector)
    // We'll store the state as a string or vector<int>.
    struct State {
        int peak;
        vector<int> x;
        bool operator>(const State &other) const {
            return peak > other.peak;
        }
    };

    // Priority queue (min-heap) based on 'peak'
    priority_queue< State, vector<State>, greater<State> > pq;

    // visited[state] will store the best known (lowest) max energy to reach 'state'.
    // We'll store states in a std::unordered_map keyed by string representation for memory reasons.
    // For bigger n, a more compact representation (bitset) might be needed.
    unordered_map<string,int> visited;
    
    // Helper to convert a vector<int> to a string key
    auto vecToString = [&](const vector<int>& v){
        string s; s.reserve(n);
        for(int b : v) s.push_back(b ? '1' : '0');
        return s;
    };

    // Start from the zero state
    vector<int> zeroState(n, 0);
    int e0 = energyOfState(H, zeroState); // Typically 0 if zeroState is a codeword
    State initState {e0, zeroState};
    pq.push(initState);
    visited[vecToString(zeroState)] = e0;

    // BFS / Dijkstra-like search
    while(!pq.empty()) {
        State curr = pq.top();
        pq.pop();

        // If we've reached c_target, curr.peak is the minimal barrier
        if(curr.x == c_target) {
            return curr.peak;
        }

        // If there's a better path to curr.x, skip
        string currKey = vecToString(curr.x);
        if(visited[currKey] < curr.peak) {
            continue;
        }

        // Explore neighbors by flipping each bit
        for(int i = 0; i < n; i++){
            vector<int> nextState = curr.x;
            nextState[i] ^= 1;  // flip bit i
            int eNext = energyOfState(H, nextState);
            int nextPeak = max(curr.peak, eNext);

            string nextKey = vecToString(nextState);
            if(!visited.count(nextKey) || visited[nextKey] > nextPeak) {
                visited[nextKey] = nextPeak;
                pq.push({nextPeak, nextState});
            }
        }
    }

    // If c_target is truly in the code, we should find it.
    // If we get here, something is off or c_target isn't actually a codeword.
    cerr << "ERROR: c_target not reachable. Is it a valid codeword?" << endl;
    return -1;
}



// ------------------- Example usage -------------------
// int main(){
//     // Example parity-check matrix H (3 checks x 5 bits)
//     // For instance:
//     //  H = [1 1 0 1 0
//     //       0 1 1 0 1
//     //       1 0 1 1 0 ]
//     vector<vector<int>> H = {
//         {1,1,0,0,0},
//         {0,1,1,0,0},
//         {0,0,1,1,0},
//         {0,0,0,1,1},
//         {1,0,0,0,1}
//     };

//     // Suppose we want the barrier to a certain target codeword c_target in {0,1}^5
//     // Make sure c_target satisfies H*c_target^T = 0
//     vector<int> c_target = {1,1,1,1,1}; // just an example; might or might not be in the code

//     // Quick check: If c_target is in the code, H*c_target^T = 0 over GF(2)
//     // Otherwise, the code might not find a path or the barrier is meaningless
//     // (We'll just proceed assuming it's valid or let the code handle it)

//     int barrier = computeEnergyBarrier(H, c_target);
//     if(barrier >= 0){
//         cout << "Energy barrier from identity to logical operator [";
//         for (size_t i = 0; i < c_target.size(); ++i) {
//             cout << c_target[i];
//             if (i < c_target.size() - 1) cout << ", ";
//         }
//         cout << "] is: " << barrier << endl;
//     } else {
//         cout << "c_target not reachable or not a valid codeword." << endl;
//     }

//     return 0;
// }


