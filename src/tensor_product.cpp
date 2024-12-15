#include <iostream>
#include <vector>
using namespace std;

/*
 * Build the tensor product parity-check matrix:
 *   H3 = [ (H1 ⊗ I_{n2}) ]   # top block
 *        [ (I_{n1} ⊗ H2)   ] # bottom block
 *
 * where H1 is m1 x n1, H2 is m2 x n2.
 * The result H3 is (m1*n2 + n1*m2) x (n1*n2).
 */
vector<vector<int>> buildTensorProductParityCheck(
    const vector<vector<int>>& H1, // m1 x n1
    const vector<vector<int>>& H2  // m2 x n2
) {
    int m1 = (int)H1.size();           
    int n1 = m1 ? (int)H1[0].size() : 0;  
    int m2 = (int)H2.size();           
    int n2 = m2 ? (int)H2[0].size() : 0;

    // Dimensions of the final matrix H3
    int rowsH3 = m1*n2 + n1*m2;  // total rows
    int colsH3 = n1*n2;         // total columns
    vector<vector<int>> H3(rowsH3, vector<int>(colsH3, 0));

    // ========== Top Block: H1 ⊗ I_{n2} ==========
    // Size: (m1*n2) x (n1*n2)
    // For each '1' in H1[i][j], place an n2 x n2 identity block at row offset (i*n2), col offset (j*n2).
    for(int i = 0; i < m1; i++){
        for(int j = 0; j < n1; j++){
            if(H1[i][j] == 1) {
                int rowOff = i*n2;
                int colOff = j*n2;
                // Place identity block I_{n2}
                for(int k = 0; k < n2; k++){
                    H3[rowOff + k][colOff + k] = 1;
                }
            }
        }
    }

    // ========== Bottom Block: I_{n1} ⊗ H2 ==========
    // Size: (n1*m2) x (n1*n2)
    // Placed in rows [m1*n2 .. m1*n2 + n1*m2 - 1].
    // For each '1' in I_{n1}[i][j] (which is true iff i==j), place H2 in the sub-block.
    int bottomBlockStart = m1*n2; // row offset for bottom block
    for(int i = 0; i < n1; i++){    // row index in I_{n1}
        for(int j = 0; j < n1; j++){ // col index in I_{n1}
            if(i == j) {
                // Place H2 at the sub-block
                int rowOff = bottomBlockStart + i*m2;
                int colOff = j*n2;
                for(int r = 0; r < m2; r++){
                    for(int c = 0; c < n2; c++){
                        if(H2[r][c] == 1){
                            H3[rowOff + r][colOff + c] = 1;
                        }
                    }
                }
            }
        }
    }

    return H3;
}


/*
 * Build the tensor product codeword from two codewords c1 and c2.
 * The result is a codeword of length n1 * n2.
 */

vector<int> buildTensorProductCodeword(const vector<int>& c1, const vector<int>& c2) {
    // Get lengths of input codewords
    int n1 = c1.size();
    int n2 = c2.size();
    
    // The tensor product codeword will have length n1 * n2
    vector<int> result(n1 * n2);
    
    // Build tensor product codeword
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            // Tensor product is the multiplication of corresponding elements
            result[i * n2 + j] = c1[i] * c2[j];
        }
    }
    
    return result;
}

// ------------------- Example usage -------------------
// int main(){
//     // Example H1: 2 x 3
//     vector<vector<int>> H1 = {
//         {1,0,1},
//         {0,1,1}
//     };
//     // Example H2: 3 x 2
//     vector<vector<int>> H2 = {
//         {1,1},
//         {0,1},
//         {1,0}
//     };

//     // Build the tensor product parity-check
//     vector<vector<int>> H3 = buildTensorProductParityCheck(H1, H2);

//     // Print the resulting matrix H3
//     cout << "H3 dimension: " 
//          << H3.size() << " x " << (H3.empty()?0:H3[0].size()) << endl;
//     for(const auto &row : H3){
//         for(int val : row) cout << val;
//         cout << endl;
//     }

//     return 0;
// }
