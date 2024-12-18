# Error Correction Code Analysis Tool

This project provides tools for analyzing error correction codes, specifically focusing on computing minimum distances and energy barriers for parity-check matrices and their tensor products.

## Features

- Computation of minimum distances for parity-check matrices
- Calculation of energy barriers for codes
- Generation of all possible codewords in GF(2)
- Tensor product construction of classical codes
- Computation of energy barrier for tensor product codes



## Usage Example

The program analyzes three parity-check matrices:
1. H1 (input matrix)
2. H2 (input matrix)
3. H3 (tensor product of H1 and H2)

For each matrix, it computes:
- The minimum distance
- The energy barrier
- The codeword achieving the minimum energy barrier


## Output Format

The program outputs:
- Matrix representations
- Minimum distances
- Energy barriers
- Achieving codewords
- A final summary comparing all three matrices

Example output:
```
=== Analysis of H1 ===
[Matrix representation]
Minimum distance of H1: X
Energy barrier of H1: Y
Achieved by codeword: Z

[... similar for H2 and H3 ...]

=== Summary ===
H1: distance = X, energy barrier = Y
H2: distance = X, energy barrier = Y
H3 (tensor product): distance = X, energy barrier = Y
```
