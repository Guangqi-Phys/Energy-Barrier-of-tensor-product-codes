#ifndef TENSOR_PRODUCT_HPP
#define TENSOR_PRODUCT_HPP

#include <vector>

/*
 * Build the tensor product parity-check matrix:
 *   H3 = [ (H1 ⊗ I_{n2}) ]   # top block
 *        [ (I_{n1} ⊗ H2)   ] # bottom block
 *
 * where H1 is m1 x n1, H2 is m2 x n2.
 * The result H3 is (m1*n2 + n1*m2) x (n1*n2).
 */
std::vector<std::vector<int>> buildTensorProductParityCheck(
    const std::vector<std::vector<int>>& H1, // m1 x n1
    const std::vector<std::vector<int>>& H2  // m2 x n2
);

#endif // TENSOR_PRODUCT_HPP
