#pragma once

/* GF(2) helpers */

#include <string>
#include <vector>

namespace spbdd {

using bits = std::vector<bool>;

bits to_bits(const std::string &s, int n);
std::string bits_to_string(const bits &v);
std::vector<std::string> bits_to_strings(const std::vector<bits> &rows);

std::vector<int> support(const bits &h); // Indices of the 1 entries.

bits swap_xz(const bits &a); // Exchange x and z of every qubit.
bits multiply(const bits &a, const bits &b); // Product of two Paulis, phase ignored.
bool symplectic_inner_product(const bits &a, const bits &b);

// Basis of { h : <m, h> = 0 for every row m of M }, <.,.> the symplectic inner product.
std::vector<bits> null_space(std::vector<bits> M, int cols);

} // namespace spbdd
