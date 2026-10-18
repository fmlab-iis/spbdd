#include "gf2.hpp"
#include "pchar.hpp"

namespace spbdd {

// bits : (x_0, z_0, x_1, z_1, ...).
bits to_bits(const std::string &s, int n) {
    bits v;
    for (pauli p : string_to_pauli(s, n)) {
        v.push_back(has_x(p));
        v.push_back(has_z(p));
    }
    return v;
}

std::string bits_to_string(const bits &v) {
    std::vector<pauli> p;
    for (std::size_t q = 0; 2 * q < v.size(); ++q)
        p.push_back(make_pauli(v[2 * q], v[2 * q + 1]));
    return pauli_to_string(p);
}

std::vector<std::string> bits_to_strings(const std::vector<bits> &rows) {
    std::vector<std::string> out;
    for (const bits &v : rows)
        out.push_back(bits_to_string(v));
    return out;
}

std::vector<int> support(const bits &h) {
    std::vector<int> vars;
    for (int i = 0; i < static_cast<int>(h.size()); ++i)
        if (h[i]) vars.push_back(i);
    return vars;
}

bits swap_xz(const bits &a) {
    bits s(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        s[i] = a[i ^ 1];  // x_q = 2q, z_q = 2q + 1, so i ^ 1 is the other one
    return s;
}

bits multiply(const bits &a, const bits &b) {
    bits out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        out[i] = a[i] ^ b[i];
    return out;
}

bool symplectic_inner_product(const bits &a, const bits &b) {
    bool s = false;
    for (std::size_t i = 0; i < a.size(); ++i)
        s ^= a[i] && b[i ^ 1];  // x pairs with z: i ^ 1 is the partner variable
    return s;
}

// Return the basis of null(M) = { h : <m, h> = 0 for every row m of M } w.r.t.
// the symplectic inner product, where M has `cols` columns.
std::vector<bits> null_space(std::vector<bits> M, int cols) {

    // <m, h> = swap_xz(m) . h, so after swapping, the condition reads M*h = 0.
    for (bits &m : M)
        m = swap_xz(m);

    /* Row operations are left-multiplication by an invertible E, so
       R = E*M has the same null space: M*h = 0  <=>  R*h = 0.

       Reduce M to RREF R by Gaussian Elimination. 
       Let p_i be the pivot column of row i of R (i.e. R[i][p_i] = 1) 
       and column p_i is 0 in every other row (i.e. R[j][p_i] = 0, j!=i). 
       
       Thanks to the pivot column, now row i of R*h = 0 reads:

           h[p_i] + ( sum_{f : free column} R[i][f] h[f] ) = 0,

       so over F_2: 

           h[p_i] = sum_{f : free column} R[i][f] h[f].

       Therefore, we now can arbitrarily decide free h[f], 
       and h[p_i] will be uniquely decided by free h[f]. 
       Here we set one free entry to 1 and the rest to 0 
       to get a basis vector per free column.
       (Here rank-nullity thm. is also clear: 
           dim(null(M)) = cols - rank(M) = # of free columns) */

    std::vector<int> pivots;
    for (int c = 0; c < cols && pivots.size() < M.size(); ++c) {
        const std::size_t row = pivots.size();
        std::size_t r = row;
        while (r < M.size() && !M[r][c])
            ++r;
        if (r == M.size())
            continue;
        std::swap(M[row], M[r]);
        for (std::size_t i = 0; i < M.size(); ++i)
            if (i != row && M[i][c])
                for (int j = 0; j < cols; ++j)
                    M[i][j] = M[i][j] ^ M[row][j];
        pivots.push_back(c);
    }
    std::vector<bool> is_pivot(cols, false);
    for (int c : pivots)
        is_pivot[c] = true;

    std::vector<bits> basis;
    for (int free = 0; free < cols; ++free) {
        if (is_pivot[free])
            continue;
        bits h(cols, false);
        h[free] = true;
        for (std::size_t i = 0; i < pivots.size(); ++i)
            h[pivots[i]] = M[i][free];
        basis.push_back(h);
    }
    return basis;
}

} // namespace spbdd
