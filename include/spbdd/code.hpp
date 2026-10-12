#pragma once

#include "spbdd/pset.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace spbdd {

class code {
public:
    // The generators need not be independent, but must commute pairwise.
    code(const context &ctx, const std::vector<std::string> &generators);

    int num_qubits() const { return n_; }
    int num_stabilizers() const { return r_; }   // rank of the generators
    int num_logical() const { return k_; }       // n - r

    // --- the symplectic basis ---------------------------------------------
    const std::vector<std::string>& stabilizers() const { return g_; }
    const std::vector<std::string>& destabilizers() const { return d_; }
    const std::vector<std::string>& logical_x() const { return lx_; }
    const std::vector<std::string>& logical_z() const { return lz_; }

    // --- the sets ---------------------------------------------------------
    pset group() const;         // S
    pset normalizer() const;    // N(S)
    int  distance() const;      // lightest element of N(S) \ S

    // --- reading one operator ---------------------------------------------
    std::vector<bool> syndrome(const std::string &e) const;           // r bits
    std::vector<bool> logical_signature(const std::string &e) const;  // 2k bits (c, f)

    // --- reading a whole class --------------------------------------------
    pset with_syndrome(const std::vector<bool> &syndrome) const;
    pset with_logical_signature(const std::vector<bool> &signature) const;

    // --- the decision problem ---------------------------------------------
    // Does `errors` contain two operators whose product lies in N(S) \ S?
    // Computed by one pass over the nodes of `errors` (see docs/dp-image.tex).
    bool has_inequivalent_pair(const pset &errors) const;

    struct pair {
        std::string       first, second;     // two members of `errors`
        std::vector<bool> syndrome;          // the one they share
        std::vector<bool> signature_first, signature_second;
    };

    std::optional<pair> find_inequivalent_pair(const pset &errors) const;

private:
    using row = std::vector<std::uint8_t>;   // a vector over GF(2), length 2n

    void require_fresh(const char *who) const;

    // Variables of the (syndrome, signature) relation: b, then c, then f.
    int b_var(int i) const { return r_ + i; }
    int c_var(int j) const { return 2 * r_ + j; }
    int f_var(int j) const { return 2 * r_ + k_ + j; }

    // G = { (syndrome(e), signature(e)) : e in f }, by the node-wise DP.
    bdd relation_by_dp(const bdd &f) const;

    int ctx_id_ = 0;
    int n_ = 0, r_ = 0, k_ = 0;

    std::vector<std::string> g_, d_, lx_, lz_;
    std::vector<row>         grow_, drow_, lxrow_, lzrow_;

    // For each of the 2n original variables, the b / c / f variables flipped
    // by the single-qubit Pauli it stands for.
    std::vector<std::vector<int>> flips_;
};

} // namespace spbdd
