#pragma once

#include "spbdd/pset.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace spbdd {

class code {
public:
    code(const context &ctx, const std::vector<std::string> &gens);

    int n() const { return n_; }        // # of physical qubit
    int k() const { return k_; }        // # of logical qubit
    int r() const { return n_ - k_; }   // # of stabilizer generators

    std::vector<std::string> stabilizer_generators() const;
    std::vector<std::string> logical_X_operators() const;
    std::vector<std::string> logical_Z_operators() const;

    pset S() const;     // S : stabilizer group
    pset N() const;     // N(S) : normalizer of stabilizer group
    pset NmS() const;   // N(S)-S : the non-trivial logical operators

    std::vector<bool> syndrome_of(const std::string &e) const;                 // r bits
    std::vector<bool> syndrome_of(const std::vector<bool> &e) const;
    std::vector<bool> logical_component_of(const std::string &e) const;        // 2k bits
    std::vector<bool> logical_component_of(const std::vector<bool> &e) const;

    pset by_syndrome(const std::vector<bool> &syndrome) const;
    pset by_logical_component(const std::vector<bool> &component) const;

    enum class fqcv_method { bruteforce, table, image, square };
    using violation = std::optional<std::pair<std::string, std::string>>;
    violation find_qec_cond_violation(const pset &errors, fqcv_method m = fqcv_method::image) const;

private:
    void check_context_id(const char *who) const;

    violation fqcv_bruteforce(const pset &errors) const;    
    violation fqcv_table(const pset &errors) const;        
    violation fqcv_image(const pset &errors) const;        
    violation fqcv_square(const pset &errors) const;       

    int ctx_id_ = 0;
    int n_ = 0, k_ = 0;

    std::vector<std::vector<bool>> g_, lx_, lz_; 

    std::vector<std::vector<int>> shift_table_; // the pattern of each xvar and zvar individually.
};

} // namespace spbdd
