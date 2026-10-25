#pragma once

#include "bdd.h"   // copy of BuDDy's bdd.h, made by make

#include <cstddef>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace spbdd {

class pset;

struct config { // for BuDDy init
    int  initial_nodes      = 1 << 20;
    int  cache_size         = 1 << 16;
    bool dynamic_reordering = false;      // sifting
};

class context {
public:
    explicit context(int num_qubits, const config &cfg = {});
    context(int num_data, int num_ancillas, const config &cfg = {});   // data qubits come first
    ~context();

    // Prevent the context object from being copied or moved.
    context(const context &)             = delete;
    context& operator= (const context &) = delete;

    // --- the live context -------------------------------------------------
    static bool     alive() { return current_ != nullptr; }
    static context& current();

    // --- set space --------------------------------------------------------
    void grow_num_qubits_to(int num_qubits);   // the new qubits are ancillas
    int num_qubits() const { return num_qubits_; }
    int num_data() const { return num_data_; }
    int num_ancillas() const { return num_qubits_ - num_data_; }
    int num_vars() const { return 2 * num_qubits_; }

    // --- variable order ---------------------------------------------------
    void set_dynamic_reordering(bool enabled);
    bool dynamic_reordering() const { return dynamic_reordering_; }
    void reorder_now();

    // --- statistics -------------------------------------------------------
    std::size_t live_nodes() const;
    std::size_t peak_nodes() const;
    std::size_t table_size() const;

    // --- primitives -------------------------------------------------------
    pset empty() const;
    pset all() const;
    pset identity() const;

    // --- explicitly listed elements ---------------------------------------
    pset from_str(const std::string &s) const;
    pset from_list(const std::vector<std::string> &list) const;
    pset from_list(std::initializer_list<std::string> list) const; // prevent compilation bug regarding input as char*

    // --- subgroups --------------------------------------------------------
    pset generated_by(const std::vector<std::string> &generators) const;
    pset coset_of(const std::string &base, const std::vector<std::string> &basis) const;
    pset coset_of(const std::string &base, const pset &set) const;

    // --- support ----------------------------------------------------------
    pset identity_on(const std::vector<int> &qubits) const;
    pset pauli_at(int qubit, char p) const;
    pset matching(const std::string &pattern) const;

    // --- weight -----------------------------------------------------------
    pset weight_exactly(int w) const;
    pset weight_at_most(int w) const;
    pset weight_between(int lo, int hi) const;
    pset weight_exactly(int w, const std::vector<int> &qubits) const;
    pset weight_at_most(int w, const std::vector<int> &qubits) const;
    pset weight_between(int lo, int hi, const std::vector<int> &qubits) const;

    // --- commutation ------------------------------------------------------
    pset commute_with(const std::string &p) const;
    pset commute_with(const std::vector<std::string> &ps) const;
    pset commute_with(std::initializer_list<std::string> ps) const;

private:
    friend class pset;
    friend class code;

    static int xvar(int q) { return 2 * q; }
    static int zvar(int q) { return 2 * q + 1; }

    static inline context *current_ = nullptr;
    static inline int      current_id_ = 0;

    int  num_qubits_;
    int  num_data_;
    bool dynamic_reordering_ = false;
};

} // namespace spbdd
