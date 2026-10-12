#pragma once

#include "spbdd/context.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace spbdd {

struct measurement_split;

class pset {
public:
    explicit pset(const context &ctx);

    pset(const pset &o);
    pset(pset &&o) noexcept;
    pset& operator= (const pset &o);
    pset& operator= (pset &&o) noexcept;
    ~pset();

    // --- context ----------------------------------------------------------
    int  num_qubits() const;

    // --- set algebra ------------------------------------------------------
    pset operator| (const pset &o) const;     // union
    pset operator& (const pset &o) const;     // intersection
    pset operator- (const pset &o) const;     // difference
    pset operator^ (const pset &o) const;     // symmetric difference (xor)
    pset operator~ () const;                  // complement

    pset& operator|= (const pset &o);
    pset& operator&= (const pset &o);
    pset& operator-= (const pset &o);
    pset& operator^= (const pset &o);

    pset  operator* (const pset &o) const;
    pset& operator*= (const pset &o);

    bool operator== (const pset &o) const;
    bool operator!= (const pset &o) const { return !(*this == o); }

    // --- predicates -------------------------------------------------------
    bool is_empty() const;
    bool contains(const std::string &p) const;
    bool subset_of(const pset &o) const;

    // --- metrics ----------------------------------------------------------
    double      size() const;                
    std::size_t node_count() const;
    int         min_weight() const;         

    // --- derived sets -----------------------------------------------------
    pset fault_inject(const std::vector<int> &qubits) const;
    pset reset(const std::vector<int> &qubits) const;

    // --- measurement ------------------------------------------------------
    measurement_split measure_z(int qubit) const;
    measurement_split measure_x(int qubit) const;
    measurement_split measure_y(int qubit) const;

    // --- Clifford gates ---------------------------------------------------
    pset x(int qubit) const;     
    pset y(int qubit) const;
    pset z(int qubit) const;

    pset h(int qubit) const;      
    pset s(int qubit) const;    
    pset sdg(int qubit) const;
    pset sx(int qubit) const;     
    pset sxdg(int qubit) const;

    pset cx(int control, int target) const;
    pset cy(int control, int target) const;
    pset cz(int a, int b) const;
    pset swap(int a, int b) const;

    // --- enumeration ------------------------------------------------------
    void for_each(const std::function<bool(const std::string &)> &fn) const;

    std::vector<std::string>   to_strings(std::size_t max_elements = 4096) const;
    std::optional<std::string> any_element() const;

private:
    friend class context;
    friend class code;

    explicit pset(const bdd &f);

    void require_fresh(const char *who) const;           

    union { bdd f_; };      
    int ctx_id_ = 0;
};

struct measurement_split {
    pset unflipped;   
    pset flipped;     
};

} // namespace spbdd
