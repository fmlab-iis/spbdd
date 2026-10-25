#pragma once

#include "context.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace spbdd {

struct meas_result;

class pset {
public:
    pset();   // the empty set, in the live context

    pset(const pset &o);
    pset& operator= (const pset &o);
    ~pset();

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

    pset  operator* (const pset &o) const;    // setwise product (or "Minkowski sum")
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
    // int         min_weight() const;

    // --- stabilizer / error correction simulation operations  -------------

    // Gates and reset act in place and return *this.

    // Gate propagation (current support only Clifford gates)
    pset& x(int qubit);
    pset& y(int qubit);
    pset& z(int qubit);
    pset& h(int qubit);
    pset& s(int qubit);
    pset& sdg(int qubit);
    pset& sx(int qubit);
    pset& sxdg(int qubit);
    pset& cx(int control, int target);
    pset& cy(int control, int target);
    pset& cz(int a, int b);
    pset& swap(int a, int b);

    // Fault injection
    pset fault_inject(const std::vector<int> &qubits) const;

    // Measure & Reset
    meas_result measure_z(int qubit) const;
    meas_result measure_x(int qubit) const;
    meas_result measure_y(int qubit) const;
    pset measure_z(int qubit, bool outcome) const;  // the part with this outcome
    pset measure_x(int qubit, bool outcome) const;
    pset measure_y(int qubit, bool outcome) const;
    pset& reset(const std::vector<int> &qubits); // reset {X,Y,Z} to I on certain qubits

    // --- enumeration ------------------------------------------------------
    // Visit every element in unspecified order; fn returns false for early stop.
    void for_each_element(const std::function<bool(const std::string &)> &fn) const;
    std::optional<std::string> any_element() const;

private:
    friend class context;
    friend class code;

    explicit pset(const bdd &f);

    bool is_old_ctx() const;
    void check_context_id(const char *who) const;
    void check_context_id(const pset &o, const char *who) const;   // both operands

    union { bdd f_; };
    int ctx_id_ = 0;
};

struct meas_result {
    pset meas_0;
    pset meas_1;
};

} // namespace spbdd
