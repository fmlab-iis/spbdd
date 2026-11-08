/*
    Flag protocol verification along a path, done the way the FTEC
    verifier's spbdd backend does it.

    A flag protocol is a control flow diagram: each node is a (flagged)
    syndrome extraction circuit, and each branch is conditioned on the
    syndrome or flag outcomes. The question checked here is the one asked at
    the end of every path: can t <= tau faults leave two errors E1, E2 that the
    decoder cannot tell apart (same measurement record, same syndrome) but
    whose product is a non-trivial logical operator?

    How the backend represents the computation:

      * A state is the set of possible errors, split by how many faults made
        them: by_t[t] = errors produced by exactly t fault locations
        (t = 0..tau). by_t[0] starts as {I}.
      * Gates conjugate every level. Every two-qubit gate is a fault location:
        by_t[t+1] |= by_t[t].fault_inject({control, target}), with all spawns
        read from the post-gate state, so one location faults at most once.
      * A measurement splits every level into the outcome-0 and outcome-1
        halves. One branch per distinct record of classical bits; branches
        that agree on all bits are the same observation and are unioned.
      * Ancillas are reset to I at the start and the end of each circuit.
      * States are interned: equal states get the same id (keyed by
        pset::hash()), which is what makes the step cache (state, circuit) and
        the check cache (state) exact.
      * At the end of a path, each level is checked on its own with
        code::find_qec_cond_violation().

    The example: the [[17,1,5]] 2D color code (d = 5, so tau = 2), measured
    WITHOUT flags. Each generator is extracted with a single syndrome
    ancilla, so one fault on that ancilla half way through spreads to the
    rest of the generator's support (a hook error), and two such faults can
    leave a pair of errors the decoder cannot tell apart. The path is

        (X-type weight-8 generator)_s      unflagged
            |
            |  [s] = [0]
            v
        (all 16 generators)_s              unflagged, one ancilla each
            |
            |  [s1..s16] = [0,...,0]
            v
           LUT

    and two faults are enough to break the QEC condition on it.
*/

#include "spbdd/spbdd.hpp"

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using namespace spbdd;

// ---- the code -------------------------------------------------------------

constexpr int n_data = 17;
constexpr int tau    = 2; // d = 5 tolerates two faults

// supports of the 8 faces (0-indexed); each face carries one X and one Z generator
const std::vector<std::vector<int>> faces = {
    {0, 1, 2, 3},   {0, 2, 4, 5},     {4, 5, 8, 9},     {6, 7, 10, 11},
    {8, 9, 12, 13}, {10, 11, 14, 15}, {7, 11, 15, 16},  {2, 3, 5, 6, 9, 10, 13, 14},
};

std::string pauli_on(char p, const std::vector<int> &support) {
    std::string g(n_data, 'I');
    for (int q : support) g[static_cast<std::size_t>(q)] = p;
    return g;
}

// ---- circuits -------------------------------------------------------------

struct instr {
    std::string gate;       // "cx", "cz", "meas_x"
    std::vector<int> qubits;
};
using circuit = std::vector<instr>;

// One generator measured with one ancilla prepared in |+> and read in the X
// basis: X-type generators by CX(ancilla -> data), Z-type by CZ. No flag.
void measure_generator(circuit &c, char p, const std::vector<int> &support, int anc) {
    for (int q : support) c.push_back({p == 'X' ? "cx" : "cz", {anc, q}});
}

// ---- the backend's state machine ------------------------------------------

using state = std::vector<pset>; // by_t[0..tau]

struct branch {
    state by_t;
    std::string bits; // classical outcomes so far, in measurement order
};

void apply_gate(state &by_t, const instr &g) {
    const int a = g.qubits[0], b = g.qubits[1];
    for (auto &level : by_t) {
        if (g.gate == "cx") level.cx(a, b);
        else if (g.gate == "cz") level.cz(a, b);
    }
    // Every two-qubit gate is a fault location; spawn from the post-gate
    // state before merging, so a location can fault at most once per lineage.
    std::vector<pset> spawn;
    for (int t = 0; t < tau; ++t) spawn.push_back(by_t[static_cast<std::size_t>(t)].fault_inject({a, b}));
    for (int t = 0; t < tau; ++t) by_t[static_cast<std::size_t>(t) + 1] |= spawn[static_cast<std::size_t>(t)];
}

std::vector<branch> split_on_measurement(const std::vector<branch> &branches, int q) {
    std::map<std::string, branch> merged; // branches with the same bits are one observation
    for (const auto &br : branches) {
        for (const bool bit : {false, true}) {
            branch next{{}, br.bits + (bit ? '1' : '0')};
            bool empty = true;
            for (const auto &level : br.by_t) {
                next.by_t.push_back(level.measure_x(q, bit));
                empty = empty && next.by_t.back().is_empty();
            }
            if (empty) continue; // this outcome cannot happen
            auto found = merged.find(next.bits);
            if (found == merged.end()) {
                merged.emplace(next.bits, std::move(next));
            } else {
                for (std::size_t t = 0; t < next.by_t.size(); ++t) found->second.by_t[t] |= next.by_t[t];
            }
        }
    }
    std::vector<branch> out;
    for (auto &[bits, br] : merged) out.push_back(std::move(br));
    return out;
}

// Run one circuit on one state: outcome bits -> resulting state.
std::map<std::string, state> step(const state &in, const circuit &c, const std::vector<int> &ancillas) {
    std::vector<branch> branches{{in, ""}};
    for (auto &level : branches.front().by_t) level.reset(ancillas);
    for (const auto &g : c) {
        if (g.gate == "meas_x") {
            branches = split_on_measurement(branches, g.qubits[0]);
        } else {
            for (auto &br : branches) apply_gate(br.by_t, g);
        }
    }
    std::map<std::string, state> out;
    for (auto &br : branches) {
        for (auto &level : br.by_t) level.reset(ancillas); // outcome recorded; ancillas are done
        out.emplace(br.bits, std::move(br.by_t));
    }
    return out;
}

// Interning: equal states -> equal ids. pset::hash() is the BDD root, so a
// bucket almost never holds more than the state itself.
struct state_table {
    std::vector<state> states;
    std::unordered_map<std::size_t, std::vector<int>> buckets;

    int intern(state s) {
        std::size_t key = 0;
        for (const auto &level : s) key = key * 1000003u ^ level.hash();
        auto &bucket = buckets[key];
        for (int id : bucket) {
            if (states[static_cast<std::size_t>(id)] == s) return id;
        }
        states.push_back(std::move(s));
        bucket.push_back(static_cast<int>(states.size()) - 1);
        return bucket.back();
    }
};

void print_sizes(const state &s) {
    for (std::size_t t = 0; t < s.size(); ++t) std::cout << "  t" << t << "=" << std::setw(8) << s[t].size();
    std::cout << "\n";
}

int main() {
    // qubits: data 0..16, then one syndrome ancilla per generator (16)
    const int n_anc = 16;
    context ctx(n_data, n_anc);

    std::vector<std::string> gens;
    for (char p : {'X', 'Z'})
        for (const auto &f : faces) gens.push_back(pauli_on(p, f));
    code c(ctx, gens);

    // node 1: the X-type weight-8 generator, unflagged, ancilla q17
    circuit node1;
    measure_generator(node1, 'X', faces[7], n_data);
    node1.push_back({"meas_x", {n_data}});
    const std::vector<int> anc1{n_data};

    // node 2: all 16 generators, unflagged, ancillas q17..q32
    circuit node2;
    std::vector<int> anc2;
    for (std::size_t i = 0; i < gens.size(); ++i) {
        const int anc = n_data + static_cast<int>(i);
        anc2.push_back(anc);
        measure_generator(node2, i < faces.size() ? 'X' : 'Z', faces[i % faces.size()], anc);
    }
    for (int anc : anc2) node2.push_back({"meas_x", {anc}});

    // the path: (node, outcome bits that take the branch)
    const std::vector<std::pair<const circuit *, std::pair<std::string, const std::vector<int> *>>> path = {
        {&node1, {"0", &anc1}},
        {&node2, {std::string(16, '0'), &anc2}},
    };

    state_table table;
    state fresh(tau + 1, ctx.empty());
    fresh[0] = ctx.identity();
    int id = table.intern(fresh);

    std::cout << "start      ";
    print_sizes(table.states[static_cast<std::size_t>(id)]);
    for (std::size_t k = 0; k < path.size(); ++k) {
        const auto &[node, choice] = path[k];
        const auto &[bits, ancillas] = choice;
        const auto outcomes = step(table.states[static_cast<std::size_t>(id)], *node, *ancillas);
        const auto taken = outcomes.find(bits);
        if (taken == outcomes.end()) {
            std::cout << "node " << k + 1 << ": outcome " << bits << " cannot happen\n";
            return 1;
        }
        id = table.intern(taken->second);
        std::cout << "node " << k + 1 << "     ";
        print_sizes(table.states[static_cast<std::size_t>(id)]);
        std::cout << "           (" << outcomes.size() << " possible outcomes, took " << bits << ")\n";
    }

    // End of the path: check each level on its own.
    const auto &final_state = table.states[static_cast<std::size_t>(id)];
    bool clean = true;
    for (std::size_t t = 0; t < final_state.size(); ++t) {
        const auto pair = c.find_qec_cond_violation(final_state[t]);
        if (!pair) continue;
        clean = false;
        const auto data = [](const std::string &p) { return p.substr(0, n_data); };
        const auto bits = [](const std::vector<bool> &v) {
            std::string s;
            for (bool x : v) s += x ? '1' : '0';
            return s;
        };
        std::cout << "violation with t = " << t << " faults:\n"
                  << "  E1 = " << data(pair->first) << "  syndrome " << bits(c.syndrome_of(pair->first))
                  << "  logical " << bits(c.logical_component_of(pair->first)) << "\n"
                  << "  E2 = " << data(pair->second) << "  syndrome " << bits(c.syndrome_of(pair->second))
                  << "  logical " << bits(c.logical_component_of(pair->second)) << "\n"
                  << "  same syndrome, different logical component: E1*E2 is a non-trivial logical\n";
        break;
    }
    if (clean) std::cout << "pass\n";
}
