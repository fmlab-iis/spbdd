/*
    Flag protocol verification along a path.

    A flag protocol is a control flow diagram: each node is a flagged
    syndrome extraction circuit, and each branch is conditioned on the
    syndrome or flag outcomes. Among the verification targets, one of the
    most important is whether the error set reaching the end of a path
    violates the quantum error-correction condition, which is what we
    show here.

    Here we choose an example from Bhatnagar et al., "Low-Depth Flag-Style
    Syndrome Extraction for Small Quantum Error-Correction Codes," Fig. 3.,
    the following path:

        (XZZXI)_f
            |
            |  [s,f] = [0,0]
            v
        (IXZZX)_f
            |
            |  [s,f] = [1,0]
            v
        (all generators)_s
            |
            |  [s1,s2,s3,s4] = [0,1,0,1]
            v
           LUT

    where the stabilizer code is the [[5,1,3]] code.
*/

#include "spbdd/spbdd.hpp"
#include <iostream>
#include <vector>
using namespace spbdd;

int main() {
    const int n_data = 5; // data qubits: q0-q4
    const int n_synd = 4; // syndrome ancilla: q5-q8
    const int n_flag = 1; // flag ancilla: q9

    context ctx(n_data, n_synd + n_flag);

    code c(ctx, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"});

    using circuit = std::vector<std::pair<std::string, std::vector<int>>>;

    // node 1: (XZZXI)_f
    const circuit node1 = {
        {"cx", {5, 0}}, {"cx", {5, 9}}, {"cz", {5, 1}}, {"cz", {5, 2}}, {"cx", {5, 9}}, {"cx", {5, 3}},
        {"meas_x", {5}}, {"meas_z", {9}},
    };

    // node 2: (IXZZX)_f
    const circuit node2 = {
        {"cx", {5, 1}}, {"cx", {5, 9}}, {"cz", {5, 2}}, {"cz", {5, 3}}, {"cx", {5, 9}}, {"cx", {5, 4}},
        {"meas_x", {5}}, {"meas_z", {9}},
    };

    // node 3: all generators, without flags
    const circuit node3 = {
        {"cx", {5, 0}}, {"cz", {5, 1}}, {"cz", {5, 2}}, {"cx", {5, 3}},
        {"cx", {6, 1}}, {"cz", {6, 2}}, {"cz", {6, 3}}, {"cx", {6, 4}},
        {"cx", {7, 0}}, {"cx", {7, 2}}, {"cz", {7, 3}}, {"cz", {7, 4}},
        {"cz", {8, 0}}, {"cx", {8, 1}}, {"cx", {8, 3}}, {"cz", {8, 4}},
        {"meas_x", {5}}, {"meas_x", {6}}, {"meas_x", {7}}, {"meas_x", {8}},
    };

    // errors[i]: errors caused by i faults
    std::vector<pset> errors(2);
    errors[0] = ctx.identity();
    errors[1] = ctx.weight_exactly(1, {0, 1, 2, 3, 4}) & ctx.identity_on({5, 6, 7, 8, 9});   // data error

    for (const auto &[g, q] : node1) {
        if (g == "meas_x" || g == "meas_z") break;
        // gate propagation
        if (g == "cx") errors[1].cx(q[0], q[1]);
        if (g == "cz") errors[1].cz(q[0], q[1]);
        // fault injection (Note: .fault_inject() returns a new pset object)
        errors[1] |= errors[0].fault_inject(q);
    }

    // measure and split the error set corresponding to the syndrome and flag outcome [s,f] = [0,0]
    errors[1] = errors[1].measure_x(5, 0).measure_z(9, 0);
    // reset and reuse the syndrome and flag qubit
    errors[1].reset({5, 9});

    for (const auto &[g, q] : node2) {
        if (g == "meas_x" || g == "meas_z") break;
        // gate propagation
        if (g == "cx") errors[1].cx(q[0], q[1]);
        if (g == "cz") errors[1].cz(q[0], q[1]);
        // fault injection
        errors[1] |= errors[0].fault_inject(q);
    }

    // measure and split the error set corresponding to [s,f] = [1,0], with no fault left
    errors[1] = errors[1].measure_x(5, 1).measure_z(9, 0);
    // reset and reuse the syndrome and flag qubit
    errors[1].reset({5, 9});

    for (const auto &[g, q] : node3) {
        if (g == "cx") errors[1].cx(q[0], q[1]);
        if (g == "cz") errors[1].cz(q[0], q[1]);
    }

    // measure and split the error set corresponding to [s1,s2,s3,s4] = [0,1,0,1]
    errors[1] = errors[1].measure_x(5, 0).measure_x(6, 1).measure_x(7, 0).measure_x(8, 1);
    errors[1].reset({5, 6, 7, 8});

    std::cout << (c.find_qec_cond_violation(errors[1]) ? "violation" : "pass") << "\n";
}
