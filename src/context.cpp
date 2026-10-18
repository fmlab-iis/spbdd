#include "spbdd/pset.hpp"
#include "gf2.hpp"
#include "pchar.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

namespace spbdd {

namespace { // helpers

bdd parity_bdd(const bits &h, bool target) {
    bdd f = target ? bddfalse : bddtrue;
    for (int v : support(swap_xz(h)))
        f ^= bdd_ithvar(v);
    return f;
}

bdd single_qubit_bdd(int q, pauli p) {
    bdd x = bdd_ithvar(2 * q), z = bdd_ithvar(2 * q + 1);
    return (has_x(p) ? x : !x) & (has_z(p) ? z : !z);
}

} // namespace

static int peak_nodes_ = 0;

// Called by BuDDy around each garbage collection; before it, the table is at its fullest.
static void track_peak(int before, bddGbcStat *s) {
    if (before)
        peak_nodes_ = std::max(peak_nodes_, s->nodes - s->freenodes);
}

static void throw_bdd_error(int e) {
    const char *msg = bdd_errstring(e);
    throw std::runtime_error(std::string("BuDDy: ") + (msg ? msg : "unknown error"));
}

context::context(int num_qubits, const config &cfg)
    : num_qubits_(num_qubits), dynamic_reordering_(cfg.dynamic_reordering) {

    if (num_qubits < 1)
        throw std::invalid_argument("context: num_qubits must be at least 1");
    if (alive())
        throw std::invalid_argument("context: another context is alive");
    if (bdd_init(cfg.initial_nodes, cfg.cache_size) < 0)
        throw std::runtime_error("context: bdd_init failed");

    bdd_error_hook(throw_bdd_error);
    peak_nodes_ = 0;
    bdd_gbc_hook(track_peak);
    try {
        bdd_setvarnum(num_vars());
        bdd_varblockall();   // reordering (automatic or reorder_now) only acts on variables in blocks
    } catch (...) {
        bdd_done();
        throw;
    }

    if (dynamic_reordering_) {
        bdd_autoreorder(BDD_REORDER_SIFT);
        bdd_enable_reorder(); // bdd_setvarnum at default switches it off
    }

    current_ = this;
    ++current_id_;
}

context::~context() {
    bdd_done();
    current_ = nullptr;
}

context& context::current() {
    if (!current_)
        throw std::logic_error("context: no live context");
    return *current_;
}

void context::grow_num_qubits_to(int num_qubits) {
    if (num_qubits < num_qubits_)
        throw std::invalid_argument("context: num_qubits cannot shrink");
    bdd_setvarnum(2 * num_qubits);
    bdd_clrvarblocks();
    bdd_varblockall();
    if (dynamic_reordering_)
        bdd_enable_reorder();
    num_qubits_ = num_qubits;
}

void context::set_dynamic_reordering(bool enabled) {
    dynamic_reordering_ = enabled;
    if (enabled) {
        bdd_autoreorder(BDD_REORDER_SIFT);
        bdd_enable_reorder();
    } else {
        bdd_disable_reorder();
    }
}

void context::reorder_now() {
    bdd_reorder(BDD_REORDER_SIFT);
}

std::size_t context::live_nodes() const {
    return bdd_getnodenum();
}

std::size_t context::peak_nodes() const {
    return std::max<std::size_t>(peak_nodes_, bdd_getnodenum());   // most nodes in use, as seen at GC time
}

std::size_t context::table_size() const {
    return bdd_getallocnum();
}

pset context::empty() const {
    return pset(bddfalse);
}

pset context::all() const {
    return pset(bddtrue);
}

pset context::identity() const {
    return from_str(std::string(num_qubits_, 'I'));
}

pset context::from_str(const std::string &s) const {
    const bits v = to_bits(s, num_qubits_);
    bdd f = bddtrue;
    for (int level = num_vars() - 1; level >= 0; --level) {   // bottom of the current order first
        const int i = bdd_level2var(level);
        f &= v[i] ? bdd_ithvar(i) : bdd_nithvar(i);
    }
    return pset(f);
}

pset context::from_list(const std::vector<std::string> &list) const {
    bdd f = bddfalse;
    for (const std::string &s : list)
        f |= from_str(s).f_;
    return pset(f);
}

pset context::from_list(std::initializer_list<std::string> list) const {
    return from_list(std::vector<std::string>(list));
}

pset context::generated_by(const std::vector<std::string> &generators) const {

    /* Let M be the matrix whose rows are the generators, so W = row(M) and
       W^⟂ = null(M). Here ⟂ is w.r.t. the symplectic inner product on F_2^N,
       which is non-degenerate, so (W^⟂)^⟂ = W. With {h_i} a basis of null(M), 
       by linearity of the inner product:

           v ∈ W  <=>  v ∈ (span{h_i})^⟂  <=>  <h_i, v> = 0  ∀ i

       Here <h_i, v> is the parity of v on supp(swap_xz(h_i)), so:

           f_W = AND_i [parity_bdd(h_i, false)]  */

    std::vector<bits> rows;
    for (const std::string &g : generators)
        rows.push_back(to_bits(g, num_qubits_));

    bdd f = bddtrue;
    for (const bits &h : null_space(rows, num_vars()))
        f &= parity_bdd(h, false);
    return pset(f);
}

pset context::coset_of(const std::string &base, const std::vector<std::string> &basis) const {
    return coset_of(base, generated_by(basis));
}

pset context::coset_of(const std::string &base, const pset &set) const {
    set.check_context_id("coset_of");
    bdd f = set.f_;
    for (int t : support(to_bits(base, num_qubits_)))
        f = bdd_compose(f, bdd_nithvar(t), t);
    return pset(f);
}

pset context::identity_on(const std::vector<int> &qubits) const {
    bdd f = bddtrue;
    for (int q : unique_qubits(qubits, num_qubits_))
        f &= single_qubit_bdd(q, pauli::I);
    return pset(f);
}

pset context::pauli_at(int qubit, char p) const {
    check_qubit_range(qubit, num_qubits_);
    return pset(single_qubit_bdd(qubit, char_to_pauli(p)));
}

pset context::matching(const std::string &pattern) const {
    if (static_cast<int>(pattern.size()) != num_qubits_)
        throw std::invalid_argument("matching: pattern length must equal num_qubits");
    bdd f = bddtrue;
    for (int q = 0; q < num_qubits_; ++q)
        if (pattern[q] != '*')
            f &= single_qubit_bdd(q, char_to_pauli(pattern[q]));
    return pset(f);
}

pset context::weight_exactly(int w) const {
    return weight_between(w, w);
}

pset context::weight_at_most(int w) const {
    return weight_between(0, w);
}

pset context::weight_between(int lo, int hi) const {
    std::vector<int> all(num_qubits_);
    std::iota(all.begin(), all.end(), 0);
    return weight_between(lo, hi, all);
}

pset context::weight_exactly(int w, const std::vector<int> &qubits) const {
    return weight_between(w, w, qubits);
}

pset context::weight_at_most(int w, const std::vector<int> &qubits) const {
    return weight_between(0, w, qubits);
}

pset context::weight_between(int lo, int hi, const std::vector<int> &qubits) const {
    const std::vector<int> qs = unique_qubits(qubits, num_qubits_);
    lo = std::max(lo, 0);
    hi = std::min(hi, static_cast<int>(qs.size()));
    if (lo > hi)
        return empty();

    /* cnt[c] is the weight c bdd as of current processing qubit 
       the update of cnt[c] on new qubit q is whether identity on q (i.e. cnt[c] & q_identity),
       or cnt[c-1] with non-identity on q (i.e. cnt[c-1] & !q_identity). */

    std::vector<bdd> cnt(hi + 1, bddfalse);
    cnt[0] = bddtrue;
    for (int q : qs) {
        const bdd q_identity = single_qubit_bdd(q, pauli::I);
        for (int c = hi; c >= 0; --c)
            cnt[c] = (cnt[c] & q_identity) | (c ? cnt[c - 1] & !q_identity : bddfalse);
    }

    bdd f = bddfalse;
    for (int c = lo; c <= hi; ++c)
        f |= cnt[c];
    return pset(f);
}

pset context::commute_with(const std::string &p) const {
    return pset(parity_bdd(to_bits(p, num_qubits_), false));
}

pset context::commute_with(const std::vector<std::string> &ps) const {
    bdd f = bddtrue;
    for (const std::string &p : ps)
        f &= commute_with(p).f_;
    return pset(f);
}

pset context::commute_with(std::initializer_list<std::string> ps) const {
    return commute_with(std::vector<std::string>(ps));
}

} // namespace spbdd
