#include "spbdd/pset.hpp"

#include "freeze_order.hpp"
#include "gf2.hpp"
#include "pchar.hpp"

#include <algorithm>
#include <climits>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

namespace spbdd {

namespace {

class sumset {
    /* Core concept: A * B = ITE(v, (A|v=1 * B|v=0) U (A|v=0 * B|v=1), (A|v=0 * B|v=0) U (A|v=1 * B|v=1))
       The reason to design sumset as a class is to build memo for dynamic programming. */
public:
    bdd run(const bdd &f, const bdd &g) {
        if (f == bddfalse || g == bddfalse)
            return bddfalse;
        if (f == bddtrue || g == bddtrue)
            return bddtrue;

        const int a = f.id(), b = g.id(); // here id is the BuDDy node id
        const std::pair<int, int> key(std::min(a, b), std::max(a, b)); // f * g = g * f
        const auto hit = memo_.find(key);
        if (hit != memo_.end())
            return hit->second;

        const int lf = level(f), lg = level(g);
        bdd r;
        if (lf < lg) {
            r = run(bdd_low(f), g) | run(bdd_high(f), g);
        } else if (lg < lf) {
            r = run(f, bdd_low(g)) | run(f, bdd_high(g));
        } else {
            const bdd same  = run(bdd_low(f), bdd_low(g))  | run(bdd_high(f), bdd_high(g));
            const bdd cross = run(bdd_low(f), bdd_high(g)) | run(bdd_high(f), bdd_low(g));
            r = bdd_ite(bdd_ithvar(bdd_var(f)), cross, same);
        }
        memo_.emplace(key, r);
        return r;
    }

private:
    static int level(const bdd &f) { return bdd_var2level(bdd_var(f)); }

    std::map<std::pair<int, int>, bdd> memo_;
};

// f with variable a replaced by (a xor b).
bdd xor_into(const bdd &f, int a, int b) {
    return bdd_compose(f, bdd_ithvar(a) ^ bdd_ithvar(b), a);
}

// f with the variables of each pair exchanged.
bdd swap_vars(const bdd &f, std::initializer_list<std::pair<int, int>> pairs) {
    std::unique_ptr<bddPair, void (*)(bddPair *)> pr(bdd_newpair(), bdd_freepair);
    for (const auto &p : pairs) {
        bdd_setpair(pr.get(), p.first, p.second);
        bdd_setpair(pr.get(), p.second, p.first);
    }
    return bdd_replace(f, pr.get());
}

// Cube of the x and z variables of the qubits, for quantification.
bdd xz_cube(const std::vector<int> &qubits) {
    bdd cube = bddtrue;
    for (int q : unique_qubits(qubits, context::current().num_qubits()))
        cube &= bdd_ithvar(2 * q) & bdd_ithvar(2 * q + 1);
    return cube;
}

void check_distinct(int a, int b, const char *who) {
    if (a == b)
        throw std::invalid_argument(std::string(who) + ": qubits must differ");
}

int num_qubits() {
    return context::current().num_qubits();
}

} // namespace

pset::pset() : pset((context::current(), bddfalse)) {}   // throws if no context is alive

pset::pset(const bdd &f) : f_(f), ctx_id_(context::current_id_) {}

pset::pset(const pset &o) : ctx_id_(o.ctx_id_) {
    if (!o.is_old_ctx())
        new (&f_) bdd(o.f_);
}

pset& pset::operator= (const pset &o) {
    if (this == &o)
        return *this;
    const bool had = !is_old_ctx(), has = !o.is_old_ctx();
    if (had && has)
        f_ = o.f_;
    else if (had)
        f_.~bdd();
    else if (has)
        new (&f_) bdd(o.f_);
    ctx_id_ = o.ctx_id_;
    return *this;
}

pset::~pset() {
    if (!is_old_ctx())
        f_.~bdd();
}

pset pset::operator| (const pset &o) const {
    check_context_id(o, "operator|");
    return pset(f_ | o.f_);
}

pset pset::operator& (const pset &o) const {
    check_context_id(o, "operator&");
    return pset(f_ & o.f_);
}

pset pset::operator- (const pset &o) const {
    check_context_id(o, "operator-");
    return pset(f_ - o.f_);
}

pset pset::operator^ (const pset &o) const {
    check_context_id(o, "operator^");
    return pset(f_ ^ o.f_);
}

pset pset::operator~ () const {
    check_context_id("operator~");
    return pset(!f_);
}

pset& pset::operator|= (const pset &o) {
    return *this = *this | o;
}

pset& pset::operator&= (const pset &o) {
    return *this = *this & o;
}

pset& pset::operator-= (const pset &o) {
    return *this = *this - o;
}

pset& pset::operator^= (const pset &o) {
    return *this = *this ^ o;
}

pset pset::operator* (const pset &o) const {
    check_context_id(o, "operator*");
    freeze_order freeze;   // sumset needs the levels fixed
    return pset(sumset().run(f_, o.f_));
}

pset& pset::operator*= (const pset &o) {
    return *this = *this * o;
}

bool pset::operator== (const pset &o) const {
    check_context_id(o, "operator==");
    return f_ == o.f_;
}

bool pset::is_empty() const {
    check_context_id("is_empty");
    return f_ == bddfalse;
}

bool pset::contains(const std::string &p) const {
    check_context_id("contains");
    return (f_ & context::current().from_str(p).f_) != bddfalse;
}

bool pset::subset_of(const pset &o) const {
    check_context_id(o, "subset_of");
    return (f_ - o.f_) == bddfalse;
}

double pset::size() const {
    check_context_id("size");
    return bdd_satcount(f_);
}

std::size_t pset::hash() const {
    check_context_id("hash");
    return std::hash<int>{}(f_.id());
}

std::size_t pset::node_count() const {
    check_context_id("node_count");
    return bdd_nodecount(f_);
}

/* Clifford gates P propagates an error E in a way that : E -> P E P† */

pset& pset::x(int qubit) {
    check_context_id("x");
    check_qubit_range(qubit, num_qubits());
    return *this;
}

pset& pset::y(int qubit) {
    check_context_id("y");
    check_qubit_range(qubit, num_qubits());
    return *this;
}

pset& pset::z(int qubit) {
    check_context_id("z");
    check_qubit_range(qubit, num_qubits());
    return *this;
}

// S: X -> Y, Z -> Z, i.e. z ^= x.
pset& pset::s(int qubit) {
    check_context_id("s");
    check_qubit_range(qubit, num_qubits());
    f_ = xor_into(f_, context::zvar(qubit), context::xvar(qubit));
    return *this;
}

pset& pset::sdg(int qubit) {
    return s(qubit);
}

// SX: Z -> Y, X -> X, i.e. x ^= z.
pset& pset::sx(int qubit) {
    check_context_id("sx");
    check_qubit_range(qubit, num_qubits());
    f_ = xor_into(f_, context::xvar(qubit), context::zvar(qubit));
    return *this;
}

pset& pset::sxdg(int qubit) {
    return sx(qubit);
}

// H: X <-> Z, i.e. exchange x and z.
pset& pset::h(int qubit) {
    check_context_id("h");
    check_qubit_range(qubit, num_qubits());
    f_ = swap_vars(f_, {{context::xvar(qubit), context::zvar(qubit)}});
    return *this;
}

// CX: x_t ^= x_c, z_c ^= z_t.
pset& pset::cx(int control, int target) {
    check_context_id("cx");
    check_qubit_range(control, num_qubits());
    check_qubit_range(target, num_qubits());
    check_distinct(control, target, "cx");
    f_ = xor_into(f_, context::xvar(target), context::xvar(control));
    f_ = xor_into(f_, context::zvar(control), context::zvar(target));
    return *this;
}

// CY = S_t CX S_t^dagger.
pset& pset::cy(int control, int target) {
    return s(target).cx(control, target).s(target);
}

// CZ: z_b ^= x_a, z_a ^= x_b.
pset& pset::cz(int a, int b) {
    check_context_id("cz");
    check_qubit_range(a, num_qubits());
    check_qubit_range(b, num_qubits());
    check_distinct(a, b, "cz");
    f_ = xor_into(f_, context::zvar(b), context::xvar(a));
    f_ = xor_into(f_, context::zvar(a), context::xvar(b));
    return *this;
}

pset& pset::swap(int a, int b) {
    check_context_id("swap");
    check_qubit_range(a, num_qubits());
    check_qubit_range(b, num_qubits());
    check_distinct(a, b, "swap");
    f_ = swap_vars(f_, {{context::xvar(a), context::xvar(b)},
                        {context::zvar(a), context::zvar(b)}});
    return *this;
}

pset pset::fault_inject(const std::vector<int> &qubits) const {
    check_context_id("fault_inject");
    return pset(bdd_exist(f_, xz_cube(qubits)));
}

pset& pset::reset(const std::vector<int> &qubits) {
    f_ = fault_inject(qubits).f_ & context::current().identity_on(qubits).f_;
    return *this;
}

meas_result pset::measure_z(int qubit) const {
    check_context_id("measure_z");
    check_qubit_range(qubit, num_qubits());
    const bdd flip = bdd_ithvar(context::xvar(qubit));
    return {pset(f_ & !flip), pset(f_ & flip)};
}

meas_result pset::measure_x(int qubit) const {
    check_context_id("measure_x");
    check_qubit_range(qubit, num_qubits());
    const bdd flip = bdd_ithvar(context::zvar(qubit));
    return {pset(f_ & !flip), pset(f_ & flip)};
}

meas_result pset::measure_y(int qubit) const {
    check_context_id("measure_y");
    check_qubit_range(qubit, num_qubits());
    const bdd flip = bdd_ithvar(context::xvar(qubit)) ^ bdd_ithvar(context::zvar(qubit));
    return {pset(f_ & !flip), pset(f_ & flip)};
}

pset pset::measure_z(int qubit, bool outcome) const {
    const meas_result r = measure_z(qubit);
    return outcome ? r.meas_1 : r.meas_0;
}

pset pset::measure_x(int qubit, bool outcome) const {
    const meas_result r = measure_x(qubit);
    return outcome ? r.meas_1 : r.meas_0;
}

pset pset::measure_y(int qubit, bool outcome) const {
    const meas_result r = measure_y(qubit);
    return outcome ? r.meas_1 : r.meas_0;
}

void pset::for_each_element(const std::function<bool(const std::string &)> &fn) const {
    check_context_id("for_each_element");
    const int n = context::current().num_vars();
    std::vector<bool> bits(n, false);

    // Visits every element of the BDD; returns false once fn asks to stop.
    std::function<bool(const bdd &, int)> walk = [&](const bdd &f, int level) {
        if (f == bddfalse)
            return true;
        if (level == n)
            return fn(bits_to_string(bits));
        const int var = bdd_level2var(level);
        for (int v = 0; v < 2; ++v) {
            bits[var] = v;
            bdd next = f;
            if (f != bddtrue && bdd_var2level(bdd_var(f)) == level)   // else don't-care: pass f down
                next = v ? bdd_high(f) : bdd_low(f);
            if (!walk(next, level + 1))
                return false;
        }
        return true;
    };
    walk(f_, 0);
}

std::optional<std::string> pset::any_element() const {
    std::optional<std::string> first;
    for_each_element([&](const std::string &s) {
        first = s;
        return false;
    });
    return first;
}

bool pset::is_old_ctx() const {
    return !context::alive() || ctx_id_ != context::current_id_;
}

void pset::check_context_id(const char *who) const {
    if (is_old_ctx())
        throw std::logic_error(std::string(who) + ": the pset's context is gone");
}

void pset::check_context_id(const pset &o, const char *who) const {
    check_context_id(who);
    o.check_context_id(who);
}

} // namespace spbdd
