#include "spbdd/code.hpp"
#include "freeze_order.hpp"
#include "gf2.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace spbdd {

namespace { // helpers

// The point whose variables first, first+1, ... take the values v.
bdd point(const bits &v, int first) {
    bdd c = bddtrue;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const int t = first + static_cast<int>(i);
        c &= v[i] ? bdd_ithvar(t) : bdd_nithvar(t);
    }
    return c;
}

// Values of the variables first .. last-1 in one satisfying assignment of f.
bits model(const bdd &f, int first, int last) {
    const bdd one = bdd_fullsatone(f);
    bits out;
    for (int t = first; t < last; ++t)
        out.push_back((one & bdd_ithvar(t)) != bddfalse);
    return out;
}

} // namespace

code::code(const context &ctx, const std::vector<std::string> &gens)
    : ctx_id_(context::current_id_), n_(ctx.num_qubits()) {

    for (const std::string &s : gens)
        g_.push_back(to_bits(s, n_));

    for (std::size_t i = 0; i < g_.size(); ++i)
        for (std::size_t j = i + 1; j < g_.size(); ++j)
            if (symplectic_inner_product(g_[i], g_[j]))
                throw std::invalid_argument("code: generators " + std::to_string(i) +
                                            " and " + std::to_string(j) + " do not commute");

    // use rank-nullity thm. to check linear independency of generators
    if (null_space(g_, 2 * n_).size() + g_.size() != static_cast<std::size_t>(2 * n_))
        throw std::invalid_argument("code: generators are not independent");
    k_ = n_ - static_cast<int>(g_.size());

    std::vector<bits> rest = null_space(g_, 2 * n_);

    while (!rest.empty()) {
        const bits u = rest.front();
        rest.erase(rest.begin());

        std::size_t p = 0;
        while (p < rest.size() && !symplectic_inner_product(u, rest[p])) ++p;
        if (p == rest.size()) continue;

        const bits w = rest[p];
        rest.erase(rest.begin() + p); 
        for (bits &v : rest) {
            const bool with_w = symplectic_inner_product(v, w), with_u = symplectic_inner_product(v, u);
            for (std::size_t i = 0; i < v.size(); ++i)
                v[i] = v[i] ^ (with_w && u[i]) ^ (with_u && w[i]);
        }
        lx_.push_back(u);
        lz_.push_back(w);
    }

    if (static_cast<int>(lx_.size()) != k_)
        throw std::logic_error("code: found " + std::to_string(lx_.size()) +
                               " logical pairs, expected " + std::to_string(k_));

    // shift_table_[v]: the bits of the label (syndrome, logical component) that variable v flips.
    for (int v = 0; v < 2 * n_; ++v) {
        bits unit(2 * n_, false);
        unit[v] = true;
        bits label = syndrome_of(unit);
        const bits logical = logical_component_of(unit);
        label.insert(label.end(), logical.begin(), logical.end());
        shift_table_.push_back(support(label));
    }
}

std::vector<std::string> code::stabilizer_generators() const {
    return bits_to_strings(g_);
}

std::vector<std::string> code::logical_X_operators() const {
    return bits_to_strings(lx_);
}

std::vector<std::string> code::logical_Z_operators() const {
    return bits_to_strings(lz_);
}

pset code::S() const {
    check_context_id("S");
    return context::current().generated_by(stabilizer_generators());
}

pset code::N() const {
    check_context_id("N");
    return context::current().commute_with(stabilizer_generators());
}

pset code::NmS() const {
    return N() - S();
}

std::vector<bool> code::syndrome_of(const std::string &e) const {
    return syndrome_of(to_bits(e, n_));
}

std::vector<bool> code::syndrome_of(const std::vector<bool> &e) const {
    if (e.size() != static_cast<std::size_t>(2 * n_))
        throw std::invalid_argument("code: expected " + std::to_string(2 * n_) + " bits");
    std::vector<bool> out;
    for (const bits &g : g_)
        out.push_back(symplectic_inner_product(g, e));
    return out;
}

std::vector<bool> code::logical_component_of(const std::string &e) const {
    return logical_component_of(to_bits(e, n_));
}

std::vector<bool> code::logical_component_of(const std::vector<bool> &e) const {
    if (e.size() != static_cast<std::size_t>(2 * n_))
        throw std::invalid_argument("code: expected " + std::to_string(2 * n_) + " bits");
    std::vector<bool> out;
    for (std::size_t j = 0; j < lx_.size(); ++j) {
        out.push_back(symplectic_inner_product(lz_[j], e));
        out.push_back(symplectic_inner_product(lx_[j], e));
    }
    return out;
}

pset code::by_syndrome(const std::vector<bool> &syndrome) const {
    check_context_id("by_syndrome");
    if (syndrome.size() != g_.size())
        throw std::invalid_argument("code: expected " + std::to_string(g_.size()) + " bits");
    const context &ctx = context::current();
    pset f = ctx.all();
    for (std::size_t i = 0; i < g_.size(); ++i) {
        const pset commute = ctx.commute_with(bits_to_string(g_[i]));
        f &= syndrome[i] ? ~commute : commute;
    }
    return f;
}

pset code::by_logical_component(const std::vector<bool> &component) const {
    check_context_id("by_logical_component");
    if (component.size() != 2 * lx_.size())
        throw std::invalid_argument("code: expected " + std::to_string(2 * lx_.size()) + " bits");
    const context &ctx = context::current();
    pset f = ctx.all();
    for (std::size_t j = 0; j < lx_.size(); ++j) {
        const pset commute_z = ctx.commute_with(bits_to_string(lz_[j]));
        const pset commute_x = ctx.commute_with(bits_to_string(lx_[j]));
        f &= component[2 * j] ? ~commute_z : commute_z;
        f &= component[2 * j + 1] ? ~commute_x : commute_x;
    }
    return f;
}

code::violation code::find_qec_cond_violation(const pset &errors, fqcv_method m) const {
    check_context_id("find_qec_cond_violation");
    errors.check_context_id("find_qec_cond_violation");

    switch (m) {
        case fqcv_method::bruteforce: return fqcv_bruteforce(errors);
        case fqcv_method::table:      return fqcv_table(errors);
        case fqcv_method::image:      return fqcv_image(errors);
        case fqcv_method::square:     return fqcv_square(errors);
    }
    throw std::invalid_argument("find_qec_cond_violation: unknown method");
}

code::violation code::fqcv_bruteforce(const pset &errors) const {
    std::vector<bits> all;
    errors.for_each_element([&](const std::string &p) {
        all.push_back(to_bits(p, n_));
        return true;
    });

    const bits zero_syndrome(r(), false), zero_component(2 * k_, false);
    for (std::size_t i = 0; i < all.size(); ++i)
        for (std::size_t j = i + 1; j < all.size(); ++j) {
            const bits product = multiply(all[i], all[j]);
            if (syndrome_of(product) == zero_syndrome && logical_component_of(product) != zero_component)
                return std::make_pair(bits_to_string(all[i]), bits_to_string(all[j]));
        }
    return std::nullopt;
}

code::violation code::fqcv_table(const pset &errors) const {
    std::unordered_map<bits, std::pair<bits, std::string>> table;   // syndrome -> (logical component, error)
    violation found;
    errors.for_each_element([&](const std::string &p) {
        const bits e = to_bits(p, n_), component = logical_component_of(e);
        const auto [it, is_new] = table.try_emplace(syndrome_of(e), component, p);
        if (is_new || it->second.first == component)
            return true;
        found = std::make_pair(it->second.second, p);
        return false;
    });
    return found;
}

code::violation code::fqcv_image(const pset &errors) const {
    const bdd &f = errors.f_;
    if (f == bddfalse)
        return std::nullopt;

    const int m = r() + 2 * k_, N = 2 * n_;
    bdd labels;

    { // for freeze_order object's lifetime
        freeze_order freeze;
        const auto level = [N](const bdd &u) { return u == bddtrue ? N : bdd_var2level(bdd_var(u)); };

        // Step 1 : record the nodes to process at each level (a node that skips a level is carried down as-is)
        std::vector<std::map<int, bdd>> todo(N + 1);
        todo[0].emplace(f.id(), f);
        for (int l = 0; l < N; ++l)
            for (const auto &[id, u] : todo[l]) {
                if (level(u) != l) {
                    todo[l + 1].emplace(id, u);
                    continue;
                }
                for (const bdd &child : {bdd_low(u), bdd_high(u)})
                    if (child != bddfalse)
                        todo[l + 1].emplace(child.id(), child);
            }

        // Step 2 : the image of the terminal is the all-zero label
        std::map<int, bdd> next;    // images from level l+1
        for (const auto &[id, u] : todo[N])
            next.emplace(id, point(bits(m, false), 0));
        // Step 3 : go up level by level, image(u) = image(low) U (shift ^ image(high))
        for (int l = N - 1; l >= 0; --l) {
            const std::vector<int> &shift = shift_table_[bdd_level2var(l)];
            const auto image = [&next](const bdd &u) { return u == bddfalse ? bddfalse : next.at(u.id()); };
            std::map<int, bdd> cur;
            for (const auto &[id, u] : todo[l]) {
                const bool tested = level(u) == l;
                bdd high = image(tested ? bdd_high(u) : u);
                for (int t : shift)
                    high = bdd_compose(high, bdd_nithvar(t), t);    // w ^ high: swap the children of the flipped variables
                cur.emplace(id, image(tested ? bdd_low(u) : u) | high);
            }
            next = std::move(cur);
        }
        labels = next.at(f.id());
    }

    // Step 4 : find a syndrome that has both values at some logical bit
    bdd logical = bddtrue;
    for (int t = r(); t < m; ++t)
        logical &= bdd_ithvar(t);

    for (int t = r(); t < m; ++t) {
        const bdd zero = labels & bdd_nithvar(t), one = labels & bdd_ithvar(t);
        const bdd ambiguous = bdd_exist(zero, logical) & bdd_exist(one, logical);   // syndromes with both values
        if (ambiguous == bddfalse) continue;

        // Step 5 : pick such a syndrome and one error for each side
        const bits s = model(ambiguous, 0, r());
        const auto error_with = [&](const bdd &side) {
            return (errors & by_syndrome(s) & by_logical_component(model(side & point(s, 0), r(), m))).any_element().value();
        };
        return std::make_pair(error_with(zero), error_with(one));
    }
    return std::nullopt;
}

code::violation code::fqcv_square(const pset &errors) const {
    const std::optional<std::string> E = ((errors * errors) & NmS()).any_element();
    if (!E) return std::nullopt;

    /* E ∈ errors * errors, so E = a * b for some a, b ∈ errors, hence b = a * E and a = b * E.
       So a ∈ errors and a ∈ errors * {E}, and the intersection is not empty.
       Any a in it works: a ∈ errors and a = c * E for some c ∈ errors, so a * E = c ∈ errors. */

    const std::string a = (errors & context::current().coset_of(*E, errors)).any_element().value();
    return std::make_pair(a, bits_to_string(multiply(to_bits(a, n_), to_bits(*E, n_))));
}

void code::check_context_id(const char *who) const {
    if (!context::alive() || ctx_id_ != context::current_id_)
        throw std::logic_error(std::string(who) + ": the code's context is gone");
}

} // namespace spbdd
