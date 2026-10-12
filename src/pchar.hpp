#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace spbdd {

// Symplectic representation : I = 00, X = 10, Z = 01, Y = 11.
enum class pauli : std::uint8_t { I = 0, X = 1, Z = 2, Y = 3 };

constexpr bool has_x(pauli p) { return static_cast<std::uint8_t>(p) & 1u; }
constexpr bool has_z(pauli p) { return static_cast<std::uint8_t>(p) & 2u; }

constexpr pauli make_pauli(bool x, bool z) {
    return static_cast<pauli>((x ? 1u : 0u) | (z ? 2u : 0u));
}

inline pauli char_to_pauli(char c) {
    switch (c) {
    case 'I': case 'i': return pauli::I;
    case 'X': case 'x': return pauli::X;
    case 'Y': case 'y': return pauli::Y;
    case 'Z': case 'z': return pauli::Z;
    }
    throw std::invalid_argument(std::string("not a Pauli letter: '") + c + "'");
}

inline std::vector<pauli> string_to_pauli(const std::string &s, int n) {
    if (static_cast<int>(s.size()) != n)
        throw std::invalid_argument("\"" + s + "\" is not a " + std::to_string(n) + "-qubit Pauli");
    std::vector<pauli> p;
    p.reserve(s.size());
    for (char c : s) p.push_back(char_to_pauli(c));
    return p;
}

inline std::string pauli_to_string(const std::vector<pauli> &p) {
    static const char letter[] = {'I', 'X', 'Z', 'Y'};
    std::string s;
    s.reserve(p.size());
    for (pauli q : p) s.push_back(letter[static_cast<int>(q)]);
    return s;
}

} // namespace spbdd
