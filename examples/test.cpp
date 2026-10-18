#include "spbdd/spbdd.hpp"

#include <iostream>

int main() {
    spbdd::context ctx(2);
    spbdd::pset s = ctx.from_list({"XI", "ZZ"});
    std::cout << (s * s).size() << std::endl;   // stdout: 2 ({II, YZ})
}