#pragma once

#include "spbdd/context.hpp"

namespace spbdd {

// Keeps the variable order fixed while object alive.
class freeze_order {
public:
    freeze_order() { bdd_disable_reorder(); }
    ~freeze_order() {
        if (context::current().dynamic_reordering())
            bdd_enable_reorder();
    }
};

} // namespace spbdd
