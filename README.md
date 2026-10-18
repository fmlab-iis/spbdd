# spbdd
A BDD-based C++ library for set operations on Pauli operators

## Usage

### 1. Clone & build

```sh
git clone https://github.com/fmlab-iis/spbdd
cd spbdd
make
```

This step produces the static library `build/libspbdd.a`. Note that this step also automatically fetches [BuDDy](https://github.com/utwente-fmt/buddy), the BDD backend used in this library.

### 2. The headers & library

To use this library in your C++ program, you only need :

1. `include/spbdd/` folder, which contains the library headers.
2. `build/libspbdd.a`, which is the static library archive.

Copy them next to your `prog.cpp`:

```
.
├── prog.cpp
├── spbdd/        <- copied from include/spbdd/
└── libspbdd.a    <- copied from build/libspbdd.a
```

### 3. Write & compile your program

Here is an example program `prog.cpp`:

```cpp
#include "spbdd/spbdd.hpp" // the only header you need
#include <iostream>

int main() {
    spbdd::context ctx(2);
    spbdd::pset s = ctx.from_list({"XI", "ZZ"});
    std::cout << (s * s).size() << std::endl;   // stdout: 2 ({II, YZ})
}
```
To compile, you only need to link `libspbdd.a`. Note that this library requires at least C++17.

```sh
g++ -std=c++17 prog.cpp libspbdd.a -o prog
```

## Examples

```sh
make examples
build/
```