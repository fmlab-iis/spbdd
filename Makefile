CXX      ?= g++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -O2

BUILD := build
BUDDY := buddy
LIB   := $(BUILD)/libspbdd.a

BUDDY_C    := bddio bddop bvec cache fdd imatrix kernel pairs prime reorder tree
BUDDY_SRC  := $(BUDDY_C:%=$(BUDDY)/src/%.c) $(BUDDY)/src/cppext.cxx
BUDDY_OBJS := $(BUDDY_C:%=$(BUILD)/buddy/%.o) $(BUILD)/buddy/cppext.o
OBJS       := $(patsubst src/%.cpp,$(BUILD)/%.o,$(wildcard src/*.cpp))
EX_BINS    := $(patsubst examples/%.cpp,$(BUILD)/%,$(wildcard examples/*.cpp))
TEST_BINS  := $(patsubst test/%.cpp,$(BUILD)/%,$(wildcard test/*.cpp))

.PHONY: all examples check clean

# From scratch: clean, clone BuDDy, build BuDDy, build spbdd.
all:
	$(MAKE) clean
	$(MAKE) $(LIB)

examples: $(EX_BINS)

check: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo $$t; ./$$t || exit 1; done

clean:
	rm -rf $(BUILD) $(BUDDY) include/spbdd/bdd.h

$(LIB): $(BUDDY_OBJS) $(OBJS)
	$(AR) rcs $@ $^

# BuDDy
$(BUDDY)/src/bdd.h:
	git clone --depth 1 https://github.com/utwente-fmt/buddy.git $(BUDDY)

$(BUDDY_SRC): $(BUDDY)/src/bdd.h

$(BUILD)/buddy/config.h: $(BUDDY)/src/bdd.h
	@mkdir -p $(@D)
	printf '#define PACKAGE "buddy"\n#define PACKAGE_VERSION "2.4"\n#define VERSION "2.4"\n#define MAJOR_VERSION 2\n#define MINOR_VERSION 4\n' > $@

$(BUILD)/buddy/%.o: $(BUDDY)/src/%.c $(BUILD)/buddy/config.h
	$(CC) -O2 -fexceptions -I$(BUILD)/buddy -c $< -o $@

$(BUILD)/buddy/cppext.o: $(BUDDY)/src/cppext.cxx $(BUILD)/buddy/config.h
	$(CXX) -O2 -I$(BUILD)/buddy -c $< -o $@

# spbdd: the public headers include a copy of BuDDy's bdd.h
include/spbdd/bdd.h: $(BUDDY)/src/bdd.h
	cp $< $@

$(BUILD)/%.o: src/%.cpp include/spbdd/bdd.h
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -Iinclude -MMD -MP -c $< -o $@

LINK = $(CXX) $(CXXFLAGS) -Iinclude $< $(LIB) -lm -o $@

$(EX_BINS): $(BUILD)/%: examples/%.cpp $(LIB)
	$(LINK)

$(TEST_BINS): $(BUILD)/%: test/%.cpp $(LIB)
	$(LINK)

-include $(OBJS:.o=.d)
