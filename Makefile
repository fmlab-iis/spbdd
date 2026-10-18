CXX      ?= g++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -O2

BUILD := build
BUDDY ?= buddy
LIB   := $(BUILD)/libspbdd.a
INC   := -Iinclude -I$(BUDDY)/src -I$(BUILD)/buddy

OBJS       := $(patsubst src/%.cpp,$(BUILD)/%.o,$(wildcard src/*.cpp))
BUDDY_C    := bddio bddop bvec cache fdd imatrix kernel pairs prime reorder tree
BUDDY_OBJS := $(BUDDY_C:%=$(BUILD)/buddy/%.o) $(BUILD)/buddy/cppext.o
EX_BINS    := $(patsubst examples/%.cpp,$(BUILD)/%,$(wildcard examples/*.cpp))
TEST_BINS  := $(patsubst test/%.cpp,$(BUILD)/%,$(wildcard test/*.cpp))

.PHONY: all examples check clean_build clean_buddy
all: $(LIB)
examples: $(EX_BINS)
check: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo $$t; ./$$t || exit 1; done

$(LIB): $(OBJS) $(BUDDY_OBJS)
	$(AR) rcs $@ $^

$(BUILD)/%.o: src/%.cpp $(BUDDY)/src/bdd.h
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC) -MMD -MP -c $< -o $@

$(EX_BINS): $(BUILD)/%: examples/%.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LIB) -lm -o $@

$(TEST_BINS): $(BUILD)/%: test/%.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LIB) -lm -o $@

$(BUILD)/buddy/%.o: $(BUDDY)/src/%.c $(BUILD)/buddy/config.h
	@mkdir -p $(@D)
	$(CC) -O2 -fexceptions $(INC) -c $< -o $@

$(BUILD)/buddy/cppext.o: $(BUDDY)/src/cppext.cxx
	@mkdir -p $(@D)
	$(CXX) -O2 $(INC) -c $< -o $@

$(BUILD)/buddy/config.h: | $(BUDDY)/src/bdd.h
	@mkdir -p $(@D)
	printf '#define PACKAGE "buddy"\n#define PACKAGE_VERSION "2.4"\n#define VERSION "2.4"\n#define MAJOR_VERSION 2\n#define MINOR_VERSION 4\n' > $@

$(addprefix $(BUDDY)/src/,$(BUDDY_C:%=%.c) cppext.cxx): $(BUDDY)/src/bdd.h

$(BUDDY)/src/bdd.h:
	git clone --depth 1 https://github.com/utwente-fmt/buddy.git $(BUDDY)

clean_build:
	rm -rf $(BUILD)

clean_buddy:
	rm -rf $(BUDDY) $(BUILD)/buddy

-include $(OBJS:.o=.d)
