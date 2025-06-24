# About
demes-c is a C library for parsing [Demes demographic models](https://popsim-consortium.github.io/demes-spec-docs/) using libyaml.

# Installation
Clone the repository and type `make`.
```
git clone --recurse-submodules https://github.com/grahamgower/demes-c.git
cd demes-c
make
```
You may need to first install libyaml and/or edit the `Makefile` to provide its location.
This will build a static library `libdemes.a` and a `resolver` binary
that can resolve Demes YAML files.

# API
See `demes.h` and `resolver.c` to understand the interface.

# Gotchas
* metadata is accepted in input files but not available in the API (issue #13)
* errors are printed to stderr instead of being returned via the API (issue #18).

# Test suite
The demes-c tests rely on the demes-spec git submodule. If you didn't clone
the repository with `--recurse-submodules` then you can get the submodule with
```
git submodule init
git submodule update
```

There are three test-related targets in the toplevel Makefile,
which are called during the continuous integration github action
that gets run when a pull request is opened.

* `make test` checks that all valid test cases in the demes-spec repo
  can be loaded and that all invalid test cases are correctly rejected.
* `make memcheck` runs the resolver under valgrind to check for memory errors.
  All error paths are checked by mocking various libc and libyaml functions
  to simulate failures. This can take a few minutes to run.
* `make pytest` compares the output of the C resolver to the
  output of the reference implementation in the demes-spec repo.
  This needs Python >= 3.10 with the packages listed in `tests/requirements.txt`.
