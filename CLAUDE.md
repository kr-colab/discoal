This project (discoal) has source code in `src/*`, and is built by running `make` which will produce a binary `build/discoal`.

It uses the library libcyaml to load a YAML. libcyaml and its dependency libyaml are built
as part of the project, with source code included in `extern/libcyaml` and `extern/libyaml`.

To test loading a yaml, run `discoal -Y example_config/new_example.yaml`. This should run
without an error, printing text output to the console.

However, this segfaults. I would like you to help me understand why it segfaults. 
I narrowed down the segfault to the function
`cyaml_load_file` which is called from the function `load_yaml_config` in
`src/core/configInterface_alt.c`.
Without modifying the source code, please help me understand where the error is coming from.

