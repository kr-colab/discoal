# this works
discoal 20 2 10000 -d 12345 10102 -t 0.1 -r 0.1 -p 3 10 5 5 -M 0.1 -N 1

# this should be equiv but has a table sorting error
build/discoal -Y config_examples/new_example.yaml

# the YAML read-in works if you set to be a single population (no migration)
# so the error seems to be with population structure
