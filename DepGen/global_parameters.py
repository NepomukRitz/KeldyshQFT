# List of allowed spin combinations
# TODO Automatize generation of this list based on input symmetries
spin_combinations = ["ssss", "sbsb", "sbbs"]

# List of diagrammatic classes. A 1 symbolizes a bare vertex on the corresponding side
diag_classes = [[1, 1],     # K1
                [0, 1],     # K2
                [1, 0],     # K2'
                [0, 0]]     # K3

# List of interaction channels
channels = ['a', 'p', 't']

# T_C prefactor change
tc_sign_change = {'1111': -1, '1112': +1, '1121': +1, '1122': -1,
                  '1211': +1, '1212': -1, '1221': -1, '1222': +1,
                  '2111': +1, '2112': -1, '2121': -1, '2122': +1,
                  '2211': -1, '2212': +1, '2221': +1, '2222': -1}

# Matsubara or Keldysh formalism
matsubara = False

# With/without frequencies
with_freqs = True

# If with frequencies, define parametrization
# Supports "bosonic": w, f, f'. Anything else defaults to standard: w, v, v'
param = "standard"

# Plot orbits or not
plot_orbits = False
