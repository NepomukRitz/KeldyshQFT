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

# Matsubara or Keldysh formalism
matsubara = False
