from typing import List
from general_purpose import combine_into_list, conjugate_keldysh, str_to_list

# List of allowed spin combinations
# TODO Automatize generation of this list based on input symmetries
spin_combinations = ["ssss", "sbsb", "sbbs"]

# List of diagrammatic classes. A 1 symbolizes a bare vertex on the corresponding side
diag_classes = [[1, 1], [0, 1], [1, 0], [0, 0]]

# List of interaction channels
channels = ['a', 'p', 't']


class Diagram:
    # TODO Implement also prefactor ?
    def __init__(self, channel: str, diag_class: List[int], indices: List[tuple], freqs: List[int], mf=False):
        """Initizalizes a diagram
            --- Params ---
                channel: either 'a', 'p' or 't'
                diag_class: 2-position list with 1's symbolizing bare vertex on corresponding side
                indices: 4-position list with 2-tuples with Keldysh and spin indices for each leg of the diagram
                MF: Matsubara formalism? (False for Keldysh formalism)
            --- Returns ---
                An initialized diagram
             """
        self.channel = channel
        self.diag_class = diag_class
        self.indices = indices
        self.mf = mf
        self.freqs = freqs

    def __str__(self):
        """ Print formatting function. Currently set to print key """
        return self.generate_key()

    def get_general_indices(self, i):
        """Returns either Keldysh indices (i=0) or spin indices (i=1) """
        indices = ''
        for index in self.indices:
            indices += index[i]
        return indices

    def get_keldysh_indices(self):
        """ Returns Keldysh indices """
        return self.get_general_indices(0)

    def get_spin_indices(self):
        """ Returns spin indices """
        return self.get_general_indices(1)

# Code commented out since it is currently not needed. Nt erased since it may be useful later
    # def set_keldysh_indices(self, new_keldysh_indices: List[str]):
    #     """ Sets Keldysh indices to input
    #     --- Params ---
    #         new_keldysh_indices: 4-position list """
    #     new_indices = []
    #     for i in range(len(self.indices)):
    #         new_indices.append((new_keldysh_indices[i], self.indices[i][1]))
    #     self.indices = new_indices

    # def set_spin_indices(self, new_spin_indices: List[str]):
    #     """ Sets spin indices to input
    #         --- Params ---
    #             new_spin_indices: 4-position list """
    #     new_indices = []
    #     for i in range(len(self.indices)):
    #         new_indices.append((self.indices[i][0], new_spin_indices[i]))
    #     self.indices = new_indices

    def generate_key(self):
        """ Generates unique key for the diagram
            --- Returns ---
                String with the unique id of the diagram """
        if self.diag_class == [1, 1]:
            dc = "1"
        elif self.diag_class == [0, 1]:
            dc = "2"
        elif self.diag_class == [1, 0]:
            dc = "2'"
        else:
            dc = "3"
        if self.mf:
            return f"K_{dc}^{self.channel} spin={self.get_spin_indices} freqs={self.freqs}"
        else:
            return f"K_{dc}^{self.channel} spin={self.get_spin_indices()} kel={self.get_keldysh_indices()} freqs={self.freqs}"

    def parity_group(self):
        """ Generates parity group respected by the diagram
            --- Returns ---
                List of ParityTrafo objects """
        completed_group = []

        if not self.mf:  # Parity group exists only in Keldysh formalism!
            L, R = False, False
            if self.diag_class[0] == 1:
                L = True
                completed_group.append(ParityTrafo(self.channel, 'L'))
            if self. diag_class[1] == 1:
                R = True
                completed_group.append(ParityTrafo(self.channel, 'R'))
            if L and R:
                completed_group.append(ParityTrafo(self.channel, '.'))
        return completed_group


# Define identity diagram. Useful to determine id of T transformations
e = Diagram('.', [-1, -1], [('a1', 's'), ('a2', 's'), ('a3', 's'), ('a4', 's')], [])


class ParityTrafo:
    def __init__(self, channel, side):
        """Initializes a ParityTrafo object
            --- Params ---
                channel: either 'a', 'p', or 't'
                side: either 'R', 'L', or '.'
            --- Returns ---
                An initialized ParityTrafo """
        self.channel = channel
        self.side = side

    def __str__(self):
        """ Print formatting. Prints Trafo as T^r_{side} """
        return f"P^{self.channel}_{self.side}"

    def T(self, diagram: Diagram):
        """Transformation operation. Uses same name as T transformations
            --- Params ---
                diagram: diagram object on which the transformation acts
            --- Returns
                New Diagram object, result of applying the symmetry """

        # Transforms the keldysh indices of the diagram into a list of strings
        keldysh_indices = str_to_list(diagram.get_keldysh_indices())

        # Determine the positions of the list of Keldysh indices to act on
        positions = []
        if self.side == '.':
            positions = [0, 1, 2, 3]
        else:
            if self.channel == 'a':
                if self.side == 'L':
                    positions = [0, 3]
                else:  # side = 'R'
                    positions = [1, 2]
            elif self.channel == 'p':
                if self.side == 'L':
                    positions = [0, 1]
                else:  # side = 'R'
                    positions = [2, 3]
            elif self.channel == 't':
                if self.side == 'L':
                    positions = [1, 3]
                else:  # side = 'R'
                    positions = [0, 2]

        # Set up list of new, conjugated, Keldysh indices
        new_keldysh_indices = "".join(conjugate_keldysh(keldysh_indices, positions))

        # Join new Keldysh indices with old, untouched, spin indices
        new_indices = combine_into_list(new_keldysh_indices, diagram.get_spin_indices())

        # Return new, transformed diagram
        return Diagram(diagram.channel, diagram.diag_class, new_indices, diagram.freqs, diagram.mf)


def generate_frequencies(diag_class: List[int]):
    """
    This function generates frequencies corresponding to input diag_class
    :param diag_class: list of integers determining the diagrammatic class
    :return: List of lists of int
    """
    # Frequencies handlers
    no_dep = [0, 0]
    pm = [+1, -1]

    freqs = [pm]    # All diag classes depend on plus/minus omega

    # Dependency on nu and nu'
    for i in range(2):
        if diag_class[i]:
            freqs.append(no_dep)
        else:
            freqs.append(pm)

    result = set()

    for w in freqs[0]:
        for v1 in freqs[1]:
            for v2 in freqs[2]:
                result.add((w, v1, v2))

    return result


def generate_diagrams(mf=False):
    """ Function that generates the complete set of diagrams
    --- Returns ---
    List of all possible diagrams """

    # Assignment of the number of Keldysh components
    if mf:
        KeldyshComponents = 1
    else:
        KeldyshComponents = 16

    # Container for all possible diagrams
    result = []

    # Start generation of diagrams
    for diag_class in diag_classes:

        # Diagrammatic class determines frequencies.
        freqs = generate_frequencies(diag_class)

        for spin in spin_combinations:
            for channel in channels:
                for iK in range(KeldyshComponents):
                    indices = []
                    # Generate adequate Keldysh and spin index combination for indices
                    for n in range(4, 0, -1):
                        indices.append((str(int(iK % (2 ** n) / (2 ** (n - 1))) + 1), spin[-n]))

                    for freq in freqs:
                        result.append(Diagram(channel, diag_class, indices, list(freq), mf))
    return result


def establish_dictionary(diagrams: List[Diagram]):
    """Establishes dictionary for storage of dependencies
    --- Params ---
        diagrams: list of Diagram objects
    --- Returns ---
        A dictionary with empty lists as values for each diagram """
    d = {}
    for diag in diagrams:
        d[diag.generate_key()] = []
    return d


def causality_enforcer(diagrams: List[Diagram]):
    for diag in diagrams:
        if diag.get_keldysh_indices() == "2222":
            return True

    return False
