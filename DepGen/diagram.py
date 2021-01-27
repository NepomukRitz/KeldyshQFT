from typing import List
import networkx as nx
from general_purpose import combine_into_list, conjugate_keldysh, str_to_list, my_sign

# List of allowed spin combinations
# TODO Automatize generation of this list based on input symmetries
# non-vanishing spin components (s stands for \sigma; b for the conjugate spin \bar{\sigma})
spin_combinations = ["sbsb", "sbbs"] # ["ssss", "sbsb", "sbbs"]       <--- ENTER indices which distinguish the particles

# List of diagrammatic classes. A 1 symbolizes a bare vertex on the corresponding side
# "read as: a 1 on the left means: both legs on the left connect to the same bare vertex (analogous for the right)"
# Hence we identify
# [1 1]: K_1 class
# [0 1]: K_2 class
# [1 0]: K_2' class
# [0 0]: K_3 class
diag_classes = [[1, 1], [0, 1], [1, 0], [0, 0]]


class Diagram:
    def __init__(self, channel: str, diag_class: List[int], indices: List[tuple], aMF = False,
                 sign_omega = 1, sign_nu = 1, sign_nup = 1, exchange_nus = False, complconj = False):
        """Initizalizes a diagram
            --- Params ---
                channel: either 'a', 'p' or 't'
                diag_class: 2-position list with 1's symbolizing bare vertex on corresponding side
                indices: 4-position list with 2-tuples with Keldysh and spin indices for each leg of the diagram
                MF: Matsubara formalism? (False for Keldysh formalism)
                sign_omega: 1 or -1
                sign_nu: 1 or -1
                sign_nup: 1 or -1
                exchange_nus: True or False (are the frequencies nu and nup to be exchanged?)
                complconj: True or False (is the diagrammatic contribution to be complex conjugated?)
            --- Returns ---
                An initialized diagram

            --- To Do's ---
                TODO Implement support of frequencies """
        self.channel = channel
        self.diag_class = diag_class
        self.indices = indices
        self.MF = aMF

        self.sign_omega = sign_omega
        self.sign_nu = sign_nu
        self.sign_nup = sign_nup
        self.exchange_nus = exchange_nus
        self.complconj = complconj

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

    def generate_key(self,with_freqs=False):
        """ Generates unique key for the diagram
            --- Params ---
            with_freqs: bool, True if the diagrammatic contributions should be distinguished by their frequencies
                            False if the diagrammatic contributions are only distinguished by channel, diagrammatic class,
                               and spin indices (in Keldysh: + by Keldysh indices)
            --- Returns ---
                String with the unique id of the diagram """
        if self.diag_class == [0, 0]:
            if self.exchange_nus==False:
                nu  = "v"
                nup = "v\'"
            else:
                nup = "v"
                nu  = "v\'"
        else:
            nu = "v"
            nup = "v"
        freqargs =  my_sign(self.sign_omega) + "w"
        if self.diag_class == [1, 1]:
            dc = "1"
            freqargs += ", , "
        elif self.diag_class == [0, 1]:
            dc = "2"
            freqargs += ","+my_sign(self.sign_nu)+ nu + ", "
        elif self.diag_class == [1, 0]:
            dc = "2\'"
            freqargs += ", ,"+my_sign(self.sign_nup)+ nup
        else:
            dc = "3"
            freqargs += ","+my_sign(self.sign_nu)+ nu +", "+my_sign(self.sign_nup)+ nup

        subscript = dc+", "+self.get_spin_indices()
        superscript = "{}".format(self.channel)
        if not self.MF:
            superscript += ",{}".format(self.get_keldysh_indices())
        key = "$K_{}^{}".format("{"+subscript+"}", "{"+superscript+"}")
        if with_freqs:
            key += " ({})".format(freqargs)
        key += "$"
        #if self.complconj:
        #    key += "^*"

        return key

    def parity_group(self):
        """ Generates parity group respected by the diagram
            --- Returns ---
                List of ParityTrafo objects """
        completed_group = []
        if not self.MF:
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
e = Diagram('a', [0, 1], [('a1', 's'), ('a2', 'b'), ('a3', 's'), ('a4', 'b')])


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
        return Diagram(diagram.channel, diagram.diag_class, new_indices,diagram.MF)

def generate_frequencies(diag_class: List[int]):
    """
    This function generates frequencies corresponding to input diag_class
    :param diag_class: list of integers determining the diagrammatic class
    :return: List of lists of int
    """
    # Frequencies handlers
    no_dep = [+1]
    pm = [+1, -1]

    freqs = [pm]    # All diag classes depend on plus/minus omega

    # Dependency on nu and nu'
    for i in range(2):
        if diag_class[i]:
            freqs.append(no_dep)
        else:
            freqs.append(pm)
    if diag_class == [0, 0]:
        freqs.append([False, True])
    else:
        freqs.append([False])



    return freqs

def generate_diagrams(MF=False, with_freqs=False):
    """ Function that generates the complete set of diagrams
    --- Params ---
        MF:         bool, True for Matsubara formalism, False for Keldysh formalism
        with_freqs: bool, True if the diagrammatic contributions should be distinguished by their frequencies
                          False if the diagrammatic contributions are only distinguished by channel, diagrammatic class,
                             and spin indices (in Keldysh: + by Keldysh indices)
    --- Returns ---
    List of all possible diagrams """
    if MF:
        KeldyshComponents = 1
    else:
        KeldyshComponents = 16
    result = []
    result_keys = []
    for diag_class in diag_classes[0:4]:
        for spin in spin_combinations:
            for channel in ['a', 'p', 't']:
                for iK in range(KeldyshComponents):
                    indices = []
                    for n in range(4, 0, -1):
                        # Generate adequate Keldysh and spin index combination for indices
                        indices.append((str(int(iK % (2 ** n) / (2 ** (n - 1))) + 1), spin[-n]))
                    freqs = generate_frequencies(diag_class)
                    for sign_omega in freqs[0]:
                        for sign_nu in freqs[1]:
                            for sign_nup in freqs[2]:
                                for exchange_nus in freqs[3]:
                                    diag = Diagram(channel, diag_class, indices, MF, sign_omega, sign_nu, sign_nup, exchange_nus)
                                    if not diag.generate_key(with_freqs) in result_keys:
                                        result.append(diag)
                                        result_keys.append(diag.generate_key(with_freqs))
    return result



def establish_dictionary(diagrams: List[Diagram], with_freqs=False):
    """Establishes dictionary for storage of dependencies
    --- Params ---
        diagrams:   list of Diagram objects
        with_freqs: bool, True if the diagrammatic contributions should be distinguished by their frequencies
                          False if the diagrammatic contributions are only distinguished by channel, diagrammatic class,
                             and spin indices (in Keldysh: + by Keldysh indices)
    --- Returns ---
        A dictionary with empty lists as values for each diagram """
    d = {}
    for diag in diagrams:
        d[diag.generate_key(with_freqs)] = {}
    return d


def causality_enforcer(diagrams: List[Diagram]):
    for diag in diagrams:
        if diag.get_keldysh_indices() == "2222":
            return True

    return False
