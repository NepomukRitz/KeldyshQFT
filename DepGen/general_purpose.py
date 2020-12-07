
def combine_into_list(keldysh_indices, spin_indices):
    """ Combines a list of Keldysh and a list of spin indices into an indices list
    --- Params ---
        keldysh_indices: 4-position list
        spin_indices: 4-position list
    --- Returns ---
        4-position list of 2-position tuples. The first entries of the tuples are the Keldysh indices
    """
    new_indices = []
    for i in range(len(keldysh_indices)):
        new_indices.append((keldysh_indices[i], spin_indices[i]))
    return new_indices


def conjugate_spin(spin_indices):
    """ Conjugates the spins
    --- Params ---
        spin_indices: string with the spin indices, each character is either 's' or 'b'
    --- Returns ---
        String with the conjugated spins
    """
    conjugated_spins = ''
    for spin in spin_indices:
        if spin == 's':
            conjugated_spins += 'b'
        elif spin == 'b':
            conjugated_spins += 's'
    return conjugated_spins


def conjugate_keldysh(indices, positions):
    """ Conjugates Keldysh indices
    --- Params ---
        indices: list of Keldysh indices, each entry is either '1' or '2'
        positions: 4-position list
    --- Returns ---
        A list with the conjugated Keldysh indices at the indicated positions """
    for pos in positions:
        if indices[pos] == '1':
            indices[pos] = '2'
        else:
            indices[pos] = '1'
    return indices


def list_to_str(list_of_values):
    """ Turns input list into a string """
    s = ""
    for num in list_of_values:
        s += str(num)
    return s


def str_to_list(word):
    """ Turns input string into a list """
    ans = []
    for char in word:
        ans.append(char)
    return ans
