from typing import List


class Diagram:
    def __init__(self, channel: str, diag_class: List[int], indices: List[tuple]):
        self.channel = channel
        self.diag_class = diag_class
        self.indices = indices

    def __str__(self):
        return self.generate_key()

    """Returns either Keldysh indices (i=0) or spin indices (i=1)"""
    def get_general_indices(self, i):
        indices = ''
        for index in self.indices:
            indices += index[i]
        return indices

    def get_keldysh_indices(self):
        return self.get_general_indices(0)

    def get_spin_indices(self):
        return self.get_general_indices(1)

    def set_keldysh_indices(self, new_keldysh_indices: List[str]):
        new_indices = []
        for i in range(len(self.indices)):
            new_indices.append((new_keldysh_indices[i], self.indices[i][1]))
        self.indices = new_indices

    def set_spin_indices(self, new_spin_indices: List[str]):
        new_indices = []
        for i in range(len(self.indices)):
            new_indices.append((self.indices[i][0], new_spin_indices[i]))
        self.indices = new_indices

    def generate_key(self):
        if self.diag_class == [1, 1]:
            dc = "1"
        elif self.diag_class == [0, 1]:
            dc = "2"
        elif self.diag_class == [1, 0]:
            dc = "2'"
        else:
            dc = "3"

        return "K_{}^{} spin_comp = {} kel_comp = {}".format(dc, self.channel, self.get_spin_indices(), self.get_keldysh_indices())

    def parity_group(self):
        completed_group = []
        L, R = False, False
        if self.diag_class[0] == 1:
            L = True
            completed_group.append(ParityTrafo(self.channel, 'L'))
        if self.diag_class[1] == 1:
            R = True
            completed_group.append(ParityTrafo(self.channel, 'R'))
        if L and R:
            completed_group.append(ParityTrafo(self.channel, '.'))
        return completed_group


e = Diagram('.', [-1, -1], [('a1', 's1'), ('a2', 's2'), ('a3', 's3'), ('a4', 's4')])


def combine_into_list(keldysh_indices, spin_indices):
    new_indices = []
    for i in range(len(keldysh_indices)):
        new_indices.append((keldysh_indices[i], spin_indices[i]))
    return new_indices


class ParityTrafo:
    def __init__(self, channel, side):
        self.channel = channel
        self.side = side

    def __str__(self):
        return "P^{}_{}".format(self.channel, self.side)

    def T(self, diagram: Diagram):
        keldysh_indices = str_to_list(diagram.get_keldysh_indices())
        spin_indices = diagram.get_spin_indices()
        positions = []
        if self.side == '.':
            positions = [0, 1, 2, 3]
        else:
            if self.channel == 'a':
                if self.side == 'L':
                    positions = [0, 3]
                else:   # side = 'R'
                    positions = [1, 2]
            elif self.channel == 'p':
                if self.side == 'L':
                    positions = [0, 1]
                else:   # side = 'R'
                    positions = [2, 3]
            elif self.channel == 't':
                if self.side == 'L':
                    positions = [1, 3]
                else:   # side = 'R'
                    positions = [0, 2]

        new_keldysh_indices = "".join(conjugate(keldysh_indices, positions))
        new_indices = combine_into_list(new_keldysh_indices, spin_indices)
        return Diagram(diagram.channel, diagram.diag_class, new_indices)


def conjugate_spin(spin_indices):
    conjugated_spins = ''
    for spin in spin_indices:
        if spin == 's':
            conjugated_spins += 'b'
        elif spin == 'b':
            conjugated_spins += 's'
    return conjugated_spins if len(conjugated_spins)>0 else spin_indices


class Trafo:
    def __init__(self, i):
        self.i = i
        self.id = self.T(e).indices

    def __str__(self):
        return "T" + str(self.i)

    def T(self, diagram: Diagram):
        new_channel = self.transform_channel(diagram.channel)
        new_diag_class = self.transform_diag_class(diagram.diag_class, diagram.channel)
        new_indices = self.transform_indices(diagram.indices)
        return Diagram(new_channel, new_diag_class, new_indices)

    def transform_channel(self, channel):
        if self.i == 1 or self.i == 2:
            if channel == 'a':
                return 't'
            elif channel == 't':
                return 'a'
        return channel

    def transform_diag_class(self, diag_class, channel):
        if self.i == 1 or self.i == 3:
            if channel == 'a' or channel == 't':
                return [diag_class[1], diag_class[0]]
        if self.i == 4:
            if channel == 'a' or channel == 'p':
                return [diag_class[1], diag_class[0]]
        return [diag_class[0], diag_class[1]]

    def transform_indices(self, indices: List[tuple]):
        new_indices = indices.copy()

        if self.i == 1 or self.i == 3:
            new_indices[-1] = indices[-2]
            new_indices[-2] = indices[-1]

        if self.i == 2 or self.i == 3:
            new_indices[1] = indices[0]
            new_indices[0] = indices[1]

        if self.i == 4:
            new_indices[-2:] = indices[:2]
            new_indices[:2] = indices[-2:]

        if self.i == 5:
            keldysh_indices = [i[0] for i in indices]
            spin_indices = [j[1] for j in indices]
            new_spin_indices = conjugate_spin(spin_indices)

            new_indices = []
            for ii in range(len(indices)):
                new_indices.append((keldysh_indices[ii], new_spin_indices[ii]))

        return new_indices


class CompositeTrafo:

    def __init__(self, trafo1: Trafo, trafo2: Trafo):
        self.trafo1 = trafo1
        self.trafo2 = trafo2
        self.id = trafo1.T(trafo2.T(e)).indices

    def __str__(self):
        return self.trafo1.__str__() + self.trafo2.__str__()

    def T(self, diagram: Diagram):
        return self.trafo1.T(self.trafo2.T(diagram))


def conjugate(indices, positions):
    for pos in positions:
        if indices[pos] == '1':
            indices[pos] = '2'
        else:
            indices[pos] = '1'
    return indices


def list_to_str(list_of_values):
    s = ""
    for num in list_of_values:
        s += str(num)
    return s


def str_to_list(word):
    ans = []
    for char in word:
        ans.append(char)
    return ans


def generate_full_group(group: list):
    ids = [trafo.id for trafo in group]

    done = False
    while not done:
        n = len(group)
        to_add = []
        for trafo1 in group:
            for trafo2 in group:
                new_trafo = CompositeTrafo(trafo1, trafo2)

                if new_trafo.id not in ids:
                    to_add.append(new_trafo)
                    ids.append(new_trafo.id)

        group.extend(to_add)

        if len(group) == n:
            done = True

    # spin_trafo = []
    # for trafo in group:
    #     if not (trafo == group[0] or trafo == group[5]):
    #         spin_trafo.append(CompositeTrafo(trafo, group[5]))      # group[5] = T_s
    #
    # group.extend(spin_trafo)

    return group
