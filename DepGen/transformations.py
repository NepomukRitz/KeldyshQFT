from typing import List
from diagram import Diagram, e, establish_dictionary, spin_combinations, causality_enforcer
import networkx as nx
from general_purpose import conjugate_spin


class Trafo:
    def __init__(self, i):
        """ Initialized a Trafo object
        --- Params ---
            i: integer. 0, 1, 2, 3, 4 and 5 and 6 supported
        --- Returns ---
            Initialized Trafo object with a calculated id """
        self.i = i
        self.channel = self.T(e).channel
        self.id = self.T(e).indices

        self.diag_class = self.T(e).diag_class
        self.sign_omega = self.T(e).sign_omega
        self.sign_nu = self.T(e).sign_nu
        self.sign_nup = self.T(e).sign_nup
        self.exchange_nus = self.T(e).exchange_nus

    def __str__(self):
        """ Print formatting. Prints Trafo as Ti """
        return "T" + str(self.i)

    def T(self, diagram: Diagram):
        """ Transformation operation
        --- Params ---
            diagram: Diagram object to be operated on
        --- Returns ---
            New Diagram object, the result of the transformation operation on the input """
        new_channel = self.transform_channel(diagram.channel)
        new_diag_class = self.transform_diag_class(diagram.diag_class, diagram.channel)
        new_indices = self.transform_indices(diagram.indices)
        new_sign_omega, new_sign_nu, new_sign_nup, new_exchange_nus = self.transform_freqs(diagram)
        return Diagram(new_channel, new_diag_class, new_indices, diagram.MF,
                       new_sign_omega, new_sign_nu, new_sign_nup, new_exchange_nus)

    def transform_channel(self, channel):
        """Transforms the channel
        --- Params ---
            channel: 'a', 'p', or 't'
        --- Retuns ---
            The channel that results when T operates """
        if self.i == 1 or self.i == 2:
            if channel == 'a':
                return 't'
            elif channel == 't':
                return 'a'
        return channel

    def transform_diag_class(self, diag_class, channel):
        """ Transforms the diagrammatic class
        --- Params ---
            diag_class: 2-position list
            channel: 'a', 'p', or 't'
        --- Returns ---
            2-position list encoding the new transformed diagrammatic class """
        if self.i == 1 or self.i == 3:
            if channel == 'a' or channel == 't':
                return [diag_class[1], diag_class[0]]
        if self.i == 4 or self.i == 6:
            if channel == 'a' or channel == 'p':
                return [diag_class[1], diag_class[0]]
        return [diag_class[0], diag_class[1]]

    def transform_indices(self, indices: List[tuple]):
        """ Transforms the indices
        --- Params ---
            indices: 4-position list of 2-position tuples
        --- Returns ---
            4-position list with the new, transformed, indices """

        # Copy old indices into new object
        new_indices = indices.copy()

        # Change indices wherever needed, according to i
        if self.i == 1 or self.i == 3:
            new_indices[-1] = indices[-2]
            new_indices[-2] = indices[-1]

        if self.i == 2 or self.i == 3:
            new_indices[1] = indices[0]
            new_indices[0] = indices[1]

        if self.i == 4 or self.i == 6:
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

    def transform_freqs(self, diagram: Diagram):
        """ Transform frequency arguments (sign)
        --- Params ---

        --- Returns ---
        """
        diag_class = diagram.diag_class
        channel = diagram.channel
        new_sign_omega = diagram.sign_omega
        new_sign_nu = diagram.sign_nu
        new_sign_nup = diagram.sign_nup
        new_exchange_nus = diagram.exchange_nus

        if self.i == 1:
            if channel == 'p':
                if new_exchange_nus:
                     new_sign_nu *= -1
                else:
                    new_sign_nup *= -1
            else:
                new_sign_omega *= -1
                new_exchange_nus = not new_exchange_nus

        if self.i == 2 and  channel == 'p':
            if new_exchange_nus:
                new_sign_nup *= -1
            else:
                new_sign_nu *= -1

        if self.i == 3:
            if channel == 'p':
                new_sign_nu *= -1
                new_sign_nup *= -1
            else:
                new_sign_omega *= -1
                new_exchange_nus = not new_exchange_nus
        if self.i == 4:
            if channel == 't':
                new_sign_omega *= -1
            else:
                new_exchange_nus = not new_exchange_nus
            if diagram.MF:
                new_sign_omega *= -1
                new_sign_nu *= -1
                new_sign_nup *= -1
        if self.i == 6:
            if channel == 'a' or channel == 'p':
                new_sign_omega *= -1
                new_sign_nu *= -1
                new_sign_nup *= -1
            else:
                new_sign_nu *= -1
                new_sign_nup *= -1
        if self.i == 7:
            new_sign_omega *= -1
            new_sign_nu *= -1
            new_sign_nup *= -1


        return new_sign_omega, new_sign_nu, new_sign_nup, new_exchange_nus


class CompositeTrafo:
    def __init__(self, trafo1: Trafo, trafo2: Trafo):
        """ Initializes a CompositeTrafo, which corresponds to a multiplication of two Trafos
        --- Params ---
            trafo1: first Trafo/CompositeTrafo object
            trafo2: second Trafo/CompositeTrafo object
        --- Returns ---
            Initialized CompositeTrafo object, with id set by the action of trafo2 and then of trafo 1 """
        self.trafo1 = trafo1
        self.trafo2 = trafo2
        self.channel = trafo1.T(trafo2.T(e)).channel
        self.diag_class = trafo1.T(trafo2.T(e)).diag_class
        self.id = trafo1.T(trafo2.T(e)).indices

        self.sign_omega = trafo1.T(trafo2.T(e)).sign_omega
        self.sign_nu = trafo1.T(trafo2.T(e)).sign_nu
        self.sign_nup = trafo1.T(trafo2.T(e)).sign_nup
        self.exchange_nus = trafo1.T(trafo2.T(e)).exchange_nus

    def __str__(self):
        """Print formatting. Prints trafos according to operation order"""
        return self.trafo1.__str__() + self.trafo2.__str__()

    def T(self, diagram: Diagram):
        """ Performs T operation
        --- Params ---
            diagram: Diagram object on which the CompositeTrafo acts
        --- Returns ---
            A new Diagram object, which is the result of applying trafo1 on the result of the application of trafo2 on
            the input
        """
        return self.trafo1.T(self.trafo2.T(diagram))


def generate_full_group(group: list):
    """ Generates a full group of symmetry transformations
    --- Params ---
        group: list of Trafo/CompositeTrafo objects
    --- Returns ---
        List of all possible, distinct combinations of symmetry transformations """

    # List of ids of transformations already in the group
    ids = [(trafo.channel, trafo.diag_class, trafo.id ,trafo.sign_omega, trafo.sign_nu, trafo.sign_nup,
            trafo.exchange_nus) for trafo in group]
    #for i in range(len(ids)):
    #    print(ids[i] )

    # Group completion step
    done = False
    while not done:
        n = len(group)
        to_add = []     # Cannot change length of iterable object while iterating -> Store in temporary list
        for trafo1 in group:
            for trafo2 in group:
                new_trafo = CompositeTrafo(trafo1, trafo2)
                #print(new_trafo.channel, new_trafo.diag_class, new_trafo.id ,new_trafo.sign_omega, new_trafo.sign_nu,
                #      new_trafo.sign_nup, new_trafo.exchange_nus)

                # If new_trafo is new, add it to to_add list and keep track of finding by adding id to list.
                if (new_trafo.channel, new_trafo.diag_class, new_trafo.id, new_trafo.sign_omega, new_trafo.sign_nu,
                        new_trafo.sign_nup, new_trafo.exchange_nus) not in ids:
                    to_add.append(new_trafo)
                    #print(new_trafo.channel, new_trafo.diag_class, new_trafo.id ,new_trafo.sign_omega, new_trafo.sign_nu,
                    #  new_trafo.sign_nup, new_trafo.exchange_nus)
                    #print('append ', new_trafo)
                    ids.append((new_trafo.channel, new_trafo.diag_class, new_trafo.id ,new_trafo.sign_omega,
                                new_trafo.sign_nu, new_trafo.sign_nup, new_trafo.exchange_nus))
                    #print(ids[-1])

        # Extend the group
        group.extend(to_add)

        # If no new transformations have been found, one is done
        if len(group) == n:
            done = True

    return group




def generate_orbit(start_diagram, all_diagrams, symmetry_group, MF, with_freqs=False, only_diff_freq_args=False):
    '''
    Generate a graph containing the diagrammatic contributions which can be obtained from start_diagram
    by application of the symmetry_group
    --- Parameters ---
    start_diagram:  Diagram
    all_diagrams:   list of all diagrams
    symmetry_group: list containing the full group of transformations
    MF:             bool, True for Matsubara formalism, False for Keldysh formalism
    with_freqs:     bool, True if the diagrammatic contributions should be distinguished by their frequencies
                          False if the diagrammatic contributions are only distinguished by channel, diagrammatic class,
                             and spin indices (in Keldysh: + by Keldysh indices)
    only_diff_freq_args: bool, False: include in the returnd graph G all diagrammatic contributions,
                               True:  only include diagrammatic contributions which differs from the start_diagram only
                                       by the frequency arguments
    --- Return ---
    G:  Graph representing the orbit of start_diagram under the action of the full group
    '''
    #if not with_freqs and only_diff_freq_args:
     #   raise ValueError("only_diff_freq_args and only be True if with_freqs is True!")
    # Create T transformation objects
    ts = Trafo(5)
    # Generate a dictionary of all diagrams to store information on dependencies
    dependencies = establish_dictionary(all_diagrams, with_freqs)
    G = nx.DiGraph()
    orbit = [start_diagram]
    orbit_keys = [start_diagram.generate_key(with_freqs)]
    indeks = 0
    leng = 1
    while indeks < leng:
            diagram = orbit[indeks]
            # Act on current diagram with whole T transformation group
            for trafo in symmetry_group[1:]:
                transformed = trafo.T(diagram)

                # Assert that spin sector is allowed i.e., spin of the transformed diagram is in the allowed spin combs
                if transformed.get_spin_indices() not in spin_combinations:
                    # False evaluation implies need for spin flip.
                    # Keep track of spin-flip transformation in redefinition of trafo
                    trafo = CompositeTrafo(trafo, ts)
                    transformed = trafo.T(diagram)

                # Check if transformed diagram has already been visited
                if not transformed.generate_key(with_freqs) in orbit_keys \
                        and not (only_diff_freq_args and transformed.generate_key(False)!=diagram.generate_key(False)):
                    # If not, mark dependency with both needed transformation and original diagram
                    dependencies[transformed.generate_key(with_freqs)] = [trafo, diagram.generate_key(with_freqs)]
                    orbit.append(transformed)
                    orbit_keys.append(transformed.generate_key(with_freqs))
                    leng += 1
                if diagram.generate_key(with_freqs)!= transformed.generate_key(with_freqs) \
                        and not (diagram.generate_key(with_freqs), transformed.generate_key(with_freqs)) in G.edges \
                        and not (transformed.generate_key(with_freqs), diagram.generate_key(with_freqs)) in G.edges \
                        and not (only_diff_freq_args and transformed.generate_key(False)!=diagram.generate_key(False)):
                    G.add_edge(diagram.generate_key(with_freqs), transformed.generate_key(with_freqs), label=str(trafo))

                # Generate parity group for the transformed diagram
                parity_group = transformed.parity_group()

                if not MF:
                    # Evaluate orbit of the parity group of the transformed diagram
                    for parity_trafo in parity_group[1::]:

                        pair_diag = parity_trafo.T(transformed)

                        # Check if parity-related diagram has already been visited
                        if not dependencies[pair_diag.generate_key(with_freqs)] \
                                and not only_diff_freq_args:
                            # If not, mark it with parity trafo and transformed diagram from mapped
                            dependencies[pair_diag.generate_key(with_freqs)] = [parity_trafo, transformed.generate_key(with_freqs)]
                            orbit.append(pair_diag)
                            orbit_keys.append(pair_diag.generate_key(with_freqs))
                            leng += 1
                        if pair_diag.generate_key(with_freqs) != transformed.generate_key(with_freqs) \
                                and not (
                                only_diff_freq_args and transformed.generate_key(False) != pair_diag.generate_key(False)):
                            G.add_edge(transformed.generate_key(with_freqs), pair_diag.generate_key(with_freqs), label='P')
            indeks +=1
    return G