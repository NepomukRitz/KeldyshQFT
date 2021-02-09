from typing import List
from diagram import Diagram, e
from general_purpose import conjugate_spin
import global_parameters as gp


class Trafo:
    def __init__(self, i):
        """ Initialized a Trafo object
        --- Params ---
            i: integer. 0, 1, 2, 3, 4 and 5 and 6 supported
        --- Returns ---
            Initialized Trafo object with a calculated id """
        self.i = i
        self.id = self.T(e).indices

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
        if gp.with_freqs:
            new_freqs = self.transform_freqs(diagram.freqs, diagram.channel, diagram.diag_class)
        else:
            new_freqs = diagram.freqs
        return Diagram(new_channel, new_diag_class, new_indices, new_freqs)

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

    def transform_freqs(self, freqs: List[int], channel: str, diag_class: List[int]):
        """ Transforms the signs of the frequencies
        --- Params ---
            freqs: List of frequencies to be transformed
            channel: str indicating the channel
            diag_class: 2-position int-list with info on the position of the bare vertices of the diagram
        --- Returns ---
            List of ints marking the transformed frequencies
        """
        new_freqs = freqs.copy()

        if len(freqs) > 0:
            if self.i == 1:
                if channel == 'a' or channel == 't':
                    new_freqs[0] *= -1
                    if gp.param == "bosonic" and diag_class == [0, 0]:
                        new_freqs[2] *= -1
                    else:
                        new_freqs[1] = freqs[2]
                        new_freqs[2] = freqs[1]
                else:  # channel p
                    if gp.param == "bosonic" and diag_class == [0, 0]:
                        new_freqs[1] = freqs[2]
                        new_freqs[2] = freqs[1]
                    else:
                        new_freqs[2] *= -1

            elif self.i == 2:
                if channel == 'p':
                    if gp.param == "bosonic" and diag_class == [0, 0]:
                        new_freqs[1] = -freqs[2]
                        new_freqs[2] = -freqs[1]
                    else:
                        new_freqs[1] *= -1

            elif self.i == 3:
                if channel == 'a' or channel == 't':
                    new_freqs[0] *= -1
                    if gp.param == "bosonic" and diag_class == [0, 0]:
                        new_freqs[2] *= -1
                    else:
                        new_freqs[1] = freqs[2]
                        new_freqs[2] = freqs[1]
                else:  # channel p      # Same effect on v, v' as on f, f'
                    new_freqs[1] *= -1
                    new_freqs[2] *= -1

            elif self.i == 4:
                if channel == 'a' or channel == 'p':
                    if gp.param == "bosonic" and diag_class == [0, 0]:
                        new_freqs[2] *= -1
                    else:
                        new_freqs[1] = freqs[2]
                        new_freqs[2] = freqs[1]
                else:  # channel t      # No effect on either v, v' or f, f'
                    new_freqs[0] *= -1

        return new_freqs


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
        self.id = trafo1.T(trafo2.T(e)).indices

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
    ids = [trafo.id for trafo in group]

    # Group completion step
    done = False
    while not done:
        n = len(group)
        to_add = []     # Cannot change length of iterable object while iterating -> Store in temporary list
        for trafo1 in group:
            for trafo2 in group:
                new_trafo = CompositeTrafo(trafo1, trafo2)

                # If new_trafo is new, add it to to_add list and keep track of finding by adding id to list.
                if new_trafo.id not in ids:
                    to_add.append(new_trafo)
                    ids.append(new_trafo.id)

        # Extend the group
        group.extend(to_add)

        # If no new transformations have been found, one is done
        if len(group) == n:
            done = True

    return group
