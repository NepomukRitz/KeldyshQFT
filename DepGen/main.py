from transformations import Trafo, generate_full_group, CompositeTrafo, generate_orbit
import networkx as nx
import matplotlib.pyplot as plt
from diagram import generate_diagrams, establish_dictionary, ParityTrafo, spin_combinations

MF = True # Matsubara formalism? (False for Keldysh formalism)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Create T transformation objects
    t0 = Trafo(0)
    t1 = Trafo(1)
    t2 = Trafo(2)
    t3 = Trafo(3)
    tc = Trafo(4)
    ts = Trafo(5)
    tph = Trafo(6)

    # Define the set of transformations to generate the group
    generators = [t0, t1, t2, t3, tc, tph]

    # Generate full group of T transformations and the full set of diagrams it will operate on
    symmetry_group = generate_full_group(generators)
    all_diagrams = generate_diagrams(MF)

    # Generate a dictionary of all diagrams to store information on dependencies
    dependencies = establish_dictionary(all_diagrams)
    num_indep = 0 # number of independent diagrammatic contributions

    # Erathostenes' Sieve-based labeling of dependencies
    # TODO Multiple visits to same diagram imply  consistency checks that can be checked with Trafo.ids equalities
    for diagram in all_diagrams:
        # If diagram has not been visited, value at current key is an empty list
        if not dependencies[diagram.generate_key()]:
            # Mark current diagram as independent i.e., related to itself through T0
            dependencies[diagram.generate_key()] = [t0]
            num_indep += 1

            # Act on current diagram with whole T transformation group
            for trafo in symmetry_group:
                transformed = trafo.T(diagram)

                # Assert that spin sector is allowed i.e., spin of the transformed diagram is in the allowed spin combs
                if transformed.get_spin_indices() not in spin_combinations:
                    # False evaluation implies need for spin flip.
                    # Keep track of spin-flip transformation in redefinition of trafo
                    trafo = CompositeTrafo(trafo, ts)
                    transformed = trafo.T(diagram)

                # Check if transformed diagram has already been visited
                if not dependencies[transformed.generate_key()]:
                    # If not, mark dependency with both needed transformation and original diagram
                    dependencies[transformed.generate_key()] = [trafo, diagram.generate_key()]
                    #G.add_edge(diagram.generate_key(), transformed.generate_key(), label=str(trafo))

                # Generate parity group for the transformed diagram
                parity_group = transformed.parity_group()

                if not MF:
                    # Evaluate orbit of the parity group of the transformed diagram
                    for parity_trafo in parity_group[1::]:

                        pair_diag = parity_trafo.T(transformed)

                        # Check if parity-related diagram has already been visited
                        if not dependencies[pair_diag.generate_key()]:
                            # If not, mark it with parity trafo and transformed diagram from mapped
                            dependencies[pair_diag.generate_key()] = [parity_trafo, transformed.generate_key()]
                            #G.add_edge(transformed.generate_key(), pair_diag.generate_key(), label='P')


    # Print out dependencies dictionary
    for key, value in dependencies.items():
        # If length of list at given key is one, it was marked with T0 and is, hence, independent
        if len(value) == 1:
            print("{} is independent!".format(key))
            continue

        # If first entry of value is a parity transformation, use the transformed diagram as key. Set value at new key
        # to a list only containing the info of the corresponding key
        if type(value[0]) == ParityTrafo:
            key_p = value[1]
            # Reassignment of value changes length to either 1 for independent key_p or 2 for a key_p that has been
            # transformed under a T trafo
            value = dependencies[key_p]

        # If length of value is 1 at this point, key and key_p are related through a parity trafo AND key_p is
        # independent. Hence, print only equality
        if len(value) == 1:
            print("{} = {}".format(key, key_p))

        # If not, key is related either through a P transformation to a diagram beforehand transformed with a T trafo OR
        # it is related through a T transformation to an independent component.
        else:
            print("{} = {} {}".format(key, value[0], value[1]))


# print number of independent components:
print("There are {} independent diagrammatic contributions out of a total of {}.".format(num_indep, len(all_diagrams)))

# Generate a graph containing the orbit of a diagram
diagram = all_diagrams[20]
print(diagram)
G = generate_orbit(diagram,all_diagrams,symmetry_group,MF)
pos = nx.spring_layout(G)
for p in pos:  # raise text positions
    pos[p][1] += 0.07
nx.draw(G, pos, font_size=16, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.isinteractive()
plt.show()