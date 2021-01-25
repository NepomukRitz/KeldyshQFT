from transformations import Trafo, generate_full_group, CompositeTrafo, generate_orbit
import networkx as nx
import matplotlib.pyplot as plt
from diagram import generate_diagrams, establish_dictionary, ParityTrafo, spin_combinations, e, causality_enforcer

MF = True           # Matsubara formalism? (False for Keldysh formalism)                         <-------ENTER bool !!!
WITH_FREQS = False  # Should the diagrammatic contributions be distinguished by their frequency arguments? (And should
                        # they be printed?)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Create T transformation objects
    t0 = Trafo(0)
    t1 = Trafo(1)       # Exchange of incoming legs
    t2 = Trafo(2)       # Exchange of outgoing legs
    t3 = Trafo(3)       # Exchange of incoming and outgoing legs
    tc = Trafo(4)       # Complex conjugation
    ts = Trafo(5)       # Spin flip
    tph = Trafo(6)      # Particle-hole symmetry
    if MF:
        tr = Trafo(7)       # Hamiltonian is real function of creation and annihilation operators

    # Define the set of transformations to generate the group
    if MF:
        generators = [t0, t1, t2, t3, tc, tr]  # , tph]
    else:
        generators = [t0, t1, t2, t3, tc] # , tph]


    # Generate full group of T transformations and the full set of diagrams it will operate on
    symmetry_group = generate_full_group(generators)
    all_diagrams = generate_diagrams(MF, with_freqs=WITH_FREQS)

    # Generate a dictionary of all diagrams to store information on dependencies
    dependencies = establish_dictionary(all_diagrams, with_freqs=WITH_FREQS)
    num_indep = 0 # number of independent diagrammatic contributions

    # Erathostenes' Sieve-based labeling of dependencies
    # TODO Multiple visits to same diagram imply  consistency checks that can be checked with Trafo.ids equalities
    for diagram in all_diagrams:
        if not dependencies[diagram.generate_key(with_freqs=WITH_FREQS)]:
            num_indep +=1
            if not MF:
                # Get first diagrams related through P-trafo's, since they're exactly equal
                parity_group = diagram.parity_group()

                parity_diagrams = [diagram]
                for p_trafo in parity_group:
                    parity_diagrams.append(p_trafo.T(diagram))

                eq_to_zero = causality_enforcer(parity_diagrams)

                if not eq_to_zero:
                    for p_trafo in parity_group:
                        dependencies[p_trafo.T(diagram).generate_key(with_freqs=WITH_FREQS)] = [p_trafo, diagram.generate_key(with_freqs=WITH_FREQS)]
            else:
                eq_to_zero = False
                parity_diagrams = [diagram]

            for p_diagram in parity_diagrams:
                # Act on current diagram with whole T transformation group
                for trafo in symmetry_group:
                    # Assert that spin sector is allowed
                    if trafo.T(p_diagram).get_spin_indices() not in spin_combinations:
                        # False evaluation implies need for spin flip.
                        trafo = CompositeTrafo(trafo, ts)

                    transformed = trafo.T(p_diagram)

                    # Check if transformed diagram has already been visited
                    if not dependencies[transformed.generate_key(with_freqs=WITH_FREQS)]:
                        # If not, mark dependency with both needed transformation and original diagram
                        if eq_to_zero:
                            dependencies[transformed.generate_key(with_freqs=WITH_FREQS)] = [0]
                        else:
                            dependencies[transformed.generate_key(with_freqs=WITH_FREQS)] = [trafo, diagram.generate_key(with_freqs=WITH_FREQS)]



    # Print out dependencies dictionary
    iter = 0
    for key, value in dependencies.items():
        if len(value) == 1:
            print(iter, f": {key} \t = 0")
            iter +=1
            continue

        if type(value[0]) == ParityTrafo:
            print(iter, f": {key} \t = P {value[1]}")

        else:
            if value[0] == t0:
                print(iter, f": {key} \t is independent!")
                iter +=1
                continue
            else:
                print(iter, f": {key} \t = {value[0]} {value[1]}")
        iter +=1

# print number of independent components:
print("There are {} independent diagrammatic contributions out of a total of {}.".format(num_indep, len(all_diagrams)))


# Generate a graph containing the orbit of a diagram
print()
DiagNumber = 0                      # number of Diagram whose orbit should be generated        <-------ENTER number !!!
plt.figure(1)
diagram = all_diagrams[DiagNumber]  # Diagram whose orbit should be generated
print("Figure 1 shows the orbit of the diagrammatic contribution ", diagram, " under the group action of full_group. "
                                                                             "(only show diagrammatic contribution which differ from diagram solely by the frequency arguments)")
G = generate_orbit(diagram,all_diagrams,symmetry_group, MF, with_freqs=True, only_diff_freq_args=True)
pos = nx.spring_layout(G)
for p in pos:  # raise text positions
    pos[p][1] += 0.07
nx.draw(G, pos, font_size=16, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.suptitle("Figure 1 shows the orbit of the diagrammatic contribution "+ str(diagram) + " under the group action of full_group. "
                                                                             +"(only show diagrammatic contribution which differ from diagram solely by the frequency arguments)")
plt.isinteractive()
#plt.show()

# Generate a graph containing the orbit of a diagram
plt.figure(2)
diagram = all_diagrams[DiagNumber]
print("Figure 2 shows the orbit of the diagrammatic contribution ", diagram, " under the group action of full_group. "
                                                                             "(shows all diagrammatic contributions in the orbit)")
G = generate_orbit(diagram,all_diagrams,symmetry_group, MF, with_freqs=True, only_diff_freq_args=False)
pos = nx.spring_layout(G)
for p in pos:  # raise text positions
    pos[p][1] += 0.07
nx.draw(G, pos, font_size=16, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.suptitle("Figure 2 shows the orbit of the diagrammatic contribution "+ str(diagram)+ " under the group action of full_group. "
                                                                             "(shows all diagrammatic contributions in the orbit)")

plt.isinteractive()
plt.show()