from transformations import Trafo, CompositeTrafo
from diagram import establish_dictionary, generate_diagrams
import global_parameters as gp
import networkx as nx


def generate_orbit(start_diagram, all_diagrams, symmetry_group, only_diff_freq_args=False):
    """
    Generate a graph containing the diagrammatic contributions which can be obtained from start_diagram
    by application of the symmetry_group
    --- Parameters ---
    start_diagram:  Diagram
    all_diagrams:   list of all diagrams
    symmetry_group: list containing the full group of transformations
    only_diff_freq_args: bool, False: include in the returned graph G all diagrammatic contributions,
                               True:  only include diagrammatic contributions which differs from the start_diagram only
                                       by the frequency arguments
    --- Return ---
    G:  Graph representing the orbit of start_diagram under the action of the full group
    """

    # if not with_freqs and only_diff_freq_args:
    #   raise ValueError("only_diff_freq_args and only be True if with_freqs is True!")
    # Create T transformation objects
    ts = Trafo(5)
    # Generate a dictionary of all diagrams to store information on dependencies
    dependencies = establish_dictionary(all_diagrams)
    G = nx.DiGraph()
    orbit = [start_diagram]
    orbit_keys = [start_diagram.generate_key()]
    indeks = 0
    leng = 1
    while indeks < leng:
        diag = orbit[indeks]
        # Act on current diagram with whole T transformation group
        for trafo in symmetry_group[1:]:
            transformed = trafo.T(diag)

            # Assert that spin sector is allowed i.e., spin of the transformed diagram is in the allowed spin combs
            if transformed.get_spin_indices() not in gp.spin_combinations:
                # False evaluation implies need for spin flip.
                # Keep track of spin-flip transformation in redefinition of trafo
                trafo = CompositeTrafo(trafo, ts)
                transformed = trafo.T(diag)

            # Check if transformed diagram has already been visited
            if not transformed.generate_key() in orbit_keys \
                    and not (only_diff_freq_args and transformed.generate_key() != diag.generate_key()):
                # If not, mark dependency with both needed transformation and original diagram
                dependencies[transformed.generate_key()] = [trafo, diag.generate_key()]
                orbit.append(transformed)
                orbit_keys.append(transformed.generate_key())
                leng += 1
            if diag.generate_key() != transformed.generate_key() \
                    and not (diag.generate_key(), transformed.generate_key()) in G.edges \
                    and not (transformed.generate_key(), diag.generate_key()) in G.edges \
                    and not (only_diff_freq_args and transformed.generate_key() != diag.generate_key()):
                G.add_edge(diag.generate_key(), transformed.generate_key(), label=str(trafo))

            # Generate parity group for the transformed diagram
            parity_group = transformed.parity_group()

            if not gp.matsubara:
                # Evaluate orbit of the parity group of the transformed diagram
                for parity_trafo in parity_group[1::]:

                    pair_diag = parity_trafo.T(transformed)

                    # Check if parity-related diagram has already been visited
                    if not dependencies[pair_diag.generate_key()] \
                            and not only_diff_freq_args:
                        # If not, mark it with parity trafo and transformed diagram from mapped
                        dependencies[pair_diag.generate_key()] = [parity_trafo, transformed.generate_key()]
                        orbit.append(pair_diag)
                        orbit_keys.append(pair_diag.generate_key())
                        leng += 1
                    if pair_diag.generate_key() != transformed.generate_key() \
                            and not (
                            only_diff_freq_args and transformed.generate_key() != pair_diag.generate_key()):
                        G.add_edge(transformed.generate_key(), pair_diag.generate_key(), label='P')
        indeks += 1
    return G
