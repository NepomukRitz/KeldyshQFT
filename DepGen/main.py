from transformations import Trafo, generate_full_group, CompositeTrafo
from diagram import generate_diagrams, establish_dictionary, ParityTrafo, spin_combinations

MF = False  # Matsubara formalism? (False for Keldysh formalism)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Create T transformation objects
    t0 = Trafo(0)
    t1 = Trafo(1)
    t2 = Trafo(2)
    t3 = Trafo(3)
    tc = Trafo(4)
    ts = Trafo(5)

    # Define the set of transformations to generate the group
    generators = [t0, t1, t2, t3, tc]

    # Generate full group of T transformations and the full set of diagrams it will operate on
    symmetry_group = generate_full_group(generators)
    all_diagrams = generate_diagrams(MF)

    # Generate a dictionary of all diagrams to store information on dependencies
    dependencies = establish_dictionary(all_diagrams)

    # Erathostenes' Sieve-based labeling of dependencies
    # TODO Improve
    for diagram in all_diagrams:

        # Act on current diagram with whole T transformation group
        for trafo in symmetry_group:
            # Assert that spin sector is allowed i.e., spin of the transformed diagram is in the allowed spin combs
            if trafo.T(diagram).get_spin_indices() not in spin_combinations:
                # False evaluation implies need for spin flip.
                # Keep track of spin-flip transformation in redefinition of trafo
                trafo = CompositeTrafo(trafo, ts)

            transformed = trafo.T(diagram)

            # Check if transformed diagram has already been visited
            if not dependencies[transformed.generate_key()]:
                # If not, mark dependency with both needed transformation and original diagram
                dependencies[transformed.generate_key()] = [trafo, diagram.generate_key()]

                # Generate parity group for the transformed diagram
                parity_group = transformed.parity_group()

                # Evaluate orbit of the parity group of the transformed diagram
                for parity_trafo in parity_group:
                    pair_diag = parity_trafo.T(transformed)
                    dependencies[pair_diag.generate_key()].append([parity_trafo, transformed.generate_key()])

    # Print out dependencies dictionary
    # TODO Revise print out method and checks
    for key, value in dependencies.items():
        # If length of list at given key is one, it was marked with T0 and is, hence, independent
        if len(value) == 1:
            print(f"{key} is independent!")
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
            print(f"{key} = {key_p}")

        # If not, key is related either through a P transformation to a diagram beforehand transformed with a T trafo OR
        # it is related through a T transformation to an independent component.
        else:
            print(f"{key} = {value[0]} {value[1]}")
