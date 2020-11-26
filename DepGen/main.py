from diagram import *


spin_combinations = ["ssss", "sbsb", "sbbs"]
diag_classes = [[1, 1], [0, 1], [1, 0], [0, 0]]


def generate_diagrams():
    result = []
    for diag_class in diag_classes:
        for spin in spin_combinations:
            for channel in ['a', 'p', 't']:
                for iK in range(16):
                    indices = []
                    for n in range(4, 0, -1):
                        indices.append((str(int(iK%(2**n)/(2**(n-1)))+1), spin[-n]))
                    diag = Diagram(channel, diag_class, indices)
                    result.append(diag)
    return result


def establish_dictionary(diagrams: List[Diagram]):
    d = {}
    for diag in diagrams:
        d[diag.generate_key()] = []
    return d


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    t0 = Trafo(0)
    t1 = Trafo(1)
    t2 = Trafo(2)
    t3 = Trafo(3)
    tc = Trafo(4)
    ts = Trafo(5)

    generators = [t0, t1, t2, t3, tc]

    symmetry_group = generate_full_group(generators)
    all_diagrams = generate_diagrams()

    dependencies = establish_dictionary(all_diagrams)

    for diagram in all_diagrams:
        if not dependencies[diagram.generate_key()]:
            dependencies[diagram.generate_key()] = [t0]
            for trafo in symmetry_group:
                transformed = trafo.T(diagram)
                if transformed.get_spin_indices() not in spin_combinations:
                    trafo = CompositeTrafo(trafo, ts)
                    transformed = trafo.T(diagram)
                if not dependencies[transformed.generate_key()]:
                    dependencies[transformed.generate_key()] = [trafo, diagram.generate_key()]
                parity_group = transformed.parity_group()
                for parity_trafo in parity_group:
                    pair_diag = parity_trafo.T(transformed)
                    if not dependencies[pair_diag.generate_key()]:
                        dependencies[pair_diag.generate_key()] = [parity_trafo, transformed.generate_key()]

    for key, value in dependencies.items():
        if len(value) == 1:
            print("{} is independent!".format(key))
            continue
        if type(value[0]) == ParityTrafo:
            key_p = value[1]
            value = dependencies[key_p]

        if len(value) == 1:
            print("{} = {}".format(key, key_p))                    
        else:
            print("{} = {} {}".format(key, value[0], value[1]))




