a.out : main.cpp
	mkdir Output
	g++ main.cpp -o a.out -fopenmp

main.cpp : parameters.h a_bubble.h p_bubble.h t_bubble.h a_vertex.h p_vertex.h t_vertex.h bubbles.h data_structures.h free_propagator_functions.h frequency_grid.h integrator.h internal_symmetries.h Keldysh_symmetries.h loop.h propagator.h selfenergy.h state.h susceptibility.h util.h vertex.h


clean : 
	rm -f -r a.out Output

