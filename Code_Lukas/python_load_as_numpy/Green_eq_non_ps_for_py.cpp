#include <boost/python.hpp>
#include <iostream>
#include "matrix.h"
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <complex> 
#include "basic.h"
#include "Syma_Matrix.h"

using namespace boost::python;

struct Green_eq_non_ps{
 	syma<complex<double> > ERet;
	syma<complex<double> > G;
	void resize(int Nges){
	 	ERet.resize(Nges);
	}
	void populate_ERet(int i, int j, complex<double> value){
	 	ERet(i,j) = value;
	}
	void compute(double f, double h, double Lambda){
	 	matrix<syma<complex<double> > > GS;
		GS=green_and_single_eq_non_ps(f, ERet, h, 1.0, Lambda);
		G = GS(0);
	}
	complex<double>  component(int i, int j){
	 	return G(i,j);
	}
};

BOOST_PYTHON_MODULE(Green_eq_non_ps_for_py)
{
//First level structures:

	class_<Green_eq_non_ps>("Green_eq_non_ps_py")
	      .def("resize", &Green_eq_non_ps::resize)
	      .def("populate_ERet", &Green_eq_non_ps::populate_ERet)
	      .def("compute", &Green_eq_non_ps::compute)
	      .def("component", &Green_eq_non_ps::component)
	;
	

	
}
