#include <boost/python.hpp>
#include <iostream>
#include "matrix.h"
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <complex> 
#include "basic.h"

using namespace boost::python;

struct Conductance_for_py{
	syma<complex<double> > ERetu_at_mu;
	syma<complex<double> > ERetd_at_mu;
	int Nges;
	double h;
	double mu;
	double Lambda;
	void set_parameters(int Nges_in,double mu_in, double h_in, double Lambda_in){
		Nges = Nges_in;
		h=h_in;
		mu=mu_in;
		Lambda=Lambda_in;
		ERetu_at_mu.resize(Nges);
		ERetd_at_mu.resize(Nges);
	}
	void populate_ERet_at_mu(int i, int j, complex<double> value_up, complex<double> value_down){
		ERetu_at_mu(i,j) = value_up;
		ERetd_at_mu(i,j) = value_down;
	}
	double compute(bool spin){
		syma<complex<double> > G;
	 	double omega;
	 	if(spin==1){
			omega = mu - 0.5*h; 	
			matrix<syma<complex<double> > > GS=green_and_single_eq_non_ps(mu, ERetu_at_mu, h, 1.0, Lambda);
			G = GS(0);
		}
		else{
			omega = mu + 0.5*h; 	
			matrix<syma<complex<double> > > GS=green_and_single_eq_non_ps(mu, ERetd_at_mu, -h, 1.0, Lambda);
			G = GS(0);
		}
		if(omega>=-2. && omega<=2.){
			return 2.0*(1.0 - (omega/2.0)*(omega/2.0))*pow(abs(G(G.dim-1,0)),2);
		}
		else{
			return 0.0;
		}
	}
};

BOOST_PYTHON_MODULE(One_particle_conductance_for_py)
{
//First level structures:

	class_<Conductance_for_py>("Conductance_for_py_py")
	      .def("populate_ERet_at_mu", &Conductance_for_py::populate_ERet_at_mu)
	      .def("set_parameters", &Conductance_for_py::set_parameters)
	      .def("compute", &Conductance_for_py::compute)
	;
	

	
}

