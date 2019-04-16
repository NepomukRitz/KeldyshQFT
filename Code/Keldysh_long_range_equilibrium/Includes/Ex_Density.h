#ifndef EX_DENSITY_21022019
#define EX_DENSITY_21022019

template <int mode> class Integrator_density{
	 public:
		bool spin;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		Integrator_density(Physics &phy_in,
		                   Numerics &num_in,
		                   Ex_Precomputation<mode> &pre_in,
		                   Substitution<mode> &sub_in);
		matrix<double> integrate(bool spin);
};

template<int mode> Integrator_density<mode>::Integrator_density(Physics &phy_in,
                                                                Numerics &num_in,
                                                                Ex_Precomputation<mode> &pre_in,
                                                                Substitution<mode> &sub_in):
                                                                phy(phy_in),
                                                                num(num_in),
                                                                pre(pre_in),
                                                                sub(sub_in){
}

template<int mode> matrix<double> Integrator_density<mode>::integrate(bool spin){
	matrix<syma<complex<double> > G; 
	if(spin==1){
		G = pre.Gu;
	}
	else{
		G = pre.Gd;
	}
	matrix<double> density(num.Nges);
	density = 0.0;
	for(int i=0; i<num.num_freq_pre-1; i++){
		double n1 = sub.weight_concatenated(freq_pre(i))*fermi(sub.resu_concatenated(freq_pre(i)),phy.mu,phy.T)/M_PI;
		double n2 = sub.weight_concatenated(freq_pre(i+1))*fermi(sub.resu_concatenated(freq_pre(i+1)),phy.mu,phy.T)/M_PI;
		for(int j=0; j<num.Nges; ++j){
			density(j) -= 0.5*(imag(G(i)(j,j))*n1 + imag(G(i+1)(j,j))*n2)*(freq_pre(i+1) - freq_pre(i));
		}
	}
	return density;
}





#endif

