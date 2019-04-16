#ifndef EX_HF_EQUATIONS_24022019
#define EX_HF_EQUATIONS_24022019

template <int mode> class Integrand_G_lesser{
	public:
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double spin;
		Integrand_G_lesser(Physics &phy_in, Numerics &num_in, Ex_Precomputation<mode> &pre_in, Substitution<mode> &sub_in, double spin_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template <int mode> Integrand_G_lesser<mode>::Integrand_G_lesser(Physics &phy_in, 
                                                                 Numerics &num_in,
                                                                 Ex_Precomputation<mode> &pre_in, 
                                                                 Substitution<mode> &sub_in,
                                                                 double spin_in): 
                                                                 phy(phy_in),
                                                                 num(num_in),
                                                                 pre(pre_in),
                                                                 sub(sub_in),
                                                                 spin(spin_in){
}

template<int mode> syma<complex<double> > Integrand_G_lesser<mode>::operator()(double internal){
	syma<complex<double> > G; 
	if(spin==1){
		G = pre.iGu(internal); 
	}
	else{
		G = pre.iGd(internal); 
	}
	syma<complex<double> > ret(num.Nges); 
	double intline = sub.resu_concatenated(internal); 
	double nf   = -sub.weight_concatenated(internal)*(fermi(intline, phy.mu, phy.T))/M_PI;
	
	for (int j=0; j<num.Nges; j++) {
		for (int i=0; i<=j; i++) {
			ret(j,i) = nf*G(j,i).imag();
		}
	}
	return ret;
}

template<int mode> matrix<double> Integrand_G_lesser<mode>::select(syma<complex<double> > &M){
	return matrix_to_vector(M);
}


template<int mode> class Ex_HF_equations{
 	public:
		static const double eps = 1e-16;
		double accuracy; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Substitution<mode> &sub;
		Barevertex &barevertex;
		Ex_Stops<mode> &stops_obj;
		matrix<double> &additional_stops;
		double Lambda;
		Ex_HF_equations(Physics &phy_in, 
		                Numerics &num_in,
		                Substitution<mode> &sub_in,
		                Barevertex &barevertex_in,
		                Ex_Stops<mode> &stops_obj_in,
		                matrix<double> &additional_stops_in,
		                double Lambda_in,
		                double accuracy_in);
		syma<complex<double> > G_lesser_integrated(Ex_Precomputation<mode> &pre, double spin);
		matrix<syma<complex<double> > > Selfenergy(Ex_Precomputation<mode> &pre);
		matrix<syma<complex<double> > > HF_iteration(int number_of_iterations);
};

template<int mode> Ex_HF_equations<mode>::Ex_HF_equations(Physics &phy_in,
                                                          Numerics &num_in,
                                                          Substitution<mode> &sub_in,
                                                          Barevertex &barevertex_in,
                                                          Ex_Stops<mode> &stops_obj_in,
                                                          matrix<double> &additional_stops_in,
                                                          double Lambda_in,
                                                          double accuracy_in):
                                                          phy(phy_in),
                                                          num(num_in),
                                                          sub(sub_in),
                                                          barevertex(barevertex_in),
                                                          stops_obj(stops_obj_in),
                                                          additional_stops(additional_stops_in),
                                                          Lambda(Lambda_in),
                                                          accuracy(accuracy_in){
}

template<int mode> syma<complex<double> > Ex_HF_equations<mode>::G_lesser_integrated(Ex_Precomputation<mode> &pre, double spin){
 	Integrand_G_lesser<mode> Int_G(phy,num,pre,sub,spin);
	matrix<double> stops = stops_obj(0.0, additional_stops);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;

	double delta = .0;
	for(int i=0; i<stops.dim_c-1; i++){
		delta = stops(i+1)-stops(i);
		if(delta>eps){
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,Int_G);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}

template<int mode> matrix<syma<complex<double> > > Ex_HF_equations<mode>::Selfenergy(Ex_Precomputation<mode> &pre){
 	matrix<syma<complex<double> > > erg(2);
	erg(0).resize(num.Nges);
	erg(1).resize(num.Nges);
 	syma<complex<double> > Gu_int = G_lesser_integrated(pre,1);
 	syma<complex<double> > Gd_int = G_lesser_integrated(pre,0);
	//Only for two particle interaction:
	for(int j=0; j<num.Nges; ++j){
	 	for(int i=0; i<num.Nges; ++i){	 
		 	erg(0)(j,j)+= barevertex(i,1,j,1,i,1,j,1)*Gu_int(i,i) + barevertex(i,0,j,1,i,0,j,1)*Gd_int(i,i); 
		 	erg(1)(j,j)+= barevertex(i,1,j,0,i,1,j,0)*Gu_int(i,i) + barevertex(i,0,j,0,i,0,j,0)*Gd_int(i,i); 
		}
	}
	for(int j=0; j<num.Nges; ++j){
	 	for(int i=0; i<j; ++i){
		 	erg(0)(j,i) = barevertex(i,1,j,1,j,1,i,1)*Gu_int(j,i) + barevertex(i,0,j,1,j,0,i,1)*Gd_int(j,i);
		 	erg(1)(j,i) = barevertex(i,1,j,0,j,1,i,0)*Gu_int(j,i) + barevertex(i,0,j,0,j,0,i,0)*Gd_int(j,i);
		}
	}
	return erg;
}

template<int mode> matrix<syma<complex<double> > > Ex_HF_equations<mode>::HF_iteration(int number_of_iterations){
 	Ex_Precomputation<mode> pre(phy,num);	
	matrix<double> wf_triv(2);
	wf_triv(0) = 0.0;
	wf_triv(1) = 1.0;
	matrix<syma<complex<double> > > Eu(2);
	Eu(0) = phy.hamiltonian;
	Eu(1) = phy.hamiltonian;
	matrix<syma<complex<double> > > Ed = Eu;
	linear_ipol_bin<syma<complex<double> > > iEu(wf_triv,Eu);
	linear_ipol_bin<syma<complex<double> > > iEd(wf_triv,Ed);
	pre.set_freq_pre(sub);
	pre.precompute(Lambda,sub,iEu,iEd);
	matrix<syma<complex<double> > > Selfenergy_hf;
	for(int i=0; i<number_of_iterations; ++i){
		Selfenergy_hf = Selfenergy(pre); 	
		Eu(0) = phy.hamiltonian + Selfenergy_hf(0);
		Eu(1) = phy.hamiltonian + Selfenergy_hf(0);
		Ed(0) = phy.hamiltonian + Selfenergy_hf(1);
		Ed(1) = phy.hamiltonian + Selfenergy_hf(1);
		linear_ipol_bin<syma<complex<double> > > iEu(wf_triv,Eu);
		linear_ipol_bin<syma<complex<double> > > iEd(wf_triv,Ed);
		pre.precompute(Lambda,sub,iEu,iEd);
	}
	Selfenergy_hf(0) += phy.hamiltonian;
	Selfenergy_hf(1) += phy.hamiltonian;
	return Selfenergy_hf;
}

#endif

