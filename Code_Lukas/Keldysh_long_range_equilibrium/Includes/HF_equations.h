#ifndef HF_EQUATIONS_02062018
#define HF_EQUATIONS_02062018

#include <integrate_new.h>
#include "Stops.h"

template <int mode> class Integrand_G_lesser_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Integrand_G_lesser_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template <int mode> Integrand_G_lesser_zero_mag<mode>::Integrand_G_lesser_zero_mag(Physics &phy_in,
                                                                                     Numerics &num_in,
                                                                                     Precomputation_zeromag<mode> &pre_in, 
																					 Substitution<mode> &sub_in): 
																					 phy(phy_in),
																					 num(num_in),
												                                     pre(pre_in),
																					 sub(sub_in){}


template<int mode> syma<complex<double> > Integrand_G_lesser_zero_mag<mode>::operator()(double internal){
	syma<complex<double> > Gu = pre.iGu(internal); 
	syma<complex<double> > ret(num.Nges); 
	double intline = sub.resu_concatenated(internal); 
	double nf   = -sub.weight_concatenated(internal)*(fermi(intline, phy.mu, phy.T))/M_PI;
	
	for (int j=0; j<num.Nges; j++) {
		for (int i=0; i<=j; i++) {
			ret(j,i) = nf*Gu(j,i).imag();
		}
	}
	return ret;
}

template<int mode> matrix<double> Integrand_G_lesser_zero_mag<mode>::select(syma<complex<double> > &M){
	int N_halfp = num.N +1;
	matrix<double> n(8*N_halfp);
	for (int i=0;i<N_halfp;i++) {
		n(i)=M(i,i).real();
		n(i+1*N_halfp)=M(i,0).real();
		n(i+2*N_halfp)=M(num.Nges-i-1,i).real();
		n(i+3*N_halfp)=M(i+1,i).real();
		n(i+4*N_halfp)=M(i,i).imag();
		n(i+5*N_halfp)=M(i,0).imag();
		n(i+6*N_halfp)=M(num.Nges-i-1,i).imag();
		n(i+7*N_halfp)=M(i+1,i).imag();
	}
	return n;
}


template<int mode> class HF_equations{
 	public:
		static const double eps = 1e-16;
		static const double accuracy=ACCURACY_HF_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Substitution<mode> &sub;
		Barevertex &barevertex;
		Stops<mode> &stops_obj;
		double Lambda;
		HF_equations(Physics &phy_in, 
		             Numerics &num_in,
					 Substitution<mode> &sub_in,
					 Barevertex &barevertex_in,
					 Stops<mode> &stops_obj_in,
					 double Lambda_in);
		syma<complex<double> > GK_integrated(Precomputation_zeromag<mode> &pre);
		syma<complex<double> > Selfenergy(Precomputation_zeromag<mode> &pre);
		syma<complex<double> > HF_iteration(int number_of_iterations);
};

template<int mode> HF_equations<mode>::HF_equations(Physics &phy_in,
                                                    Numerics &num_in,
				                                    Substitution<mode> &sub_in,
					                                Barevertex &barevertex_in,
													Stops<mode> &stops_obj_in,
													double Lambda_in):
													phy(phy_in),
													num(num_in),
													sub(sub_in),
													barevertex(barevertex_in),
													stops_obj(stops_obj_in),
													Lambda(Lambda_in){};


template<int mode> syma<complex<double> > HF_equations<mode>::GK_integrated(Precomputation_zeromag<mode> &pre){
 	Integrand_G_lesser_zero_mag<mode> Int_G(phy,num,pre,sub);
	matrix<double> stops = stops_obj.Self_energy_stops(0.0);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,Int_G);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}

template<int mode> syma<complex<double> > HF_equations<mode>::Selfenergy(Precomputation_zeromag<mode> &pre){
 	syma<complex<double> > erg(num.Nges);
 	syma<complex<double> > GK_int = GK_integrated(pre);
	//Only for two particle interaction:
	for(int j=0; j<num.Nges; ++j){
	 	for(int i=0; i<num.Nges; ++i){	 
		 	erg(j,j)+= (barevertex(i,1,j,1,i,1,j,1) +barevertex(i,0,j,1,i,0,j,1)) *GK_int(i,i); 
		}
	}
	for(int j=0; j<num.Nges; ++j){
	 	for(int i=0; i<j; ++i){
		 	erg(j,i) = (barevertex(i,1,j,1,j,1,i,1) + barevertex(i,0,j,1,j,0,i,1))*GK_int(j,i);
		}
	}
	return erg;
	
}

template<int mode> syma<complex<double> > HF_equations<mode>::HF_iteration(int number_of_iterations){
 	Precomputation_zeromag<mode> pre(phy,num);	
	matrix<double> zero_one(2);
	zero_one(0) = 0.0;
	zero_one(1) = 1.0;
	matrix<syma<complex<double> > > Eu_trivial(2);
	Eu_trivial(0) = phy.hamiltonian;
	Eu_trivial(1) = phy.hamiltonian;
	linear_ipol_bin<syma<complex<double> > > iEu(zero_one,Eu_trivial);
	pre.precompute(Lambda,sub,iEu);
	for(int i=0; i<number_of_iterations; ++i){
		syma<complex<double> > Selfenergy_hf = Selfenergy(pre); 	
		Eu_trivial(0) = phy.hamiltonian + Selfenergy_hf;
		Eu_trivial(1) = phy.hamiltonian + Selfenergy_hf;
		linear_ipol_bin<syma<complex<double> > > iEu(zero_one,Eu_trivial);
		pre.precompute(Lambda,sub,iEu);
	}
	return Selfenergy(pre);
}







#endif
