#ifndef EX_CONDUCTANCE_21022019
#define EX_CONDUCTANCE_21022019

#include "integrate_new.h"
#include "Ex_Precomputation.h"
#include "Ex_stops.h"
 

//This class is not meant for speed!
template <int mode> class Integrand_cond_vertex_cont_naive{
	 public:
		double external_freq;
		bool spin;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		Ex_Vertex<mode> &gamma;
		Integrand_cond_vertex_cont_naive(double external_freq_in,
		                                 bool spin_in,
		                                 Physics &phy_in,
		                                 Numerics &num_in,
		                                 Ex_Precomputation<mode> &pre_in,
		                                 Substitution<mode> &sub_in,
		                                 Ex_Vertex<mode> &gamma_in);
		matrix<double> select(matrix<complex<double> >  &M);
		matrix<complex<double> > operator()(double internal); 
};

template<int mode> Integrand_cond_vertex_cont_naive<mode>::Integrand_cond_vertex_cont_naive(double external_freq_in,
                                                                                            bool spin_in,
                                                                                            Physics &phy_in,
                                                                                            Numerics &num_in,
                                                                                            Ex_Precomputation<mode> &pre_in,
                                                                                            Substitution<mode> &sub_in,
                                                                                            Ex_Vertex<mode> &gamma_in):
                                                                                            external_freq(external_freq_in),
                                                                                            spin(spin_in),
                                                                                            phy(phy_in),
                                                                                            num(num_in),
                                                                                            pre(pre_in),
                                                                                            sub(sub_in),
                                                                                            gamma(gamma_in){
}

template<int mode> matrix<double> Integrand_cond_vertex_cont_naive<mode>::select(matrix<complex<double> >  &M){
	return matrix_to_vector(M);	
}

template<int mode> matrix<complex<double> > Integrand_cond_vertex_cont_naive<mode>::operator()(double internal){
	//cout<<"internal="<<internal<<endl;
	double intline = sub.resu_concatenated(internal);
	double prefactor = sub.weight_concatenated(internal)*2.*sqrt(1. - (intline/2.0)*(intline/2.0) )/M_PI;
	double fp = 1./tanh(((intline + external_freq)/2. - phy.mu)/phy.T) - 1. + 2.*fermi(intline, phy.mu, phy.T);
	double fx = 1./tanh(((intline - external_freq)/2.         )/phy.T) - 1. + 2.*fermi(intline, phy.mu, phy.T);
	//cout<<"phy.T="<<phy.T<<endl;
	//cout<<"fp="<<fp<<endl;
	//cout<<"fx="<<fx<<endl;

	syma<complex<double> > Gu = pre.iGu(internal);
	syma<complex<double> > Gd = pre.iGd(internal);
	matrix<complex<double> > Kuu(num.Nges,num.Nges);
	matrix<complex<double> > Kdd(num.Nges,num.Nges);
	matrix<complex<double> > Kud(num.Nges,num.Nges);
	matrix<complex<double> > Kdu(num.Nges,num.Nges);
	
	matrix<matrix<complex<double> > > APud = gamma.aPud.complete_ipol(external_freq + intline, num.pos_NfbP_2mu,num.Lp_bounds); 
	APud = fp*imaginary_part_of_str(APud); 
	

	matrix<matrix<complex<double> > > APss;
	matrix<matrix<complex<double> > > ADss;
	matrix<matrix<complex<double> > > AXud;
	
	if(spin==1){
		APss = gamma.aPuu.complete_ipol(external_freq + intline, num.pos_NfbP_2mu,num.Lp_bounds); 
		APss = fp*imaginary_part_of_str(APss); 
		ADss = gamma.aDuu.complete_ipol(intline - external_freq, num.pos_NfbX_0, num.Lx_bounds); 
		ADss = fx*imaginary_part_of_str(ADss); 
		AXud = gamma.aXud.complete_ipol(intline - external_freq, num.pos_NfbX_0, num.Lx_bounds); 
		AXud = fx*imaginary_part_of_str(AXud); 
	}
	else{
		APss = gamma.aPdd.complete_ipol(external_freq + intline, num.pos_NfbP_2mu,num.Lp_bounds); 
		APss = fp*imaginary_part_of_str(APss); 
		ADss = gamma.aDdd.complete_ipol(intline - external_freq, num.pos_NfbX_0, num.Lx_bounds); 
		ADss = fx*imaginary_part_of_str(ADss); 
		AXud = gamma.aXud.complete_ipol(external_freq - intline, num.pos_NfbX_0, num.Lx_bounds); 
		AXud = fx*imaginary_part_of_str(AXud); 
	}

	int Lp = (APud.dim_r-1)/2;
	int Lx = (AXud.dim_r-1)/2;
	matrix<complex<double> > ret(num.Nges,num.Nges);
	ret = (complex<double>) 0.0;
	for(int j2 = 0; j2<num.Nges; ++j2){ 
		for(int i2 = 0; i2<num.Nges; ++i2){ 
			for(int j1 = 0; j1<num.Nges; ++j1){ 
				for(int i1 = 0; i1<num.Nges; ++i1){ 
					if(spin==1){
						complex<double> tmp= (complex<double>) 0.0;
						if(inrange(APss,j2-j1,i2-i1,j1,i1)) tmp += APss(j2-j1+Lp,i2-i1+Lp)(j1-max(0,-(j2-j1)),i1-max(0,-(i2-i1)));
						if(inrange(ADss,i2-j1,j2-i1,j1,i1)) tmp -= ADss(i2-j1+Lx,j2-i1+Lx)(j1-max(0,-(i2-j1)),i1-max(0,-(j2-i1)));
						//ret(j2,i2) += prefactor*conj(Gu(i1,0))*Gu(j1,0)*tmp; 
						ret(j2,i2) += prefactor*conj(Gu(num.Nges-1,i1))*Gu(num.Nges-1,j1)*tmp; 

						tmp = (complex<double>) 0.0;
						if(inrange(APud,j1-j2,i1-i2,j2,i2)) tmp += APud(j1-j2+Lp,i1-i2+Lp)(j2-max(0,-(j1-j2)),i2-max(0,-(i1-i2)));
						if(inrange(AXud,j1-i2,i1-j2,i2,j2)) tmp -= AXud(j1-i2+Lx,i1-j2+Lx)(i2-max(0,-(j1-i2)),j2-max(0,-(i1-j2))); 
						//ret(j2,i2) += prefactor*conj(Gd(i1,0))*Gd(j1,0)*tmp; 
						ret(j2,i2) += prefactor*conj(Gd(num.Nges-1,i1))*Gd(num.Nges-1,j1)*tmp; 
					}
					else{
						complex<double> tmp= (complex<double>) 0.0;
						if(inrange(APss,j2-j1,i2-i1,j1,i1)) tmp += APss(j2-j1+Lp,i2-i1+Lp)(j1-max(0,-(j2-j1)),i1-max(0,-(i2-i1)));
						if(inrange(ADss,i2-j1,j2-i1,j1,i1)) tmp -= ADss(i2-j1+Lx,j2-i1+Lx)(j1-max(0,-(i2-j1)),i1-max(0,-(j2-i1)));
						//ret(j2,i2) += prefactor*conj(Gd(i1,0))*Gd(j1,0)*tmp; 
						ret(j2,i2) += prefactor*conj(Gd(num.Nges-1,i1))*Gd(num.Nges-1,j1)*tmp; 

						tmp = (complex<double>) 0.0;
						if(inrange(APud,j2-j1,i2-i1,j1,i1)) tmp += APud(j2-j1+Lp,i2-i1+Lp)(j1-max(0,-(j2-j1)),i1-max(0,-(i2-i1)));
						if(inrange(AXud,i2-j1,j2-i1,j1,i1)) tmp += AXud(i2-j1+Lx,j2-i1+Lx)(j1-max(0,-(i2-j1)),i1-max(0,-(j2-i1)));
						//ret(j2,i2) += prefactor*conj(Gu(i1,0))*Gu(j1,0)*tmp; 
						ret(j2,i2) += prefactor*conj(Gu(num.Nges-1,i1))*Gu(num.Nges-1,j1)*tmp; 
					}
				}
			}
		}
	}
	return ret;
}

template <int mode> class Integrand_cond_vertex_cont_dyn{
	 public:
		double external_freq;
		bool spin;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		Ex_Vertex<mode> &gamma;
		Integrand_cond_vertex_cont_dyn(double external_freq_in,
		                               bool spin_in,
		                               Physics &phy_in,
		                               Numerics &num_in,
		                               Ex_Precomputation<mode> &pre_in,
		                               Substitution<mode> &sub_in,
		                               Ex_Vertex<mode> &gamma_in);
		matrix<double> select(matrix<complex<double> >  &M);
		matrix<complex<double> > operator()(double internal); 
};

template<int mode> Integrand_cond_vertex_cont_dyn<mode>::Integrand_cond_vertex_cont_dyn(double external_freq_in,
                                                                                        bool spin_in,
                                                                                        Physics &phy_in,
                                                                                        Numerics &num_in,
                                                                                        Ex_Precomputation<mode> &pre_in,
                                                                                        Substitution<mode> &sub_in,
                                                                                        Ex_Vertex<mode> &gamma_in):
                                                                                        external_freq(external_freq_in),
                                                                                        spin(spin_in),
                                                                                        phy(phy_in),
                                                                                        num(num_in),
                                                                                        pre(pre_in),
                                                                                        sub(sub_in),
                                                                                        gamma(gamma_in){
}

template<int mode> matrix<double> Integrand_cond_vertex_cont_dyn<mode>::select(matrix<complex<double> > &M){
	return matrix_to_vector(M);	
}

template<int mode> matrix<complex<double> > Integrand_cond_vertex_cont_dyn<mode>::operator()(double internal){
	//cout<<"internal="<<internal<<endl;
	double intline = sub.resu_concatenated(internal);
	double prefactor = sub.weight_concatenated(internal)*2.*sqrt(1. - (intline/2.0)*(intline/2.0) )/M_PI;
	double fp = 1./tanh(((intline + external_freq)/2. - phy.mu)/phy.T) - 1. + 2.*fermi(intline, phy.mu, phy.T);
	double fx = 1./tanh(((intline - external_freq)/2.         )/phy.T) - 1. + 2.*fermi(intline, phy.mu, phy.T);
	//cout<<"phy.T="<<phy.T<<endl;
	//cout<<"fp="<<fp<<endl;
	//cout<<"fx="<<fx<<endl;

	syma<complex<double> > Gu = pre.iGu(internal);
	syma<complex<double> > Gd = pre.iGd(internal);
	
	matrix<matrix<complex<double> > > APud = gamma.aPud.dyn_ipol(external_freq + intline, num.pos_NfbP_2mu); 
	APud = fp*imaginary_part_of_str(APud); 
	

	matrix<matrix<complex<double> > > APss;
	matrix<matrix<complex<double> > > ADss;
	matrix<matrix<complex<double> > > AXud;
	
	if(spin==1){
		APss = gamma.aPuu.dyn_ipol(external_freq + intline, num.pos_NfbP_2mu); 
		APss = fp*imaginary_part_of_str(APss); 
		ADss = gamma.aDuu.dyn_ipol(intline - external_freq, num.pos_NfbX_0); 
		ADss = fx*imaginary_part_of_str(ADss); 
		AXud = gamma.aXud.dyn_ipol(intline - external_freq, num.pos_NfbX_0); 
		AXud = fx*imaginary_part_of_str(AXud); 
	}
	else{
		APss = gamma.aPdd.dyn_ipol(external_freq + intline, num.pos_NfbP_2mu); 
		APss = fp*imaginary_part_of_str(APss); 
		ADss = gamma.aDdd.dyn_ipol(intline - external_freq, num.pos_NfbX_0); 
		ADss = fx*imaginary_part_of_str(ADss); 
		AXud = gamma.aXud.dyn_ipol(external_freq - intline, num.pos_NfbX_0); 
		AXud = fx*imaginary_part_of_str(AXud); 
	}

	int Lp = (APud.dim_r-1)/2;
	int Lx = (AXud.dim_r-1)/2;
	matrix<complex<double> > ret(num.Nges,num.Nges);
	ret = (complex<double>) 0.0;
	//P-Contribution:
	for(int l = -Lp; l<=Lp; ++l){ 
		for(int k = -Lp; k<=Lp; ++k){ 
			for(int jc = 0, j=max(-l,0); jc<num.Nges-abs(l); ++jc,++j){ 
				for(int ic = 0, i=max(-k,0); ic<num.Nges-abs(k); ++ic,++i){ 
					if(spin==1){
						ret(j+l,i+k) += prefactor*conj(Gu(num.Nges-1,i))*Gu(num.Nges-1,j)*APss(l+Lp,k+Lp)(jc,ic); 
						
						ret(j,i) += prefactor*conj(Gd(num.Nges-1,i+k))*Gd(num.Nges-1,j+l)*APud(l+Lp,k+Lp)(jc,ic); 
					}
					else{
						ret(j+l,i+k) += prefactor*conj(Gd(num.Nges-1,i))*Gd(num.Nges-1,j)*APss(l+Lp,k+Lp)(jc,ic); 
						
						ret(j+l,i+k) += prefactor*conj(Gu(num.Nges-1,i))*Gu(num.Nges-1,j)*APud(l+Lp,k+Lp)(jc,ic); 
					}
				}
			}
		}
	}
	//X-Contribution:
	for(int l = -Lx; l<=Lx; ++l){ 
		for(int k = -Lx; k<=Lx; ++k){ 
			for(int jc = 0, j=max(-l,0); jc<num.Nges-abs(l); ++jc,++j){ 
				for(int ic = 0, i=max(-k,0); ic<num.Nges-abs(k); ++ic,++i){ 
					if(spin==1){
						ret(i+k,j+l) -= prefactor*conj(Gu(num.Nges-1,i))*Gu(num.Nges-1,j)*ADss(l+Lp,k+Lp)(jc,ic); 
						
						ret(i,j) -= prefactor*conj(Gd(num.Nges-1,i+k))*Gd(num.Nges-1,j+l)*AXud(l+Lp,k+Lp)(jc,ic); 
					}
					else{
						ret(i+k,j+l) -= prefactor*conj(Gd(num.Nges-1,i))*Gd(num.Nges-1,j)*ADss(l+Lp,k+Lp)(jc,ic); 
						
						ret(i+k,j+l) += prefactor*conj(Gu(num.Nges-1,i))*Gu(num.Nges-1,j)*AXud(l+Lp,k+Lp)(jc,ic); 
					}
				}
			}
		}
	}

	return ret;
}



//Is not needed for #LONG_RANGE_EXTRAPOLATION==0
//template<int mode> class Precomputation_cond_vertex_cont{
// 	public:
//		static const double eps = 1e-10; 
//		static const double delta = 1e-10; 
//		Physics &phy;	
//		Numerics &num;
//		Ex_Precomputation<mode> &pre;
//		Substitution<mode> &sub;
//		int N_cond_pre;
//		matrix<double> freq_pre;
//		Precomputation_cond_vertex_cont(Physics &phy_in,
//		                                Numerics &num_in,
//		                                Ex_Precomputation<mode> &pre_in,
//		                                Substitution<mode> &sub_in,
//		                                int N_cond_pre);
//		std::pair<matrix<matrix<complex<double> > >, matrix<matrix<complex<double> > > > preintegrate(double external_freq, bool spin);
//};
//
//template<int mode> Precomputation_cond_vertex_cont<mode>::Precomputation_cond_vertex_cont(Physics &phy_in,
//                                                                                          Numerics &num_in,
//                                                                                          Ex_Precomputation<mode> &pre_in,
//                                                                                          Substitution<mode> &sub_in,
//                                                                                          int N_cond_pre_in):
//                                                                                          phy(phy_in),
//                                                                                          num(num_in),
//                                                                                          pre(pre_in),
//                                                                                          sub(sub_in),
//                                                                                          N_cond_pre(N_cond_pre_in),
//                                                                                          freq_pre(N_cond_pre_in){
//	freq_pre = linspace(N_cond_pre,-2.+eps,2.-eps);
//}
//
//template<int mode> std::pair<matrix<matrix<complex<double> > >, matrix<matrix<complex<double> > > > Precomputation_cond_vertex_cont<mode>::preintegrate(double external_freq, bool spin){
//	
//	matrix<matrix<complex<double> > > gp(N_cond_pre);
//	matrix<matrix<complex<double> > > gx(N_cond_pre);
//	matrix<matrix<complex<double> > > retp(N_cond_pre);
//	matrix<matrix<complex<double> > > retx(N_cond_pre);
//	syma<complex<double> > G; 
//	//Set integrand values:
//	for(int i=0; i<N_cond_pre; ++i){
//		double intline   = freq_pre(i);
//		double internal  = sub.subst_concatenated(intline);
//		if(spin==1){
//			G = pre.iGu(internal);
//		}
//		else{
//			G = pre.iGd(internal);
//		}
//		double fp = 1./tanh( ((intline + external_freq)/2 - phy.mu)/phy.T ) -1. + 2.*fermi(intline, phy.mu, phy.T);
//		double fx = 1./tanh( ((intline - external_freq)/2         )/phy.T ) -1. + 2.*fermi(intline, phy.mu, phy.T);
//		double prefactor = 2.*sqrt(1. - (intline/2.0)*(intline/2.0) )/M_PI;
//		gp(i).resize(num.Nges,num.Nges);
//		gx(i).resize(num.Nges,num.Nges);
//		for(int j1=0; j1<num.Nges; ++j1){
//			for(int j2=0; j2<num.Nges; ++j2){
//				gp(i)(j1,j2) = prefactor*conj(G(j1,0))*G(j2,0)*fp;
//				gx(i)(j1,j2) = prefactor*conj(G(j1,0))*G(j2,0)*fx;
//			}
//		}
//	}
//	//Preintegrate
//	matrix<complex<double> > sum_p(num.Nges,num.Nges);
//	matrix<complex<double> > sum_x(num.Nges,num.Nges);
//	sum_p = (complex<double>) 0.0;
//	sum_x = (complex<double>) 0.0;
//	for(int i=0; i<N_cond_pre-1; ++i){
//		double freq_diff = freq_pre(i+1) - freq_pre(i);
//		if(freq_diff>delta){
//			sum_p += 0.5*(gp(i+1)+gp(i))*freq_diff;
//			sum_x += 0.5*(gx(i+1)+gx(i))*freq_diff;
//		}
//		retp(i) = sum_p;
//		retx(i) = sum_x;
//	}
//	std::pair<matrix<matrix<complex<double> > >, matrix<matrix<complex<double> > > > ret = make_pair(retp,retx);
//	return ret;
//}



template<int mode, typename T> class Integrator_cond_vertex_cont{
 	public:
		static const double eps = 1e-10; 
		Physics &phy;	
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double accuracy;
		matrix<double> &additional_stops;
		Ex_Stops<mode> &stops_obj;
		Ex_Vertex<mode> &gamma;
		matrix<complex<double> > operator()(double external_freq, bool spin);
		Integrator_cond_vertex_cont(Physics &phy_in,
		                            Numerics &num_in,
		                            Ex_Precomputation<mode> &pre_in,
		                            Substitution<mode> &sub_in,
		                            double accuracy_in,
		                            matrix<double> &additional_stops_in,
		                            Ex_Stops<mode> &stops_obj_in,
		                            Ex_Vertex<mode> &gamma_in);
};

template <int mode, typename T> Integrator_cond_vertex_cont<mode,T>::Integrator_cond_vertex_cont(Physics &phy_in,
                                                                                   Numerics &num_in,
                                                                                   Ex_Precomputation<mode> &pre_in,
                                                                                   Substitution<mode> &sub_in,
                                                                                   double accuracy_in,
                                                                                   matrix<double> &additional_stops_in,
                                                                                   Ex_Stops<mode> &stops_obj_in,
                                                                                   Ex_Vertex<mode> &gamma_in):
                                                                                   phy(phy_in),
                                                                                   num(num_in),
                                                                                   pre(pre_in),
                                                                                   sub(sub_in),
                                                                                   accuracy(accuracy_in),
                                                                                   additional_stops(additional_stops_in),
                                                                                   stops_obj(stops_obj_in),
                                                                                   gamma(gamma_in){
}

template<int mode, typename T> matrix<complex<double> > Integrator_cond_vertex_cont<mode,T>::operator()(double external_freq, bool spin){
	T integrand(external_freq, spin, phy, num, pre, sub, gamma);
	matrix<double> stops = stops_obj(external_freq, additional_stops); 
	matrix<complex<double> > ret;
	ret.resize(num.Nges,num.Nges);
	ret = (complex<double>) 0.0;
	for (int i=0; i<stops.dim_c-1; i++){
		if(stops(i+1)-stops(i) > eps){
			intgk(ret,stops(i),stops(i+1),accuracy,1e-4,1e-14,integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return ret;
}


template<int mode> class Compute_Conductance{
	public:
		static const double eps=1e-7; 
		static const double delta=1e-7; 
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		Ex_Vertex<mode> &gamma;
		Compute_Conductance(Physics &phy_in,
		                    Numerics &num_in,
		                    Ex_Precomputation<mode> &pre_in,
		                    Substitution<mode> &sub_in,
		                    Ex_Vertex<mode> &gamma_in);
		double one_particle_contribution_zero_temp(bool spin, double Lambda);
		std::pair<complex<double> ,complex<double> > full_conductance(bool spin, double Lambda, double accuracy,int number_of_freq);
		double check_ward_identity(bool spin, double Lambda, double accuracy,int number_of_freq);
};

template<int mode> Compute_Conductance<mode>::Compute_Conductance(Physics &phy_in,
                                                                  Numerics &num_in,
                                                                  Ex_Precomputation<mode> &pre_in,
                                                                  Substitution<mode> &sub_in,
                                                                  Ex_Vertex<mode> &gamma_in):
                                                                  phy(phy_in),
                                                                  num(num_in),
                                                                  pre(pre_in),
                                                                  sub(sub_in),
                                                                  gamma(gamma_in){
}

template<int mode> double Compute_Conductance<mode>::one_particle_contribution_zero_temp(bool spin, double Lambda){
	double omega;
	syma<complex<double> > G;
	if(spin==1){
		omega = phy.mu - 0.5*phy.h; 	
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = gamma.ERetu_ipol(phy.mu);
		//syma<complex<double> > E = phy.hamiltonian;
		GS = green_and_single_eq_non_ps(phy.mu, E, phy.h, 1.0, Lambda);
		G=GS(0);
	}
	else{
		omega = phy.mu + 0.5*phy.h; 	
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = gamma.ERetd_ipol(phy.mu);
		GS = green_and_single_eq_non_ps(phy.mu, E, -phy.h, 1.0, Lambda);
		G=GS(0);
	}
	if(omega>=-2. && omega<=2.){
		return 2.0*(1.0 - (omega/2.0)*(omega/2.0))*pow(abs(G(num.Nges-1,0)),2);
	}
	else{
	      return 0.0;
	}
}

template<int mode> double Compute_Conductance<mode>::check_ward_identity(bool spin, double Lambda, double accuracy, int number_of_freq){
	complex<double> I(0.0,1.0);
 	//Precompute Vertex part:
	matrix<double> additional_stops(0);
	Vertex_Conductance_Stops<mode> stop_obj(phy,sub,Lambda);
	//Integrator_cond_vertex_cont<mode,Integrand_cond_vertex_cont_naive<mode> > integrator_vertex(phy,num,pre,sub,accuracy,additional_stops,stop_obj,gamma);
	Integrator_cond_vertex_cont<mode,Integrand_cond_vertex_cont_dyn<mode> > integrator_vertex(phy,num,pre,sub,accuracy,additional_stops,stop_obj,gamma); //Takes only #LONG_RANGE_EXTRAPOLATION==0 into account
	matrix<double> frequencies = linspace(number_of_freq,-2.0,2.0); 
	matrix<matrix<complex<double> > > vertex_contribution_for_ward(number_of_freq);
	//#pragma omp parallel for
	for(int i=0; i<number_of_freq; ++i){
		cout<<"i="<<i<<endl;
		double intline = frequencies(i);
		double internal = sub.subst_concatenated(intline);
		matrix<complex<double> > vertex_contribution =  integrator_vertex(intline,spin);
		////Begin Reproduce errors of old code
		//matrix<complex<double> > tmp2 =  integrator_vertex(intline,spin);
		//syma<complex<double> > vertex_contribution =  0.5*matrix_to_syma(tmp2);
		////End Reproduce errors of old code
		vertex_contribution_for_ward(i) =  vertex_contribution;
	}
	//Check ward identity:
	double diff_ward=0.0;
	for(int i=0; i<number_of_freq;++i){
		double intline = frequencies(i);
		syma<complex<double> > E;
		if(spin==1){
			E = gamma.ERetu_ipol(intline);
		}
		else{
			E = gamma.ERetd_ipol(intline);
		}
		syma<complex<double> > A = I*(E - E.conj());
		matrix<complex<double> > Am = syma_to_matrix(A); 
		matrix<complex<double> > B = vertex_contribution_for_ward(i);
		matrix<complex<double> > C =B;
		for(int i=0; i<C.dim_r; ++i){
			for(int j=0; j<C.dim_r; ++j){
				C(i,j)+=B(C.dim_r-1-i,C.dim_r-1-j);
			}
		}
		cout<<"abs(Am-C)="<<abs(Am-C)<<endl;
		diff_ward = max(diff_ward,abs(Am-C));
		cout<<"abs(Am)="<<abs(Am)<<endl;
		cout<<"abs(C-C.transp())="<<abs(C-C.transp())<<endl;
		cout<<"abs(B-B.transp())="<<abs(B-B.transp())<<endl;
	}
	cout<<"diff_ward="<<diff_ward<<endl;
	return diff_ward;	
}

template<int mode> std::pair<complex<double> ,complex<double> > Compute_Conductance<mode>::full_conductance(bool spin, double Lambda, double accuracy, int number_of_freq){
	complex<double> I(0.0,1.0);
 	//Precompute Vertex part:
	matrix<double> additional_stops(0);
	Vertex_Conductance_Stops<mode> stop_obj(phy,sub,Lambda);
	//Integrator_cond_vertex_cont<mode,Integrand_cond_vertex_cont_naive<mode> > integrator_vertex(phy,num,pre,sub,accuracy,additional_stops,stop_obj,gamma);
	Integrator_cond_vertex_cont<mode,Integrand_cond_vertex_cont_dyn<mode> > integrator_vertex(phy,num,pre,sub,accuracy,additional_stops,stop_obj,gamma); //Takes only #LONG_RANGE_EXTRAPOLATION==0 into account
	matrix<double> frequencies = linspace(number_of_freq,max(-2.0+eps,phy.mu-15.*phy.T),min(2.0-eps,phy.mu+15.*phy.T)); 
	matrix<matrix<complex<double> > > vertex_contribution_for_ward(number_of_freq);
	matrix<complex<double> > vertex_part(number_of_freq);
	matrix<complex<double> > vertex_part_2(number_of_freq);
	matrix<complex<double> > one_particle_part(number_of_freq);
	matrix<double> norm_part(number_of_freq);
	#pragma omp parallel for
	for(int i=0; i<number_of_freq; ++i){
		cout<<"i="<<i<<endl;
		double intline = frequencies(i);
		double internal = sub.subst_concatenated(intline);
		double df = -1./phy.T*pow(fermi(intline,phy.mu,phy.T),2)*exp((intline-phy.mu)/phy.T); 
		matrix<complex<double> > vertex_contribution =  integrator_vertex(intline,spin);
		////Begin Reproduce errors of old code
		//cout<<"Old code vertex correction"<<endl;
		//matrix<complex<double> > tmp9 =  integrator_vertex(intline,spin);
		//syma<complex<double> > vertex_contribution =  0.5*matrix_to_syma(tmp9);
		////End Reproduce errors of old code
		vertex_contribution_for_ward(i) =  vertex_contribution;
		matrix<complex<double> > G;
		double omega;
		if(spin==1){
			omega = phy.mu - 0.5*phy.h; 	
			G = pre.iGu(internal); 
		}
		else{
			omega = phy.mu + 0.5*phy.h; 	
			G = pre.iGd(internal); 
		}
		double g = 2.*sqrt(1-(intline/2.)*(intline/2.));
		matrix<complex<double> > tmp = G.conj()*vertex_contribution*G; 
		vertex_part(i) = df*g*tmp(0,0);
		//Alternative two particle contribution via ward identity:
		syma<complex<double> > E;
		if(spin==1){
			E = gamma.ERetu_ipol(intline);
		}
		else{
			E = gamma.ERetd_ipol(intline);
		}
		syma<complex<double> > tmp2 = I*(E - E.conj());
		matrix<complex<double> > tmp3 = G*tmp2*G.conj(); 
		vertex_part_2(i) = df*g*(tmp3(num.Nges-1,num.Nges-1) - tmp(num.Nges-1,num.Nges-1) );
		one_particle_part(i) = df*g*g*pow(abs(G(num.Nges-1,0)),2);
		norm_part(i) = df;
	}
	//Check ward identity:
	double diff_ward=0.0;
	for(int i=0; i<number_of_freq;++i){
		double intline = frequencies(i);
		syma<complex<double> > E;
		if(spin==1){
			E = gamma.ERetu_ipol(intline);
		}
		else{
			E = gamma.ERetd_ipol(intline);
		}
		syma<complex<double> > A = I*(E - E.conj());
		matrix<complex<double> > Am = syma_to_matrix(A); 
		matrix<complex<double> > B = vertex_contribution_for_ward(i);
		matrix<complex<double> > C =B;
		for(int i=0; i<C.dim_r; ++i){
			for(int j=0; j<C.dim_r; ++j){
				C(i,j)+=B(C.dim_r-1-i,C.dim_r-1-j);
			}
		}
		cout<<"abs(Am-C)="<<abs(Am-C)<<endl;
		diff_ward = max(diff_ward,abs(Am-C));
		cout<<"abs(Am)="<<abs(Am)<<endl;
		cout<<"abs(C-C.transp())="<<abs(C-C.transp())<<endl;
		cout<<"abs(B-B.transp())="<<abs(B-B.transp())<<endl;
	}
	cout<<"diff_ward="<<diff_ward<<endl;
	//End check ward identity
	//Integration:
	complex<double> cond_op=0.0;
	complex<double> cond_tp=0.0;
	complex<double> cond_tp_2=0.0;
	double norm=0;
	for(int i=0; i<number_of_freq-1; ++i){
		double freq_diff = frequencies(i+1)-frequencies(i);
		if(freq_diff>delta){
			cond_op -= 0.5*(one_particle_part(i+1) + one_particle_part(i))*freq_diff;
			cond_tp -= 0.5*(vertex_part(i+1) + vertex_part(i))*freq_diff;
			cond_tp_2 -= 0.5*(vertex_part_2(i+1) + vertex_part_2(i))*freq_diff;
			norm -= 0.5*(norm_part(i+1) + norm_part(i))*freq_diff;
		}
	}
	cout<<"conductance: norm="<<norm<<endl;
	cout<<"cond_op="<<cond_op<<endl;
	cout<<"cond_tp="<<cond_tp<<endl;
	cout<<"cond_tp_2="<<cond_tp_2<<endl;
	std::pair<complex<double> ,complex<double> > ret = make_pair(cond_op,cond_tp);
	return ret;
}

#endif
