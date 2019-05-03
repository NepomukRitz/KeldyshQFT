#ifndef Ex_BUBBLE_11122018
#define Ex_BUBBLE_11122018

#include <integrate_new.h>
#include "Ex_stops.h"
#include "Ex_Precomputation.h"
#include "Ex_functions.h"



//Memory allocation in Integrand_bubble_base could possibly be improved!


template <int mode, typename T> class Integrand_bubble{
	 public:
		//static int number_of_eval_p;
		//static int number_of_eval_x;
		double external_freq;
		int l;
		int k;
		bool spin1;
		bool spin2;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		int Ngesl;
		int Ngesk;
		syma<double> G1;
		syma<double> S1;
		syma<complex<double> > G2;
		syma<complex<double> > S2;
		Integrand_bubble(double external_freq_in,
		                 int l_in,
		                 int k_in,
		                 bool spin1_in,
		                 bool spin2_in,
		                 Physics &phy_in,
		                 Numerics &num_in,
		                 Ex_Precomputation<mode> &pre_in,
		                 Substitution<mode> &sub_in,
		                 double measure_flow_in);
		matrix<double> select(T &M);	
		void set_prop_p(double internal);
		void set_prop_x(double internal);
		matrix<complex<double> > compute_general(double internal);
		matrix<double> compute_feedback(double internal);
		syma<complex<double> > compute_complex_diag(double internal);
		virtual T operator()(double internal)=0;
};

//template<int mode, typename T> int Integrand_bubble<mode,T>::number_of_eval_p=0; 
//template<int mode, typename T> int Integrand_bubble<mode,T>::number_of_eval_x=0; 

template<int mode, typename T> Integrand_bubble<mode,T>::Integrand_bubble(double external_freq_in,
                                                                          int l_in,
                                                                          int k_in,
                                                                          bool spin1_in,
                                                                          bool spin2_in,
                                                                          Physics &phy_in,
                                                                          Numerics &num_in,
                                                                          Ex_Precomputation<mode> &pre_in,
                                                                          Substitution<mode> &sub_in,
                                                                          double measure_flow_in):
                                                                          external_freq(external_freq_in),
                                                                          l(l_in),
                                                                          k(k_in),
                                                                          spin1(spin1_in),
                                                                          spin2(spin2_in),
                                                                          phy(phy_in),
                                                                          num(num_in),
                                                                          pre(pre_in),
                                                                          sub(sub_in),
                                                                          measure_flow(measure_flow_in),
                                                                          Ngesl(num.Nges-abs(l)),
                                                                          Ngesk(num.Nges-abs(k)),
                                                                          G1(num.Nges),
                                                                          S1(num.Nges),
                                                                          G2(num.Nges),
                                                                          S2(num.Nges){
} 

template<int mode, typename T> matrix<double> Integrand_bubble<mode,T>::select(T &M){
	return matrix_to_vector(M);	
}

template<int mode, typename T> void Integrand_bubble<mode,T>::set_prop_p(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff = external_freq - intline;
	double diff_subst = sub.subst_concatenated(diff);
	double nf  = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nfm = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	#if(H_EQUAL_ZERO==0) //This could still be more optimized by taking the spin cases out of the integrand structure, and using minimal interpolation possible!
		if(spin1==spin2){
			if(spin1==1){
				G2 = pre.iGu(diff_subst);
				G1 = nfm*G2.imag();
				#if(COMPUTE_RPA_BUBBLE==0)
					S2 = pre.iSu(internal  );
				#else
					S2 = 0.5*pre.iGu(internal  )*sub.weight_concatenated(internal);
				#endif
				S1 = nf *S2.imag();
			}
			else{
				G2 = pre.iGd(diff_subst);
				G1 = nfm*G2.imag();
				#if(COMPUTE_RPA_BUBBLE==0)
					S2 = pre.iSd(internal  );
				#else
					S2 = 0.5*pre.iGd(internal  )*sub.weight_concatenated(internal);
				#endif
				S1 = nf *S2.imag();
			}
		}
		else{
			if(spin1==1){
				G2 = pre.iGd(diff_subst);
				G1 = nfm*pre.iGu(diff_subst).imag();
				#if(COMPUTE_RPA_BUBBLE==0)
					S2 = pre.iSd(internal  );
					S1 = nf *pre.iSu(internal  ).imag();
				#else
					S2 = 0.5*pre.iGd(internal  )*sub.weight_concatenated(internal);
					S1 = 0.5*nf*pre.iGu(internal  ).imag()*sub.weight_concatenated(internal);
				#endif
			}
			else{
				G2 = pre.iGu(diff_subst);
				G1 = nfm*pre.iGd(diff_subst).imag();
				#if(COMPUTE_RPA_BUBBLE==0)
					S2 = pre.iSu(internal  );
					S1 = nf *pre.iSd(internal  ).imag();
				#else
					S2 = 0.5*pre.iGu(internal  )*sub.weight_concatenated(internal);
					S1 = 0.5*nf *pre.iGd(internal  ).imag()*sub.weight_concatenated(internal);
				#endif
			}
		}
	#else
		G2 = pre.iGu(diff_subst);
		G1 = nfm*G2.imag();
		#if(COMPUTE_RPA_BUBBLE==0)
			S2 = pre.iSu(internal  );
		#else
			S2 = 0.5*pre.iGu(internal  )*sub.weight_concatenated(internal);
		#endif
		S1 = nf* S2.imag();
	#endif
}

template<int mode, typename T> void Integrand_bubble<mode,T>::set_prop_x(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff = intline - external_freq;
	double diff_subst = sub.subst_concatenated(diff);
	double sum = intline + external_freq;
	double sum_subst = sub.subst_concatenated(sum);
	double nf  = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nfm = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	double nfs = -measure_flow*(1.-2.*fermi(sum    , phy.mu, phy.T))/M_PI; 
	#if(H_EQUAL_ZERO==0) //This could still be more optimized by taking the spin cases out of the integrand structure, and using minimal interpolation possible!
		if(this->spin1==1){
			G1 = nfm*pre.iGu(diff_subst).imag();
			#if(COMPUTE_RPA_BUBBLE==0)
				S1 = nf *pre.iSu(internal  ).imag();
			#else
				S1 = 0.5*nf *pre.iGu(internal  ).imag()*sub.weight_concatenated(internal);
			#endif
		}
		else{
			G1 = nfm*pre.iGd(diff_subst).imag();
			#if(COMPUTE_RPA_BUBBLE==0)
				S1 = nf *pre.iSd(internal  ).imag();
			#else
				S1 = 0.5*nf *pre.iGd(internal  ).imag()*sub.weight_concatenated(internal);
			#endif
		}
		if(spin2==1){
			G2 = pre.iGu(sum_subst).conj();
			#if(COMPUTE_RPA_BUBBLE==0)
				S2 = pre.iSu(internal ).conj();
			#else
				S2 = 0.5*pre.iGu(internal ).conj()*sub.weight_concatenated(internal);
			#endif
		}
		else{
			G2 = pre.iGd(sum_subst).conj();
			#if(COMPUTE_RPA_BUBBLE==0)
				S2 = pre.iSd(internal ).conj();
			#else
				S2 = 0.5*pre.iGd(internal ).conj()*sub.weight_concatenated(internal);
			#endif
		}
	#else
		G2 = pre.iGu(sum_subst).conj();
		#if(COMPUTE_RPA_BUBBLE==0)
			S2 = pre.iSu(internal ).conj();
		#else
			S2 = 0.5*pre.iGu(internal ).conj()*sub.weight_concatenated(internal);
		#endif
		G1 = nfm*pre.iGu(diff_subst).imag();
		S1 = -nf*S2.imag();
	#endif
}

template<int mode, typename T> matrix<complex<double> > Integrand_bubble<mode,T>::compute_general(double internal){
	matrix<complex<double> > ret(Ngesl,Ngesk);
	for(int jc=0, j=max(0,-l), jsum=j+l; jc<Ngesl; ++jc, ++j, ++jsum){  
		int ic=0;
		int i=max(0,-k);
		int isum=i+k;
		for( ; ic<Ngesk && i<j; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(j,i)*S2(jsum,isum) + G2(jsum,isum)*S1(j,i); 
		}
		for( ; ic<Ngesk && i<jsum-k; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(i,j)*S2(jsum,isum) + G2(jsum,isum)*S1(i,j); 
		}
		for( ; ic<Ngesk; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(i,j)*S2(isum,jsum) + G2(isum,jsum)*S1(i,j); 
		}
	}
	return ret;

} 

template<int mode, typename T> matrix<double> Integrand_bubble<mode,T>::compute_feedback(double internal){
	matrix<double> ret(Ngesl,Ngesk);
	for(int jc=0, j=max(0,-l), jsum=j+l; jc<Ngesl; ++jc, ++j, ++jsum){  
		int ic=0;
		int i=max(0,-k);
		int isum=i+k;
		for( ; ic<Ngesk && i<j; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(j,i)*real(S2(jsum,isum)) + real(G2(jsum,isum))*S1(j,i); 
		}
		for( ; ic<Ngesk && i<jsum-k; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(i,j)*real(S2(jsum,isum)) + real(G2(jsum,isum))*S1(i,j); 
		}
		for( ; ic<Ngesk; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(i,j)*real(S2(isum,jsum)) + real(G2(isum,jsum))*S1(i,j); 
		}
	}
	return ret;

} 

template<int mode, typename T> syma<complex<double> > Integrand_bubble<mode,T>::compute_complex_diag(double internal){
	syma<complex<double> > ret(Ngesl);
	for(int jc=0, j=max(0,-l), jsum=j+l; jc<Ngesl; ++jc, ++j, ++jsum){  
		int ic=0;
		int i=max(0,-k);
		int isum=i+k;
		for( ; i<=j; ++ic, ++i, ++isum){ 
			ret(jc,ic) = G1(j,i)*S2(jsum,isum) + G2(jsum,isum)*S1(j,i); 
		}
	}
	return ret;
}



template<int mode> class Integrand_P_bubble_general: public Integrand_bubble<mode,matrix<complex<double> > >{
	public:
                using Integrand_bubble<mode,matrix<complex<double> > >::Integrand_bubble; 
		matrix<complex<double> > operator()(double internal){this->set_prop_p(internal); return this->compute_general(internal);}
};

template<int mode> class Integrand_P_bubble_feedback: public Integrand_bubble<mode,matrix<double> >{
	public:
                using Integrand_bubble<mode,matrix<double> >::Integrand_bubble; 
		matrix<double> operator()(double internal){this->set_prop_p(internal); return this->compute_feedback(internal);}
};

template<int mode> class Integrand_P_bubble_complex_diag: public Integrand_bubble<mode,syma<complex<double> > >{
	public:
                using Integrand_bubble<mode,syma<complex<double> > >::Integrand_bubble; 
		syma<complex<double> > operator()(double internal){this->set_prop_p(internal); return this->compute_complex_diag(internal);}
};

template<int mode> class Integrand_X_bubble_general: public Integrand_bubble<mode,matrix<complex<double> > >{
	public:
                using Integrand_bubble<mode,matrix<complex<double> > >::Integrand_bubble; 
		matrix<complex<double> > operator()(double internal){this->set_prop_x(internal); return this->compute_general(internal);}
};

template<int mode> class Integrand_X_bubble_feedback: public Integrand_bubble<mode,matrix<double> >{
	public:
                using Integrand_bubble<mode,matrix<double> >::Integrand_bubble; 
		matrix<double> operator()(double internal){this->set_prop_x(internal); return this->compute_feedback(internal);}
};

template<int mode> class Integrand_X_bubble_complex_diag: public Integrand_bubble<mode,syma<complex<double> > >{
	public:
                using Integrand_bubble<mode,syma<complex<double> > >::Integrand_bubble; 
		syma<complex<double> > operator()(double internal){this->set_prop_x(internal); return this->compute_complex_diag(internal);}
};




template<int mode, typename T> class Integrator_bubble{
	public:
		static const double eps = 1e-10; 
		bool spin1;
		bool spin2;
		Physics &phy;	
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		double accuracy;
		matrix<double> &additional_stops;
		Ex_Stops<mode> &stops_obj;
		Integrator_bubble(bool spin1_in,
		                  bool spin2_in, 
		                  Physics &phy_in,
		                  Numerics &num_in,
		                  Ex_Precomputation<mode> &pre_in,
		                  Substitution<mode> &sub_in,
		                  double measure_flow_in,
		                  double accuracy_in,
		                  matrix<double> &additional_stops_in,
		                  Ex_Stops<mode> &stops_obj_in);
		typename std::result_of<T(double)>::type operator()(double external_freq, int l, int k);
		typename std::result_of<T(double)>::type operator()(matrix<double> &job);
		int dim_r(matrix<double> &job);	
		int dim_c(matrix<double> &job);	
		int volume(matrix<double> &job);	
}; 

template <int mode, typename T> Integrator_bubble<mode,T>::Integrator_bubble(bool spin1_in,
                                                                                           bool spin2_in, 
                                                                                           Physics &phy_in,
                                                                                           Numerics &num_in,
                                                                                           Ex_Precomputation<mode> &pre_in,
                                                                                           Substitution<mode> &sub_in,
                                                                                           double measure_flow_in,
                                                                                           double accuracy_in,
                                                                                           matrix<double> &additional_stops_in,
                                                                                           Ex_Stops<mode> &stops_obj_in):
                                                                                           spin1(spin1_in),
                                                                                           spin2(spin2_in),
                                                                                           phy(phy_in),
                                                                                           num(num_in),
                                                                                           pre(pre_in),
                                                                                           sub(sub_in),
                                                                                           measure_flow(measure_flow_in),
                                                                                           accuracy(accuracy_in),
                                                                                           additional_stops(additional_stops_in),
                                                                                           stops_obj(stops_obj_in){
}

template<int mode, typename T> typename std::result_of<T(double)>::type Integrator_bubble<mode,T>::operator()(double external_freq, int l, int k){
	T Integrand(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);
	matrix<double> stops = stops_obj(external_freq, additional_stops); 
	typename std::result_of<T(double)>::type ret;
	ret.resize(num.Nges-abs(l),num.Nges-abs(k));
	auto w = ret.p[0];
	using W = decltype(w);
	ret = (W) 0.0;
	for (int i=0; i<stops.dim_c-1; i++){
		if(stops(i+1)-stops(i) > eps){
			intgk(ret,stops(i),stops(i+1),accuracy,1e-4,1e-14,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return ret;
}

template<int mode, typename T> typename std::result_of<T(double)>::type Integrator_bubble<mode,T>::operator()(matrix<double> &job){
	return operator()(job(0), (int) job(1), (int) job(2));
}

template <int mode, typename T> int Integrator_bubble<mode,T>::dim_r(matrix<double> &job){
	int l = (int) job(1);
	return num.Nges - abs(l);
}
                                                         
template <int mode, typename T> int Integrator_bubble<mode,T>::dim_c(matrix<double> &job){
	int k = (int) job(2);
	return num.Nges - abs(k);
}

template <int mode, typename T> int Integrator_bubble<mode,T>::volume(matrix<double> &job){
	int l = (int) job(1);
	using V = typename std::result_of<T(double)>::type;
	if(typeid(V)==typeid(syma<complex<double> >) || typeid(V)==typeid(syma<double>)){
		return ((num.Nges-abs(l)+1)*(num.Nges-abs(l)))/2;
	}
	else{
		return dim_r(job)*dim_c(job);
	}
}

#endif
