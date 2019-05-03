#ifndef Ex_SELF_11122018
#define Ex_SELF_11122018

#include <integrate_new.h>
#include "Ex_stops.h"
#include "Ex_Precomputation.h"
#include "Ex_functions.h"
#include "Ex_Vertex.h"
#include "Barevertex.h"

using namespace std;

template <int mode> class Integrand_S{
	 public:
		bool spin;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		Integrand_S(bool spin_in,
		            Physics &phy_in,
		            Numerics &num_in,
		            Ex_Precomputation<mode> &pre_in,
		            Substitution<mode> &sub_in,
		            double measure_flow_in);
		matrix<double> select(syma<complex<double> > &M);	
		syma<complex<double> > compute_ret(double internal);
		syma<complex<double> > compute_kel(double internal);
		virtual syma<complex<double> >  operator()(double internal)=0;
};

template<int mode> Integrand_S<mode>::Integrand_S(bool spin_in,
                                                  Physics &phy_in,
                                                  Numerics &num_in,
                                                  Ex_Precomputation<mode> &pre_in,
                                                  Substitution<mode> &sub_in,
                                                  double measure_flow_in):
                                                  spin(spin_in),
                                                  phy(phy_in),
                                                  num(num_in),
                                                  pre(pre_in),
                                                  sub(sub_in),
                                                  measure_flow(measure_flow_in){
}

template<int mode> matrix<double> Integrand_S<mode>::select(syma<complex<double> > &M){
	return matrix_to_vector(M);	
}

template<int mode> syma<complex<double> > Integrand_S<mode>::compute_ret(double internal){
	double f = measure_flow/M_PI; 
	#if(H_EQUAL_ZERO==0)
		if(spin==1){
			return f*pre.iSu(internal); 
		}
		else{
			return f*pre.iSd(internal); 
		}
	#else
		return f*pre.iSu(internal); 
	#endif
}

template<int mode> syma<complex<double> > Integrand_S<mode>::compute_kel(double internal){
	double intline = sub.resu_concatenated(internal);
	double f = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI; 
	#if(H_EQUAL_ZERO==0)
		if(spin==1){
			return f*pre.iSu(internal).imag(); 
		}
		else{
			return f*pre.iSd(internal).imag(); 
		}
	#else
		return f*pre.iSu(internal).imag(); 
	#endif
}

template<int mode> class Integrand_S_ret : public Integrand_S<mode>{
	public:
		using Integrand_S<mode>::Integrand_S; 
		syma<complex<double> > operator()(double internal){return this->compute_ret(internal);}
};

template<int mode> class Integrand_S_kel : public Integrand_S<mode>{
	public:
		using Integrand_S<mode>::Integrand_S; 
		syma<complex<double> > operator()(double internal){return this->compute_kel(internal);}
};



template<int mode, typename T> class Compute_S_semi_static{
	public:
		static const double eps = 1e-10; 
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		double accuracy;
		matrix<double> &additional_stops;
		Ex_Stops<mode> &stops_obj;
		Compute_S_semi_static(Physics &phy_in,
		                      Numerics &num_in,
		                      Ex_Precomputation<mode> &pre_in,
		                      Substitution<mode> &sub_in,
		                      double measure_flow_in,
		                      double accuracy_in,
		                      matrix<double> &additional_stops_in,
		                      Ex_Stops<mode> &stops_obj_in);
		matrix<syma<complex<double> > > operator()(bool spin, matrix<double> &freq_steps);
		matrix<syma<complex<double> > > get_ret_preintegrated(bool spin, matrix<double> &freq_steps);
		matrix<syma<complex<double> > > get_kel_preintegrated(bool spin, matrix<double> &freq_steps);
};

template<int mode, typename T> Compute_S_semi_static<mode,T>::Compute_S_semi_static(Physics &phy_in,
                                                                                    Numerics &num_in,
                                                                                    Ex_Precomputation<mode> &pre_in,
                                                                                    Substitution<mode> &sub_in,
                                                                                    double measure_flow_in,
                                                                                    double accuracy_in,
                                                                                    matrix<double> &additional_stops_in,
                                                                                    Ex_Stops<mode> &stops_obj_in):
                                                                                    phy(phy_in),
                                                                                    num(num_in),
                                                                                    pre(pre_in),
                                                                                    sub(sub_in),
                                                                                    measure_flow(measure_flow_in),
                                                                                    accuracy(accuracy_in),
                                                                                    additional_stops(additional_stops_in),
                                                                                    stops_obj(stops_obj_in){
}

template<int mode, typename T> matrix<syma<complex<double> > > Compute_S_semi_static<mode,T>::operator()(bool spin, matrix<double> &freq_steps){
	int N_steps = freq_steps.dim_c;
	matrix<syma<complex<double> > > ret(N_steps-1);
	syma<complex<double> > sum(num.Nges);
	sum = (complex<double>) 0.0;
	T Integrand(spin, phy, num, pre, sub, measure_flow);
	matrix<double> stops = stops_obj(0.0, additional_stops); 
	for(int j=0; j<N_steps-1; ++j){
		matrix<double> stops_interval = get_stops_in_interval(stops, min(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))),  max(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))));
		for (int i=0; i<stops_interval.dim_c-1; i++){
			if(stops_interval(i+1)-stops_interval(i) > eps){
				intgk(sum,stops_interval(i),stops_interval(i+1),accuracy,1e-4,1e-14,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
			}
		}
		ret(j)=sum;
	}
	return ret;
}

template<int mode, typename T> matrix<syma<complex<double> > > Compute_S_semi_static<mode,T>::get_ret_preintegrated(bool spin, matrix<double> &freq_steps){
	int N_steps = freq_steps.dim_c;
	matrix<syma<complex<double> > > ret(N_steps-1);
	syma<complex<double> > sum(num.Nges);
	sum = (complex<double>) 0.0;
	matrix<double> stops = stops_obj(0.0, additional_stops); 
	if(spin==1){
		for(int j=0; j<N_steps-1; ++j){
			matrix<double> stops_interval = get_stops_in_interval(stops, min(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))),  max(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))));
			for (int i=0; i<stops_interval.dim_c-1; i++){
				if(stops_interval(i+1)-stops_interval(i) > eps){
					sum += pre.iSu_ret_integrated(stops_interval(i+1)) - pre.iSu_ret_integrated(stops_interval(i));
				}
			}
			ret(j)=(measure_flow/M_PI)*sum;
			//ret(j)=sum;
		}
	}
	else{
		for(int j=0; j<N_steps-1; ++j){
			matrix<double> stops_interval = get_stops_in_interval(stops, min(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))),  max(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))));
			for (int i=0; i<stops_interval.dim_c-1; i++){
				if(stops_interval(i+1)-stops_interval(i) > eps){
					sum += pre.iSd_ret_integrated(stops_interval(i+1)) - pre.iSd_ret_integrated(stops_interval(i));
				}
			}
			ret(j)=(measure_flow/M_PI)*sum;
			//ret(j)=sum;
		}
	}
	return ret;
}

template<int mode, typename T> matrix<syma<complex<double> > > Compute_S_semi_static<mode,T>::get_kel_preintegrated(bool spin, matrix<double> &freq_steps){
	int N_steps = freq_steps.dim_c;
	matrix<syma<complex<double> > > ret(N_steps-1);
	syma<complex<double> > sum(num.Nges);
	sum = (complex<double>) 0.0;
	matrix<double> stops = stops_obj(0.0, additional_stops); 
	if(spin==1){
		for(int j=0; j<N_steps-1; ++j){
			matrix<double> stops_interval = get_stops_in_interval(stops, min(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))),  max(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))));
			for (int i=0; i<stops_interval.dim_c-1; i++){
				if(stops_interval(i+1)-stops_interval(i) > eps){
					sum += pre.iSu_kel_integrated(stops_interval(i+1)) - pre.iSu_kel_integrated(stops_interval(i));
				}
			}
			ret(j)=(measure_flow/M_PI)*sum;
			//ret(j)=sum;
		}
	}
	else{
		for(int j=0; j<N_steps-1; ++j){
			matrix<double> stops_interval = get_stops_in_interval(stops, min(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))),  max(sub.subst_concatenated(freq_steps(j)),sub.subst_concatenated(freq_steps(j+1))));
			for (int i=0; i<stops_interval.dim_c-1; i++){
				if(stops_interval(i+1)-stops_interval(i) > eps){
					sum += pre.iSd_kel_integrated(stops_interval(i+1)) - pre.iSd_kel_integrated(stops_interval(i));
				}
			}
			ret(j)=(measure_flow/M_PI)*sum;
			//ret(j)=sum;
		}
	}
	return ret;
}

template<int mode> class Compute_self_semi_static{
	public:
		static const double eps_at_inf = 1e-7;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		double accuracy;
		matrix<double> &additional_stops;
		Ex_Stops<mode> &stops_obj;
		Ex_Vertex<mode> &gamma;
		Barevertex &bare;
		matrix<matrix<int> > L_steps_p;
		matrix<matrix<int> > L_steps_x;
		Compute_self_semi_static(Physics &phy_in,
		                         Numerics &num_in,
		                         Ex_Precomputation<mode> &pre_in,
		                         Substitution<mode> &sub_in,
		                         double measure_flow_in,
		                         double accuracy_in,
		                         matrix<double> &additional_stops_in,
		                         Ex_Stops<mode> &stops_obj_in,
		                         Ex_Vertex<mode> &gamma_in,
		                         Barevertex &bare_in);
		matrix<syma<complex<double> > > full_static();
		matrix<syma<complex<double> > > semi_static(double external_freq);
};

template<int mode> Compute_self_semi_static<mode>::Compute_self_semi_static(Physics &phy_in,
                                                                            Numerics &num_in,
                                                                            Ex_Precomputation<mode> &pre_in,
                                                                            Substitution<mode> &sub_in,
                                                                            double measure_flow_in,
                                                                            double accuracy_in,
                                                                            matrix<double> &additional_stops_in,
                                                                            Ex_Stops<mode> &stops_obj_in,
                                                                            Ex_Vertex<mode> &gamma_in,
                                                                            Barevertex &bare_in):
                                                                            phy(phy_in),
                                                                            num(num_in),
                                                                            pre(pre_in),
                                                                            sub(sub_in),
                                                                            measure_flow(measure_flow_in),
                                                                            accuracy(accuracy_in),
                                                                            additional_stops(additional_stops_in),
                                                                            stops_obj(stops_obj_in),
                                                                            gamma(gamma_in),
                                                                            bare(bare_in){
	L_steps_p = determine_L_steps(num.L, num.Lp_structure, num.pos_NfbP_2mu);
	L_steps_x = determine_L_steps(num.L, num.Lx_structure, num.pos_NfbX_0);
}


template<int mode> matrix<syma<complex<double> > > Compute_self_semi_static<mode>::semi_static(double external_freq){
	matrix<matrix<double> > freq_steps_p = determine_freq_steps_shifted(external_freq, num.L, num.Lp_bounds, L_steps_p, num.wbP, 2*phy.mu);
	matrix<matrix<double> > freq_steps_x = determine_freq_steps_shifted(-external_freq, num.L, num.Lx_bounds, L_steps_x, num.wbX, 0.0);
	//Integrate different S:
	Compute_S_semi_static<mode,Integrand_S_ret<mode> > Int_ret(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj);
	Compute_S_semi_static<mode,Integrand_S_kel<mode> > Int_kel(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj);
	matrix<matrix<syma<complex<double> > > > Su_p1(2); 
	matrix<matrix<syma<complex<double> > > > Su_p2(2); 
	matrix<matrix<syma<complex<double> > > > Su_x1(2); 
	matrix<matrix<syma<complex<double> > > > Su_x2(2); 
        #if(H_EQUAL_ZERO==0)
		matrix<matrix<syma<complex<double> > > > Sd_p1(2); 
		matrix<matrix<syma<complex<double> > > > Sd_p2(2); 
		matrix<matrix<syma<complex<double> > > > Sd_x1(2); 
		matrix<matrix<syma<complex<double> > > > Sd_x2(2); 
	#else
		matrix<matrix<syma<complex<double> > > > &Sd_p1 = Su_p1; 
		matrix<matrix<syma<complex<double> > > > &Sd_p2 = Su_p2; 
		matrix<matrix<syma<complex<double> > > > &Sd_x1 = Su_x1; 
		matrix<matrix<syma<complex<double> > > > &Sd_x2 = Su_x2; 
	#endif
	for(int i=0; i<2; ++i){
		//Su_p1(i) = Int_ret(1,freq_steps_p(i)); 
		//Su_p2(i) = Int_kel(1,freq_steps_p(i)); 
		//Su_x1(i) = Int_ret(1,freq_steps_x(i)); 
		//Su_x2(i) = Int_kel(1,freq_steps_x(i)); 
		Su_p1(i) = Int_ret.get_ret_preintegrated(1,freq_steps_p(i)); 
		Su_p2(i) = Int_kel.get_kel_preintegrated(1,freq_steps_p(i)); 
		Su_x1(i) = Int_ret.get_ret_preintegrated(1,freq_steps_x(i)); 
		Su_x2(i) = Int_kel.get_kel_preintegrated(1,freq_steps_x(i)); 
		//{
		//	int rank, error;
		//	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		//	if(rank==0 && i==0){
		//		cout<<"freq_steps_p(0) ="<<freq_steps_p(i)<<endl;
		//		cout<<"measure_flow="<<measure_flow<<endl;
		//		cout<<"diff Su_p1="<<abs(Su_p1(i) - Int_ret.get_ret_preintegrated(1,freq_steps_p(i)))<<endl; 
		//		cout<<"diff Su_p2="<<abs(Su_p2(i) - Int_kel.get_kel_preintegrated(1,freq_steps_p(i)))<<endl; 
		//		cout<<"diff Su_x1="<<abs(Su_x1(i) - Int_ret.get_ret_preintegrated(1,freq_steps_x(i)))<<endl; 
		//		cout<<"diff Su_x2="<<abs(Su_x2(i) - Int_kel.get_kel_preintegrated(1,freq_steps_x(i)))<<endl; 
		//	}
		//}
        	#if(H_EQUAL_ZERO==0)
			//Sd_p1(i) = Int_ret(0,freq_steps_p(i)); 
			//Sd_p2(i) = Int_kel(0,freq_steps_p(i)); 
			//Sd_x1(i) = Int_ret(0,freq_steps_x(i)); 
			//Sd_x2(i) = Int_kel(0,freq_steps_x(i)); 
			Sd_p1(i) = Int_ret.get_ret_preintegrated(0,freq_steps_p(i)); 
			Sd_p2(i) = Int_kel.get_kel_preintegrated(0,freq_steps_p(i)); 
			Sd_x1(i) = Int_ret.get_ret_preintegrated(0,freq_steps_x(i)); 
			Sd_x2(i) = Int_kel.get_kel_preintegrated(0,freq_steps_x(i)); 
		#endif
	}
	matrix<syma<complex<double> > > ret(2);
	//Spin up component:
	matrix<complex<double> > self_up(num.Nges,num.Nges); 
	self_up = (complex<double>) 0.0; 

	matrix<double> b_pre_p =  b_prefactor_p(num.wbP,phy.T,phy.mu);
	matrix<double> b_pre_x =  b_prefactor_x(num.wbX,phy.T);
	matrix<double> b_pre_d =  b_prefactor_d(num.wbX,phy.T);
	matrix<matrix<matrix<complex<double> > > > tmp;
	#if(LONG_RANGE_EXTRAPOLATION==1)
		tmp = gamma.aDuu.static_extensions_b(num.Lx_bounds, b_pre_d); 
		tmp = flip_vector(tmp);
		self_up -= compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Su_x1);

		tmp = gamma.aXud.static_extensions_b(num.Lx_bounds, b_pre_x); 
		self_up += compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Sd_x1);
		
		tmp = gamma.aPuu.static_extensions_b(num.Lp_bounds, b_pre_p); 
		self_up += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Su_p1);
		
		tmp = gamma.aPud.static_extensions_b(num.Lp_bounds, b_pre_p); 
		self_up += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Sd_p1);
	#endif
	
	tmp = gamma.aPuu.static_extensions_a(num.Lp_bounds); 
	self_up += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Su_p2);
	
	tmp = gamma.aDuu.static_extensions_a(num.Lx_bounds); 
	tmp = flip_vector(tmp);
	self_up -= compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Su_x2);
	
	tmp = gamma.aPud.static_extensions_a(num.Lp_bounds); 
	self_up += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Sd_p2);
	
	tmp = gamma.aXud.static_extensions_a(num.Lx_bounds); 
	self_up += compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Sd_x2);
	

	ret(0) = matrix_to_syma(self_up);
	
	//Spin down component:
	#if(H_EQUAL_ZERO==0)
		matrix<complex<double> > self_down(num.Nges,num.Nges); 
		self_down = (complex<double>) 0.0; 

		tmp = gamma.aDdd.static_extensions_b(num.Lx_bounds, b_pre_d);
		tmp = flip_vector(tmp);
		self_down -= compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Sd_x1);
	
		tmp = gamma.aXud.static_extensions_b(num.Lx_bounds, b_pre_x); 
		tmp = flip_vector(tmp);
		tmp(0) = mirror_lr_str(tmp(0));
		tmp(1) = mirror_lr_str(tmp(1));
		self_down += compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Su_x1);
	
		tmp = gamma.aPdd.static_extensions_b(num.Lp_bounds, b_pre_p); 
		self_down += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Sd_p1);
	
		tmp = gamma.aPud.static_extensions_b(num.Lp_bounds, b_pre_p); 
		tmp(0) = mirror_lr_str(tmp(0));
		tmp(1) = mirror_lr_str(tmp(1));
		self_down += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Su_p1);
	
		tmp = gamma.aPdd.static_extensions_a(num.Lp_bounds); 
		self_down += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Sd_p2);
	
		tmp = gamma.aDdd.static_extensions_a(num.Lx_bounds); 
		tmp = flip_vector(tmp);
		self_down -= compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Sd_x2);
	
		tmp = gamma.aPud.static_extensions_a(num.Lp_bounds); 
		tmp(0) = mirror_lr_str(tmp(0));
		tmp(1) = mirror_lr_str(tmp(1));
		self_down += compute_ring_contributions(num.L, num.N, L_steps_p, tmp, Su_p2);
	
		tmp = gamma.aXud.static_extensions_a(num.Lx_bounds); 
		tmp = flip_vector(tmp);
		tmp(0) = conjugate_mirror_lr_str(tmp(0));
		tmp(1) = conjugate_mirror_lr_str(tmp(1));
		self_down += compute_ring_contributions(num.L, num.N, L_steps_x, tmp, Su_x2);
	
		ret(1) = matrix_to_syma(self_down);
	#else
		ret(1)=ret(0);
	#endif
	return ret;
}

template<int mode> matrix<syma<complex<double> > > Compute_self_semi_static<mode>::full_static(){
	//Compute_S_semi_static<mode,Integrand_S_kel<mode> > Int_kel(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj);
	//matrix<double> freq_steps(2);
	//freq_steps(0) = sub.resu_concatenated(-7.+eps_at_inf); 
	//freq_steps(1) = sub.resu_concatenated(+7.-eps_at_inf); 
	//matrix<syma<complex<double> > > tmp1 = Int_kel(1, freq_steps); 
	//matrix<syma<complex<double> > > tmp2 = Int_kel(0, freq_steps); 
	//syma<complex<double> > S_up = tmp1(0);
	//syma<complex<double> > S_down = tmp2(0);
	syma<complex<double> > S_up = (measure_flow/M_PI)*(pre.iSu_kel_integrated(+7.-eps_at_inf) - pre.iSu_kel_integrated(-7.+eps_at_inf));
	syma<complex<double> > S_down = (measure_flow/M_PI)*(pre.iSd_kel_integrated(+7.-eps_at_inf) -pre.iSd_kel_integrated(-7.+eps_at_inf)) ;

	matrix<syma<complex<double> > > ret(2);
	////unoptimized full_access version
	//matrix<complex<double> > self_up(num.Nges,num.Nges);
	//self_up = (complex<double>) 0.0;
	//#if(H_EQUAL_ZERO==0)
	//	matrix<complex<double> > self_down(num.Nges,num.Nges);
	//	self_down = (complex<double>) 0.0;
	//#endif
	//for(int l=-num.L; l<=num.L; ++l){
	//	for(int k=-num.L; k<=num.L; ++k){
	//		for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l); ++jc,++j){
	//			for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k); ++ic,++i){
	//				self_up(j,j+l) += (0.5*bare(j,1,i,1,j+l,1,i+k,1) + gamma.aDuu.static_str(l+num.L,-k+num.L)(jc,ic))*S_up.full_access(i+k,i);
	//				self_up(j,j+l) += (0.5*bare(j,1,i,0,j+l,1,i+k,0) + gamma.aDud.static_str(l+num.L,-k+num.L)(jc,ic))*S_down.full_access(i+k,i);

	//				#if(H_EQUAL_ZERO==0)
	//					self_down(j,j+l) += (0.5*bare(j,0,i,0,j+l,0,i+k,0) + gamma.aDdd.static_str(l+num.L,-k+num.L)(jc,ic))*S_down.full_access(i+k,i);
	//					self_down(j,j+l) += (0.5*bare(i,1,j,0,i+k,1,j+l,0) + gamma.aDud.static_str(-k+num.L,l+num.L)(ic,jc))*S_up.full_access(i+k,i);
	//				#endif
	//			}
	//		}
	//	}
	//}

	//optimized version:
	ret(0).resize(num.Nges); 
	ret(0) = (complex<double>) 0.0;
	#if(H_EQUAL_ZERO==0)
		ret(1).resize(num.Nges); 
		ret(1) = (complex<double>) 0.0;
	#endif
	for(int l=-num.L; l<=0; ++l){
		for(int k=-num.L; k<0; ++k){
			for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l); ++jc,++j){
				for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k); ++ic,++i){
					ret(0)(j,j+l) += (0.5*bare(j,1,i,1,j+l,1,i+k,1) + gamma.aDuu.static_str(l+num.L,-k+num.L)(jc,ic))*S_up(i,i+k);
					ret(0)(j,j+l) += (0.5*bare(j,1,i,0,j+l,1,i+k,0) + gamma.aDud.static_str(l+num.L,-k+num.L)(jc,ic))*S_down(i,i+k);
					#if(H_EQUAL_ZERO==0)
						ret(1)(j,j+l) += (0.5*bare(j,0,i,0,j+l,0,i+k,0) + gamma.aDdd.static_str(l+num.L,-k+num.L)(jc,ic))*S_down(i,i+k);
						ret(1)(j,j+l) += (0.5*bare(i,1,j,0,i+k,1,j+l,0) + gamma.aDud.static_str(-k+num.L,l+num.L)(ic,jc))*S_up(i,i+k);
					#endif
				}
			}
		}
	}
	for(int l=-num.L; l<=0; ++l){
		for(int k=0; k<=num.L; ++k){
			for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l); ++jc,++j){
				for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k); ++ic,++i){
					ret(0)(j,j+l) += (0.5*bare(j,1,i,1,j+l,1,i+k,1) + gamma.aDuu.static_str(l+num.L,-k+num.L)(jc,ic))*S_up(i+k,i);
					ret(0)(j,j+l) += (0.5*bare(j,1,i,0,j+l,1,i+k,0) + gamma.aDud.static_str(l+num.L,-k+num.L)(jc,ic))*S_down(i+k,i);
					#if(H_EQUAL_ZERO==0)
						ret(1)(j,j+l) += (0.5*bare(j,0,i,0,j+l,0,i+k,0) + gamma.aDdd.static_str(l+num.L,-k+num.L)(jc,ic))*S_down(i+k,i);
						ret(1)(j,j+l) += (0.5*bare(i,1,j,0,i+k,1,j+l,0) + gamma.aDud.static_str(-k+num.L,l+num.L)(ic,jc))*S_up(i+k,i);
					#endif
				}
			}
		}
	}

	//ret(0) = matrix_to_syma(self_up);
	#if(H_EQUAL_ZERO==0)
		//ret(1) = matrix_to_syma(self_down);
	#else
		ret(1) =ret(0);
	#endif
	return ret;
}

template<int mode> class Integrator_self_semi_static: public Compute_self_semi_static<mode>{
	public:
		using Compute_self_semi_static<mode>::Compute_self_semi_static;
		matrix<complex<double> > operator()(double external_freq);
		int dim_r(double &job){return this->num.Nges;}
		int dim_c(double &job){return this->num.Nges+1;}
		int volume(double &job){return (this->num.Nges)*(this->num.Nges+1);}
};

template<int mode> matrix<complex<double> > Integrator_self_semi_static<mode>::operator()(double external_freq){
	matrix<syma<complex<double> > > A = this->semi_static(external_freq);
	matrix<complex<double> > ret = symas_combine_to_matrix(A(0),A(1));
	return ret; 
}

template <int mode> class Integrand_self{
	 public:
		double external_freq;
		bool spin;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		Ex_Vertex<mode> &gamma;
		syma<complex<double> > S;
		matrix<matrix<complex<double> > > a1;
		matrix<matrix<complex<double> > > a2;
		double f1;
		double f2; 
		double nf;
		Integrand_self(double external_freq_in,
		               bool spin_in,
		               Physics &phy_in,
                               Numerics &num_in,
                               Ex_Precomputation<mode> &pre_in,
                               Substitution<mode> &sub_in,
                               double measure_flow_in,
		               Ex_Vertex<mode> &gamma_in);
		matrix<double> select(syma<complex<double> > &M);	
		void set_same_spin(double internal);
		void set_opposite_spin(double internal);
		syma<complex<double> >  compute(double internal);
		virtual syma<complex<double> >  operator()(double internal)=0;
};

template<int mode> Integrand_self<mode>::Integrand_self(double external_freq_in,
                                                        bool spin_in,
                                                        Physics &phy_in,
                                                        Numerics &num_in,
                                                        Ex_Precomputation<mode> &pre_in,
                                                        Substitution<mode> &sub_in,
                                                        double measure_flow_in,
                                                        Ex_Vertex<mode> &gamma_in):
                                                        external_freq(external_freq_in),
                                                        spin(spin_in),
                                                        phy(phy_in),
                                                        num(num_in),
                                                        pre(pre_in),
                                                        sub(sub_in),
                                                        measure_flow(measure_flow_in),
                                                        gamma(gamma_in),
                                                        S(num.Nges){
}

template<int mode> matrix<double> Integrand_self<mode>::select(syma<complex<double> > &M){
	return matrix_to_vector(M);	
}

template<int mode> void Integrand_self<mode>::set_same_spin(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff = external_freq - intline;
	double sum  = external_freq + intline;
	f1 = measure_flow/(M_PI*tanh((diff/2.         )/phy.T)); 
	f2 = measure_flow/(M_PI*tanh((sum/2.  - phy.mu)/phy.T)); 
	nf  = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  

	if(spin==1){
		S = pre.iSu(internal);
		a1 = -gamma.aDuu.dyn_ipol(diff, num.pos_NfbX_0); 
		a2 =  gamma.aPuu.dyn_ipol(sum, num.pos_NfbP_2mu); 
	}
	else{
		S = pre.iSd(internal);
		a1 = -gamma.aDdd.dyn_ipol(diff, num.pos_NfbX_0); 
		a2 =  gamma.aPdd.dyn_ipol(sum, num.pos_NfbP_2mu); 
	}
}

template<int mode> void Integrand_self<mode>::set_opposite_spin(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff = intline - external_freq;
	double sum  = external_freq + intline;
	f1 = -measure_flow/(M_PI*tanh((diff/2.         )/phy.T)); 
	f2 = measure_flow/(M_PI*tanh((sum/2.  - phy.mu)/phy.T)); 
	nf = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  

	if(spin==1){
		S = pre.iSd(internal);
		a1 = gamma.aXud.dyn_ipol(diff, num.pos_NfbX_0); 
		a2 = gamma.aPud.dyn_ipol(sum, num.pos_NfbP_2mu); 
	}
	else{
		S = pre.iSu(internal);
		matrix<matrix<complex<double> > > tmp = gamma.aXud.dyn_ipol(-diff, num.pos_NfbX_0);
		a1 = conjugate_mirror_lr_str(tmp); 
		tmp = gamma.aPud.dyn_ipol(sum, num.pos_NfbP_2mu);
		a2 = mirror_lr_str(tmp); 
	}
}

////Full access can still be optimized away!
//template<int mode> syma<complex<double> > Integrand_self<mode>::compute(double internal){
//	syma<complex<double> > ret(num.Nges); 
//	for(int j=0; j<num.Nges; ++j){
//		for(int i=0; i<=j; ++i){
//			ret(j,i) = (complex<double>) 0.0;
//			int L1ges = a1.dim_c;
//			int L1 = (L1ges-1)/2;
//			for(int l=max(-L1,-j); l<=min(L1,num.Nges-j-1); ++l){ 
//				for(int k=max(-L1,-i); k<=min(L1,num.Nges-i-1); ++k){ 
//					int jc = j - max(0,-l); 
//					int ic = i - max(0,-k); 
//					ret(j,i) += f1*imag(a1(l+L1,k+L1)(jc,ic))*S.full_access(j+l,i+k) 
//					           +a1(l+L1,k+L1)(jc,ic)*nf*imag(S.full_access(j+l,i+k));
//					//ret(j,i) += //f1*imag(a1(l+L1,k+L1)(jc,ic))*S.full_access(j+l,i+k);
//					//           +a1(l+L1,k+L1)(jc,ic)*nf*imag(S.full_access(j+l,i+k));
//				}
//			}
//			int L2ges = a2.dim_c;
//			int L2 = (L2ges-1)/2;
//			for(int l=max(-L2,-j); l<=min(L2,num.Nges-j-1); ++l){ 
//				for(int k=max(-L2,-i); k<=min(L2,num.Nges-i-1); ++k){ 
//					int jc = j - max(0,-l); 
//					int ic = i - max(0,-k); 
//					ret(j,i) += f2*imag(a2(l+L2,k+L2)(jc,ic))*conj(S.full_access(j+l,i+k)) 
//					           +a2(l+L2,k+L2)(jc,ic)*nf*imag(S.full_access(j+l,i+k));
//					//ret(j,i) +=  f2*imag(a2(l+L2,k+L2)(jc,ic))*conj(S.full_access(j+l,i+k));                                                          
//					//          //+a2(l+L2,k+L2)(jc,ic)*nf*imag(S.full_access(j+l,i+k));
//				}
//			}
//		}
//	}
//	return ret;
//}

//Optimized version:
template<int mode> syma<complex<double> > Integrand_self<mode>::compute(double internal){
	syma<complex<double> > ret(num.Nges); 
	ret = (complex<double>) 0.0;  
	int L1ges = a1.dim_c;
	int L1 = (L1ges-1)/2;
	int L2ges = a2.dim_c;
	int L2 = (L2ges-1)/2;
	for(int l=-L1; l<=L1; ++l){
		int j_min = max(-l,0); 
		int j_max = min(num.Nges, num.Nges-l); 
		for(int k=-L1; k<=L1; ++k){
		 	int i_min = max(-k,0); 
			for(int j=j_min, jc=0; j<j_max; ++j, ++jc){
		 		int i_max = min(min(num.Nges,num.Nges-k),j+1); 
			 	int i_inter = max(i_min,min(i_max,j+l-k));
				int ic=0;
				for(int i=i_min; i<i_inter; ++i,++ic){
					ret(j,i) += f1*imag(a1(l+L1,k+L1)(jc,ic))*S(j+l,i+k) 
					           +a1(l+L1,k+L1)(jc,ic)*nf*imag(S(j+l,i+k));
				}
				for(int i=i_inter; i<i_max; ++i, ++ic){
					ret(j,i) += f1*imag(a1(l+L1,k+L1)(jc,ic))*S(i+k,j+l) 
					           +a1(l+L1,k+L1)(jc,ic)*nf*imag(S(i+k,j+l));
				}
			}
		}
	}
	for(int l=-L2; l<=L2; ++l){
		int j_min = max(-l,0); 
		int j_max = min(num.Nges, num.Nges-l); 
		for(int k=-L2; k<=L2; ++k){
		 	int i_min = max(-k,0); 
			for(int j=j_min, jc=0; j<j_max; ++j, ++jc){
		 		int i_max = min(min(num.Nges,num.Nges-k),j+1); 
			 	int i_inter = max(i_min,min(i_max,j+l-k));
				int ic=0;
				for(int i=i_min; i<i_inter; ++i,++ic){
					ret(j,i) += f2*imag(a2(l+L2,k+L2)(jc,ic))*conj(S(j+l,i+k)) 
					           +a2(l+L2,k+L2)(jc,ic)*nf*imag(S(j+l,i+k));
				}
				for(int i=i_inter; i<i_max; ++i, ++ic){
					ret(j,i) += f2*imag(a2(l+L2,k+L2)(jc,ic))*conj(S(i+k,j+l)) 
					           +a2(l+L2,k+L2)(jc,ic)*nf*imag(S(i+k,j+l));
				}
			}
		}
	}
	return ret;
}



template<int mode> class Integrand_self_same_spin: public Integrand_self<mode>{
	public:
		using Integrand_self<mode>::Integrand_self;
		syma<complex<double> > operator()(double internal){this->set_same_spin(internal); return this->compute(internal);} 
};

template<int mode> class Integrand_self_opposite_spin: public Integrand_self<mode>{
	public:
		using Integrand_self<mode>::Integrand_self;
		syma<complex<double> > operator()(double internal){this->set_opposite_spin(internal); return this->compute(internal);} 
};






template<int mode> class Integrand_self_naive: public Integrand_self<mode>{
	public:
		using Integrand_self<mode>::external_freq; 
		using Integrand_self<mode>::spin; 
		using Integrand_self<mode>::phy; 
		using Integrand_self<mode>::num; 
		using Integrand_self<mode>::pre; 
		using Integrand_self<mode>::sub; 
		using Integrand_self<mode>::measure_flow; 
		using Integrand_self<mode>::gamma; 
		matrix<double> b_pre_p;	
		matrix<double> b_pre_x;	
		matrix<double> b_pre_d;	
		Barevertex &barevertex;
		Integrand_self_naive(double external_freq_in,
		                     bool spin_in,
		                     Physics &phy_in,
                                     Numerics &num_in,
                                     Ex_Precomputation<mode> &pre_in,
                                     Substitution<mode> &sub_in,
                                     double measure_flow_in,
		                     Ex_Vertex<mode> &gamma_in,
		                     Barevertex &barevertex_in);
	syma<complex<double> > operator()(double internal);
};

template<int mode> Integrand_self_naive<mode>::Integrand_self_naive(double external_freq_in,
                                                                    bool spin_in,
                                                                    Physics &phy_in,
                                                                    Numerics &num_in,
                                                                    Ex_Precomputation<mode> &pre_in,
                                                                    Substitution<mode> &sub_in,
                                                                    double measure_flow_in,
                                                                    Ex_Vertex<mode> &gamma_in,
                                                                    Barevertex &barevertex_in):
                                                                    Integrand_self<mode>::Integrand_self(external_freq_in,
                                                                                   spin_in,
                                                                                   phy_in,
                                                                                   num_in,
                                                                                   pre_in,
                                                                                   sub_in,
                                                                                   measure_flow_in,
                                                                                   gamma_in),
                                                                    barevertex(barevertex_in){
	b_pre_p = b_prefactor_p(num.wbP,phy.T,phy.mu);	
	b_pre_x = b_prefactor_x(num.wbX,phy.T);	
	b_pre_d = b_prefactor_d(num.wbX,phy.T);	
}


template<int mode> syma<complex<double> > Integrand_self_naive<mode>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff_d = external_freq - intline; 
	double diff_x = intline - external_freq; 
	double sum = external_freq + intline;
	double f1 = measure_flow/M_PI;
	double f2 = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI; 
	syma<complex<double> > Sru = pre.iSu(internal); 
	syma<complex<double> > Srd = pre.iSd(internal); 
	syma<complex<double> > Sku = f2*Sru.imag(); 
	syma<complex<double> > Skd = f2*Srd.imag(); 

	matrix<complex<double> > self(num.Nges,num.Nges);	
	self = (complex<double>) 0.0;
	
	if(spin==1){	
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
				for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l); ++jc,++j){
					for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k); ++ic,++i){
						self(j,i) -= f1*gamma.aDuu.complete_ipol_b(diff_d, num.pos_NfbX_0,num.Lx_bounds,b_pre_d)(l+num.L,k+num.L)(jc,ic)*Sru.full_access(j+l,i+k); 
						self(j,i) += f1*gamma.aXud.complete_ipol_b(diff_x, num.pos_NfbX_0,num.Lx_bounds,b_pre_x)(l+num.L,k+num.L)(jc,ic)*Srd.full_access(j+l,i+k); 
						self(j,i) += f1*gamma.aPuu.complete_ipol_b(sum, num.pos_NfbP_2mu,num.Lp_bounds,b_pre_p)(l+num.L,k+num.L)(jc,ic)*conj(Sru.full_access(j+l,i+k)); 
						self(j,i) += f1*gamma.aPud.complete_ipol_b(sum, num.pos_NfbP_2mu,num.Lp_bounds,b_pre_p)(l+num.L,k+num.L)(jc,ic)*conj(Srd.full_access(j+l,i+k)); 

						self(j,i) += gamma.aPuu.complete_ipol(sum, num.pos_NfbP_2mu,num.Lp_bounds)(l+num.L,k+num.L)(jc,ic)*Sku.full_access(j+l,i+k); 
						self(j,i) -= gamma.aDuu.complete_ipol(diff_d, num.pos_NfbX_0,num.Lx_bounds)(l+num.L,k+num.L)(jc,ic)*Sku.full_access(j+l,i+k); 
						
						self(j,i) += gamma.aPud.complete_ipol(sum, num.pos_NfbP_2mu,num.Lp_bounds)(l+num.L,k+num.L)(jc,ic)*Skd.full_access(j+l,i+k); 
						self(j,i) += gamma.aXud.complete_ipol(diff_x, num.pos_NfbX_0,num.Lx_bounds)(l+num.L,k+num.L)(jc,ic)*Skd.full_access(j+l,i+k); 
						
						self(j,j+l) += (0.5*barevertex(j,1,i,1,j+l,1,i+k,1) + gamma.aDuu.static_str(l+num.L,-k+num.L)(jc,ic))*Sku.full_access(i+k,i);
						self(j,j+l) += (0.5*barevertex(j,1,i,0,j+l,1,i+k,0) + gamma.aDud.static_str(l+num.L,-k+num.L)(jc,ic))*Skd.full_access(i+k,i);

					}
				}
			}
		}
		return matrix_to_syma(self);
	}
	else{
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
				for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l); ++jc,++j){
					for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k); ++ic,++i){
						self(j,i) -= f1*gamma.aDdd.complete_ipol_b(diff_d, num.pos_NfbX_0,num.Lx_bounds,b_pre_d)(l+num.L,k+num.L)(jc,ic)*Srd.full_access(j+l,i+k); 
						self(j,i) += f1*gamma.aXud.complete_ipol_b(diff_d, num.pos_NfbX_0,num.Lx_bounds,b_pre_x)(-l+num.L,-k+num.L)(jc,ic)*Sru.full_access(j+l,i+k); 
						self(j,i) += f1*gamma.aPdd.complete_ipol_b(sum, num.pos_NfbP_2mu,num.Lp_bounds,b_pre_p)(l+num.L,k+num.L)(jc,ic)*conj(Srd.full_access(j+l,i+k)); 
						self(j,i) += f1*gamma.aPud.complete_ipol_b(sum, num.pos_NfbP_2mu,num.Lp_bounds,b_pre_p)(-l+num.L,-k+num.L)(jc,ic)*conj(Sru.full_access(j+l,i+k)); 

						self(j,i) += gamma.aPdd.complete_ipol(sum, num.pos_NfbP_2mu,num.Lp_bounds)(l+num.L,k+num.L)(jc,ic)*Skd.full_access(j+l,i+k); 
						self(j,i) -= gamma.aDdd.complete_ipol(diff_d, num.pos_NfbX_0,num.Lx_bounds)(l+num.L,k+num.L)(jc,ic)*Skd.full_access(j+l,i+k); 
						
						self(j,i) += gamma.aPud.complete_ipol(sum, num.pos_NfbP_2mu,num.Lp_bounds)(-l+num.L,-k+num.L)(jc,ic)*Sku.full_access(j+l,i+k); 
						self(j,i) += conj(gamma.aXud.complete_ipol(-diff_x, num.pos_NfbX_0,num.Lx_bounds)(-l+num.L,-k+num.L)(jc,ic))*Sku.full_access(j+l,i+k); 
						
						self(j,j+l) += (0.5*barevertex(j,0,i,0,j+l,0,i+k,0) + gamma.aDdd.static_str(l+num.L,-k+num.L)(jc,ic))*Skd.full_access(i+k,i);
						self(j,j+l) += (0.5*barevertex(i,1,j,0,i+k,1,j+l,0) + gamma.aDud.static_str(-k+num.L,l+num.L)(ic,jc))*Sku.full_access(i+k,i);
					}
				}
			}
		}
		return matrix_to_syma(self);
	}
}

template<int mode, typename T> class Integrator_self{
	public:
		static const double eps = 1e-6; 
		bool spin;
		Physics &phy;	
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		Ex_Vertex<mode> &gamma;
		double accuracy;
		matrix<double> &additional_stops;
		Ex_Stops<mode> &stops_obj;
		Integrator_self(bool spin_in,
		                Physics &phy_in,
		                Numerics &num_in,
		                Ex_Precomputation<mode> &pre_in,
		                Substitution<mode> &sub_in,
		                double measure_flow_in,
		                Ex_Vertex<mode> &gamma_in,
		                double accuracy_in,
		                matrix<double> &additional_stops_in,
		                Ex_Stops<mode> &stops_obj_in);
		syma<complex<double> > operator()(double external_freq);
		int dim_r(double external_freq);	
		int dim_c(double external_freq);	
		int volume(double external_freq);	
}; 

template <int mode, typename T> Integrator_self<mode,T>::Integrator_self(bool spin_in,
                                                                         Physics &phy_in,
                                                                         Numerics &num_in,
                                                                         Ex_Precomputation<mode> &pre_in,
                                                                         Substitution<mode> &sub_in,
                                                                         double measure_flow_in,
                                                                         Ex_Vertex<mode> &gamma_in,
                                                                         double accuracy_in,
                                                                         matrix<double> &additional_stops_in,
                                                                         Ex_Stops<mode> &stops_obj_in):
                                                                         spin(spin_in),
                                                                         phy(phy_in),
                                                                         num(num_in),
                                                                         pre(pre_in),
                                                                         sub(sub_in),
                                                                         measure_flow(measure_flow_in),
                                                                         gamma(gamma_in),
                                                                         accuracy(accuracy_in),
                                                                         additional_stops(additional_stops_in),
                                                                         stops_obj(stops_obj_in){
}

template<int mode, typename T> syma<complex<double> > Integrator_self<mode,T>::operator()(double external_freq){
	matrix<double> additional_stops_ext = determine_additional_stops_dyn(external_freq, phy, num);
	matrix<double> additional_stops_tmp = concatenate_vectors(additional_stops,additional_stops_ext);
	T Integrand(external_freq, spin, phy, num, pre, sub, measure_flow, gamma);
	matrix<double> stops = stops_obj(external_freq, additional_stops_tmp); 
	//cout<<"stops="<<stops<<endl;
	syma<complex<double> > ret(num.Nges);
	ret = (complex<double>) 0.0;
	for (int i=0; i<stops.dim_c-1; i++){
		if(stops(i+1)-stops(i) > eps){
			intgk(ret,stops(i),stops(i+1),accuracy,1e-4,1e-14,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return ret;
}

template <int mode, typename T> int Integrator_self<mode,T>::dim_r(double external_freq){
	return num.Nges;
}

template <int mode, typename T> int Integrator_self<mode,T>::dim_c(double external_freq){
	return num.Nges;
}

template <int mode, typename T> int Integrator_self<mode,T>::volume(double external_freq){
	return ((num.Nges+1)*num.Nges)/2;
}



template<int mode, typename T> class Integrator_self_naive: public Integrator_self<mode,T>{
	public:
		using Integrator_self<mode,T>::eps;
		using Integrator_self<mode,T>::spin;
		using Integrator_self<mode,T>::phy;
		using Integrator_self<mode,T>::num;
		using Integrator_self<mode,T>::pre;
		using Integrator_self<mode,T>::sub;
		using Integrator_self<mode,T>::measure_flow;
		using Integrator_self<mode,T>::gamma;
		using Integrator_self<mode,T>::accuracy;
		using Integrator_self<mode,T>::additional_stops;
		using Integrator_self<mode,T>::Integrator_self;
		syma<complex<double> > operator()(double external_freq, Barevertex &barevertex);
};

template<int mode, typename T> syma<complex<double> > Integrator_self_naive<mode,T>::operator()(double external_freq, Barevertex &barevertex){
	T Integrand(external_freq, spin, phy, num, pre, sub, measure_flow, gamma, barevertex);
	matrix<double> stops = stops_obj(external_freq, additional_stops); 
	syma<complex<double> > ret(num.Nges);
	ret = (complex<double>) 0.0;
	for (int i=0; i<stops.dim_c-1; i++){
		if(stops(i+1)-stops(i) > eps){
			intgk(ret,stops(i),stops(i+1),accuracy,1e-4,1e-14,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return ret;
}

#endif
