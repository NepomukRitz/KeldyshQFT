#ifndef EX_PRECOMPUTATION_11122018
#define EX_PRECOMPUTATION_11122018

#ifndef H_EQUAL_ZERO
	#define H_EQUAL_ZERO 0
#endif

#ifndef PARITY_SYMMETRY
	#define PARITY_SYMMETRY 0
#endif

#ifndef USE_MPI_FOR_PRECOMPUTATION 
	#define USE_MPI_FOR_PRECOMPUTATION 0
#endif

#include <approxxpp.h> 
#include <basic.h>

#include "Numerics.h"
#include "Ex_mpi.h"

template<int mode> class Green_mpi_eval{
	public:
		double h;
		int Nges;
		linear_ipol_bin<syma<complex<double> > > &iE;
		double Lambda;
		Substitution<mode> &sub;
		Green_mpi_eval(double h_in,
		               int Nges_in,
		               linear_ipol_bin<syma<complex<double> > > &iE_in,
		               double Lambda_in,
		               Substitution<mode> &sub_in);
		matrix<complex<double> > operator()(double fs);
		int dim_r(double job);	
		int dim_c(double job);	
		int volume(double job);	
};

template<int mode> Green_mpi_eval<mode>::Green_mpi_eval(double h_in,
                                                        int Nges_in,
                                                        linear_ipol_bin<syma<complex<double> > > &iE_in,
                                                        double Lambda_in,
                                                        Substitution<mode> &sub_in):
                                                        h(h_in),
                                                        Nges(Nges_in),
                                                        iE(iE_in),
                                                        Lambda(Lambda_in),
                                                        sub(sub_in){
}

template<int mode> matrix<complex<double> > Green_mpi_eval<mode>::operator()(double fs){
	matrix<syma<complex<double> > > GS(2);
	syma<complex<double> > E = iE(fs);
	#if(PARITY_SYMMETRY==0)
		GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, h, 1.0, Lambda);
	#else
		GS = green_and_single_Req_ps(sub.resu_concatenated(fs), E, h, 1.0, Lambda);
	#endif
	GS(1) = GS(1)*sub.weight_concatenated(fs);
	return symas_combine_to_matrix(GS(0), GS(1));
}

template<int mode> int Green_mpi_eval<mode>::dim_r(double job){
	return Nges;
}

template<int mode> int Green_mpi_eval<mode>::dim_c(double job){
	return Nges+1;
}

template<int mode> int Green_mpi_eval<mode>::volume(double job){
	return dim_r(job)*dim_c(job);
}

void pre_sort(matrix<matrix<complex<double> > > &res_mpi, matrix<syma<complex<double> > > &G, matrix<syma<complex<double> > > &S){
	for(int i=0; i<res_mpi.dim_c; ++i){
		auto pair = matrix_split_to_symas(res_mpi(i));
		G(i) = pair.first;
		S(i) = pair.second;
	}
}

template <int mode> class Ex_Precomputation{
	public:
		static const int N_intervals = 10; //should be a divisor of num_freq_pre, otherwise preintegration will be slow. This feature could be improved upon;
		static const double delta_around_dangerous_frequencies = 1e-7;
		static const double delta_around_dangerous_frequencies_at_pm7 = 1e-12;
		static const double delta_around_specific_frequencies = 1e-2;
		Physics &phy;
		Numerics &num;
		double Lambda;
		matrix<double> freq_pre;
		matrix<syma<complex<double> > > Gu;
		matrix<syma<complex<double> > > Su;
		matrix<syma<complex<double> > > Su_ret_integrated;
		matrix<syma<complex<double> > > Su_kel_integrated;
		#if(H_EQUAL_ZERO==0)
			matrix<syma<complex<double> > > Gd;
			matrix<syma<complex<double> > > Sd;
			matrix<syma<complex<double> > > Sd_ret_integrated;
			matrix<syma<complex<double> > > Sd_kel_integrated;
		#else
			matrix<syma<complex<double> > > &Gd = Gu;
			matrix<syma<complex<double> > > &Sd = Su;
			matrix<syma<complex<double> > > &Sd_ret_integrated = Su_ret_integrated;
			matrix<syma<complex<double> > > &Sd_kel_integrated = Su_kel_integrated;
		#endif
		linear_ipol_bin<syma<complex<double> > > iGu; 
		linear_ipol_bin<syma<complex<double> > > iSu; 
		linear_ipol_bin<syma<complex<double> > > iGd; 
		linear_ipol_bin<syma<complex<double> > > iSd; 
		linear_ipol_bin<syma<complex<double> > > iSu_ret_integrated; 
		linear_ipol_bin<syma<complex<double> > > iSu_kel_integrated; 
		linear_ipol_bin<syma<complex<double> > > iSd_ret_integrated; 
		linear_ipol_bin<syma<complex<double> > > iSd_kel_integrated; 
		Ex_Precomputation(Physics &phy_in,
		                  Numerics &num_in);
		void set_freq_pre(Substitution<mode> &sub, matrix<double> spec_freq);
		void set_freq_pre(Substitution<mode> &sub);
		void precompute(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu, linear_ipol_bin<syma<complex<double> > > &iEd);
		void preintegrate(Substitution<mode> &sub);
		void preintegrate_fast(Substitution<mode> &sub);
		void add_katanin(Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &idEu, linear_ipol_bin<syma<complex<double> > > &idEd);
};

template <int mode> Ex_Precomputation<mode>::Ex_Precomputation(Physics &phy_in,
                                                               Numerics &num_in):
                                                               phy(phy_in),
                                                               num(num_in), 
                                                               freq_pre(num.num_freq_pre),
                                                               Gu(num.num_freq_pre),
                                                               Su(num.num_freq_pre),
                                                               Su_ret_integrated(num.num_freq_pre),
                                                               Su_kel_integrated(num.num_freq_pre),
                                                               #if(H_EQUAL_ZERO==0)
                                                               	Gd(num.num_freq_pre),
                                                               	Sd(num.num_freq_pre),
                                                               	Sd_ret_integrated(num.num_freq_pre),
                                                               	Sd_kel_integrated(num.num_freq_pre),
                                                               #endif
                                                               iGu(freq_pre, Gu),
                                                               iSu(freq_pre, Su),
                                                               iGd(freq_pre, Gd),
                                                               iSd(freq_pre, Sd),
                                                               iSu_ret_integrated(freq_pre, Su_ret_integrated),
                                                               iSu_kel_integrated(freq_pre, Su_kel_integrated),
                                                               iSd_ret_integrated(freq_pre, Sd_ret_integrated),
                                                               iSd_kel_integrated(freq_pre, Sd_kel_integrated){

}

template<int mode> void Ex_Precomputation<mode>::set_freq_pre(Substitution<mode> &sub, matrix<double> spec_freq){
	int N_spec_freq = spec_freq.dim_c;
	int num_freq_pre = num.num_freq_pre;
	int num_freq_pre_eff = num_freq_pre - 14 - 3*N_spec_freq; 
	for(int i=0; i<num_freq_pre_eff; ++i){
		freq_pre(i)= -7.+14.*(double)(i+1)/(double)(num_freq_pre_eff+1);
	}
	for(int i=0; i<N_spec_freq; i=i+3){
		freq_pre(num_freq_pre-15-i) = sub.subst_concatenated(spec_freq(i));
		freq_pre(num_freq_pre-15-i-1) = sub.subst_concatenated(spec_freq(i)) + delta_around_specific_frequencies;
		freq_pre(num_freq_pre-15-i-2) = sub.subst_concatenated(spec_freq(i)) - delta_around_specific_frequencies;
	}
	freq_pre(num_freq_pre-14) = -7. + delta_around_dangerous_frequencies_at_pm7;
	freq_pre(num_freq_pre-13) =  7. - delta_around_dangerous_frequencies_at_pm7;
	freq_pre(num_freq_pre-12) = -6. + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre-11) = -6. - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre-10) = -2. + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 9) = -2. - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 8) =  .0 + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 7) =  .0 - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 6) =  2. + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 5) =  2. - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 4) =  6. + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 3) =  6. - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
	freq_pre.sort();
}

template<int mode> void Ex_Precomputation<mode>::set_freq_pre(Substitution<mode> &sub){
	matrix<double> spec_freq(0);
	set_freq_pre(sub,spec_freq);
}

template <int mode> void Ex_Precomputation<mode>::precompute(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu, linear_ipol_bin<syma<complex<double> > > &iEd){
	Lambda=Lambda_in;
	#if(USE_MPI_FOR_PRECOMPUTATION==0)
		#pragma omp parallel for
		for (int i=0; i<num.num_freq_pre; i++) {
			double fs = freq_pre(i);
			matrix<syma<complex<double> > > GS(2);
			syma<complex<double> > E = iEu(fs);
			#if(PARITY_SYMMETRY==0)
				GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
			#else
				GS = green_and_single_Req_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
			#endif
			Gu(i) = GS(0);
			Su(i) = GS(1)*sub.weight_concatenated(fs);
			#if(H_EQUAL_ZERO==0)
				E = iEd(fs);
				#if(PARITY_SYMMETRY==0)
					GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, -phy.h, 1.0, Lambda);
				#else
					GS = green_and_single_Req_ps(sub.resu_concatenated(fs), E, -phy.h, 1.0, Lambda);
				#endif
				Gd(i) = GS(0);
				Sd(i) = GS(1)*sub.weight_concatenated(fs);
			#endif	
		}
	#else
		Green_mpi_eval<mode> green_mpi_up(phy.h, num.Nges, iEu, Lambda, sub);
		matrix<matrix<complex<double> > > Green_list = ex_mpi_computation<complex<double>, double, matrix<complex<double> >, Green_mpi_eval<mode> >(freq_pre, green_mpi_up);
		pre_sort(Green_list,Gu,Su);
		#if(H_EQUAL_ZERO==0)
			Green_mpi_eval<mode> green_mpi_down(-phy.h, num.Nges, iEd, Lambda, sub);
			Green_list = ex_mpi_computation<complex<double>, double, matrix<complex<double> >, Green_mpi_eval<mode> >(freq_pre, green_mpi_down);
			pre_sort(Green_list,Gd,Sd);
		#endif
	#endif
}


template <int mode> void Ex_Precomputation<mode>::add_katanin(Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &idEu, linear_ipol_bin<syma<complex<double> > > &idEd){
	#pragma omp parallel for
	for (int i=0; i<num.num_freq_pre; i++) {
		double fs = freq_pre(i);
		syma<complex<double> > dEu = idEu(fs);
		syma<complex<double> > dEd = idEd(fs);
		syma<complex<double> > tmp = Su(i);
		Su(i) = sub.weight_concatenated(fs)*Gu(i)*dEu*Gu(i); 
		Su(i) += tmp; 
		tmp = Sd(i);
		Sd(i) = sub.weight_concatenated(fs)*Gd(i)*dEd*Gd(i); 
		Sd(i) += tmp; 
	}
}

template <int mode> void Ex_Precomputation<mode>::preintegrate(Substitution<mode> &sub){
	if(num.num_freq_pre%N_intervals==0){
		preintegrate_fast(sub);
	}
	else{
		syma<complex<double> > tmp_ret_u(num.Nges);
		syma<complex<double> > tmp_kel_u(num.Nges);
		tmp_ret_u = (complex<double>) 0.0;
		tmp_kel_u = (complex<double>) 0.0;
		Su_ret_integrated(0) = tmp_ret_u;
		Su_kel_integrated(0) = tmp_kel_u;
		#if(H_EQUAL_ZERO==0)
			syma<complex<double> > tmp_ret_d(num.Nges);
			syma<complex<double> > tmp_kel_d(num.Nges);
			tmp_ret_d = (complex<double>) 0.0;
			tmp_kel_d = (complex<double>) 0.0;
			Sd_ret_integrated(0) = tmp_ret_u;
			Sd_kel_integrated(0) = tmp_kel_u;
		#endif
		for (int i=0; i<num.num_freq_pre-1; i++) {
			double n1 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(i)),phy.mu,phy.T);
			double n2 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(i+1)),phy.mu,phy.T);
			tmp_ret_u += 0.5*(Su(i) + Su(i+1))*(freq_pre(i+1) - freq_pre(i));
			tmp_kel_u += 0.5*(Su(i).imag()*n1 + Su(i+1).imag()*n2)*(freq_pre(i+1) - freq_pre(i));
			Su_ret_integrated(i+1) = tmp_ret_u;
			Su_kel_integrated(i+1) = tmp_kel_u;
			#if(H_EQUAL_ZERO==0)
				tmp_ret_d += 0.5*(Sd(i) + Sd(i+1))*(freq_pre(i+1) - freq_pre(i));
				tmp_kel_d += 0.5*(Sd(i).imag()*n1 + Sd(i+1).imag()*n2)*(freq_pre(i+1) - freq_pre(i));
				Sd_ret_integrated(i+1) = tmp_ret_d;
				Sd_kel_integrated(i+1) = tmp_kel_d;
			#endif
		}
	}
}

template <int mode> void Ex_Precomputation<mode>::preintegrate_fast(Substitution<mode> &sub){
	int nf_per_core = num.num_freq_pre/N_intervals; 
	
	#pragma omp parallel for
	for(int n=0; n<N_intervals; ++n){
		int index = n*nf_per_core;
		syma<complex<double> > tmp_ret_u(num.Nges);
		syma<complex<double> > tmp_kel_u(num.Nges);
		tmp_ret_u = (complex<double>) 0.0;
		tmp_kel_u = (complex<double>) 0.0;
		Su_ret_integrated(index) = tmp_ret_u;
		Su_kel_integrated(index) = tmp_kel_u;
		#if(H_EQUAL_ZERO==0)
			syma<complex<double> > tmp_ret_d(num.Nges);
			syma<complex<double> > tmp_kel_d(num.Nges);
			tmp_ret_d = (complex<double>) 0.0;
			tmp_kel_d = (complex<double>) 0.0;
			Sd_ret_integrated(index) = tmp_ret_d;
			Sd_kel_integrated(index) = tmp_kel_d;
		#endif
		for (int i=0; i<nf_per_core-1; ++i, ++index) {
			double n1 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(index)),phy.mu,phy.T);
			double n2 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(index+1)),phy.mu,phy.T);
			tmp_ret_u += 0.5*(Su(index) + Su(index+1))*(freq_pre(index+1) - freq_pre(index));
			tmp_kel_u += 0.5*(Su(index).imag()*n1 + Su(index+1).imag()*n2)*(freq_pre(index+1) - freq_pre(index));
			Su_ret_integrated(index+1) = tmp_ret_u;
			Su_kel_integrated(index+1) = tmp_kel_u;
			#if(H_EQUAL_ZERO==0)
				tmp_ret_d += 0.5*(Sd(index) + Sd(index+1))*(freq_pre(index+1) - freq_pre(index));
				tmp_kel_d += 0.5*(Sd(index).imag()*n1 + Sd(index+1).imag()*n2)*(freq_pre(index+1) - freq_pre(index));
				Sd_ret_integrated(index+1) = tmp_ret_d;
				Sd_kel_integrated(index+1) = tmp_kel_d;
			#endif
		}
	}
	for(int n=1; n<N_intervals; ++n){
		int ni = n*nf_per_core;
		#pragma omp parallel for
		for(int k=0; k<nf_per_core; ++k){
			double n1 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(ni-1)),phy.mu,phy.T);
			double n2 = 1. - 2.*fermi(sub.resu_concatenated(freq_pre(ni)),phy.mu,phy.T);
			Su_ret_integrated(ni + k) += Su_ret_integrated(ni-1)
			                          + 0.5*(Su(ni-1) + Su(ni))*(freq_pre(ni) - freq_pre(ni-1)); 
			Su_kel_integrated(ni + k) += Su_kel_integrated(ni-1)
			                          + 0.5*(Su(ni-1).imag()*n1 + Su(ni).imag()*n2)*(freq_pre(ni) - freq_pre(ni-1)); 
			#if(H_EQUAL_ZERO==0)
				Sd_ret_integrated(ni + k) += Sd_ret_integrated(ni-1)
				                          + 0.5*(Sd(ni-1) + Sd(ni))*(freq_pre(ni) - freq_pre(ni-1)); 
				Sd_kel_integrated(ni + k) += Sd_kel_integrated(ni-1)
				                          + 0.5*(Sd(ni-1).imag()*n1 + Sd(ni).imag()*n2)*(freq_pre(ni) - freq_pre(ni-1)); 
			#endif
		}
	}
}

#endif
