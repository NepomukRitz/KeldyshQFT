#ifndef EX_COMPUTE_BUBBLE_30102018
#define EX_COMPUTE_BUBBLE_30102018

#include <iostream> 
#include <stdio.h>
#include <string.h>
#include <typeinfo>
#include <chrono>

#ifndef SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 
	#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 0
#endif

#include "matrix.h" 
#include "Numerics.h"
#include "Ex_freq_str.h"
#include "Ex_mpi.h"
#include "Ex_bubble.h"
#include "Ex_Diagnostics.h"

using namespace std;

int ex_determine_number_of_offdiag_matrix_components(matrix<int> &L_structure, bool without_cross){
	int N_total=0;
	int Nf = L_structure.dim_c;
	for(int i=0; i<Nf; ++i){
		int L_dyn = L_structure(i);
		int L_dyn_ges = 2*L_dyn+1;
		N_total += (L_dyn_ges*(L_dyn_ges-1))/2;
		if(without_cross){
			N_total -= 2*L_structure(i);
		}
	}
	return N_total;
}

int ex_determine_number_of_diag_matrix_components(matrix<int> &L_structure, bool without_cross){
	int N_total=0;
	int Nf = L_structure.dim_c;
	for(int i=0; i<Nf; ++i){
		N_total += 2*L_structure(i)+1;
	}
	if(without_cross){
		N_total-=L_structure.dim_c;
	}
	return N_total;
}

matrix<matrix<double> > ex_job_list_for_dyn_offdiag_bubble(matrix<int> &L_structure, matrix<double> &wb, bool without_cross){
	int N_total = ex_determine_number_of_offdiag_matrix_components(L_structure, without_cross);
	#if(SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE==0)
		N_total+= ex_determine_number_of_diag_matrix_components(L_structure, without_cross);
	#endif
	int Nf = L_structure.dim_c;
	matrix<matrix<double> > matrix_list(N_total);
	for(int i=0, z=0; i<Nf; ++i){
		int L_dyn = L_structure(i);
		int L_dyn_ges = 2*L_dyn+1;
		for(int l=-L_dyn; l<=L_dyn; ++l){
			#if(SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE==0)
				for(int k=-L_dyn; k<=l; ++k){
			#else
				for(int k=-L_dyn; k<l; ++k){
			#endif
					if( (!without_cross) || (l!=0 && k!=0) ){ 
						matrix<double> tmp(4);
						tmp(0) = wb(i);
						tmp(1) = (double) l;
						tmp(2) = (double) k;
						tmp(3) = (double) i;
						matrix_list(z) = tmp; 
						++z;
					}
				}
		}
	}
	return matrix_list;
}

matrix<matrix<double> > ex_job_list_for_dyn_diag_bubble(matrix<int> &L_structure, matrix<double> &wb, bool without_cross){
	int N_total = ex_determine_number_of_diag_matrix_components(L_structure, without_cross);
	int Nf = L_structure.dim_c;
	matrix<matrix<double> > matrix_list(N_total);
	for(int i=0, z=0; i<Nf; ++i){
		int L_dyn = L_structure(i);
		int L_dyn_ges = 2*L_dyn+1;
		for(int l=-L_dyn; l<=L_dyn; ++l){
			if( (!without_cross) || (l!=0) ){
				matrix<double> tmp(4);
				tmp(0) = wb(i);
				tmp(1) = (double) l;
				tmp(2) = (double) l;
				tmp(3) = (double) i;
				matrix_list(z) = tmp; 
				++z;
			}
			
		}
	}
	return matrix_list;
}

matrix<matrix<double> > ex_job_list_for_static_bubble(int L, double static_freq){
	int Lges = 2*L+1;
	int N_total = (Lges*(Lges+1))/2;
	matrix<matrix<double> > matrix_list(N_total);
	for(int l=-L, z=0; l<=L; ++l){
		for(int k=-L; k<=l; ++k, ++z){
			matrix<double> tmp(3);
			tmp(0) = static_freq;
			tmp(1) = l;
			tmp(2) = k;
			matrix_list(z) = tmp; 
		}
	}
	return matrix_list;
}


template <typename T1, typename T2> void ex_compute_dyn_bubble(Ex_freq_str &bubble, T1 &Bubble_dyn_offdiag, T2 &Bubble_dyn_diag, bool without_cross){

	//Compute (off)diagonal part:
	matrix<matrix<double> > job_list=ex_job_list_for_dyn_offdiag_bubble(bubble.L_structure, bubble.wb, without_cross);
	matrix<matrix<complex<double> > > bubble_list_offdiag = ex_mpi_computation<complex<double>, matrix<double>, matrix<complex<double> >, T1 >(job_list, Bubble_dyn_offdiag); 	
	for(int z=0; z<job_list.dim_c; ++z){
		int l = (int) job_list(z)(1);
		int k = (int) job_list(z)(2);
		int i = (int) job_list(z)(3);
		int Li = bubble.L_structure(i);
		bubble.dynamic_str(i)(Li+l, Li+k) = bubble_list_offdiag(z); 
		bubble.dynamic_str(i)(Li+k, Li+l) = bubble_list_offdiag(z).transp(); 
	}
	//Set cross to zero:
	if(without_cross){
		for(int i=0; i<bubble.N_freq; ++i){
			int Li = bubble.L_structure(i);
			int Liges = 2*Li+1;
			for(int l=0; l<Liges; ++l){
				bubble.dynamic_str(i)(l,Li)= (complex<double>) 0.0;
				bubble.dynamic_str(i)(Li,l)= (complex<double>) 0.0;
			} 
		}
	}
	#if(SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE==1)
		job_list = ex_job_list_for_dyn_diag_bubble(bubble.L_structure, bubble.wb, without_cross);
		matrix<syma<complex<double> > > bubble_list_diag = ex_mpi_computation<complex<double>, matrix<double>, syma<complex<double> >, T2>(job_list, Bubble_dyn_diag);
		for(int z=0; z<job_list.dim_c; ++z){
			int l = (int) job_list(z)(1);
			int i = (int) job_list(z)(3);
			int Li = bubble.L_structure(i);
			bubble.dynamic_str(i)(Li+l,Li+l) = bubble_list_diag(z);
		}
	#endif
}


template <typename T> void ex_compute_static_bubble(double static_freq, Ex_freq_str &bubble, T &Bubble_stat){

	//Compute static part:
	int L = bubble.L;
	matrix<matrix<double> > job_list=ex_job_list_for_static_bubble(L, static_freq);
	matrix<matrix<double> > bubble_list = ex_mpi_computation<double, matrix<double>, matrix<double>, T >(job_list, Bubble_stat); 	
	for(int z=0; z<job_list.dim_c; ++z){
		int l = (int) job_list(z)(1);
		int k = (int) job_list(z)(2);
		bubble.static_str(l+L, k+L) = bubble_list(z); 
		bubble.static_str(k+L, l+L) = bubble_list(z).transp(); 
	}
}


void generate_Pud_bubble(Ex_freq_str &Pud_bubble, Ex_freq_str &bubble_epud, Ex_freq_str &bubble_epdu){
	Pud_bubble.dynamic_str = bubble_epud.dynamic_str + mirror_lr_str(bubble_epdu.dynamic_str); 
	Pud_bubble.static_str = bubble_epud.static_str + mirror_lr_str(bubble_epdu.static_str); 
}

void generate_Xud_bubble(Ex_freq_str &Xud_bubble, Ex_freq_str &bubble_exud, Ex_freq_str &bubble_exdu){
	Xud_bubble.dynamic_str = bubble_exud.dynamic_str + conjugate_mirror_freq_lr_str(bubble_exdu.dynamic_str); 
	Xud_bubble.static_str = bubble_exud.static_str + mirror_lr_str(bubble_exdu.static_str); 
}

void generate_Dss_bubble(Ex_freq_str &Dss_bubble, Ex_freq_str &bubble_exss){
	Dss_bubble.dynamic_str = mirror_freq_lr_str(bubble_exss.dynamic_str) + conjugate_mirror_lr_str(bubble_exss.dynamic_str); 
	Dss_bubble.static_str = bubble_exss.static_str + mirror_lr_str(bubble_exss.static_str); 
}

template<int mode> class Ex_Compute_bubble{
	public:
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		double Lambda;
		Ex_Diagnostics &diagnostics;
		Ex_Compute_bubble(Physics &phy_in,
		                    Numerics &num_in,
		                    Ex_Precomputation<mode> &pre_in,
		                    Substitution<mode> &sub_in,
		                    double measure_flow_in,
		                    double Lambda_in,
		                    Ex_Diagnostics &diagnostics_in);
		void compute_ePst_bubble(bool spin1, bool spin2, matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_eXst_bubble(bool spin1, bool spin2, matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Puu_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Pdd_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Pud_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Xud_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Duu_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
		void compute_Ddd_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy);
};

template<int mode> Ex_Compute_bubble<mode>::Ex_Compute_bubble(Physics &phy_in,
                                                                  Numerics &num_in,
                                                                  Ex_Precomputation<mode> &pre_in,
                                                                  Substitution<mode> &sub_in,
                                                                  double measure_flow_in,
                                                                  double Lambda_in,
                                                                  Ex_Diagnostics &diagnostics_in):
                                                                  phy(phy_in),
                                                                  num(num_in),
                                                                  pre(pre_in),
                                                                  sub(sub_in),
                                                                  measure_flow(measure_flow_in),
                                                                  Lambda(Lambda_in),
                                                                  diagnostics(diagnostics_in){
}

template<int mode> void Ex_Compute_bubble<mode>::compute_ePst_bubble(bool spin1, bool spin2, matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	auto start_time= std::chrono::system_clock::now();

	P_Stops<mode> stops(phy, sub, Lambda);
	Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	Integrator_bubble<mode,Integrand_P_bubble_feedback<mode> > Bubble_feedback(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	#if(ONLY_ONSITE_INTERACTION==0)
		bool without_cross=0;
	#else
		bool without_cross=(spin1==spin2);
	#endif
	ex_compute_dyn_bubble(Bubble, Bubble_dyn_off, Bubble_dyn_on, without_cross);
	ex_compute_static_bubble(2.*phy.mu,Bubble, Bubble_feedback);
	
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.feed_P_bubble(spin1, spin2, time.count());
	if(rank==root){
		cout<<"Time for eP-bubble, spin1="<<spin1<<", spin2="<<spin2<<", time="<<time.count()<<endl;
	}
}

template<int mode> void Ex_Compute_bubble<mode>::compute_eXst_bubble(bool spin1, bool spin2, matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	auto start_time= std::chrono::system_clock::now();

	X_Stops<mode> stops(phy, sub, Lambda);
	Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	Integrator_bubble<mode,Integrand_X_bubble_feedback<mode> > Bubble_feedback(spin1,spin2,phy,num,pre,sub,measure_flow,accuracy, additional_stops, stops);
	bool without_cross=0;
	ex_compute_dyn_bubble(Bubble, Bubble_dyn_off, Bubble_dyn_on, without_cross);
	ex_compute_static_bubble(0.0,Bubble, Bubble_feedback);
	
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.feed_X_bubble(spin1, spin2, time.count());
	if(rank==root){
		cout<<"Time for eX-bubble, spin1="<<spin1<<", spin2="<<spin2<<", time="<<time.count()<<endl;
	}
}
		
template<int mode> void Ex_Compute_bubble<mode>::compute_Puu_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	compute_ePst_bubble(1,1,additional_stops,Bubble, accuracy);
}

template<int mode> void Ex_Compute_bubble<mode>::compute_Pdd_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	compute_ePst_bubble(0,0,additional_stops,Bubble, accuracy);
}

template<int mode> void Ex_Compute_bubble<mode>::compute_Pud_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	matrix<matrix<matrix<complex<double> > > >  bubble_dyn_ud(num.NfbP);
	matrix<matrix<double> >  bubble_stat_ud;
	Ex_freq_str Bubble_ud(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_ud, bubble_stat_ud);
	Bubble_ud.resize();
	compute_ePst_bubble(1,0,additional_stops,Bubble_ud, accuracy);
	
	#if(H_EQUAL_ZERO==0)
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn_du(num.NfbP);
		matrix<matrix<double> >  bubble_stat_du;
		Ex_freq_str Bubble_du(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_du, bubble_stat_du);
		Bubble_du.resize();
		compute_ePst_bubble(0,1,additional_stops,Bubble_du, accuracy);
		generate_Pud_bubble(Bubble,Bubble_ud,Bubble_du);
	#else
		generate_Pud_bubble(Bubble,Bubble_ud,Bubble_ud);
	#endif	
	
}

template<int mode> void Ex_Compute_bubble<mode>::compute_Xud_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	matrix<matrix<matrix<complex<double> > > >  bubble_dyn_ud(num.NfbX);
	matrix<matrix<double> >  bubble_stat_ud;
	Ex_freq_str Bubble_ud(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_ud, bubble_stat_ud);
	Bubble_ud.resize();
	compute_eXst_bubble(1,0,additional_stops,Bubble_ud, accuracy);
	
	#if(H_EQUAL_ZERO==0)
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn_du(num.NfbX);
		matrix<matrix<double> >  bubble_stat_du;
		Ex_freq_str Bubble_du(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_du, bubble_stat_du);
		Bubble_du.resize();
		compute_eXst_bubble(0,1,additional_stops,Bubble_du, accuracy);
		generate_Xud_bubble(Bubble,Bubble_ud,Bubble_du);
	#else
		generate_Xud_bubble(Bubble,Bubble_ud,Bubble_ud);
	#endif
	
}

template<int mode> void Ex_Compute_bubble<mode>::compute_Duu_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	matrix<matrix<matrix<complex<double> > > >  bubble_e_dyn(num.NfbX);
	matrix<matrix<double> >  bubble_e_stat;
	Ex_freq_str Bubble_e(num.L, num.N, num.Lx_structure, num.wbX, bubble_e_dyn, bubble_e_stat);
	Bubble_e.resize();
	compute_eXst_bubble(1,1,additional_stops,Bubble_e, accuracy);
	
	generate_Dss_bubble(Bubble,Bubble_e);
}

template<int mode> void Ex_Compute_bubble<mode>::compute_Ddd_bubble(matrix<double> &additional_stops, Ex_freq_str &Bubble, double accuracy){
	matrix<matrix<matrix<complex<double> > > >  bubble_e_dyn(num.NfbX);
	matrix<matrix<double> >  bubble_e_stat;
	Ex_freq_str Bubble_e(num.L, num.N, num.Lx_structure, num.wbX, bubble_e_dyn, bubble_e_stat);
	Bubble_e.resize();
	compute_eXst_bubble(0,0,additional_stops,Bubble_e, accuracy);
	
	generate_Dss_bubble(Bubble,Bubble_e);
}

#endif
