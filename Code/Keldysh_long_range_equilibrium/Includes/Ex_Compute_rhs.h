#ifndef EX_COMPUTE_RHS_06112018
#define EX_COMPUTE_RHS_06112018

#include <iostream> 
#include <stdio.h>
#include <string.h>
#include <chrono>

#ifndef ONLY_ONSITE_INTERACTION 
	#define ONLY_ONSITE_INTERACTION 0
#endif

#include "matrix.h" 
#include "Numerics.h"
#include "Ex_freq_str.h"
#include "Ex_Diagnostics.h"
#include "Ex_Compute_bubble.h"
#include "Ex_Compute_self.h"
#include "Ex_multiplication.h"

using namespace std;

template<int mode> class Ex_Flow_static_background{
	public:
		Numerics &num;
		Barevertex &barevertex;
		Ex_Vertex<mode> &gamma;
		matrix<matrix<double> > Puu;
		matrix<matrix<double> > Pdd;
		matrix<matrix<double> > Pud;
		matrix<matrix<double> > Xud;
		matrix<matrix<double> > Dud;
		matrix<matrix<double> > Duu;
		matrix<matrix<double> > Ddd;
		Ex_Flow_static_background(Numerics &num_in,
		                          Barevertex &barevertex_in,
		                          Ex_Vertex<mode> &gamma_in);
		void assemble();
};




template<int mode> Ex_Flow_static_background<mode>::Ex_Flow_static_background(Numerics &num_in,
                                                                              Barevertex &barevertex_in,
                                                                              Ex_Vertex<mode> &gamma_in):
                                                                              num(num_in),
                                                                              barevertex(barevertex_in),
                                                                              gamma(gamma_in){
}

template<int mode> void Ex_Flow_static_background<mode>::assemble(){
	int L=num.L;
	int Lges=2*L+1;
	resize_str(Puu,L,num.N);
	resize_str(Pdd,L,num.N);
	resize_str(Pud,L,num.N);
	resize_str(Xud,L,num.N);
	resize_str(Duu,L,num.N);
	resize_str(Ddd,L,num.N);
	resize_str(Dud,L,num.N);
	#pragma omp parallel for collapse(2)
	for(int lc=0, l=-L; lc<Lges; ++lc, ++l){
		for(int kc=0, k=-L; kc<Lges; ++kc, ++k){
			int Ngesl = num.Nges-abs(l);
			int Ngesk = num.Nges-abs(k);
			for(int jc=0, j=max(0,-l); jc<Ngesl; ++jc, ++j){
				for(int ic=0, i=max(0,-k); ic<Ngesk; ++ic, ++i){
					Puu(lc,kc)(jc,ic) = 0.5*barevertex(j,1,j+l,1,i,1,i+k,1);

					Pdd(lc,kc)(jc,ic) = 0.5*barevertex(j,0,j+l,0,i,0,i+k,0);

					Pud(lc,kc)(jc,ic) = 0.5*barevertex(j,1,j+l,0,i,1,i+k,0);

					Xud(lc,kc)(jc,ic) = 0.5*barevertex(j,1,i+k,0,i,1,j+l,0);
					
					Duu(lc,kc)(jc,ic) = 0.5*barevertex(j,1,i+k,1,j+l,1,i,1);
					
					Ddd(lc,kc)(jc,ic) = 0.5*barevertex(j,0,i+k,0,j+l,0,i,0);
					
					Dud(lc,kc)(jc,ic) = 0.5*barevertex(j,1,i+k,0,j+l,1,i,0);
					#if(RPA_MODE==0)
						Puu(lc,kc)(jc,ic) -= gamma.aDuu.get_stat(i+k-j,j+l-i,j,i);
						Puu(lc,kc)(jc,ic) += gamma.aDuu.get_stat(i-j,j+l-i-k,j,i+k);
						
						Pdd(lc,kc)(jc,ic) -= gamma.aDdd.get_stat(i+k-j,j+l-i,j,i);
						Pdd(lc,kc)(jc,ic) += gamma.aDdd.get_stat(i-j,j+l-i-k,j,i+k);
						
						Pud(lc,kc)(jc,ic) += gamma.aXud.get_stat(i+k-j,j+l-i,j,i);
						Pud(lc,kc)(jc,ic) += gamma.aDud.get_stat(i-j,j+l-i-k,j,i+k);
						
						Xud(lc,kc)(jc,ic) += gamma.aPud.get_stat(i+k-j,j+l-i,j,i);
						Xud(lc,kc)(jc,ic) += gamma.aDud.get_stat(i-j,i+k-j-l,j,j+l);
						
						Duu(lc,kc)(jc,ic) += gamma.aPuu.get_stat(i+k-j,i-j-l,j,j+l);
						Duu(lc,kc)(jc,ic) -= gamma.aDuu.get_stat(i-j,i+k-j-l,j,j+l);
						
						Ddd(lc,kc)(jc,ic) += gamma.aPdd.get_stat(i+k-j,i-j-l,j,j+l);
						Ddd(lc,kc)(jc,ic) -= gamma.aDdd.get_stat(i-j,i+k-j-l,j,j+l);
						
						Dud(lc,kc)(jc,ic) += gamma.aPud.get_stat(i+k-j,i-j-l,j,j+l);
						Dud(lc,kc)(jc,ic) += gamma.aXud.get_stat(i-j,i+k-j-l,j,j+l);
					#endif
				}
			}
		}
	}
}



template<int mode> class Ex_Compute_rhs{
	public:
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		double Lambda;
		Ex_Vertex<mode> &gamma;
		Ex_Vertex<mode> &dgamma;
		Ex_Flow_static_background<mode> &background;
		Ex_Diagnostics &diagnostics;
		Ex_Compute_bubble<mode> bubble_computer;
		Ex_Compute_rhs(Physics &phy_in,
		               Numerics &num_in,
		               Ex_Precomputation<mode> &pre_in,
		               Substitution<mode> &sub_in,
		               double measure_flow_in,
		               double Lambda_in,
		               Ex_Vertex<mode> &gamma_in,
		               Ex_Vertex<mode> &dgamma_in,
		               Ex_Flow_static_background<mode> &background_in,
		               Ex_Diagnostics &diagnostics_in);
		void p_channel(double accuracy, matrix<double> &additional_stops);
		void xd_channel(double accuracy, matrix<double> &additional_stops);
		void self_energy(double accuracy, matrix<double> &additional_stops);
};

template<int mode> Ex_Compute_rhs<mode>::Ex_Compute_rhs(Physics &phy_in,
                                                        Numerics &num_in,
                                                        Ex_Precomputation<mode> &pre_in,
                                                        Substitution<mode> &sub_in,
                                                        double measure_flow_in,
                                                        double Lambda_in,
                                                        Ex_Vertex<mode> &gamma_in,
                                                        Ex_Vertex<mode> &dgamma_in,
                                                        Ex_Flow_static_background<mode> &background_in,
                                                        Ex_Diagnostics &diagnostics_in):
                                                        phy(phy_in),
                                                        num(num_in),
                                                        pre(pre_in),
                                                        sub(sub_in),
                                                        measure_flow(measure_flow_in),
                                                        Lambda(Lambda_in),
                                                        gamma(gamma_in),
                                                        dgamma(dgamma_in),
                                                        background(background_in),
                                                        diagnostics(diagnostics_in),
                                                        bubble_computer(phy,num,pre,sub,measure_flow,Lambda,diagnostics){
}


template<int mode> void Ex_Compute_rhs<mode>::p_channel(double accuracy, matrix<double> &additional_stops){ 
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	auto start_time= std::chrono::system_clock::now();
	#if(H_EQUAL_ZERO==1)
		matrix<matrix<matrix<complex<double> > > >  bubble_e_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_e_stat;
		Ex_freq_str Bubble_e(num.L, num.N, num.Lp_structure, num.wbP, bubble_e_dyn, bubble_e_stat);
		Bubble_e.resize();
		bubble_computer.compute_ePst_bubble(1,1,additional_stops,Bubble_e,accuracy);
	#endif
	#if(WITHOUT_PUU_PDD_DUD==0)
		{//Puu:
			matrix<matrix<matrix<complex<double> > > > A_dyn = add_static_term(gamma.aPuu.dynamic_str,background.Puu); 
			matrix<matrix<double> > A_stat = gamma.aPuu.static_str + background.Puu;
			#if(H_EQUAL_ZERO==0)
				matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
				matrix<matrix<double> >  bubble_stat;
				Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
				Bubble.resize();
				bubble_computer.compute_Puu_bubble(additional_stops,Bubble,accuracy);
				auto pair = total_mult(A_dyn, A_stat, Bubble.dynamic_str, Bubble.static_str, A_dyn, A_stat, num.Lp_structure, num.Lp_bounds, num.pos_NfbP_2mu);
			#else
				auto pair = total_mult(A_dyn, A_stat, Bubble_e.dynamic_str, Bubble_e.static_str, A_dyn, A_stat, num.Lp_structure, num.Lp_bounds, num.pos_NfbP_2mu);
			#endif
			dgamma.aPuu.dynamic_str = pair.first;
			dgamma.aPuu.static_str = pair.second;
			#if(RPA_BUBBLE_ONLY==1)
				dgamma.aPuu.dynamic_str = Bubble.dynamic_str; 
				dgamma.aPuu.static_str  = Bubble.static_str;
			#endif
		}
		{//Pdd:
			#if(H_EQUAL_ZERO==0)
				matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
				matrix<matrix<double> >  bubble_stat;
				Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
				Bubble.resize();
				bubble_computer.compute_Pdd_bubble(additional_stops,Bubble,accuracy);

				matrix<matrix<matrix<complex<double> > > > A_dyn = add_static_term(gamma.aPdd.dynamic_str,background.Pdd); 
				matrix<matrix<double> > A_stat = gamma.aPdd.static_str + background.Pdd;
				auto pair = total_mult(A_dyn, A_stat, Bubble.dynamic_str, Bubble.static_str, A_dyn, A_stat, num.Lp_structure, num.Lp_bounds, num.pos_NfbP_2mu);
				dgamma.aPdd.dynamic_str = pair.first;
				dgamma.aPdd.static_str = pair.second;
				#if(RPA_BUBBLE_ONLY==1)
					dgamma.aPdd.dynamic_str = Bubble.dynamic_str; 
					dgamma.aPdd.static_str  = Bubble.static_str;
				#endif
			#else
				dgamma.aPdd.dynamic_str = dgamma.aPuu.dynamic_str;
				dgamma.aPdd.static_str = dgamma.aPuu.static_str;
			#endif
		}
	#endif
	{//Pud:
		int rank;
		int root=0;
		int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto start_time= std::chrono::system_clock::now();
		matrix<matrix<matrix<complex<double> > > > A_dyn = add_static_term(gamma.aPud.dynamic_str,background.Pud); 
		matrix<matrix<double> > A_stat = gamma.aPud.static_str + background.Pud;
		auto end_time= std::chrono::system_clock::now();
		std::chrono::duration<double> time = end_time - start_time;
		if(rank==root) cout<<"Internal Pud time: assemble static factors="<<time.count()<<endl;;

		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		#if(H_EQUAL_ZERO==0)
			bubble_computer.compute_Pud_bubble(additional_stops,Bubble,accuracy);
		#else
			start_time= std::chrono::system_clock::now();
			generate_Pud_bubble(Bubble,Bubble_e,Bubble_e);	
			end_time= std::chrono::system_clock::now();
			time = end_time - start_time;
			if(rank==root) cout<<"Internal Pud time: generate_Pud_bubble="<<time.count()<<endl;;
		#endif
		start_time= std::chrono::system_clock::now();
		auto pair = total_mult(A_dyn, A_stat, Bubble.dynamic_str, Bubble.static_str, A_dyn, A_stat, num.Lp_structure, num.Lp_bounds, num.pos_NfbP_2mu);
		end_time= std::chrono::system_clock::now();
		time = end_time - start_time;
		if(rank==root) cout<<"Internal Pud time: total_mult="<<time.count()<<endl;;

		start_time= std::chrono::system_clock::now();
		dgamma.aPud.dynamic_str = pair.first;
		dgamma.aPud.static_str = pair.second;
		end_time= std::chrono::system_clock::now();
		time = end_time - start_time;
		if(rank==root) cout<<"Internal Pud time: write results="<<time.count()<<endl;;
		#if(RPA_BUBBLE_ONLY==1)
			dgamma.aPud.dynamic_str = Bubble.dynamic_str; 
			dgamma.aPud.static_str  = Bubble.static_str;
		#endif
	}
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.P_rhs_time.push_back(time.count());
	if(rank==root){
		cout<<"Time for rhs-P-channel="<<time.count()<<endl;
	}
}


template<int mode> void Ex_Compute_rhs<mode>::xd_channel(double accuracy, matrix<double> &additional_stops){ 
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	auto start_time= std::chrono::system_clock::now();
	#if(H_EQUAL_ZERO==1)
		matrix<matrix<matrix<complex<double> > > >  bubble_e_dyn(num.NfbX);
		matrix<matrix<double> >  bubble_e_stat;
		Ex_freq_str Bubble_e(num.L, num.N, num.Lx_structure, num.wbX, bubble_e_dyn, bubble_e_stat);
		Bubble_e.resize();
		bubble_computer.compute_eXst_bubble(1,1,additional_stops,Bubble_e,accuracy);
	#endif
	{//Xud:
		matrix<matrix<matrix<complex<double> > > > A_dyn = add_static_term(gamma.aXud.dynamic_str,background.Xud); 
		matrix<matrix<double> > A_stat = gamma.aXud.static_str + background.Xud;

		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbX);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn, bubble_stat);
		#if(H_EQUAL_ZERO==0)
			bubble_computer.compute_Xud_bubble(additional_stops,Bubble,accuracy);
		#else
			generate_Xud_bubble(Bubble,Bubble_e,Bubble_e);
		#endif

		auto pair = total_mult(A_dyn, A_stat, Bubble.dynamic_str, Bubble.static_str, A_dyn, A_stat, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
		dgamma.aXud.dynamic_str = pair.first;
		dgamma.aXud.static_str = pair.second;
		#if(RPA_BUBBLE_ONLY==1)
			dgamma.aXud.dynamic_str = Bubble.dynamic_str; 
			dgamma.aXud.static_str  = Bubble.static_str;
		#endif
	}
	{// Duu, Ddd, and Dud:
		matrix<matrix<matrix<complex<double> > > > A_dyn_uu = add_static_term(gamma.aDuu.dynamic_str,background.Duu); 
		matrix<matrix<double> > A_stat_uu = gamma.aDuu.static_str + background.Duu;
		
		#if(H_EQUAL_ZERO==0)
			matrix<matrix<matrix<complex<double> > > > A_dyn_dd = add_static_term(gamma.aDdd.dynamic_str,background.Ddd); 
			matrix<matrix<double> > A_stat_dd = gamma.aDdd.static_str + background.Ddd;
		#else
			matrix<matrix<matrix<complex<double> > > > &A_dyn_dd = A_dyn_uu; 
			matrix<matrix<double> > &A_stat_dd = A_stat_uu;
		#endif
		
		matrix<matrix<matrix<complex<double> > > > A_dyn_ud = add_static_term(gamma.aDud.dynamic_str,background.Dud); 
		matrix<matrix<double> > A_stat_ud = gamma.aDud.static_str + background.Dud;
		
		matrix<matrix<matrix<complex<double> > > > A_dyn_ud_transp = transpose_lr_str(A_dyn_ud); 
		matrix<matrix<double> > A_stat_ud_transp = transpose_lr_str(A_stat_ud);

		{ //Bubble_up_up contribution:
			
			matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbX);
			matrix<matrix<double> >  bubble_stat;
			Ex_freq_str Bubble(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn, bubble_stat);
			#if(H_EQUAL_ZERO==0)
				bubble_computer.compute_Duu_bubble(additional_stops,Bubble,accuracy);
			#else
				generate_Dss_bubble(Bubble,Bubble_e);
			#endif

			auto pair = total_mult(A_dyn_uu, A_stat_uu, Bubble.dynamic_str, Bubble.static_str, A_dyn_uu, A_stat_uu, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
			dgamma.aDuu.dynamic_str = -pair.first;
			dgamma.aDuu.static_str =  -pair.second;
			     
			#if(WITHOUT_PUU_PDD_DUD==0)
				     pair = total_mult(A_dyn_uu, A_stat_uu, Bubble.dynamic_str, Bubble.static_str, A_dyn_ud, A_stat_ud, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
				dgamma.aDud.dynamic_str = -pair.first;
				dgamma.aDud.static_str =  -pair.second;
			#endif

			#if(H_EQUAL_ZERO==0)
				     pair = total_mult(A_dyn_ud_transp, A_stat_ud_transp, Bubble.dynamic_str, Bubble.static_str, A_dyn_ud, A_stat_ud, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
				dgamma.aDdd.dynamic_str = -pair.first;
				dgamma.aDdd.static_str =  -pair.second;
			#else
				     pair = total_mult(A_dyn_ud, A_stat_ud, Bubble.dynamic_str, Bubble.static_str, A_dyn_ud_transp, A_stat_ud_transp, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
				dgamma.aDuu.dynamic_str -= pair.first;
				dgamma.aDuu.static_str  -= pair.second;
				     
				#if(WITHOUT_PUU_PDD_DUD==0)
					     pair = total_mult(A_dyn_ud, A_stat_ud, Bubble.dynamic_str, Bubble.static_str, A_dyn_dd, A_stat_dd, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
					dgamma.aDud.dynamic_str -= pair.first;
					dgamma.aDud.static_str  -= pair.second;
				#endif

				dgamma.aDdd.dynamic_str = dgamma.aDuu.dynamic_str;
				dgamma.aDdd.static_str  = dgamma.aDuu.static_str;
			#endif
		}
		{ //Bubble_down_down contribution:
			#if(H_EQUAL_ZERO==0)
				matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbX);
				matrix<matrix<double> >  bubble_stat;
				Ex_freq_str Bubble(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn, bubble_stat);
				#if(H_EQUAL_ZERO==0)
					bubble_computer.compute_Ddd_bubble(additional_stops,Bubble,accuracy);
				#else
					generate_Dss_bubble(Bubble,Bubble_e);
				#endif

				auto pair = total_mult(A_dyn_ud, A_stat_ud, Bubble.dynamic_str, Bubble.static_str, A_dyn_ud_transp, A_stat_ud_transp, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
				dgamma.aDuu.dynamic_str -= pair.first;
				dgamma.aDuu.static_str  -= pair.second;
				     
				#if(WITHOUT_PUU_PDD_DUD==0)
					     pair = total_mult(A_dyn_ud, A_stat_ud, Bubble.dynamic_str, Bubble.static_str, A_dyn_dd, A_stat_dd, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
					dgamma.aDud.dynamic_str -= pair.first;
					dgamma.aDud.static_str  -= pair.second;
				#endif

				     pair = total_mult(A_dyn_dd, A_stat_dd, Bubble.dynamic_str, Bubble.static_str, A_dyn_dd, A_stat_dd, num.Lx_structure, num.Lx_bounds, num.pos_NfbX_0);
				dgamma.aDdd.dynamic_str -= pair.first;
				dgamma.aDdd.static_str  -= pair.second;
			#endif
		}
	}
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.XD_rhs_time.push_back(time.count());
	if(rank==root){
		cout<<"Time for rhs-XD-channel="<<time.count()<<endl;
	}
	
	#if(ONLY_D_CHANNEL==1)
		dgamma.aPuu.initialize(0.0);
		dgamma.aPdd.initialize(0.0);
		dgamma.aPud.initialize(0.0);
		dgamma.aXud.initialize(0.0);
	#endif
	

}

template<int mode> void Ex_Compute_rhs<mode>::self_energy(double accuracy, matrix<double> &additional_stops){ 
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	auto start_time= std::chrono::system_clock::now();
	matrix<matrix<syma<complex<double> > > > self_energy = ex_compute_self_complete(phy, num, pre, sub, Lambda, measure_flow, gamma, background.barevertex, accuracy, additional_stops);
	dgamma.ERetu = self_energy(0);
	dgamma.ERetd = self_energy(1);
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.Self_rhs_time.push_back(time.count());
	if(rank==root){
		cout<<"Time for rhs-Selfenergy="<<time.count()<<endl;
	}
}


#endif
