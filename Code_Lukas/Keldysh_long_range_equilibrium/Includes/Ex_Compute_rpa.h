#ifndef EX_COMPUTE_RPA_10022019
#define EX_COMPUTE_RPA_10022019

#define COMPUTE_RPA_BUBBLE 1

#include "Ex_Precomputation.h"
#include "Ex_Diagnostics.h"
#include "Ex_Compute_bubble.h"
#include "Ex_multiplication.h"

template<int mode> class Ex_Compute_rpa{
	public:
		static const double measure_flow = 1.0;
		double Lambda;
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution<mode> &sub;
		Barevertex &barevertex;
		Ex_Vertex<mode> &gamma_rpa;
		Ex_Diagnostics &diagnostics;
		Ex_Compute_bubble<mode> bubble_computer;
		Ex_Compute_rpa(double Lambda_in,
		               Physics &phy_in,
		               Numerics &num_in,
		               Ex_Precomputation<mode> &pre_in,
		               Substitution<mode> &sub_in,
		               Barevertex &barevertex_in,
		               Ex_Vertex<mode> &gamma_rpa_in,
		               Ex_Diagnostics &diagnostics_in);
		void p_vertex(double accuracy, matrix<double> &additional_stops);
		void xd_vertex(double accuracy, matrix<double> &additional_stops);
};

template<int mode> Ex_Compute_rpa<mode>::Ex_Compute_rpa(double Lambda_in,
                                                        Physics &phy_in,
                                                        Numerics &num_in,
                                                        Ex_Precomputation<mode> &pre_in,
                                                        Substitution<mode> &sub_in,
                                                        Barevertex &barevertex_in,
                                                        Ex_Vertex<mode> &gamma_rpa_in,
                                                        Ex_Diagnostics &diagnostics_in):
                                                        Lambda(Lambda_in),
                                                        phy(phy_in),
                                                        num(num_in),
                                                        pre(pre_in),
                                                        sub(sub_in),
                                                        barevertex(barevertex_in),
                                                        gamma_rpa(gamma_rpa_in),
                                                        diagnostics(diagnostics_in),
                                                        bubble_computer(phy,num,pre,sub,measure_flow,Lambda,diagnostics){
}

template<int mode> void Ex_Compute_rpa<mode>::p_vertex(double accuracy, matrix<double> &additional_stops){
	//Puu contribution:
	{
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		Bubble.resize();
		bubble_computer.compute_Puu_bubble(additional_stops,Bubble,accuracy);
		int Lges = 2*num.L+1;
		matrix<matrix<complex<double> > > str_bare = 0.5*barevertex_to_p_str(1,1,num.L,num.N,barevertex);
		matrix<matrix<complex<double> > > str_unit = unit_str<complex<double> >(num.L, num.N);
		for(int i=0; i<num.NfbP; ++i){
			double freq = num.wbP(i);
			matrix<matrix<complex<double> > > str_bubble = Bubble.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds);
			matrix<matrix<complex<double> > > tmp = str_unit - str_bare*str_bubble;  
			tmp = invert_str(tmp);
			tmp = tmp*str_bare -str_bare;
			gamma_rpa.aPuu.dynamic_str(i) = block_core(tmp,num.Lp_structure(i)); 
			if(i==num.pos_NfbP_2mu){
				gamma_rpa.aPuu.static_str = real_real_part_of_str(tmp);
			}
		}
		#if(RPA_BUBBLE_ONLY==1)	
			gamma_rpa.aPuu.dynamic_str = Bubble.dynamic_str; 
			gamma_rpa.aPuu.static_str = Bubble.static_str; 
	
		#endif
	}
	//Pdd contribution:
	{
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		Bubble.resize();
		bubble_computer.compute_Pdd_bubble(additional_stops,Bubble,accuracy);
		int Lges = 2*num.L+1;
		matrix<matrix<complex<double> > > str_bare = 0.5*barevertex_to_p_str(0,0,num.L,num.N,barevertex);
		matrix<matrix<complex<double> > > str_unit = unit_str<complex<double> >(num.L, num.N);
		for(int i=0; i<num.NfbP; ++i){
			double freq = num.wbP(i);
			matrix<matrix<complex<double> > > str_bubble = Bubble.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds);
			matrix<matrix<complex<double> > > tmp = str_unit - str_bare*str_bubble;  
			tmp = invert_str(tmp);
			tmp = tmp*str_bare-str_bare;
			gamma_rpa.aPdd.dynamic_str(i) = block_core(tmp,num.Lp_structure(i)); 
			if(i==num.pos_NfbP_2mu){
				gamma_rpa.aPdd.static_str = real_real_part_of_str(tmp);
			}
		}
		#if(RPA_BUBBLE_ONLY==1)	
			gamma_rpa.aPdd.dynamic_str = Bubble.dynamic_str; 
			gamma_rpa.aPdd.static_str = Bubble.static_str; 
	
		#endif
	}
	//Pud contribution:
	{
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		Bubble.resize();
		bubble_computer.compute_Pud_bubble(additional_stops,Bubble,accuracy);
		int Lges = 2*num.L+1;
		matrix<matrix<complex<double> > > str_bare = 0.5*barevertex_to_p_str(1,0,num.L,num.N,barevertex);
		matrix<matrix<complex<double> > > str_unit = unit_str<complex<double> >(num.L, num.N);
		for(int i=0; i<num.NfbP; ++i){
			double freq = num.wbP(i);
			matrix<matrix<complex<double> > > str_bubble = Bubble.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds);
			matrix<matrix<complex<double> > > tmp = str_unit - str_bare*str_bubble;  
			tmp = invert_str(tmp);
			tmp = tmp*str_bare-str_bare;
			gamma_rpa.aPud.dynamic_str(i) = block_core(tmp,num.Lp_structure(i)); 
			if(i==num.pos_NfbP_2mu){
				gamma_rpa.aPud.static_str = real_real_part_of_str(tmp);
			}
		}
		#if(RPA_BUBBLE_ONLY==1)	
			gamma_rpa.aPud.dynamic_str = Bubble.dynamic_str; 
			gamma_rpa.aPud.static_str = Bubble.static_str; 
	
		#endif
		//Testing of various relations (here a contains the barevertex):
		#if(TEST_RPA_RELATIONS==1)
			cout<<"TEST_RPA_RELATIONS"<<endl;
			//RPA:
			matrix<matrix<matrix<complex<double> > > > a_semi_static_rpa(num.NfbP); 
			for(int i=0; i<num.NfbP; ++i){
			 	resize_str(a_semi_static_rpa(i),num.L,num.N);
				init(a_semi_static_rpa(i),(complex<double>)0.0);
			}
			for(int i=0; i<num.NfbP; ++i){
			 	matrix<matrix<complex<double> > > a_semi_static;
				double freq = num.wbP(i);
				matrix<matrix<complex<double> > > Bubble_complete_str = Bubble.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds);
				a_semi_static = str_unit - str_bare*Bubble_complete_str; 
				a_semi_static = invert_str(a_semi_static);
				a_semi_static = a_semi_static * str_bare;
				a_semi_static_rpa(i) = a_semi_static;
			}

			//Alternative:
			cout<<"Alternative"<<endl;
			matrix<matrix<matrix<complex<double> > > > A_dyn = add_static_term(gamma_rpa.aPud.dynamic_str,str_bare); 
			matrix<matrix<double> > A_stat = gamma_rpa.aPud.static_str + real_real_part_of_str(str_bare);
			std::pair<matrix<matrix<matrix<complex<double> > > >,matrix<matrix<matrix<complex<double> > > > > a_semi_static_alt = rpa_semi_static(A_dyn,A_stat,A_dyn,A_stat,num.Lp_structure,num.pos_NfbP_2mu); 

			//matrix<matrix<matrix<complex<double> > > > a_semi_static_alt(num.NfbP); 
			//for(int i=0; i<num.NfbP; ++i){
			//	resize_str(a_semi_static_alt(i),num.L,num.N);
			//	init(a_semi_static_alt(i),(complex<double>)0.0);
			//}
			//for(int n=0; n<2*num.NfbP; ++n){
			// 	//Compute symmetrically around feedback freq:
			// 	int i;
			// 	if(n%2==0){
			//		i= num.pos_NfbP_2mu + (n+1)/2;
			//	}
			//	else{
			//		i= num.pos_NfbP_2mu - (n+1)/2;
			//	}
			//	if(i<0 || i>=num.NfbP) continue;

			// 	double freq=num.wbP(i);
			//	int L_inner = num.Lp_structure(i);
			//	matrix<matrix<complex<double> > > a_str = gamma_rpa.aPud.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds) + str_bare;
			//	matrix<matrix<complex<double> > > a_str_small = block_core(a_str,L_inner);
			//	//get prev index:
			//	int i_prev=get_index_of_prev_larger_L(i,num.L,num.Lp_structure,num.pos_NfbP_2mu);
			//	cout<<"i="<<i<<", i_prev="<<i_prev<<endl;
			//	if(i_prev==i){
			//		a_semi_static_alt(i) = a_str;
			//	}
			//	else{
			//		matrix<matrix<complex<double> > > str_unit_small = block_core(str_unit,L_inner);
			//		matrix<matrix<complex<double> > > Bubble_complete_str = Bubble.complete_ipol(freq,num.pos_NfbP_2mu,num.Lp_bounds);
			//		matrix<matrix<complex<double> > > Bubble_complete_str_small = block_core(Bubble_complete_str,L_inner);
			//	 	matrix<matrix<complex<double> > > Bubble_prev_str = Bubble.complete_ipol(num.wbP(i_prev),num.pos_NfbP_2mu,num.Lp_bounds);
			//	 	matrix<matrix<complex<double> > > Bubble_prev_str_small = block_core(Bubble_prev_str,L_inner);
			//		matrix<matrix<complex<double> > > a_prev_str = gamma_rpa.aPud.complete_ipol(num.wbP(i_prev),num.pos_NfbP_2mu,num.Lp_bounds) + str_bare;
			//		matrix<matrix<complex<double> > > a_prev_str_small = block_core(a_prev_str,L_inner);

			//		matrix<matrix<complex<double> > > tmp1 = str_unit_small + Bubble_prev_str_small*a_prev_str_small; 
			//		tmp1=invert_str(tmp1);
			//		matrix<matrix<complex<double> > > tmp2 = str_unit_small + Bubble_complete_str_small*a_str_small; 
			//		matrix<matrix<complex<double> > > tmp3 = tmp1*tmp2; 
			//		for(int l=-num.L; l<=num.L; ++l){
			//			for(int k=-L_inner; k<=L_inner; ++k){
			//				for(int k1=-L_inner; k1<=L_inner; ++k1){
			//				a_semi_static_alt(i)(l+num.L,k+num.L) += a_semi_static_alt(i_prev)(l+num.L,k1+num.L)*tmp3(k1+L_inner,k+L_inner);
			//				}
			//			}
			//		}
			//	}
			//}

			//Comparison:
			matrix<double> differences(num.NfbP);
			differences=0.0;
			for(int i=0; i<num.NfbP; ++i){
				int L_inner = num.Lp_structure(i);
				double diff=0.0;
				for(int l=-num.L; l<=num.L; ++l){
					for(int k=-L_inner; k<=L_inner; ++k){
						diff = max(diff, abs(a_semi_static_rpa(i)(k+num.L,l+num.L)-a_semi_static_alt.first(i)(k+num.L,l+num.L)) );
						diff = max(diff, abs(a_semi_static_rpa(i)(l+num.L,k+num.L)-a_semi_static_alt.second(i)(l+num.L,k+num.L)) );
					}
				}
				differences(i)=diff;
				cout<<"i="<<i<<", diff="<<diff<<endl;
			}
			cout<<"abs(differences)="<<abs(differences)<<endl;
		#endif
	}
}

template<int mode> void Ex_Compute_rpa<mode>::xd_vertex(double accuracy, matrix<double> &additional_stops){
	//Xud contribution:
	{
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbX);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn, bubble_stat);
		Bubble.resize();
		bubble_computer.compute_Xud_bubble(additional_stops,Bubble,accuracy);
		int Lges = 2*num.L+1;
		matrix<matrix<complex<double> > > str_bare = 0.5*barevertex_to_x_str(1,0,num.L,num.N,barevertex);
		matrix<matrix<complex<double> > > str_unit = unit_str<complex<double> >(num.L, num.N);
		for(int i=0; i<num.NfbX; ++i){
			double freq = num.wbX(i);
			matrix<matrix<complex<double> > > str_bubble = Bubble.complete_ipol(freq,num.pos_NfbX_0,num.Lx_bounds);
			matrix<matrix<complex<double> > > tmp = str_unit - str_bare*str_bubble;  
			tmp = invert_str(tmp);
			tmp = tmp*str_bare-str_bare;
			gamma_rpa.aXud.dynamic_str(i) = block_core(tmp,num.Lx_structure(i)); 
			if(i==num.pos_NfbX_0){
				gamma_rpa.aXud.static_str = real_real_part_of_str(tmp);
			}
		}
		#if(RPA_BUBBLE_ONLY==1)	
			gamma_rpa.aXud.dynamic_str = Bubble.dynamic_str; 
			gamma_rpa.aXud.static_str = Bubble.static_str; 
	
		#endif
	}
	//D contribution
	{
		matrix<matrix<matrix<complex<double> > > > block_unit_str(2,2); 		
		resize_str(block_unit_str(0,1),num.L,num.N);
		init(block_unit_str(0,1),(complex<double>) 0.0);
		resize_str(block_unit_str(1,0),num.L,num.N);
		init(block_unit_str(1,0),(complex<double>) 0.0);
		block_unit_str(0,0) = unit_str<complex<double> >(num.L, num.N);
		block_unit_str(1,1) = unit_str<complex<double> >(num.L, num.N);
		

		matrix<matrix<matrix<complex<double> > > > str_block_bare(2,2);
		str_block_bare(0,0) = 0.5*barevertex_to_d_str(1,1,num.L,num.N,barevertex);
		str_block_bare(0,1) = 0.5*barevertex_to_d_str(1,0,num.L,num.N,barevertex);
		str_block_bare(1,0) = 0.5*barevertex_to_d_str(0,1,num.L,num.N,barevertex);
		str_block_bare(1,1) = 0.5*barevertex_to_d_str(0,0,num.L,num.N,barevertex);


		//up-up contribution:
		matrix<matrix<matrix<complex<double> > > >  bubble_uu_dyn(num.NfbX);
		matrix<matrix<double> >  bubble_uu_stat;
		Ex_freq_str Bubble_uu(num.L, num.N, num.Lx_structure, num.wbX, bubble_uu_dyn, bubble_uu_stat);
		Bubble_uu.resize();
		bubble_computer.compute_Duu_bubble(additional_stops,Bubble_uu,accuracy);

		//down-down contribution:
		matrix<matrix<matrix<complex<double> > > >  bubble_dd_dyn(num.NfbX);
		matrix<matrix<double> >  bubble_dd_stat;
		Ex_freq_str Bubble_dd(num.L, num.N, num.Lx_structure, num.wbX, bubble_dd_dyn, bubble_dd_stat);
		Bubble_dd.resize();
		bubble_computer.compute_Ddd_bubble(additional_stops,Bubble_dd,accuracy);

		for(int i=0; i<num.NfbX; ++i){
			double freq = num.wbX(i);
			matrix<matrix<matrix<complex<double> > > > str_block_bubble(2,2);
			resize_str(str_block_bubble(0,1),num.L,num.N);
			init(str_block_bubble(0,1),(complex<double>) 0.0);
			resize_str(str_block_bubble(1,0),num.L,num.N);
			init(str_block_bubble(1,0),(complex<double>) 0.0);
			str_block_bubble(0,0) = Bubble_uu.complete_ipol(freq,num.pos_NfbX_0,num.Lx_bounds);
			str_block_bubble(1,1) = Bubble_dd.complete_ipol(freq,num.pos_NfbX_0,num.Lx_bounds);
		
			matrix<matrix<matrix<complex<double> > > > tmp = block_unit_str + str_block_bare*str_block_bubble; 
			//Invert tmp as blockmatrix:	
			matrix<matrix<matrix<complex<double> > > > tmp_inv(2,2);  
			matrix<matrix<complex<double> > > tmp2 = tmp(1,1);
			tmp2 = invert_str(tmp2);
			tmp_inv(0,0) = tmp(0,0) - tmp(0,1)*tmp2*tmp(1,0);
			tmp_inv(0,0) = invert_str(tmp_inv(0,0));
			tmp_inv(0,1) = -tmp_inv(0,0)*tmp(0,1)*tmp2;
			tmp_inv(1,0) = -tmp2*tmp(1,0)*tmp_inv(0,0);
			tmp_inv(1,1) = tmp2 + tmp2*tmp(1,0)*tmp_inv(0,0)*tmp(0,1)*tmp2;

			tmp = tmp_inv*str_block_bare - str_block_bare;
			gamma_rpa.aDuu.dynamic_str(i) = block_core(tmp(0,0),num.Lx_structure(i)); 
			gamma_rpa.aDdd.dynamic_str(i) = block_core(tmp(1,1),num.Lx_structure(i)); 
			gamma_rpa.aDud.dynamic_str(i) = block_core(tmp(0,1),num.Lx_structure(i)); 
			if(i==num.pos_NfbX_0){
				gamma_rpa.aDuu.static_str = real_real_part_of_str(tmp(0,0));
				gamma_rpa.aDdd.static_str = real_real_part_of_str(tmp(1,1));
				gamma_rpa.aDud.static_str = real_real_part_of_str(tmp(0,1));
			}
		}
	}
}











#endif

