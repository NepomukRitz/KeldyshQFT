#ifndef EX_FREQ_STR_30102018 
#define EX_FREQ_STR_30102018

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "Numerics.h"
#include "Ex_functions.h"

using namespace std;

class Ex_freq_str{
	public:
		int L;
		int N;
		matrix<int> &L_structure;
		matrix<double> &wb;
		matrix<matrix<matrix<complex<double> > > > &dynamic_str;
		matrix<matrix<double> > &static_str;
		int N_freq;
		int Nges;
		Ex_freq_str(int L_in,
		            int N_in,
		            matrix<int> &L_structure_in,
		            matrix<double> &wb_in,
		            matrix<matrix<matrix<complex<double> > > > &dynamic_str_in,
		            matrix<matrix<double> > &static_str_in
		            );
		void resize();
		void initialize(double init);
		void initialize_random(int D);
		bool inrange_stat(int l, int k, int j, int i); 
		double get_stat(int l, int k, int j, int i); 
		matrix<matrix<complex<double> > > dyn_ipol(double freq, int pos_feedback);
		matrix<matrix<complex<double> > > dyn_ipol_b(double freq, int pos_feedback, matrix<double> &b_pre);
		matrix<matrix<complex<double> > > complete_ipol(double freq, int pos_feedback, matrix<matrix<int> > &L_bounds);
		matrix<matrix<complex<double> > > complete_ipol_b(double freq, int pos_feedback, matrix<matrix<int> > &L_bounds, matrix<double> &b_pre);
		matrix<matrix<matrix<complex<double> > > > static_extensions_a(matrix<matrix<int> > &L_bounds);
		matrix<matrix<matrix<complex<double> > > > static_extensions_b(matrix<matrix<int> > &L_bounds, matrix<double> &prefactors);
		void save(string filename, string variable);
		void load(string filename, string variable);
};

Ex_freq_str::Ex_freq_str(int L_in,
                         int N_in,
                         matrix<int> &L_structure_in,
                         matrix<double> &wb_in,
                         matrix<matrix<matrix<complex<double> > > > &dynamic_str_in,
                         matrix<matrix<double> > &static_str_in):
                         L(L_in),
                         N(N_in),
                         L_structure(L_structure_in),
                         wb(wb_in),
                         dynamic_str(dynamic_str_in),
                         static_str(static_str_in),
                         N_freq(L_structure.dim_c),
                         Nges(2*N+1){
}

void Ex_freq_str::resize(){
	resize_str(dynamic_str,L_structure,N);
	resize_str(static_str,L,N);
}

void Ex_freq_str::initialize(double d){
	init(dynamic_str, (complex<double>) d);
	init(static_str, d);
}

void Ex_freq_str::initialize_random(int D){
	init_random(dynamic_str, D);
	init_random(static_str, D);
}

bool Ex_freq_str::inrange_stat(int l, int k, int j, int i){
	return (-L<=l && l<=L && -L<=k && k<=L && max(0,-l)<=j && j<min(Nges,Nges-l) && max(0,-k)<=i && i<min(Nges,Nges-k));
}

double Ex_freq_str::get_stat(int l, int k, int j, int i){
	if(inrange_stat(l,k,j,i)){
		return static_str(l+L,k+L)(j-max(0,-l),i-max(0,-k));
	}
	else{
		return 0.0;
	}
}
		
matrix<matrix<complex<double> > > Ex_freq_str::dyn_ipol(double freq, int pos_feedback){
	int index = geq_freq_index(wb, freq); 
	if(index==0){
		return dynamic_str(0);
	}
	else{
		double f1 = wb(index-1); 
		double f2 = wb(index); 
		double diff = f2-f1;
		double w1 = (f2 - freq)/diff; 
		double w2 = 1.0 - w1;
		matrix<matrix<complex<double> > > tmp1;
		matrix<matrix<complex<double> > > tmp2;
		if((index-1) !=pos_feedback){
			tmp1 = dynamic_str(index-1);
		}
		else{
			cast(tmp1,static_str); 
		}
		if(index != pos_feedback){
			tmp2 = dynamic_str(index);
		}
		else{
			cast(tmp2,static_str); 
		}
		return str_weighted_sum(w1,tmp1, w2, tmp2);
	}
}

matrix<matrix<complex<double> > > Ex_freq_str::dyn_ipol_b(double freq, int pos_feedback, matrix<double> &b_pre){
	int index = geq_freq_index(wb, freq); 
	if(index==0){
		double pre = b_pre(index); 
		return pre*imaginary_part_of_str(dynamic_str(index));
	}
	else{
		double f1 = wb(index-1); 
		double f2 = wb(index); 
		double diff = f2-f1;
		double w1 = (f2 - freq)/diff; 
		double w2 = 1.0 - w1;
		double pre1 = b_pre(index-1); 
		double pre2 = b_pre(index); 
		matrix<matrix<complex<double> > > tmp1 = pre1*dynamic_str(index-1);
		matrix<matrix<complex<double> > > tmp2 = pre2*dynamic_str(index);
		if((index-1) !=pos_feedback){
			tmp1 = pre1*imaginary_part_of_str(dynamic_str(index-1));
		}
		else{
			resize_str(tmp1,L,N);
			init(tmp1,(complex<double>) 0.0);
		}
		if(index != pos_feedback){
			tmp2 = pre2*imaginary_part_of_str(dynamic_str(index));
		}
		else{
			resize_str(tmp2,L,N);
			init(tmp2,(complex<double>) 0.0);
		}
		return str_weighted_sum(w1,tmp1, w2, tmp2);
	}
}

//This function is not optimized for speed!	
matrix<matrix<complex<double> > > Ex_freq_str::complete_ipol(double freq, int pos_feedback, matrix<matrix<int> > &L_bounds){
	int Lges = 2*L+1;
	double freq_feedback = wb(pos_feedback);
	matrix<matrix<complex<double> > > ret(Lges,Lges);
	matrix<matrix<complex<double> > > ret_dyn = dyn_ipol(freq, pos_feedback);
	int Ldynges = ret_dyn.dim_c;
	int Ldyn= (Ldynges - 1)/2;
	for(int l=-L; l<=L; ++l){	
		for(int k=-L; k<=L; ++k){	
			if(abs(l)<=Ldyn && abs(k)<=Ldyn){
				ret(l+L,k+L) = ret_dyn(l+Ldyn,k+Ldyn);
			}
			else{
				#if(LONG_RANGE_EXTRAPOLATION==1)
					//ret(l+L,k+L).resize(Nges-abs(l),Nges-abs(k)); //Debug
					//ret(l+L,k+L) = (complex<double>) 0.0; //Debug
					int index;
					if(freq<freq_feedback){
						index = L_bounds(0)(l+L,k+L); 
					}
					else{
						index = L_bounds(1)(l+L,k+L); 
					}
					if(index==-1){
						ret(l+L,k+L) = static_str(l+L,k+L); 
					}
					else{
						int Ls = L_structure(index);
						ret(l+L,k+L) = dynamic_str(index)(l+Ls,k+Ls); 
					}
				#else
					ret(l+L,k+L) = static_str(l+L,k+L);
				#endif
			}
		}
	}
	return ret;
}

matrix<matrix<complex<double> > > Ex_freq_str::complete_ipol_b(double freq, int pos_feedback, matrix<matrix<int> > &L_bounds, matrix<double> &b_pre){
	int Lges = 2*L+1;
	double freq_feedback = wb(pos_feedback);
	matrix<matrix<complex<double> > > ret(Lges,Lges);
	matrix<matrix<complex<double> > > ret_dyn = dyn_ipol_b(freq, pos_feedback, b_pre);
	int Ldynges = ret_dyn.dim_c;
	int Ldyn= (Ldynges - 1)/2;
	for(int l=-L; l<=L; ++l){	
		for(int k=-L; k<=L; ++k){	
			if(abs(l)<=Ldyn && abs(k)<=Ldyn){
				ret(l+L,k+L) = ret_dyn(l+Ldyn,k+Ldyn);
			}
			else{
				#if(LONG_RANGE_EXTRAPOLATION==1)
					//ret(l+L,k+L).resize(Nges-abs(l),Nges-abs(k)); //Debug
					//ret(l+L,k+L) = (complex<double>) 0.0; //Debug
					int index;
					if(freq<freq_feedback){
						index = L_bounds(0)(l+L,k+L); 
					}
					else{
						index = L_bounds(1)(l+L,k+L); 
					}
					if(index==-1){
						ret(l+L,k+L).resize(Nges-abs(l),Nges-abs(k)); 
						ret(l+L,k+L) = (complex<double>) 0.0; 
					}
					else{
						int Ls = L_structure(index);
						ret(l+L,k+L) = b_pre(index)*dynamic_str(index)(l+Ls,k+Ls).imag(); 
					}
				#else
					ret(l+L,k+L).resize(Nges-abs(l),Nges-abs(k)); 
					ret(l+L,k+L) = (complex<double>) 0.0; 
				#endif
			}
		}
	}
	return ret;
}

matrix<matrix<matrix<complex<double> > > > Ex_freq_str::static_extensions_a(matrix<matrix<int> > &L_bounds){
	matrix<matrix<matrix<complex<double> > > > ret(2);
	#if(LONG_RANGE_EXTRAPOLATION==1)
		for(int s=0; s<2; ++s){
			ret(s).resize(2*L+1,2*L+1);
			for(int l=-L; l<=L; ++l){
				for(int k=-L; k<=L; ++k){
					int i = L_bounds(s)(l+L,k+L);
					if(i==-1){
						ret(s)(l+L,k+L) = static_str(l+L,k+L);
					}
					else{
						int Li = L_structure(i);
						ret(s)(l+L,k+L) = dynamic_str(i)(l+Li,k+Li);
					}
				}
			}
		}
	#else
		cast(ret(0),static_str); 
		cast(ret(1),static_str); 
	#endif
	return ret;
}

matrix<matrix<matrix<complex<double> > > > Ex_freq_str::static_extensions_b(matrix<matrix<int> > &L_bounds, matrix<double> &prefactors){
	matrix<matrix<matrix<complex<double> > > > ret(2);
	for(int s=0; s<2; ++s){
		ret(s).resize(2*L+1,2*L+1);
		for(int l=-L; l<=L; ++l){
			for(int k=-L; k<=L; ++k){
				int i = L_bounds(s)(l+L,k+L);
				if(i==-1){
					ret(s)(l+L,k+L).resize(Nges-abs(l),Nges-abs(k));
					ret(s)(l+L,k+L) = (complex<double>) 0.0;
				}
				else{
					int Li = L_structure(i);
					ret(s)(l+L,k+L) = prefactors(i)*dynamic_str(i)(l+Li,k+Li).imag();
				}
			}
		}
	}
	return ret;
}

void Ex_freq_str::save(string filename, string variable){
	matrix<double> tmp(1);
	tmp(0) = L;
	tmp.save(filename.c_str(),(variable+"_L").c_str());
	tmp(0) = N;
	tmp.save(filename.c_str(),(variable+"_N").c_str());
	L_structure.save(filename.c_str(),(variable+"_L_structure").c_str());
	wb.save(filename.c_str(),(variable+"_wb").c_str());
	save_str(dynamic_str,filename,variable+"_dynamic_str");
	save_str(static_str,filename,variable+"_static_str");
}

void Ex_freq_str::load(string filename, string variable){
	matrix<double> tmp;
	tmp.load(filename.c_str(),(variable+"_L").c_str());
	L = tmp(0);
	tmp.load(filename.c_str(),(variable+"_N").c_str());
	N = tmp(0);
	L_structure.load(filename.c_str(),(variable+"_L_structure").c_str());	
	wb.load(filename.c_str(),(variable+"_wb").c_str());
	load_str(dynamic_str,filename,variable+"_dynamic_str");
	load_str(static_str,filename,variable+"_static_str");
	N_freq = L_structure.dim_c;
	Nges = 2*N+1;
}

#endif
