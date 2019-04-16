#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=0;
	srand(seed);
	int Nff=5;
	int L=2;
	int N=1;
	int D=100;
	int pos_feedback = 2;
	int Nges=2*N+1;
	int Lges=2*L+1;
	matrix<int> L_structure(Nff);
	//init_monoton_L_structure(L,Nff,pos_feedback,L_structure);	
	L_structure(0) = 1;
	L_structure(1) = 1;
	L_structure(2) = 0;
	L_structure(3) = 1;
	L_structure(4) = 0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<int> > L_steps = determine_L_steps(L, L_structure, pos_feedback);
		cout<<"L_steps="<<endl;
		cout<<L_steps<<endl;
		matrix<matrix<syma<complex<double> > > > A(2);
		for(int s=0; s<2; ++s){
			A(s).resize(L_steps(s).dim_c-1);
			A(s)(0).resize(Nges);
			init_random(A(s)(0),D);
			for(int i=1; i<L_steps(s).dim_c-1; ++i){
				syma<complex<double> > tmp(Nges);
				init_random(tmp,D);
				A(s)(i) = A(s)(i-1) + tmp;
			}
		}
		matrix<matrix<complex<double> > > B; 
		resize_str(B,L,N);
		complex<double> I(0.0,1.0);
		init(B, 1.0+0.*I);
		//for(int l=-L; l<=L; ++l){
		//	for(int k=-L; k<=L; ++k){
		//		if(in_ring(l,k,1,2)){
		//			B(l+L,k+L) = (complex<double>) 0.0;
		//		}
		//	}
		//}
		//Test:
		cout<<"A(0).dim_c="<<A(0).dim_c<<endl;
		matrix<complex<double> > C = compute_ring_contributions(L, N, L_steps, B, A);
		cout<<"abs(C)="<<abs(C)<<endl;
		cout<<"A(0)(0)="<<endl;
		cout<<A(0)(0)<<endl;
		cout<<"A(1)(0)="<<endl;
		cout<<A(1)(0)<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;
		
		matrix<complex<double> > D(Nges,Nges);
		D=(complex<double>) 0.0;
		//lower contribution:
		for(int k=-L; k<=L; ++k){
			int l=-2;
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
					D(jc,ic) += B(l+L,k+L)(jc,ic)*A(0)(0).full_access(j+l,i+k);
					D(ic,jc) += B(k+L,l+L)(ic,jc)*A(0)(0).full_access(i+k,j+l);
				}
			}
			l=+2;
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
					D(jc,ic) += B(l+L,k+L)(jc,ic)*A(0)(0).full_access(j+l,i+k);
					D(ic,jc) += B(k+L,l+L)(ic,jc)*A(0)(0).full_access(i+k,j+l);
				}
			}
		}
		int l=-2, k=-2;
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
			for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
				D(jc,ic) -= B(l+L,k+L)(jc,ic)*A(0)(0).full_access(j+l,i+k);
			}
		}
		l=+2, k=-2;
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
			for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
				D(jc,ic) -= B(l+L,k+L)(jc,ic)*A(0)(0).full_access(j+l,i+k);
				D(ic,jc) -= B(k+L,l+L)(ic,jc)*A(0)(0).full_access(i+k,j+l);
			}
		}
		l=+2, k=+2;
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
			for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
				D(jc,ic) -= B(l+L,k+L)(jc,ic)*A(0)(0).full_access(j+l,i+k);
			}
		}


		//higher contribution:
		for(int l=-1; l<=1; ++l){
			for(int k=-1; k<=1; ++k){
				for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
					for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
						D(jc,ic)+= B(l+L,k+L)(jc,ic)*A(1)(0).full_access(j+l,i+k);
					}
				}
			}
		}
		l=0,k=0;
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
			for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
				D(jc,ic)-= B(l+L,k+L)(jc,ic)*A(1)(0).full_access(j+l,i+k);
			}
		}
		for(int l=-2; l<=2; ++l){
			for(int k=-2; k<=2; ++k){
				for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
					for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
						D(jc,ic)+= B(l+L,k+L)(jc,ic)*A(1)(1).full_access(j+l,i+k);
					}
				}
			}
		}
		for(int l=-1; l<=1; ++l){
			for(int k=-1; k<=1; ++k){
				for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
					for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
						D(jc,ic)-= B(l+L,k+L)(jc,ic)*A(1)(1).full_access(j+l,i+k);
					}
				}
			}
		}
		cout<<"abs(C-D)="<<abs(C-D)<<endl;
		cout<<"C="<<endl;
		cout<<C<<endl;
		cout<<"D="<<endl;
		cout<<D<<endl;
		
		
	}
	return 0;
}

