#ifndef FLOW_ZERO_MAG_28072017
#define FLOW_ZERO_MAG_28072017

#include <omp.h>
#include <ctime>
#include <integrate_new.h>
#include "Precomputation.h"
#include "Generalmatrix.h"
#include "Substitution_flow.h"
#include "Barevertex.h"
#include "Vertex.h"
#include "P_bubble_feedback_zero_mag.h"
#include "P_bubble_central_zero_mag.h"
#include "X_bubble_feedback_zero_mag.h"
#include "X_bubble_central_zero_mag.h"
#include "Self_energy_central_zero_mag.h"

//Reduce memory allocation!

template <int mode> class Flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution_flow sub_flow;
		Barevertex &barevertex;
		Flow_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Barevertex &barevertex_in);
		void operator()(double x, Generalmatrix &y, Generalmatrix &dy);
};

template <int mode> Flow_zero_mag<mode>::Flow_zero_mag(Physics &phy_in,
                                                       Numerics &num_in,
										               Precomputation_zeromag<mode> &pre_in,
													   Barevertex &barevertex_in): 
													   phy(phy_in),
													   num(num_in),
													   pre(pre_in), 
													   barevertex(barevertex_in){
}

template <int mode> void Flow_zero_mag<mode>::operator()(double x, Generalmatrix &y, Generalmatrix &dy){
	time_t t1, t2;
	double Lambda = sub_flow.resu(x);
	double measure_flow = sub_flow.weight(x);
 	cout<<"Lambda="<<Lambda<<", x="<<x<<endl;
	Substitution<mode> sub(Lambda);
	Vertex<mode> gamma(num,sub,y);
	dy.resize(num);
	dy.initialize(0.0);
	Vertex<mode> dgamma(num,sub,dy);
 	pre.precompute(Lambda,sub,gamma.ERetu_ipol_subst);


	/*P_channel:*/


	P_bubble_feedback_zero_mag<mode> P_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow); 


	/*First the static contributions to aP:*/

	/*static P_bubble:*/
	
	matrix<matrix<double> > Bubble_data;
	Blockmatrix<double> Bubble(num.L, num.N, Bubble_data);
	Bubble.resize(num.L, num.N);
	time(&t1);
	omp_set_num_threads(16);
	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=l; ++k){
		 	Bubble_data(l+num.L,k+num.L) = P_bubble_feedback(l,k);
		 	Bubble_data(k+num.L,l+num.L) = Bubble_data(l+num.L,k+num.L).transp();
		}
	}
	time(&t2);
	cout<<"Time for static P_bubble="<<t2 - t1<<endl;
	//Bubble_data.save("operator.mat","Bubble_data");


	time(&t1);
	matrix<double> static_subtraction_uu(num.Nges); //Not elegant yet
	matrix<double> static_subtraction_ud(num.Nges); //Not elegant yet
	matrix<double> static_subtraction_ud_2(num.Nges); //Not elegant yet
	
	/*daPuu:*/
	{	
		matrix<matrix<double> > P_data;
		Blockmatrix<double> P(num.L, num.N, P_data);
		P.resize(num.L, num.N);

		omp_set_num_threads(16); //Checken, ob omp hier sinnvoll ist
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
				for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
					for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	    P(l,k,j,i)  = 0.5*barevertex(j,1,j+l,1,i,1,i+k,1)
						                 +gamma.aPuu_feedback(l,k,j,i);
						if(gamma.aDuu_feedback.inrange(i+k-j,j+l-i,j,i)){
						 	P(l,k,j,i) -= gamma.aDuu_feedback(i+k-j,j+l-i,j,i);
						}
						if(gamma.aDuu_feedback.inrange(i-j,j+l-i-k,j,i+k)){
						 	P(l,k,j,i) += gamma.aDuu_feedback(i-j,j+l-i-k,j,i+k);
						}
						 	
					}
				}
			}
		}


		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){
			 	for(int q=-num.L; q<=num.L; ++q){
				 	for(int p=-num.L; p<=num.L; ++p){
					 	dgamma.aPuu_feedback_data(l+num.L, k+num.L) = P_data(l+num.L,q+num.L)*Bubble_data(q+num.L,p+num.L)*P_data(p+num.L,k+num.L);	
					}
				}
			}
		}
	 	static_subtraction_uu = P_data(num.L, num.L)*Bubble_data(num.L,num.L)*P_data(num.L, num.L); 
		
	}

	/*daPdd:*/
	dgamma.aPdd_feedback_data = dgamma.aPuu_feedback_data;

	
	/*daPud:*/
	{
		matrix<matrix<double> > P_data;
		Blockmatrix<double> P(num.L, num.N, P_data);
		P.resize(num.L, num.N);

		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
				for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
					for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	    P(l,k,j,i)  = 0.5*barevertex(j,1,j+l,0,i,1,i+k,0)
						                 +gamma.aPud_feedback(l,k,j,i);
						if(gamma.aXud_feedback.inrange(i+k-j,j+l-i,j,i)){
						 	P(l,k,j,i) += gamma.aXud_feedback(i+k-j,j+l-i,j,i);
						}
						if(gamma.aDud_feedback.inrange(i-j,j+l-i-k,j,i+k)){
						 	P(l,k,j,i) += gamma.aDud_feedback(i-j,j+l-i-k,j,i+k);
						}
						 	
					}
				}
			}
		}
		//P_data.save("operator.mat","P_data");


		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){
			 	for(int q=-num.L; q<=num.L; ++q){
				 	for(int p=-num.L; p<=num.L; ++p){
					 	dgamma.aPud_feedback_data(l+num.L, k+num.L) += P_data(l+num.L,q+num.L)*(Bubble_data(q+num.L,p+num.L)+Bubble_data(-q+num.L,-p+num.L))*P_data(p+num.L,k+num.L);	
					}
				}
			}
		}
	 	static_subtraction_ud = P_data(num.L, num.L)*Bubble_data(num.L,num.L)*P_data(num.L, num.L); 
	}
	time(&t2);
	cout<<"Time for static P_flow multiplication="<<t2 - t1<<endl;

	/*Dynamic Contribution to aP:*/
	
	time(&t1);
	{

		P_bubble_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda, measure_flow); 
	
		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
		 	syma<complex<double> > Puu(num.Nges);
		 	syma<complex<double> > Pud(num.Nges);
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
				 	        Puu(j1,j2)  = 0.5*barevertex(j1,1,j1,1,j2,1,j2,1)
					                     +gamma.aPuu_central(i)(j1,j2);
				 	        Pud(j1,j2)  = 0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
					                     +gamma.aPud_central(i)(j1,j2);
						if(j1==j2){
						 	Pud(j1,j2)  += gamma.aXud_feedback_data(num.L,num.L)(j1,j2) 
							              +gamma.aDud_feedback_data(num.L,num.L)(j1,j2); 
						}
				}
			}
		
			syma<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Bubble(num.wbP(i));
			dgamma.aPuu_central(i) = Puu*Bubble_at_freq*Puu;
			dgamma.aPud_central(i) = Pud*Bubble_at_freq*Pud;
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<j1; ++j2){
		 			dgamma.aPuu_central(i)(j1,j2) += (complex<double>)( dgamma.aPuu_feedback_data(num.L,num.L)(j1,j2) 
										                               -static_subtraction_uu(j1,j2)); 
		 			dgamma.aPud_central(i)(j1,j2) += (complex<double>)( dgamma.aPud_feedback_data(num.L,num.L)(j1,j2) 
										                               -static_subtraction_ud(j1,j2)); 
				}
			}
		}





	}
	time(&t2);
	cout<<"Time for dynamic P_flow="<<t2 - t1<<endl;

	/*X_channel and D_channel*/
	time(&t1);

	X_bubble_feedback_zero_mag<mode> X_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow);
	
	/*First the static contributions to aX and aD:*/

	/*static X_bubble and D_bubble:*/
	
	omp_set_num_threads(16);
	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=l; ++k){
		 	Bubble_data(l+num.L,k+num.L) = X_bubble_feedback(l,k);
		 	Bubble_data(k+num.L,l+num.L) = Bubble_data(l+num.L,k+num.L).transp();
		}
	}
	//Bubble_data.save("operator.mat","Bubble_data");
	time(&t2);
	cout<<"Time for static X_bubble="<<t2 - t1<<endl;
	
	
	/*daXud:*/
	time(&t1);
	{
		matrix<matrix<double> > X_data;
		Blockmatrix<double> X(num.L, num.N, X_data);
		X.resize(num.L, num.N);

		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
				for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
					for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	    X(l,k,j,i)  = 0.5*barevertex(j,1,i+k,0,i,1,j+l,0)
						                 +gamma.aXud_feedback(l,k,j,i);
						if(gamma.aPud_feedback.inrange(i+k-j,j+l-i,j,i)){
						 	X(l,k,j,i) += gamma.aPud_feedback(i+k-j,j+l-i,j,i);
						}
						if(gamma.aDud_feedback.inrange(i-j,i+k-j-l,j,j+l)){
						 	X(l,k,j,i) += gamma.aDud_feedback(i-j,i+k-j-l,j,j+l);
						}
						 	
					}
				}
			}
		}
		//X_data.save("operator.mat","X_data");


		omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){
			 	for(int q=-num.L; q<=num.L; ++q){
				 	for(int p=-num.L; p<=num.L; ++p){
					 	dgamma.aXud_feedback_data(l+num.L, k+num.L) += X_data(l+num.L,q+num.L)*Bubble_data(q+num.L,p+num.L)*X_data(p+num.L,k+num.L);	
					}
				}
			}
		}
	 	static_subtraction_ud = X_data(num.L, num.L)*Bubble_data(num.L,num.L)*X_data(num.L, num.L); 
	}
	

	/*daDud, daDuu and daDdd:*/ //Minuszeichen beachten!
	{
		matrix<matrix<double> > D1_data;
		Blockmatrix<double> D1(num.L, num.N, D1_data);
		D1.resize(num.L, num.N);

		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
				for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
					for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	    D1(l,k,j,i)  = 0.5*barevertex(j,1,i+k,0,j+l,1,i,0)
						                 +gamma.aDud_feedback(l,k,j,i);
						if(gamma.aPud_feedback.inrange(i+k-j,i-j-l,j,j+l)){
						 	D1(l,k,j,i) += gamma.aPud_feedback(i+k-j,i-j-l,j,j+l);
						}
						if(gamma.aXud_feedback.inrange(i-j,i+k-j-l,j,j+l)){
						 	D1(l,k,j,i) += gamma.aXud_feedback(i-j,i+k-j-l,j,j+l);
						}
						 	
					}
				}
			}
		}
		
		matrix<matrix<double> > D2_data;
		Blockmatrix<double> D2(num.L, num.N, D2_data);
		D2.resize(num.L, num.N);

		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
				for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
					for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	    D2(l,k,j,i)  = 0.5*barevertex(j,0,i+k,0,j+l,0,i,0)
						                 +gamma.aDdd_feedback(l,k,j,i);
						if(gamma.aPdd_feedback.inrange(i+k-j,i-j-l,j,j+l)){
						 	D2(l,k,j,i) += gamma.aPdd_feedback(i+k-j,i-j-l,j,j+l);
						}
						if(gamma.aDdd_feedback.inrange(i-j,i+k-j-l,j,j+l)){
						 	D2(l,k,j,i) -= gamma.aDdd_feedback(i-j,i+k-j-l,j,j+l);
						}
						 	
					}
				}
			}
		}

		//D1_data.save("operator.mat","D1_data");
		//D2_data.save("operator.mat","D2_data");


		omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=num.L; ++k){
			 	for(int q=-num.L; q<=num.L; ++q){
				 	for(int p=-num.L; p<=num.L; ++p){
					 	dgamma.aDud_feedback_data(l+num.L, k+num.L) -= D1_data(l+num.L,q+num.L)*Bubble_data(-q+num.L,-p+num.L)*D2_data(p+num.L,k+num.L);	
						dgamma.aDuu_feedback_data(l+num.L, k+num.L) -= D2_data(l+num.L,q+num.L)*Bubble_data(-q+num.L,-p+num.L)*D2_data(p+num.L,k+num.L)
						                                              +D1_data(l+num.L,q+num.L)*Bubble_data(-q+num.L,-p+num.L)*D1_data(k+num.L,p+num.L).transp();
					}
				}
			}
		}
		omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
		#pragma omp parallel for
		for(int l=-num.L; l<=num.L; ++l){
		 	for(int k=-num.L; k<=l; ++k){
			 	dgamma.aDud_feedback_data(l+num.L, k+num.L) += dgamma.aDud_feedback_data(k+num.L, l+num.L).transp();	
				dgamma.aDud_feedback_data(k+num.L, l+num.L) = dgamma.aDud_feedback_data(l+num.L, k+num.L).transp();
			}
		}
		dgamma.aDdd_feedback_data = dgamma.aDuu_feedback_data;



	 	static_subtraction_ud_2 = -D1_data(num.L, num.L)*Bubble_data(num.L,num.L)*D2_data(num.L, num.L); 
		static_subtraction_ud_2 += static_subtraction_ud_2.transp();

		static_subtraction_uu = -D2_data(num.L, num.L)*Bubble_data(num.L,num.L)*D2_data(num.L, num.L)
		                        -D1_data(num.L, num.L)*Bubble_data(num.L,num.L)*D1_data(num.L, num.L).transp();

	}
	time(&t2);
	cout<<"Time for static XD_flow multiplication="<<t2 - t1<<endl;
	
	
	
	
	
	
	
	/*Dynamic Contribution to aX and aD:*/
	
	time(&t1);
	{

		X_bubble_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda, measure_flow); 
	
		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
			syma<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Bubble(num.wbX(i));
			{
			 	syma<complex<double> > Xud(num.Nges);
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
					 	        Xud(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j2,1,j1,0)
						                     +gamma.aXud_central(i)(j1,j2);
							if(j1==j2){
							 	Xud(j1,j2) += gamma.aPud_feedback_data(num.L,num.L)(j1,j1) 
								             +gamma.aDud_feedback_data(num.L,num.L)(j1,j1); 
							}
					}
				}
			
				//dgamma.aXud_central(i) = Xud*Bubble_at_freq*Xud;
				dgamma.aXud_central(i) = Bubble_at_freq;
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
	//		 			dgamma.aXud_central(i)(j1,j2) += (complex<double>)( dgamma.aXud_feedback_data(num.L,num.L)(j1,j2) 
	//										                               -static_subtraction_ud(j1,j2)); 
					}
				}
			}
			
			{
			 	syma<complex<double> > Dud_1(num.Nges);
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
					 	        Dud_1(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j1,1,j2,0)
						                       +gamma.aDud_central(i)(j1,j2);
							if(j1==j2){
							 	Dud_1(j1,j2) += gamma.aPud_feedback_data(num.L,num.L)(j1,j1) 
								               +gamma.aXud_feedback_data(num.L,num.L)(j1,j1); 
							}
					}
				}
			 	syma<complex<double> > Dud_2(num.Nges);
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
					 	        Dud_2(j1,j2)  = 0.5*barevertex(j1,0,j2,0,j1,0,j2,0)
						                       +gamma.aDdd_central(i)(j1,j2);
							if(j1==j2){
							 	Dud_2(j1,j2) += gamma.aPdd_feedback_data(num.L,num.L)(j1,j1) 
								               -gamma.aDdd_feedback_data(num.L,num.L)(j1,j1); 
							}
					}
				}
			
				dgamma.aDud_central(i) = -Dud_1*Bubble_at_freq.conj()*Dud_2;
				dgamma.aDuu_central(i) = -Dud_2*Bubble_at_freq.conj()*Dud_2 - Dud_1*Bubble_at_freq.conj()*Dud_1;
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
			 			//dgamma.aDud_central(i)(j1,j2) += (complex<double>)( dgamma.aDud_feedback_data(num.L,num.L)(j1,j2) 
						//					                               -static_subtraction_ud_2(j1,j2)); 
			 			//dgamma.aDuu_central(i)(j1,j2) += (complex<double>)( dgamma.aDuu_feedback_data(num.L,num.L)(j1,j2) 
						//					                               -static_subtraction_uu(j1,j2)); 
					}
				}
				dgamma.aDdd_central(i) = dgamma.aDuu_central(i);
			}



		}


	}
	time(&t2);
	cout<<"Time for dynamic XD_flow ="<<t2 - t1<<endl;
//	
//	/*Static Contribution to selfenergy:*/
//	Self_energy_static_zero_mag<mode> self_stat(phy,num,pre,sub,Lambda,measure_flow,gamma,barevertex);
//	cout<<"measure_flow="<<measure_flow<<endl;
//	time(&t1);
//	syma<complex<double> > Self_stat = self_stat();
//	time(&t2);
//	cout<<"Time for static Self_energy part ="<<t2 - t1<<endl;
//
//	/*Dynamic Contribution to selfenergy:*/
//	time(&t1);
//	Self_energy_dynamic_zero_mag<mode> self_dyn(phy,num,pre,sub,Lambda,measure_flow,gamma);
//	omp_set_num_threads(16);
//	#pragma omp parallel for
//	for(int i=0; i<num.Nff; ++i){
//	 	dgamma.ERetu(i) =  self_dyn(num.wf(i))+ Self_stat;
//	 	dgamma.ERetd(i) =  dgamma.ERetu(i);
//	}
//	time(&t2);
	cout<<"Time for dynamic Self_energy part ="<<t2 - t1<<endl;
	 	
}
	







#endif
