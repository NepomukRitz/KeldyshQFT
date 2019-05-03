#ifndef P_FLOW_ZERO_MAG_DEBUG_23082017
#define P_FLOW_ZERO_MAG_DEBUG_23082017

#include "Vertex.h"
#include "P_bubble_central_zero_mag.h"
#include "P_bubble_feedback_zero_mag.h"

//Get rid of the dynamic aPuu, aPdd!

template<int mode> class P_flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_flow_zero_mag(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
		              double measure_flow,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma,
		              Vertex<mode> &dgamma
					  );
};


template <int mode> P_flow_zero_mag<mode>::P_flow_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void P_flow_zero_mag<mode>::operator()(double Lambda,
                                                      double measure_flow,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma,
													  Vertex<mode> &dgamma){
 	#if RPA_MODE==1
	cout<<"Achtung: RPA Mode an "<<endl;
	#endif
 	time_t t1, t2;
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


	time(&t1);
	matrix<double> static_addition_ud(num.Nges,num.Nges); //Not elegant yet
	matrix<double> dynamic_addition_half_1(num.Nges,num.Nges); //Not elegant yet
	matrix<double> dynamic_addition_half_2(num.Nges,num.Nges); //Not elegant yet
	static_addition_ud = 0.0;
	dynamic_addition_half_1 = 0.0;
	dynamic_addition_half_2 = 0.0;


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
						#if RPA_MODE==0
						if(gamma.aXud_feedback.inrange(i+k-j,j+l-i,j,i)){
						 	P(l,k,j,i) += gamma.aXud_feedback(i+k-j,j+l-i,j,i);
						}
						if(gamma.aDud_feedback.inrange(i-j,j+l-i-k,j,i+k)){
						 	P(l,k,j,i) += gamma.aDud_feedback(i-j,j+l-i-k,j,i+k);
						}
						#endif
						 	
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
					 	dgamma.aPud_feedback_data(l+num.L, k+num.L) += P_data(l+num.L,q+num.L)*(Bubble_data(q+num.L,p+num.L)+Bubble_data(-q+num.L,-p+num.L))*P_data(p+num.L,k+num.L);	
					}
				}
			}
		}
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
			 	if( (k!=0) && (l!=0)){
				 	static_addition_ud += P_data(num.L,l+num.L)*(Bubble_data(l+num.L,k+num.L)+Bubble_data(-l+num.L,-k+num.L))*P_data(k+num.L,num.L);  
				}
			}
		}
		
		for(int l=-num.L; l<=num.L; ++l){
		 	if( l!=0){
			 	dynamic_addition_half_1 += P_data(num.L,l+num.L)*(Bubble_data(l+num.L,num.L)+Bubble_data(-l+num.L,num.L));
			 	dynamic_addition_half_2 += (Bubble_data(num.L,l+num.L)+Bubble_data(num.L,-l+num.L))*P_data(l+num.L,num.L);
			}
		}
		

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
				 	        Pud(j1,j2)  = 0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
					                     +gamma.aPud_central(i)(j1,j2);
						//if(j1==j2){
						// 	Pud(j1,j2)  += gamma.aXud_feedback_data(num.L,num.L)(j1,j2) 
						//	              +gamma.aDud_feedback_data(num.L,num.L)(j1,j2); 
						//				  /*Hier fehlt nicht das Feedback in Puu, da es sich cancelt!*/
						//}
						#if RPA_MODE==0
						if(gamma.aXud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
						 	Pud(j1,j2) += gamma.aXud_feedback(j2-j1,j1-j2,j1,j2);
						}
						if(gamma.aDud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
						 	Pud(j1,j2) += gamma.aDud_feedback(j2-j1,j1-j2,j1,j2);
						}
						#endif
				}
			}
		
			syma<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Bubble(num.wbP(i));
			dgamma.aPud_central(i) = Pud*Bubble_at_freq*Pud;
			matrix<complex<double> > tmp = dynamic_addition_half_1*Pud + Pud*dynamic_addition_half_2;

			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
		 			dgamma.aPud_central(i)(j1,j2) += (complex<double>)( static_addition_ud(j1,j2) + tmp(j1,j2)
																	  ); 
				}
			}
		}





	}
	time(&t2);
	cout<<"Time for dynamic P_flow="<<t2 - t1<<endl;



}


#endif
