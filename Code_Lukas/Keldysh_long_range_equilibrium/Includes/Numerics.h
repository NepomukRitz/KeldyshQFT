#ifndef NUMERICS_14032017
#define NUMERICS_14032017

#include <iostream> 
#include <stdio.h>
#include <string>
#include <cstring>

#include "matrix.h" 
#include "approxxpp.h"
#include "Physics.h"
#include "Norm.h"
#include "Ex_functions.h"

using namespace std;

class Numerics{
	public:
		int L;
		int N;
		int Nges;
		int twoN; //Teste diese Variable noch und fuege sie in compare hinzu!
		int Nff;
		int NfbP;
		int NfbX;
		int Nff_input;
		int NfbP_input;
		int NfbX_input;
		int pos_Nff_mu;
		int pos_NfbX_0;
		int pos_NfbP_2mu;
		matrix<double> wf;
		matrix<double> wbP;
		matrix<double> wbX;
		int num_freq_pre;
		Physics phy;
		matrix<int> Lp_structure;
		matrix<int> Lx_structure;
		matrix<matrix<int> > Lp_bounds;
		matrix<matrix<int> > Lx_bounds;
		
		Numerics(){};
		Numerics(int L_in, int N_in, int Nff_in, int NfbP_in, int NfbX_in, int num_freq_pre_in, Physics &phy_in);
		Numerics(int L_in, int N_in, matrix<double> & wf_in, matrix<double> &wbP_in, matrix<double> wbX_in, int num_freq_pre_in, Physics &phy_in);
		
		Numerics(int L_in, int N_in, int Nff_in, int num_freq_pre_in, Physics &phy_in); //Old frequency discretization from the Keldysh short range programs 
		
		void save(char *filename);
		void save(string filename);
		
		double compare(Numerics num);
		int binary_search_for_interpolation(double x, int imin, int imax, matrix<double> &xi);
		int interpol_wf(double freq);
		int interpol_P(double freq);
		int interpol_X(double freq);
		void initialize_long_str();
		void initialize_long_str(int N_L_full_around_feedback);
		void initialize_long_range_bounds();
};


Numerics::Numerics(int L_in, int N_in, int Nff_in, int NfbP_in, int NfbX_in, int num_freq_pre_in, Physics &phy_in): L(L_in), N(N_in), Nges(2*N+1), twoN(2*N), Nff(Nff_in), NfbP(NfbP_in), NfbX(NfbX_in), num_freq_pre(num_freq_pre_in), phy(phy_in), Nff_input(Nff_in), NfbP_input(NfbP_in), NfbX_input(NfbX_in){
 
	//Zusaetzliche Frequenzen
	int N_exp; 
	int N_Tsym=50;
	#if(ADD_FINITE_TEMP_FREQ_IN_S==1)
		if(phy.T==0.0){
			N_exp = 13;
		}
		else{
			N_exp = 13+2*N_Tsym+1; 
		}
	#else
		N_exp = 13;
	#endif


	//Setze die Nff Frequenzen
	int Nff_sec=max( (Nff-N_exp)/6, 1);
	int Nff_dist= 6*Nff_sec;
	int Nff_real= 6*Nff_sec + N_exp;

	wf.resize(Nff_real);
	wf(Nff_real-1) = 1e06;
	wf(Nff_real-2) =-1e06;
	wf(Nff_real-3) = .0;
	wf(Nff_real-4) = 2.;
	wf(Nff_real-5) =-2.;
	wf(Nff_real-6) = 4.;
	wf(Nff_real-7) =-4.;
	wf(Nff_real-8) = 6.;
	wf(Nff_real-9) =-6.;
	wf(Nff_real-10) = phy.mu;
	wf(Nff_real-11) = -phy.mu;
	wf(Nff_real-12) = 2.* phy.mu;
	wf(Nff_real-13) = -2.* phy.mu;
	#if(ADD_FINITE_TEMP_FREQ_IN_S==1)
		if(phy.T!=0){
			for(int i=-N_Tsym, z=0; i<=N_Tsym; ++i, ++z){
				wf(Nff_real-14-z)= phy.mu + double(i)/N_Tsym*5.*phy.T;
			}
		}
	#endif
	
	for(int i=0; i<Nff_sec; ++i){
		double x=(double)(i)/Nff_sec;
		wf(i)= -4 - pow(1e6,x)+ 1.; //hier wurde frueher e^12 statt 1e6 benutzt. 
		wf(5*Nff_sec+i) = 4 + pow(1e6,x) - 1.;
	}
	
	for(int i=0; i< 4*Nff_sec; ++i){
		double x=(double)(i)/(4*Nff_sec);
		wf(Nff_sec+i)=(-4.)*(1-x) + 4.*x;
	}
	wf.sort();
	Nff=Nff_real;
	
	//Werfe doppelte Frequenzen weg:
	anfang_Nff:
	{
	 	for(int i=0; i<Nff-1; ++i){
		    if(wf(i) ==wf(i+1)){
		   		Nff = Nff-1;
		   		matrix<double> wf_tmp(Nff);	
		   		for(int j=0; j<i; ++j){
		   			wf_tmp(j) = wf(j);
		   		}
		   		for(int j=i;j<Nff;++j){
		   		 	wf_tmp(j) = wf(j+1);
		   		}
		   		wf = wf_tmp;
		   		goto anfang_Nff;
		   	}
	       }
	}
	       		 

	 
	#if(ADD_FINITE_TEMP_FREQ_IN_PX==1)
		if(phy.T==0.0){
			N_exp = 13;
		}
		else{
			N_exp = 13+2*N_Tsym+1; 
		}
	#else
		N_exp = 13;
	#endif

	//Setze die NfbP Frequenzen
	int NfbP_sec=max( (NfbP-N_exp)/6, 1);
	int NfbP_dist= 6*NfbP_sec;
	int NfbP_real= 6*NfbP_sec + N_exp;

	wbP.resize(NfbP_real);
	wbP(NfbP_real-1) = 1e06;
	wbP(NfbP_real-2) =-1e06;
	wbP(NfbP_real-3) = .0;
	wbP(NfbP_real-4) = 2.;
	wbP(NfbP_real-5) =-2.;
	wbP(NfbP_real-6) = 4.;
	wbP(NfbP_real-7) =-4.;
	wbP(NfbP_real-8) = 6.;
	wbP(NfbP_real-9) =-6.;
	wbP(NfbP_real-10) = phy.mu;
	wbP(NfbP_real-11) = -phy.mu;
	wbP(NfbP_real-12) = 2.* phy.mu;
	wbP(NfbP_real-13) = -2.* phy.mu;
	#if(ADD_FINITE_TEMP_FREQ_IN_PX==1)
		if(phy.T!=0){
			for(int i=-N_Tsym, z=0; i<=N_Tsym; ++i, ++z){
				wbP(NfbP_real-14-z)= 2.*phy.mu + double(i)/N_Tsym*5.*phy.T;
			}
		}
	#endif
	
	for(int i=0; i<NfbP_sec; ++i){
		double x=(double)(i)/NfbP_sec;
		wbP(i)= -4 - pow(1e6,x) + 1.; 
		wbP(5*NfbP_sec+i) = 4 + pow(1e6,x) - 1.;
	}
	
	for(int i=0; i< 4*NfbP_sec; ++i){
		double x=(double)(i)/(4*NfbP_sec);
		wbP(NfbP_sec+i)=(-4.)*(1-x) + 4.*x;
	}
	wbP.sort();
	NfbP=NfbP_real;
	
	//Werfe doppelte Frequenzen weg:
	anfang_NfbP:
	{
	 	for(int i=0; i<NfbP-1; ++i){
		    if(wbP(i) ==wbP(i+1)){
		   		NfbP = NfbP-1;
		   		matrix<double> wbP_tmp(NfbP);	
		   		for(int j=0; j<i; ++j){
		   			wbP_tmp(j) = wbP(j);
		   		}
		   		for(int j=i;j<NfbP;++j){
		   		 	wbP_tmp(j) = wbP(j+1);
		   		}
		   		wbP = wbP_tmp;
		   		goto anfang_NfbP;
		   	}
	       }
	}

	//Setze die NfbX Frequenzen
	int NfbX_sec=max( (NfbX-N_exp)/6, 1);
	int NfbX_dist= 6*NfbX_sec;
	int NfbX_real= 6*NfbX_sec + N_exp;

	wbX.resize(NfbX_real);
	wbX(NfbX_real-1) = 1e06;
	wbX(NfbX_real-2) =-1e06;
	wbX(NfbX_real-3) = .0;
	wbX(NfbX_real-4) = 2.;
	wbX(NfbX_real-5) =-2.;
	wbX(NfbX_real-6) = 4.;
	wbX(NfbX_real-7) =-4.;
	wbX(NfbX_real-8) = 6.;
	wbX(NfbX_real-9) =-6.;
	wbX(NfbX_real-10) = phy.mu;
	wbX(NfbX_real-11) = -phy.mu;
	wbX(NfbX_real-12) = 2.* phy.mu;
	wbX(NfbX_real-13) = -2.* phy.mu;
	#if(ADD_FINITE_TEMP_FREQ_IN_PX==1)
		if(phy.T!=0){
			for(int i=-N_Tsym, z=0; i<=N_Tsym; ++i, ++z){
				wbX(NfbX_real-14-z)= double(i)/N_Tsym*5.*phy.T;
			}
		}
	#endif
	
	for(int i=0; i<NfbX_sec; ++i){
		double x=(double)(i)/NfbX_sec;
		wbX(i)= -4 -  pow(1e6,x) + 1.; 
		wbX(5*NfbX_sec+i) = 4 + pow(1e6,x) - 1.;
	}
	
	for(int i=0; i< 4*NfbX_sec; ++i){
		double x=(double)(i)/(4*NfbX_sec);
		wbX(NfbX_sec+i)=(-4.)*(1-x) + 4.*x;
	}
	wbX.sort();
	NfbX=NfbX_real;
	
	//Werfe doppelte Frequenzen weg:
	anfang_NfbX:
	{
	 	for(int i=0; i<NfbX-1; ++i){
		    if(wbX(i) ==wbX(i+1)){
		   		NfbX = NfbX-1;
		   		matrix<double> wbX_tmp(NfbX);	
		   		for(int j=0; j<i; ++j){
		   			wbX_tmp(j) = wbX(j);
		   		}
		   		for(int j=i;j<NfbX;++j){
		   		 	wbX_tmp(j) = wbX(j+1);
		   		}
		   		wbX = wbX_tmp;
		   		goto anfang_NfbX;
		   	}
	       }
	}

	//Setze Feedback Positionen
	for(int i=0;i<Nff; ++i){
		if(wf(i) == phy.mu){
			pos_Nff_mu=i;
			break;
		}
	}
	for(int i=0;i<NfbP; ++i){
		if(wbP(i) == 2.*phy.mu){
			pos_NfbP_2mu=i;
			break;
		}
	}
	for(int i=0;i<NfbX; ++i){
		if(wbX(i) == 0.0){
			pos_NfbX_0=i;
			break;
		}
	}
	//Symmetrisiere explizit in X:
	{
		NfbX = pos_NfbX_0*2+1;
		matrix<double> wbX_sym(NfbX); 
		for(int i=0; i<=pos_NfbX_0; ++i){
			wbX_sym(i)=wbX(i);
			wbX_sym(NfbX-1-i) = -wbX(i);
		}
		wbX = wbX_sym;
	}
	//Set Lp,Lx structure to zero by default:
	initialize_long_str();	
}


Numerics::Numerics(int L_in, int N_in, matrix<double> & wf_in, matrix<double> &wbP_in, matrix<double> wbX_in, int num_freq_pre_in, Physics &phy_in): L(L_in), N(N_in), Nges(2*N+1), twoN(2*N), Nff(wf_in.dim_c), NfbP(wbP_in.dim_c), NfbX(wbX_in.dim_c), wf(wf_in), wbP(wbP_in), wbX(wbX_in), num_freq_pre(num_freq_pre_in), phy(phy_in){
 //Setze Feedback Positionen
	pos_Nff_mu=-1;
	pos_NfbP_2mu=-1;
	pos_NfbX_0=-1;
	for(int i=0;i<Nff; ++i){
		if(wf(i) == phy.mu){
			pos_Nff_mu=i;
			break;
		}
	}
	for(int i=0;i<NfbP; ++i){
		if(wbP(i) == 2.*phy.mu){
			pos_NfbP_2mu=i;
			break;
		}
	}
	for(int i=0;i<NfbP; ++i){
		if(wbP(i) == 0.0){
			pos_NfbX_0=i;
			break;
		}
	}


}


void Numerics::save(char *filename){
	matrix<double> tmp(1);
	tmp(0)=L;
	tmp.save(filename,"L");
	tmp(0)=N;
	tmp.save(filename,"N");
	tmp(0)=Nges;
	tmp.save(filename,"Nges");
	tmp(0)=Nff;
	tmp.save(filename,"Nff");
	tmp(0)=NfbP;
	tmp.save(filename,"NfbP");
	tmp(0)=NfbX;
	tmp.save(filename,"NfbX");
	tmp(0)=pos_Nff_mu;
	tmp.save(filename,"pos_Nff_mu");
	tmp(0)=pos_NfbP_2mu;
	tmp.save(filename,"pos_NfbP_2mu");
	tmp(0)=pos_NfbX_0;
	tmp.save(filename,"pos_NfbX_0");
	wf.save(filename,"wf");
	wbP.save(filename,"wbP");
	wbX.save(filename,"wbX");
	tmp(0)=num_freq_pre;
	tmp.save(filename,"num_freq_pre");
	tmp(0)=Nff_input;
	tmp.save(filename,"Nff_input");
	tmp(0)=NfbP_input;
	tmp.save(filename,"NfbP_input");
	tmp(0)=NfbX_input;
	tmp.save(filename,"NfbX_input");
	Lp_structure.save(filename,"Lp_structure");
	Lx_structure.save(filename,"Lx_structure");
	Lp_bounds.save(filename,"Lp_bounds");
	Lx_bounds.save(filename,"Lx_bounds");

}

void Numerics::save(string filename){
	char * cstr = new char [filename.length()+1];
  	std::strcpy (cstr, filename.c_str());
	save(cstr);
}








































Numerics::Numerics(int L_in, int N_in, int Nff_in, int num_freq_pre_in, Physics &phy_in): L(L_in), N(N_in), Nff(Nff_in), NfbP(Nff), NfbX(Nff), Nges(2*N+1), twoN(2*N), phy(phy_in), wf(Nff), wbP(Nff), wbX(Nff), num_freq_pre(num_freq_pre_in){ 

 //Setze Positionen
 pos_Nff_mu=0;
 pos_NfbX_0=0;
 pos_NfbP_2mu=0;
 int Nfb=Nff;
 double taul=1.0;
 double mu=phy.mu;
 double T=phy.T;
 double Vg=phy.Vg;
 double h=phy.h;

 //Frequenzdiskretisierung
 
 //int Nfb_reduced = Nfb-10;   //Frequencies are guaranteed to lie at 0, 2taul, -2taul, 4taul, -4taul, 6taul, -6taul, 1e06(infinity), -1e06(-infinity), mu
 int Nfb_reduced = Nfb-11;   //Frequencies are guaranteed to lie at 0, 2taul, -2taul, 4taul, -4taul, 6taul, -6taul, 1e06(infinity), -1e06(-infinity), mu, 2mu
 int cut = Nfb_reduced/6;
 int cut2= 5*cut;
 double offset = -4.*taul;
 double offset2=  4.*taul;
 
 for (int i = 0; i < Nfb_reduced; i++) {
  if (i<cut)
   wbP(i) = offset-exp(12.*(double)i/(double)cut)+1.;
  else if (i<cut2)
   wbP(i) = (offset*(double)(i-cut2+1)-offset2*(double)(i-cut))/(double)(cut-cut2+1);
  else
   wbP(i) = offset2+exp(12.*(double)(i-cut2)/(double)(Nfb_reduced-cut2))-1.;
 }
 wbP(Nfb-1) = 1e06;
 wbP(Nfb-2) =-1e06;
 wbP(Nfb-3) = .0;
 wbP(Nfb-4) = 2.*taul;
 wbP(Nfb-5) =-2.*taul;
 wbP(Nfb-6) = 4.*taul;
 wbP(Nfb-7) =-4.*taul;
 wbP(Nfb-8) = 6.*taul;
 wbP(Nfb-9) =-6.*taul;
 wbP(Nfb-10) = mu;
 wbP(Nfb-11) = 2.* mu;
 wbP.sort();
 matrix<double> helper(Nfb); //Schmeisse doppelte Frequenzen weg
 for (int i=0; i<Nfb-1; i++){
  if (wbP(i)==wbP(i+1)) {
   helper.resize(wbP.dim_c-1);
   for (int j=0; j<i; j++) {
    helper(j) = wbP(j);
   }
   for (int j=i; j<helper.dim_c; j++) {
    helper(j) = wbP(j+1);
   }
  }
 }
 matrix<double> helper2(helper.dim_c);
 int add=100;
 if (T!=0) {
  helper2.resize(helper.dim_c+add);
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c-1+i) = mu+(double)(i-add/2)/(double)add*2*T;
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c+add/2-1+i) = 2.*mu+(double)(i-add/2)/(double)add*2*T;
 }
 else {
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
 }
 
 
 Nff = helper2.dim_c;
 Nfb = helper2.dim_c;
 wbP.resize(helper2.dim_c);
 wbX.resize(helper2.dim_c);
 wf.resize(helper2.dim_c);
 wbP = helper2;
 
 wbP.sort();
 
 wbX = wbP;
 wf = wbP;

}
 
double Numerics::compare(Numerics num){
 double diff=0.0;
 diff=max(diff,abs((double)(L -num.L)));
 diff=max(diff,abs((double)(N -num.N))); 
 diff=max(diff,abs((double)(Nges -num.Nges)));
 diff=max(diff,abs((double)(Nff -num.Nff))); 
 diff=max(diff,abs((double)(NfbP -num.NfbP)));
 diff=max(diff,abs((double)(NfbX -num.NfbX))); 
 diff=max(diff,abs((double)(pos_Nff_mu -num.pos_Nff_mu)));
 diff=max(diff,abs((double)(pos_NfbX_0 -num.pos_NfbX_0))); 
 diff=max(diff,abs((double)(pos_NfbP_2mu -num.pos_NfbP_2mu)));
 if( (wf.dim_r != num.wf.dim_r) || (wf.dim_c != num.wf.dim_c) ){
  diff=1e16;
 }
 else{
  diff=max(diff,maximumsnorm(wf - num.wf));
 }
 if( (wbP.dim_r != num.wbP.dim_r) || (wbP.dim_c != num.wbP.dim_c) ){
  diff=1e16;
 }
 else{
  diff=max(diff,maximumsnorm(wbP - num.wbP));
 }
 if( (wbX.dim_r != num.wbX.dim_r) || (wbX.dim_c != num.wbX.dim_c) ){
  diff=1e16;
 }
 else{
  diff=max(diff,maximumsnorm(wbX - num.wbX));
 }
 diff=max(diff,abs((double)(num_freq_pre - num.num_freq_pre)));
 diff=max(diff,phy.compare(num.phy));

 return diff;
}
 
int Numerics::binary_search_for_interpolation(double x, int imin, int imax, matrix<double> &xi){
 if ((imax-imin)==1 || (imax-imin)==0)
  return imin;
 if (imax < imin){
  myerror("this should not happen! Search did not converge",__FILE__,__LINE__);
  return 999999999;
 }
 else {
  int imid = imin+(imax-imin)/2;
  if (x>xi(imid))
   return (binary_search_for_interpolation(x, imid, imax, xi));
  else if (x<xi(imid))
   return (binary_search_for_interpolation(x, imin, imid, xi));
  else
   return imid;
 }
};

int Numerics::interpol_wf(double freq){
	return binary_search_for_interpolation(freq,0,Nff-1,wf)+1;
};  
int Numerics::interpol_P(double freq){
	return binary_search_for_interpolation(freq,0,NfbP-1,wbP)+1;
};  
int Numerics::interpol_X(double freq){
	return binary_search_for_interpolation(freq,0,NfbX-1,wbX)+1;
};  
 
void Numerics::initialize_long_str(){
	Lp_structure.resize(NfbP); 
	Lx_structure.resize(NfbX); 
	Lp_structure = 0;
	Lx_structure = 0;
	initialize_long_range_bounds();	
}

void Numerics::initialize_long_str(int N_L_full_around_feedback){
	Lp_structure.resize(NfbP); 
	Lx_structure.resize(NfbX); 
	Lp_structure = 0;
	Lx_structure = 0;
	for(int i=1; i<=N_L_full_around_feedback; ++i){
		if(in_range(pos_NfbP_2mu+i,NfbP)) Lp_structure(pos_NfbP_2mu+i) = L;
		if(in_range(pos_NfbP_2mu-i,NfbP)) Lp_structure(pos_NfbP_2mu-i) = L;
		if(in_range(pos_NfbX_0+i,NfbX)) Lx_structure(pos_NfbX_0+i) = L;
		if(in_range(pos_NfbX_0-i,NfbX)) Lx_structure(pos_NfbX_0-i) = L;
	}
	initialize_long_range_bounds();	
}

void Numerics::initialize_long_range_bounds(){
	Lp_bounds = init_long_range_bounds(L,Lp_structure,pos_NfbP_2mu);
	Lx_bounds = init_long_range_bounds(L,Lx_structure,pos_NfbX_0);
}


#endif
