#ifndef LONG_RANGE_04112016
#define LONG_RANGE_04112016
 
#include <iostream> 
#include "basic.h"
#include <stdio.h>
#include <string.h>
#include "matrix.h" 
#include "approxxpp.h"

using namespace std;

template <class T> double maximumsnorm(matrix<T> A);

template <class T> double maximumsnorm(syma<T> A);

double parity_check(syma<complex<double> > A);

double parity_check(matrix<double> A);

double parity_check(matrix<complex<double> > A);

double parity_check(matrix<matrix<double> > A);

double parity_check(matrix<matrix<complex<double> > > A);

double symmetry_check(matrix<double> A);



//Aus Barevertex.h:


class Barevertex_new{
 public:
 matrix<double> U;
 int Lu;
 Barevertex_new(int Lu_in, matrix<double>  U_in){
  if(U_in.dim_c!=U_in.dim_r){
   cout<<"Fehler Barevertex: Wechselwirkung keine Quadratmatrix"<<endl;
   throw 1;
  }
  if((U_in.dim_c-1)%2!=0){
   cout<<"Fehler Barevertex: Wechselwirkung hat gerade Anzahl an Sites"<<endl;
   throw 1;
  }
  Lu=Lu_in;
  U=U_in;
 }
 double operator()(int j1, bool s1, int j2, bool s2, int j3, bool s3, int j4, bool
 s4){
   double value=0.;
   if(j1==j2 && j2==j3 && j3==j4 && s1!=s2 && s3!=s4){
    if(s1==s3){
       value=U(j1,j1);
    }
    else{
       value=-U(j1,j1);
    }
   }
   if(j1!=j2 && (abs(j1-j2)<=Lu)){
    if(j1==j3 && j2==j4 && s1==s3 && s2==s4){
 	  value=U(j1,j2);
    }
    if(j1==j4 && j2==j3 && s1==s4 && s2==s3){
     value=-U(j1,j2);
    }
   }
   return value;
 }
 void save(char *filename){
  matrix<double> sv(1);
  sv(0)=Lu;
  sv.save(filename,"Lu");
  U.save(filename,"U");
 }
};


//hamilton_zero von Dennis

inline syma<complex<double> > hamilton_zero(int N,short pot_type,double Vg,double Vsg, int pot_width, double taul){
// if (pot_type==0) {
  if (N%2==0) {
   syma<complex<double> > H(N);
   H = (complex<double>) .0;
   double x = .0;
   double v = .0;
   double tau = .0;
   for (int j=-(N/2-1); j<N/2; j++) {
    x = (double) j/(double)(N/2-1);
    v = Vg*exp(-x*x/(1.-x*x));
    tau = -1. + v;
    int i = j+N/2-1;
    H(i+1, i) = tau;
   }
   return H;
  }
  else {
   syma<complex<double> > H0(N);
   H0 = (complex<double>) .0;
   double N2 = (double)(N-1)/2.;
   double v=.0;
   int i=0; 
   for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
    v=-1.+(2.*Vg)/2.*exp(-j*j/(1.-j*j));
    H0(i+1, i)=v;
    i++;
   }
   return H0;
  }
//  }
//	//TODO: This needs to be debugged/optimized (insert good stopping points for interior integrators).
//	 else if (pot_type==1)
//	  return(QD_hamiltonian(N, pot_width,Vsg,Vg));
//	 else if (pot_type==2)
//	  return(Wire_hamiltonian(N, pot_width,Vg));
//	 else if (pot_type==3) {
//	  if (N%2==0) {
//	   syma<complex<double> > H(N);
//	   H = (complex<double>) .0;
//	   double x = .0;
//	   double v = .0;
//	   double tau = .0;
//	   for (int j=-(N/2-1); j<N/2; j++) {
//	    int i = j+N/2-1;
//	    if((j<-pot_width/2. && abs(j+pot_width/2.)<2.3) || (j>pot_width/2. && abs(j-pot_width/2.)<2.3))
//	     tau = -1.+Vg;
//	    else
//	     tau = -1.;
//	    H(i+1, i) = tau;
//	   }
//	   return H;
//	  }
//	 }
//	 else if (pot_type==4){
//	  syma<complex<double> > H0(N);
//	  H0 = (complex<double>) .0;
//	  double N2 = (double)(N-1)/2.;
//	  double v=.0;
//	//   int i=0; 
//	//   for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
//	//    int i = j+N/2-1;
//	//    if(j<-pot_width/2. || j>pot_width/2.)
//	//     v = -1.+Vg;
//	//    else
//	//     v = -1.;
//	//    H0(i+1, i)=v;
//	//    i++;
//	//   }
//	  int bar=N-pot_width;
//	  for (int i=0; i<N-1; i++) {
//	   H0(i+1, i)=-1.;
//	  }
//	  for (int i=(int)((double)bar/2.)-2; i<(double)bar/2.; i++) {
//	   H0(i+1, i)=-1.+Vg;
//	   H0(N-i-1, N-i-2)=-1.+Vg;
//	  }
//	  return H0;
//	 }
//	 else if (pot_type==5) {
//	  syma<complex<double> > H(N);
//	  H=(complex<double>)0.;
//	  for (int i=0;i<N;i++)
//	   H(i,i)=Vg/(1.+2.*.55)*exp(-pow(((i-.5*N+.5)/(N-1.)*2.),2)/(1.-pow(((i-.5*N+.5)/(N-1)*2.),2)));
//	 
//	  for (int i=0;i<N-1;i++)
//	   H(i+1,i)=-taul+.5*.55*(H(i,i)+H(i+1,i+1));
//	  return H;
//	 }
}


// declaration Blockmatrix

int komp(int L, int N, int l, int j);

int dim_all(int L, int N);


template <class T> void resize_block(int L, int N, matrix<matrix<T> > &A);

template <class T> void initialize_block(int L, int N, matrix<matrix<T> > &A, T x);

bool inrange_block(int L, int N, int l, int k, int j, int i);
bool inrange_block_phy(int L, int N, int l, int k, int j, int i);

template <class T> class Blockmatrix{
 public:
 int L;
 int N;
 matrix<matrix<T> > &A;
 Blockmatrix(int L_in, int N_in, matrix<matrix<T> > &A_in);
 T & operator()(int l, int k, int j, int i);
 T & fastaccess(int ltilde, int ktilde, int jtilde, int itilde);
 T & phyaccess(int l, int k, int j_asymmetric, int i_asymmetric);
 matrix<T> convert_to_matrix();
 void resize(int L, int N);
 void initialize(int L, int N, T x);
 void save(char *filename, char *variable);
 bool inrange(int l, int k, int j, int i);
};

template <class T> matrix<matrix<T> > convert_to_blockmatrix(matrix<T> A, int L, int N);

template <int mode> class Substitution{
 public:
 double taul;
 double Lambda;
 double h;
 Substitution(double taul_in, double Lambda_in, double h_in): taul(taul_in), Lambda(Lambda_in), h(h_in){};
 double subst_concatenated (double omega);
 double resu_concatenated (double y);
 double weight_concatenated (double y);
};

class Substitution_flow{
 public:
 double subst_flow(double x){
  return exp(x)/(1.-exp(x));
 }
 double resu_flow(double x){
  return log(x/(1.+x));
 }
 double weight_flow(double x){
  return exp(x)/(1.-exp(x))/(1.-exp(x));
 }
};


class Numerics{
 public:
 const int number_of_frequencies_at_which_g_and_s_are_precomputed;
 const double delta_around_dangerous_frequencies;
 const double delta_around_dangerous_frequencies_at_pm7;
 int L;
 int N;
 int Nges;
 int Nff;
 int Nfb;
 int pos_Nff_0;
 int pos_Nfb_0;
 int pos_Nfb_2mu;
 matrix<double> &wf;
 matrix<double> &wbP;
 matrix<double> &wbX;
 matrix<double> frequencies_of_precomputation;
 double Vg;
 double mu;
 double T;
 Numerics(const int number_of_frequencies_at_which_g_and_s_are_precomputed_in, const double delta_around_dangerous_frequencies_in,const double delta_around_dangerous_frequencies_at_pm7_in, int L_in, int N_in, int Nff_in, int Nfb_in, matrix<double> & wf_in, matrix<double> & wbP_in, matrix<double> & wbX_in, double Vg_in, double mu_in, double T_in);

 Numerics(const int number_of_frequencies_at_which_g_and_s_are_precomputed_in, const double delta_around_dangerous_frequencies_in, const double delta_around_dangerous_frequencies_at_pm7_in, int L_in, int N_in, matrix<double> & wf_in, double Vg_in, double mu_in, double T_in):
 number_of_frequencies_at_which_g_and_s_are_precomputed(number_of_frequencies_at_which_g_and_s_are_precomputed_in),
 delta_around_dangerous_frequencies(delta_around_dangerous_frequencies_in),
 delta_around_dangerous_frequencies_at_pm7(delta_around_dangerous_frequencies_at_pm7_in),
 L(L_in),
 N(N_in),
 Nges(2*N+1),
 Nff(wf_in.dim_c),
 Nfb(Nff),
 pos_Nff_0(0),
 pos_Nfb_0(0),
 pos_Nfb_2mu(0),
 wf(wf_in),
 wbP(wf),
 wbX(wf),
 Vg(Vg_in),
 mu(mu_in),
 T(T_in) {

//Frequencies of precomputation
double taul=1.0;
double Lambda_fake=0.0; // Solange mu im Band liegt ist dieser Wert egal.
frequencies_of_precomputation.resize(number_of_frequencies_at_which_g_and_s_are_precomputed+22);

for (int i=0; i<number_of_frequencies_at_which_g_and_s_are_precomputed; i++)
 frequencies_of_precomputation(i) = -7.+14.*(double)(i+1)/(double)(number_of_frequencies_at_which_g_and_s_are_precomputed+1);
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-22) = -7.+delta_around_dangerous_frequencies_at_pm7;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-21) =  7.-delta_around_dangerous_frequencies_at_pm7;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-20) = -6.*taul+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-19) = -6.*taul-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-18) = -6.*taul+6.*Vg+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-17) = -6.*taul+6.*Vg-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-16) = -2.*taul+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-15) = -2.*taul-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-14) = -2.*taul+2.*Vg+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-13) = -2.*taul+2.*Vg-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-12) =  6.*taul+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-11) =  6.*taul-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-10) =  6.*taul-6.*Vg+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-9)  =  6.*taul-6.*Vg-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-8)  =  2.*taul+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-7)  =  2.*taul-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-6)  =  2.*taul-2.*Vg+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-5)  =  2.*taul-2.*Vg-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-4)  =  subst_concatenated(mu, taul, Lambda_fake)+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-3)  =  subst_concatenated(mu, taul, Lambda_fake)-delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-2)  =  .0+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-1)  =  .0-delta_around_dangerous_frequencies;
frequencies_of_precomputation.sort();
}

 void save(char *filename);
 };

//Generalmatrix for short and long structure of vertices

class Generalmatrix{
public:
int L;
int N;
int wf;
matrix < matrix < syma < complex < double > > > > short_str;
matrix < matrix < matrix < double > > > long_str;
Generalmatrix(){};
Generalmatrix(Numerics num): L(num.L), N(num.N), wf(num.Nff), short_str(9), long_str(7){
 for(int i=0;i<9;++i){
  short_str(i).resize(wf);
  for(int w=0;w<wf;++w){
   short_str(i)(w).resize(2*N+1);
  }
 }
 for(int i=0;i<7;++i){
  long_str(i).resize(2*L+1,2*L+1);
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
    long_str(i)(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
   }
  }
 }
} 

Generalmatrix(int L_in, int N_in, int Nff_in): L(L_in), N(N_in), wf(Nff_in), short_str(9), long_str(7){
 for(int i=0;i<9;++i){
  short_str(i).resize(wf);
  for(int w=0;w<wf;++w){
   short_str(i)(w).resize(2*N+1);
  }
 }
 for(int i=0;i<7;++i){
  long_str(i).resize(2*L+1,2*L+1);
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
    long_str(i)(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
   }
  }
 }
} 
void initialize_random(){
 for(int c=0;c<9;++c){
  for(int w=0;w<wf;++w){
   for(int i=0;i<2*N+1;++i){
	for(int j=0;j<=i;++j){
	 short_str(c)(w)(i,j)=complex<double>((double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX);
    }
   }
  }
 }
 for(int c=0;c<7;++c){
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
	for(int j=max(-N,-N-l);j<=min(N,N-l);++j){
	 for(int i=max(-N,-N-k);i<=min(N,N-k);++i){
	  long_str(c)(l+L,k+L)(j-max(-N,-N-l),i-max(-N,-N-k))=(double)rand()/(double)RAND_MAX;
	 }
	}
   }
  }
 }
}
void initialize(double init){
 for(int c=0;c<9;++c){
  for(int w=0;w<wf;++w){
   for(int i=0;i<2*N+1;++i){
	for(int j=0;j<=i;++j){
	 short_str(c)(w)(i,j)=(complex<double>)init;
    }
   }
  }
 }
 for(int c=0;c<7;++c){
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
	for(int j=max(-N,-N-l);j<=min(N,N-l);++j){
	 for(int i=max(-N,-N-k);i<=min(N,N-k);++i){
	  long_str(c)(l+L,k+L)(j-max(-N,-N-l),i-max(-N,-N-k))=init;
	 }
	}
   }
  }
 }
}


Generalmatrix operator+ (const Generalmatrix & gm2){
 cout<<"operator+: L="<<L<<", N="<<N<<", wf="<<wf<<endl;
 Generalmatrix temp(L, N, wf);
 temp.short_str=short_str+gm2.short_str;
 temp.long_str=long_str+gm2.long_str;
 return temp;
}

void resize(Generalmatrix &gm){
 L=gm.L;
 N=gm.N;
 wf=gm.wf;
 cout<<"in resize: L="<<L<<endl;
 cout<<"in resize: N="<<N<<endl;
 short_str.resize(9);
 long_str.resize(7);
 for(int i=0;i<9;++i){
  short_str(i).resize(wf);
  for(int w=0;w<wf;++w){
   short_str(i)(w).resize(2*N+1);
  }
 }
 for(int i=0;i<7;++i){
  long_str(i).resize(2*L+1,2*L+1);
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
    long_str(i)(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
   }
  }
 }
}

void resize(int L_in, int N_in, int Nff_in){
 L=L_in;
 N=N_in;
 wf=Nff_in;
 cout<<"in resize: L="<<L<<endl;
 cout<<"in resize: N="<<N<<endl;
 short_str.resize(9);
 long_str.resize(7);
 for(int i=0;i<9;++i){
  short_str(i).resize(wf);
  for(int w=0;w<wf;++w){
   short_str(i)(w).resize(2*N+1);
  }
 }
 for(int i=0;i<7;++i){
  long_str(i).resize(2*L+1,2*L+1);
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
    long_str(i)(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
   }
  }
 }
}
  


void save(char *filename, const char *variable_in){
 char variable[255]="";
 char variable2[255]="";
 strcat(variable,variable_in);
 strcpy(variable2, variable);
 strcat(variable,"_short_str");
 short_str.save(filename,variable);
 strcat(variable2,"_long_str");
 long_str.save(filename,variable2);
}

double errnorm(double atol,double rtol, Generalmatrix &y1, Generalmatrix &y2) {
 double err, err_short, err_long;
 err_short=short_str.errnorm(atol,rtol,y1.short_str,y2.short_str);
 err_long=long_str.errnorm(atol,rtol,y1.long_str,y2.long_str);
 err=sqrt((err_short*err_short+err_long*err_long)/2.0);
 return err;
}


};

Generalmatrix operator*(const double &a, const Generalmatrix &gm);


class Physics{
 public:
 Numerics &num;
 double Vg;
 double h;
 double mu;
 double T;
 syma<complex<double> > hamiltonian;
 Physics(Numerics &num_in, double Vg_in, double h_in, double mu_in, double T_in): num(num_in), Vg(Vg_in), h(h_in), mu(mu_in), T(T_in) {
 hamiltonian=hamilton_zero(num.Nges,0,Vg,0.0, 0, 1.0);
 };
 void save(char *filename);
};

class Channel{
 public:
 Numerics & numerics;
 Blockmatrix<double> &freq_z_str;
 matrix<syma<complex<double> > > &freq_str;
 Channel(Numerics & numerics_in, Blockmatrix<double> & freq_z_str_in, matrix<syma<complex<double> > > & freq_str_in );
  
};

class Selfenergy{
 public:
 Numerics & numerics;
 matrix<syma<complex<double> > > & freq_str;
 Selfenergy(Numerics & numerics_in, matrix<syma<complex<double> > > & freq_str_in): numerics(numerics_in), freq_str(freq_str_in) {};
};

// declaration Vertex
// under construction

template <int mode>
class Vertex{
public:
Numerics & num;
Physics & phy;
Generalmatrix & data;
double Lambda;
Substitution<mode> sub;
matrix<double> wf_sub;
matrix<double> wbP_sub;
matrix<double> wbX_sub;
linear_ipol_bin<syma<complex<double> > > aPuu; 
linear_ipol_bin<syma<complex<double> > > aPdd; 
linear_ipol_bin<syma<complex<double> > > aPud; 
linear_ipol_bin<syma<complex<double> > > aXud; 
linear_ipol_bin<syma<complex<double> > > aDuu; 
linear_ipol_bin<syma<complex<double> > > aDdd; 
linear_ipol_bin<syma<complex<double> > > aDud; 
linear_ipol_bin<syma<complex<double> > > ERetu; 
linear_ipol_bin<syma<complex<double> > > ERetd; 
matrix<matrix< double > > &a0Puu_data;
matrix<matrix< double > > &a0Pdd_data;
matrix<matrix< double > > &a0Pud_data;
matrix<matrix< double > > &a0Xud_data;
matrix<matrix< double > > &a0Duu_data;
matrix<matrix< double > > &a0Ddd_data;
matrix<matrix< double > > &a0Dud_data;
Blockmatrix<double> a0Puu;
Blockmatrix<double> a0Pdd;
Blockmatrix<double> a0Pud;
Blockmatrix<double> a0Xud;
Blockmatrix<double> a0Duu;
Blockmatrix<double> a0Ddd;
Blockmatrix<double> a0Dud;
Barevertex_new &bare;
 
Vertex(Numerics & num_in, Physics & phy_in, Generalmatrix &data_in, double Lambda_in, Barevertex_new &bare_in): num(num_in), phy(phy_in), data(data_in), Lambda(Lambda_in), sub(1.0, Lambda, phy.h), wf_sub(num.Nff), wbP_sub(num.Nff), wbX_sub(num.Nff), 

aPuu(wbP_sub,data_in.short_str(0)), aPdd(wbP_sub,data_in.short_str(1)), aPud(wbP_sub,data_in.short_str(2)), aXud(wbX_sub,data_in.short_str(3)), aDuu(wbX_sub,data_in.short_str(4)), aDdd(wbX_sub,data_in.short_str(5)), aDud(wbX_sub,data_in.short_str(6)), ERetu(wf_sub, data_in.short_str(7)), ERetd(wf_sub, data_in.short_str(8))

,a0Puu_data(data_in.long_str(0)), a0Pdd_data(data_in.long_str(1)), a0Pud_data(data_in.long_str(2)), a0Xud_data(data_in.long_str(3)), a0Duu_data(data_in.long_str(4)), a0Ddd_data(data_in.long_str(5)), a0Dud_data(data_in.long_str(6))

,a0Puu(num.L, num.N, a0Puu_data), a0Pdd(num.L, num.N, a0Pdd_data), a0Pud(num.L, num.N, a0Pud_data), a0Xud(num.L, num.N, a0Xud_data), a0Duu(num.L, num.N, a0Duu_data), a0Ddd(num.L, num.N, a0Ddd_data), a0Dud(num.L, num.N, a0Dud_data) 

,bare(bare_in)
{

 for(int i=0;i<num.Nff;++i){
  wf_sub(i)=sub.subst_concatenated(num.wf(i)); 
  wbP_sub(i)=sub.subst_concatenated(num.wf(i)); 
  wbX_sub(i)=sub.subst_concatenated(num.wf(i)); 
 }
}
void set_aPuu(int freq_number, syma<complex<double> > A){
 data.short_str(0)(freq_number)=A;
}
void set_aPdd(int freq_number, syma<complex<double> > A){
 data.short_str(1)(freq_number)=A;
}
void set_aPud(int freq_number, syma<complex<double> > A){
 data.short_str(2)(freq_number)=A;
}
void set_aXud(int freq_number, syma<complex<double> > A){
 data.short_str(3)(freq_number)=A;
}
void set_aDuu(int freq_number, syma<complex<double> > A){
 data.short_str(4)(freq_number)=A;
}
void set_aDdd(int freq_number, syma<complex<double> > A){
 data.short_str(5)(freq_number)=A;
}
void set_aDud(int freq_number, syma<complex<double> > A){
 data.short_str(6)(freq_number)=A;
}
void set_ERetu(int freq_number, syma<complex<double> > A){
 data.short_str(7)(freq_number)=A;
}
void set_ERetd(int freq_number, syma<complex<double> > A){
 data.short_str(8)(freq_number)=A;
}

syma<complex<double> > get_aPuu(int freq_number){
 return data.short_str(0)(freq_number);
}
syma<complex<double> > get_aPdd(int freq_number){
 return data.short_str(1)(freq_number);
}
syma<complex<double> > get_aPud(int freq_number){
 return data.short_str(2)(freq_number);
}
syma<complex<double> > get_aXud(int freq_number){
 return data.short_str(3)(freq_number);
}
syma<complex<double> > get_aDuu(int freq_number){
 return data.short_str(4)(freq_number);
}
syma<complex<double> > get_aDdd(int freq_number){
 return data.short_str(5)(freq_number);
}
syma<complex<double> > get_aDud(int freq_number){
 return data.short_str(6)(freq_number);
}

};


#endif
