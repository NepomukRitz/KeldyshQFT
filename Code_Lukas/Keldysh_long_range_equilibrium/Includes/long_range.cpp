#include "long_range.h" 
 
using namespace std;

template <class T> double maximumsnorm(matrix<T> A){
 double norm=0.0;
 for(int j=0;j<A.dim_r;++j){
  for(int i=0;i<A.dim_c;++i){
   norm=max(norm,abs(A(j,i)));
  }
 }
 return norm;
}

template <class T> double maximumsnorm(syma<T> A){
 double norm=0.0;
 for(int j=0;j<A.dim;++j){
  for(int i=0;i<=j;++i){
   norm=max(norm,abs(A(j,i)));
  }
 }
 return norm;
}

double parity_check(matrix<double> A){
 double diff=0.0;
 int Nges=A.dim_r;
 for(int j=0;j<Nges;++j){
  for(int i=0;i<Nges;++i){
   diff=max(diff, abs( A(j,i) - A(Nges-j-1,Nges-i-1) )); 
  }
 }
 return diff;
}

double parity_check(matrix<complex<double> > A){
 double diff=0.0;
 int Nges=A.dim_r;
 for(int j=0;j<Nges;++j){
  for(int i=0;i<Nges;++i){
   diff=max(diff, abs( A(j,i) - A(Nges-j-1,Nges-i-1) )); 
  }
 }
 return diff;
}

double parity_check(syma<complex<double> > A){
 double diff=0.0;
 int Nges=A.dim;
 for(int j=0;j<Nges;++j){
  for(int i=0;i<=j;++i){
   diff=max(diff, abs( A(j,i) - A(Nges-i-1, Nges-j-1) )); 
  }
 }
 return diff;
}

double parity_check(matrix<matrix<double> > A){
 double diff=0.0;
 int L=(A.dim_r -1)/2;
 int N= (A(L,L).dim_r -1)/2;
 int twoN=2*N;
 Blockmatrix<double> Block_A(L,N,A);
 for(int jmin, jmax, l=-L;l<=L;++l){ jmin=max(0,-l); jmax=min(twoN,twoN-l);
  for(int imin, imax, k=-L;k<=L;++k){ imin=max(0,-k); imax=min(twoN,twoN-k);
   for(int j=jmin;j<=jmax;++j){
    for(int i=imin;i<=imax;++i){
	 diff=max(diff, abs(Block_A.phyaccess(l,k,j,i) - Block_A.phyaccess(-l,-k,twoN-j,twoN-i))); 
	}
   }
  }
 }
 return diff;
}

double parity_check(matrix<matrix<complex<double> > > A){
 double diff=0.0;
 int L=(A.dim_r -1)/2;
 int N= (A(L,L).dim_r -1)/2;
 int twoN=2*N;
 Blockmatrix<complex<double> > Block_A(L,N,A);
 for(int jmin, jmax, l=-L;l<=L;++l){ jmin=max(0,-l); jmax=min(twoN,twoN-l);
  for(int imin, imax, k=-L;k<=L;++k){ imin=max(0,-k); imax=min(twoN,twoN-k);
   for(int j=jmin;j<=jmax;++j){
    for(int i=imin;i<=imax;++i){
	 diff=max(diff, abs(Block_A.phyaccess(l,k,j,i) - Block_A.phyaccess(-l,-k,twoN-j,twoN-i))); 
	}
   }
  }
 }
 return diff;
}

double symmetry_check(matrix<double> A){
 double diff=0.0;
 int Nges=A.dim_r;
 for(int j=0;j<Nges;++j){
  for(int i=0;i<Nges;++i){
   diff=max(diff,abs(A(j,i)-A(i,j)));
  }
 }
 return diff;
}

double symmetry_check(matrix<matrix<double> > A){
 double diff=0.0;
 int L=(A.dim_r -1)/2;
 int N= (A(L,L).dim_r -1)/2;
 int twoN=2*N;
 Blockmatrix<double> Block_A(L,N,A);
 for(int jmin, jmax, l=-L;l<=L;++l){ jmin=max(0,-l); jmax=min(twoN,twoN-l);
  for(int imin, imax, k=-L;k<=L;++k){ imin=max(0,-k); imax=min(twoN,twoN-k);
   for(int j=jmin;j<=jmax;++j){
    for(int i=imin;i<=imax;++i){
	 diff=max(diff, abs(Block_A.phyaccess(l,k,j,i) - Block_A.phyaccess(k,l,i,j))); 
	}
   }
  }
 }
 return diff;
}



 

// definition blockmatrix

Generalmatrix operator*(const double &a, const Generalmatrix &gm){
 Generalmatrix temp(gm.L, gm.N, gm.wf);
 temp.short_str=a*gm.short_str;
 temp.long_str=a*gm.long_str;
 return temp;
}

int komp(int L, int N, int l, int j){
 int Ng=2*N+1;
 if(l<=1){
  return (Ng-L-1)*(L+l)+((l+L)*(l+L+1))/2+j-max(-N,-N-l);
 }
 else{
  return Ng*(L+l)-(L*(L+1)+(l-1)*l)/2+j-max(-N,-N-l);
 }
}

int dim_all(int L, int N){
 return (2*N+1)*(1+2*L)-L*(L+1);
}

template <class T> void resize_block(int L, int N, matrix<matrix<T> > &A){
 A.resize(2*L+1,2*L+1);
 for(int l=-L;l<=L;++l){
  for(int k=-L;k<=L;++k){
   A(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
  }
 }
}

template <class T> void initialize_block(int L, int N, matrix<matrix<T> > &A, T x){
 A.resize(2*L+1,2*L+1);
 for(int l=-L;l<=L;++l){
  for(int k=-L;k<=L;++k){
   A(l+L,k+L).resize(2*N+1-abs(l),2*N+1-abs(k));
   A(l+L,k+L)=x; 
  }
 }
}

bool inrange_block(int L, int N, int l, int k, int j, int i){
 return (-L<=l && l<=L && -L<=k && k<=L && max(-N,-N-l)<=j && j<=min(N,N-l) && max(-N,-N-k)<=i && i<=min(N,N-k));
}

bool inrange_block_phy(int L, int N, int l, int k, int j, int i){
 return (-L<=l && l<=L && -L<=k && k<=L && max(0,-l)<=j && j<=min(2*N,2*N-l) && max(0,-k)<=i && i<=min(2*N,2*N-k));
}

// constructor and memberfunctions of Blockmatrix

 
template <class T> Blockmatrix<T>::Blockmatrix(int L_in, int N_in, matrix<matrix<T> > &A_in): L(L_in), N(N_in), A(A_in) {};

template <class T> T & Blockmatrix<T>::operator()(int l, int k, int j, int i){
  return A(l+L,k+L)(j-max(-N,-N-l),i-max(-N,-N-k));
}

template <class T> T & Blockmatrix<T>::fastaccess(int ltilde, int ktilde, int jtilde, int itilde){
 return A(ltilde,ktilde)(jtilde,itilde);
}

template <class T> T & Blockmatrix<T>::phyaccess(int l, int k, int j_asymmetric, int i_asymmetric){
 return A(l+L,k+L)(j_asymmetric-max(0,-l),i_asymmetric-max(0,-k));
}

template <class T> matrix<T> Blockmatrix<T>::convert_to_matrix(){
  int Ng=dim_all(L,N);
  matrix<T> erg(Ng,Ng);
  for(int l=-L;l<=L;++l){
   for(int k=-L;k<=L;++k){
	for(int j=max(-N,-N-l);j<=min(N,N-l);++j){
	 for(int i=max(-N,-N-k);i<=min(N,N-k);++i){
	  erg(komp(L,N,l,j),komp(L,N,k,i))=(*this)(l,k,j,i);
	 }
	}
   }
  }
 return erg;
}
 
template <class T> void Blockmatrix<T>::resize(int L, int N){
  resize_block(L,N,A);
 }
 
template <class T> void Blockmatrix<T>::initialize(int L, int N, T x){
  initialize_block(L,N,A,x);
 }
 
template <class T> void Blockmatrix<T>::save(char *filename, char *variable){
  matrix<T> B=(*this).convert_to_matrix();
  B.save(filename,variable);
 }
 
template <class T> bool Blockmatrix<T>::inrange(int l, int k, int j, int i){
  return inrange_block(L,N,l,k,j,i);
 }

//functions that use Blockmatrix

template <class T> matrix<matrix<T> > convert_to_blockmatrix(matrix<T> A, int L, int N){
 int Nmitte=2*N+1;
 matrix<matrix<T> > erg(2*L+1,2*L+1);
 for(int l=-L;l<=L;++l){
  for(int k=-L;k<=L;++k){
   erg(l+L,k+L).resize(Nmitte-abs(l),Nmitte-abs(k));
  }
 }
 for(int l=-L;l<=L;++l){
  for(int k=-L;k<=L;++k){
   for(int j=max(-N,-N-l);j<=min(N,N-l);++j){
	for(int i=max(-N,-N-k);i<=min(N,N-k);++i){
	 erg(l+L,k+L)(j-max(-N,-N-l),i-max(-N,-N-k))=A(komp(L,N,l,j),komp(L,N,k,i));
	}
   }
  }
 }
 return erg;
}

//line to finite line segments
template <>
double Substitution<0>::subst_concatenated (double omega)
{
	if (omega<-6.*taul)
	{
		double oshifted = (omega/taul+6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)-6.*taul;
	}
	else if (omega==-6.*taul)
	{
		return -6.*taul;
	}
	else if (omega<-2.*taul)
	{
		return -2.*sqrt(taul)*sqrt(-omega-2.*taul)-2.*taul;
	}
	else if (omega==-2.*taul)
	{
		return -2.*taul;
	}
	else if (omega<2.*taul)
	{
		double x = omega/taul;
		return taul*(sqrt(2.+x)-sqrt(2.-x));
	}
	else if (omega==2.*taul)
	{
		return 2.*taul;
	}
	else if (omega<6.*taul)
	{
		return 2.*sqrt(taul)*sqrt(omega-2.*taul)+2.*taul;
	}
	else if (omega==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double oshifted = (omega/taul-6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)+6.*taul;
	}
}

//finite line segments to line
template <>
double Substitution<0>::resu_concatenated (double y)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		return (-2.*taul*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)-6.*taul;
	}
	else if (y==-6.*taul)
	{
		return -6.*taul;
	}
	else if (y<-2.*taul)
	{
		return -2.*taul-(y+2.*taul)*(y+2.*taul)/4./taul;
	}
	else if (y==-2.*taul)
	{
		return -2.*taul;
	}
	else if (y<2.*taul)
	{
		if (y==.0)
			return .0;
		else
			return y*sqrt(4.*taul*taul/y/y-(y*y/taul/taul-4.)*(y*y/taul/taul-4.)*taul*taul/4./y/y);
	}
	else if (y==2.*taul)
	{
		return 2.*taul;
	}
	else if (y<6.*taul)
	{
		return 2.*taul+(y-2.*taul)*(y-2.*taul)/4./taul;
	}
	else if (y==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double yshifted = (y-6.*taul)/(7.-6.*taul);
		return (-2.*taul*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)+6.*taul;
	}
}

//measure on finite line segments
template <>
double Substitution<0>::weight_concatenated (double y)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul)*(1.+Lambda);
	}
	else if (y==-6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-6taul. This should not happen" << endl;
		return 1e20;
	}
	else if (y<-2.*taul)
	{
		return -1.-y/2./taul;
	}
	else if (y==-2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-2taul. This should not happen" << endl;
		return 1e20;
	}
	else if (y<2.*taul)
	{
		double x = Substitution<0>::resu_concatenated(y);
		return 2./((1./sqrt(2.+x/taul)+1./sqrt(2.-x/taul)));
	}
	else if (y==2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=2taul. This should not happen" << endl;
		return 1e20;
	}
	else if (y<6.*taul)
	{
		return -1.+y/2./taul;
	}
	else if (y==6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=6taul. This should not happen" << endl;
		return 1e20;
	}
	else
	{
		double yshifted = (y-6.)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul)*(1.+Lambda);
	}
}


template <>
double Substitution<1>::subst_concatenated (double omega){
 matrix<double> breaks(4);
 breaks(0)=-2.;
 breaks(1)=-2.;
 breaks(2)=2.;
 breaks(3)=2.;
 return subst_concatenated_fast_decay(omega,breaks,Lambda); 
}

template <>
double Substitution<1>::resu_concatenated (double x){
 matrix<double> breaks(4);
 breaks(0)=-2.;
 breaks(1)=-2.;
 breaks(2)=2.;
 breaks(3)=2.;
 return resu_concatenated_fast_decay(x,breaks,Lambda); 
}

template <>
double Substitution<1>::weight_concatenated (double x){
 matrix<double> breaks(4);
 breaks(0)=-2.;
 breaks(1)=-2.;
 breaks(2)=2.;
 breaks(3)=2.;
 return weight_concatenated_fast_decay(x,breaks,Lambda); 
}


//memberfunctions of Physics:

void Physics::save(char *filename){
 num.save(filename);
 matrix<double> tmp(1);
 tmp(0)=Vg;
 tmp.save(filename,"Vg");
 tmp(0)=h;
 tmp.save(filename,"h");
 tmp(0)=mu;
 tmp.save(filename,"mu");
 tmp(0)=T;
 tmp.save(filename,"T");
 hamiltonian.save(filename,"hamiltonian");
}



//constructor and memberfunctions of Numerics

Numerics::Numerics(const int number_of_frequencies_at_which_g_and_s_are_precomputed_in, const double delta_around_dangerous_frequencies_in, const double delta_around_dangerous_frequencies_at_pm7_in, int L_in, int N_in, int Nff_in, int Nfb_in, matrix<double> & wf_in,  matrix<double> & wbP_in, matrix<double> & wbX_in, double Vg_in, double mu_in, double T_in):
number_of_frequencies_at_which_g_and_s_are_precomputed(number_of_frequencies_at_which_g_and_s_are_precomputed_in),
delta_around_dangerous_frequencies(delta_around_dangerous_frequencies_in),
delta_around_dangerous_frequencies_at_pm7(delta_around_dangerous_frequencies_at_pm7_in),
L(L_in),
N(N_in),
Nff(Nff_in),
Nfb(Nfb_in),
wf(wf_in),
wbP(wbP_in),
wbX(wbX_in),
Vg(Vg_in),
mu(mu_in),
T(T_in) {
Nges=2*N_in+1; 
wf.resize(Nff);
wbP.resize(Nfb);
wbX.resize(Nfb);
double taul=1.0;

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
matrix<double> helper(Nfb);
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

//Symmetrisiere die Frequenzen:
matrix<double> wf_tmp;
int z=0;
for(int fj=0; wf(fj)<0; ++fj, ++z){};
wf_tmp.resize(2*z+1);
for(int fj=0; wf(fj)<0; ++fj){
 wf_tmp(2*fj)=wf(fj);
 wf_tmp(2*fj+1)=-wf(fj);
}
wf_tmp(2*z)=0;
wf_tmp.sort();
wf=wf_tmp;
wbP=wf;
wbX=wf;
Nff =wf.dim_c;
Nfb =wbP.dim_c;

//Bestimme Positionen von 0 und 2mu:
for(int fj=0; fj<Nff; ++fj){
 if(wf(fj)==0){ 
  pos_Nff_0=fj;
  break; 
 }
}  
for(int fj=0; fj<Nff; ++fj){
 if(wbP(fj)==2*mu){ 
  pos_Nfb_2mu=fj;
  break; 
 }
}  
for(int fj=0; fj<Nff; ++fj){
 if(wbX(fj)==0){ 
  pos_Nfb_0=fj;
  break; 
 }
}  
 



//Frequencies of precomputation

double Lambda_fake=0.0; // Solange mu im Band liegt ist dieser Wert egal.
frequencies_of_precomputation.resize(number_of_frequencies_at_which_g_and_s_are_precomputed+22);

for (int i=0; i<number_of_frequencies_at_which_g_and_s_are_precomputed; i++)
 frequencies_of_precomputation(i) = -7.+14.*(double)(i+1)/(double)(number_of_frequencies_at_which_g_and_s_are_precomputed+1);
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-22) = -7.+delta_around_dangerous_frequencies;
frequencies_of_precomputation(frequencies_of_precomputation.dim_c-21) =  7.-delta_around_dangerous_frequencies;
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

 };

void Numerics::save(char *filename){
matrix<double> tmp(1);
tmp(0)=number_of_frequencies_at_which_g_and_s_are_precomputed;
tmp.save(filename,"number_of_frequencies_at_which_g_and_s_are_precomputed");
tmp(0)=delta_around_dangerous_frequencies;
tmp.save(filename,"delta_around_dangerous_frequencies");
tmp(0)=delta_around_dangerous_frequencies_at_pm7;
tmp.save(filename,"delta_around_dangerous_frequencies_at_pm7");
tmp(0)=L;
tmp.save(filename,"L");
tmp(0)=N;
tmp.save(filename,"N");
tmp(0)=Nges;
tmp.save(filename,"Nges");
tmp(0)=Nff;
tmp.save(filename,"Nff");
tmp(0)=Nfb;
tmp.save(filename,"Nfb");
tmp(0)=pos_Nff_0;
tmp.save(filename,"pos_Nff_0");
tmp(0)=pos_Nfb_0;
tmp.save(filename,"pos_Nfb_0");
tmp(0)=pos_Nfb_2mu;
tmp.save(filename,"pos_Nfb_2mu");
wf.save(filename,"wf");
wbP.save(filename,"wbP");
wbX.save(filename,"wbX");
frequencies_of_precomputation.save(filename,"frequencies_of_precomputation");
tmp(0)=Vg;
tmp.save(filename,"Vg");
tmp(0)=mu;
tmp.save(filename,"mu");
tmp(0)=T;
tmp.save(filename,"T");
};


//constructor and memberfunctions of Channel

Channel::Channel(Numerics & numerics_in, Blockmatrix<double> & freq_z_str_in, matrix<syma<complex<double> > > & freq_str_in ): numerics(numerics_in), freq_z_str(freq_z_str_in), freq_str(freq_str_in) {}; 


template class Blockmatrix<double>; 
template class Blockmatrix<complex<double> >; 
template void resize_block(int L, int N, matrix<matrix<double> > &A);
template void initialize_block(int L, int N, matrix<matrix<double> > &A, double x);
template matrix<matrix<double> > convert_to_blockmatrix(matrix<double> A, int L, int N);
template matrix<matrix<complex<double> > > convert_to_blockmatrix(matrix<complex<double> > A, int L, int N);
template double maximumsnorm(matrix<double> A);
template double maximumsnorm(syma<complex<double> > A);
