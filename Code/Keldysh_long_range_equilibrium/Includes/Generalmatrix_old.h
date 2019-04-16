#ifndef GENERALMATRIX_17032017
#define GENERALMATRIX_17032017

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "Blockmatrix.h"
#include "Numerics.h"

using namespace std;


class Generalmatrix{
public:
 Numerics num;
 matrix < matrix < syma < complex < double > > > > short_str;
 matrix < matrix < matrix < double > > > long_str;
 Generalmatrix(){};
 Generalmatrix(Numerics &num_in);

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


Generalmatrix::Generalmatrix(Numerics &num_in): num(num_in), short_str(9), long_str(7){
 short_str(0).resize(num.NfbP);
 short_str(1).resize(num.NfbP);
 short_str(2).resize(num.NfbP);
 for(int i=0; i<num.NfbP; ++i){
  short_str(0)(i).resize(num.Nges);
  short_str(1)(i).resize(num.Nges);
  short_str(2)(i).resize(num.Nges);
 }

 short_str(3).resize(num.NfbX);
 short_str(4).resize(num.NfbX);
 short_str(5).resize(num.NfbX);
 short_str(6).resize(num.NfbX);
 for(int i=0; i<num.NfbX; ++i){
  short_str(3)(i).resize(num.Nges);
  short_str(4)(i).resize(num.Nges);
  short_str(5)(i).resize(num.Nges);
  short_str(6)(i).resize(num.Nges);
 }

 short_str(7).resize(num.Nff);
 short_str(8).resize(num.Nff);
 for(int i=0; i<num.Nff; ++i){
  short_str(7)(i).resize(num.Nges);
  short_str(8)(i).resize(num.Nges);
 }

 for(int i=0;i<7;++i){
  Blockmatrix<double> block(num.L, num.N, long_str(i));
  block.resize(num.L, num.N);
 }

}

void Generalmatrix::initialize_random(){
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

 
Generalmatrix operator*(const double &a, const Generalmatrix &gm);







#endif
 
