#ifndef _kagome_facil_hpp_
#define _kagome_facil_hpp_

#include <complex>
#include <vector>
#include "H5Cpp.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include <mpi.h>
#include <fstream>
#include <stdlib.h>

#ifndef MPI
#define MPI 1
#endif

#ifndef temp
#define temp 0
#endif

#ifndef sym
#define sym 2
#endif
#ifndef reg
#define reg 2
#endif
#ifndef grid //grid 1= linear, grid2 = log+lin
#define grid 2
#endif

using namespace std;



//const H5std_string	DATASET_R("R");
//const H5std_string	DATASET_K1("K1");
//const H5std_string	DATASET_K2("K2");
//const H5std_string	DATASET_K2b("K2b");
//const H5std_string	DATASET_irred("irred");
//const H5std_string	DATASET_sus("sus");
//const H5std_string	FERM_FREQS_LIST("ferm_freqslist");
//const H5std_string	BOS_FREQS_LIST("bos_freqslist");
//const H5std_string	SELF_LIST("selflist");
//const H5std_string	LAMBDA_LIST("lambdas");
//const H5std_string      MEMBER1( "spin_re" );
////const H5std_string      MEMBER2( "spin_im" );
//const H5std_string      MEMBER2( "dens_re" );
////const H5std_string      MEMBER4( "dens_im" );
////const H5std_string      MEMBER5( "re" );
//const H5std_string      MEMBER6( "im" );



const H5std_string	DATASET_irred("irred");

const H5std_string	DATASET_R_s("R_s");
const H5std_string	DATASET_K1_s("K1_s");
const H5std_string	DATASET_K2_s("K2_s");
const H5std_string	DATASET_K2b_s("K2b_s");


const H5std_string	DATASET_R_t("R_t");
const H5std_string	DATASET_K1_t("K1_t");
const H5std_string	DATASET_K2_t("K2_t");
const H5std_string	DATASET_K2b_t("K2b_t");



const H5std_string	DATASET_R_u("R_u");
const H5std_string	DATASET_K1_u("K1_u");
const H5std_string	DATASET_K2_u("K2_u");
const H5std_string	DATASET_K2b_u("K2b_u");



const H5std_string	DATASET_sus("sus");
const H5std_string	FERM_FREQS_LIST("ferm_freqslist");
const H5std_string	BOS_FREQS_LIST("bos_freqslist");
const H5std_string	SELF_LIST("selflist");
const H5std_string	LAMBDA_LIST("lambdas");
const H5std_string      MEMBER1( "spin_re" );
//const H5std_string      MEMBER2( "spin_im" );
const H5std_string      MEMBER2( "dens_re" );
//const H5std_string      MEMBER4( "dens_im" );
//const H5std_string      MEMBER5( "re" );
const H5std_string      MEMBER6( "im" );





//IMPORTANT: COMPILE WITH C++11 !!
//the first two of the three spacial lattice indices on the kagome lattice specify the unit cell and the third one denotes the site within this UC (nuc,nuc,3) components. the first two can take the values [-(nuc-1)/2,(nuc-1)/2] and the third can take the values 1,2,3
// For different lattice structures that can be described by three independent indices, one only needs to change to conversion from input to lattice values in definition of vertex classes and the site-conversion in the t-bubble.
extern const double pi;
extern const complex<double> ci; //complex i
//extern const int reg;
//extern const int sym;//with all symmetries if sym = 2, with only diagrammatic symmetries if sym=1 and without if sym =0.
extern const double sharp; // governs the sharpness of the second regularizer
extern const double wt ;//transition frequency between log and lin grid.
extern const double w0;
extern const double wmax;
extern const int nlog;
extern const double k;
extern const double delw;

extern const double T;

extern  int nw;
extern const int nw3;//number of lattice points in for the rest function
extern const int nw1;//number of lattice points in for the K1 - channel where only one freuency is saved
extern const int nw2;//number of lattice points in for the K2 - channel
extern const int nw1_q, nw2_q,nw2_w1,nw3_q,nw3_w1,nw3_w2;//number of frequency points that are effectively saved due to symmetries. Need to be adjusted in "kagome.cpp"
extern const int nuc_eff;//this value must be adjusted in "kagome.cpp" if nuc >5 due to the lattice site parametrization



extern const int nlin;//total number of point in linear part of grid
extern vector<double> ffreqs;//fermionic matsubara freqs
extern vector<double> bfreqs;//bosnonic matsubara freqs, NOTE: in the case T=0 (which is implemented here), these two grids are equivalent
extern vector<double> Lambdas;//grid for RG flow parameter Lambda
extern const int nuc ; //number of unit cells in one direction on the kagome lattice that are included in the flow. Should be odd such that geometry is symmetric
extern const double d_c ;//cutoff distance
extern const int     NX;                     // dataset dimensions
extern const int     NY;
extern const int     NZ;                     // dataset dimensions
extern const int     RANK_R;
extern const int     RANK_K1;
extern const int     RANK_K2;
extern const int     RANK_K2b;
extern const int     RANK_irreducible;
extern const int     RANK_sus;
extern const int     RANK_self;


extern const int K3_dim6;
extern const int K3_dim5;
extern const int K3_dim4;
extern const int K3_dim3;
extern const int K3_dim2;
extern const int K3_dim1;


//K1:


extern const int K1_dim4 ;
extern const int K1_dim3;
extern const int K1_dim2;
extern const int K1_dim1;


//K2:

extern const int K2_dim5;
extern const int K2_dim4;
extern const int K2_dim3;
extern const int K2_dim2;
extern const int K2_dim1;

extern const int irred_dim3;
extern const int irred_dim2;
extern const int irred_dim1;


extern const complex<double> zero;


#if temp==0
const double wlimit = 1e10;//frequency used to perform the numerical limits for the asymptotic functions
int fconv(double w); //function that converts physical freqs to the lattice site on the freq mesh that corresponds to the closest grid frequency below
int fconv_n(double w, int n); //function that converts physical freqs to the lattice site on a mesh of size n that corresponds to the closest grid frequency below
#elif temp==1
const int wlimit = 1e5;//frequency used to perform the numerical limits for the asymptotic functions
int fconv(int w); //function that converts physical freqs to the lattice site on the freq mesh that corresponds to the closest grid frequency below
int fconv_n(int w, int n); //function that converts physical freqs to the lattice site on a mesh of size n that corresponds to the closest grid frequency below
#endif
inline bool exists_test (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}
//int bconv(double q); //function that converts physical freqs to the lattice site on the freq mesh that corresponds to the closest grid frequency below
//int bconv_n(double q, int n); //function that converts physical freqs to the lattice site on a mesh of size n that corresponds to the closest grid frequency below

inline bool comp(double a,double b)
{
    return (a < b);
}





struct site{
    int a;
    int b;
    int c;
public:
    site(int x,int y,int z){a = x; b = y; c = z;}
    void set(int x,int y,int z){a = x; b = y; c = z;}
    site(){}
};
site site_switch(int a, int b, int c);
site site_project(int a, int b, int c);//projects any site to its equivalent site in the upper half of the lattice
template <typename T0> int sgn(T0 val) {//signum function
    return (T0(0) < val) - (val < T0(0));
}
double distance(int a, int b, int c);//yields the distance of an arbitrary site to the reference site at the origin.



/*******************************CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX*********************************/

class svert{



    //K3:

    vector<double> K3 = vector<double>(K3_dim1);

    //K1:

    vector<double> K1 = vector<double>(K1_dim1);

    //K2:

    vector<double> K2 = vector<double>(K2_dim1);


#if sym==0

    //K2b:


    vector<double> K2b = vector<double>(K2_dim1);//same dimensionality as K2-clas
#endif

public:
    template<typename T0>
    double vvalsmooth(int, int, int,  T0,T0,T0,  char);
   template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
template<typename T0>
    double vvalsmooth(int,int,int, int, int, T0,T0,T0, char, int,char);//first two arguments: int red_side, int map
   template<typename T0>
    double vvalsmooth(int, int, int,  T0,T0,T0);

    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int, int, double);
        void R_direct_set( int,  double);
        void K1_direct_set( int,  double);
        void K2_direct_set( int, double);
    void K1_copy(vector<double>);
    vector<double> allout_K1(){return K1;}
    void K2_copy(vector<double>);
    vector<double> allout_K2(){return K2;}
#if sym==0
    void K2b_copy(vector<double>);
    vector<double> allout_K2b(){return K2b;}
#endif
    void K3_copy(vector<double>);
    vector<double> allout_K3(){return K3;}
#if sym==0
    void K2b_setvert( int, double);
    void K2b_direct_set( intdouble);
#endif
    double R_vval(int, int, int, int,int,int);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int,int);
    double K2b_vval(int, int, int, int,int);
    double R_vvalsmooth(int, int, int, double, double, double);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double,  double);
    double K2b_vvalsmooth(int, int, int, double,  double);
    friend svert operator*(double alpha, const svert& vertex);
    friend svert operator*(const svert& vertex, double alpha);
    friend svert operator+(const svert& vertex1, const svert& vertex2);
    friend svert abs_sum_tiny(const svert& vertex1, const svert& vertex2, double tiny);
    double K1_acc(int i){return K1[i];}
    double K2_acc(int i){return K2[i];}
#if sym ==0
    double K2b_acc(int i){return K2b[i];}
#endif
    double K3_acc(int i){return K3[i];}
};


//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template<typename T0>
double svert::vvalsmooth(int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel){

    if(distance(a,b,c) <= d_c){//cutoff distance

        T0 s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2.;
            w2_s = (w1-w2+q)/2.;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2.;
            w2_s = (-w1+w2-q)/2.;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2.;
            w2_s = w2-(q+w1)/2.;};
        double value=0;

#if temp==0
        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};

        value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations
#elif temp==1
            //K1:
            site x(a,b,c);
            int i1 = fconv_n(s,nw1);
            if(i1<nw1/2 && i1>0){
                x = site_switch(x.a,x.b,x.c);
                i1 = nw1-i1;};
            value +=  K1_vval(x.a,x.b,x.c,i1);

            //K2
            x.set(a,b,c);
            int i2 = fconv_n(s,nw2);
            int j2 = fconv_n(w1_s,nw2);
            if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
            if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
            else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            //K2b
            x=site_switch(a,b,c);
            int i2b = fconv_n(s,nw2);
            int j2b = fconv_n(w1_s,nw2);
            if(s%4==0 ){j2b = nw2-1-j2b;}
            else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
            if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
            if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
            else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            //K3
           x.set(a,b,c);
            int i3 =fconv_n(s,nw3);
            int j3 = fconv_n(w1_s,nw3);
            int k3 = fconv_n(w2_s,nw3);
            int i3_eff,j3_eff,k3_eff;
            if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
            if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
            if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
            else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
            value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);


#endif

        return value;}
    else{return 0;}

}
////additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class



//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
template<typename T0>
double svert::vvalsmooth(int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel, int p, char f){

    if(distance(a,b,c) <= d_c){//cutoff distance

        T0  s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2.;
            w2_s = (w1-w2+q)/2.;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2.;
            w2_s = (-w1+w2-q)/2.;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2.;
            w2_s = w2-(q+w1)/2.;};
        double value=0;



   #if temp==0
        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};

        if(p==1){
            if(channel=='s'){
                if(f == 'R' || f == 'M'){value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2b_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){value = K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'M'){value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        }

        else if(p==2){
            if(channel=='s'){
                if(f == 'R' || f == 'L'){value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to differentbare vertex
                else if(f == 'K' || f == 'M'){value = K1_vvalsmooth(a,b,c,s) + K2b_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'L'){value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        };

#elif temp==1
        if(p==1){
            if(channel=='s'){
                if(f == 'R' || f == 'M'){

                    //K2b
                    site x=site_switch(a,b,c);
                    int i2b = fconv_n(s,nw2);
                    int j2b = fconv_n(w1_s,nw2);
                    if(s%4==0 ){j2b = nw2-1-j2b;}
                    else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
                    if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                    //K3
                      x.set(a,b,c);
                    int i3 =fconv_n(s,nw3);
                    int j3 = fconv_n(w1_s,nw3);
                    int k3 = fconv_n(w2_s,nw3);
                    int i3_eff,j3_eff,k3_eff;
                    if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
                    if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                    if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                    value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
                }//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){
                    //K1:
                    site x(a,b,c);
                    int i1 = fconv_n(s,nw1);
                    if(i1<nw1/2 && i1>0){
                        x = site_switch(x.a,x.b,x.c);
                        i1 = nw1-i1;};
                    value +=  K1_vval(x.a,x.b,x.c,i1);

                    //K2
                      x.set(a,b,c);
                    int i2 = fconv_n(s,nw2);
                    int j2 = fconv_n(w1_s,nw2);
                    if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2,j2) ;
                };//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'M'){
                    //K1:
                    site x(a,b,c);
                    int i1 = fconv_n(s,nw1);
                    if(i1<nw1/2 && i1>0){
                        x = site_switch(x.a,x.b,x.c);
                        i1 = nw1-i1;};
                    value +=  K1_vval(x.a,x.b,x.c,i1);

                    //K2
                      x.set(a,b,c);
                    int i2 = fconv_n(s,nw2);
                    int j2 = fconv_n(w1_s,nw2);
                    if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                    //K2b
                    x=site_switch(a,b,c);
                    int i2b = fconv_n(s,nw2);
                    int j2b = fconv_n(w1_s,nw2);
                    if(s%4==0 ){j2b = nw2-1-j2b;}
                    else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
                    if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                    //K3
                    x.set(a,b,c);
                    int i3 =fconv_n(s,nw3);
                    int j3 = fconv_n(w1_s,nw3);
                    int k3 = fconv_n(w2_s,nw3);
                    int i3_eff,j3_eff,k3_eff;
                    if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
                    if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                    if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                    value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
                                    };
            };
        }

        else if(p==2){
            if(channel=='s'){
                if(f == 'R' || f == 'L'){

                    //K2
                    site x(a,b,c);
                    int i2 = fconv_n(s,nw2);
                    int j2 = fconv_n(w1_s,nw2);
                    if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                    //K3
                      x.set(a,b,c);
                    int i3 =fconv_n(s,nw3);
                    int j3 = fconv_n(w1_s,nw3);
                    int k3 = fconv_n(w2_s,nw3);
                    int i3_eff,j3_eff,k3_eff;
                    if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
                    if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                    if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                    value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

                }//if outer legs are conntected to differentbare vertex
                else if(f == 'K' || f == 'M'){

                    //K1:
                    site x(a,b,c);
                    int i1 = fconv_n(s,nw1);
                    if(i1<nw1/2 && i1>0){
                        x = site_switch(x.a,x.b,x.c);
                        i1 = nw1-i1;};
                    value +=  K1_vval(x.a,x.b,x.c,i1);


                    //K2b
                    x=site_switch(a,b,c);
                    int i2b = fconv_n(s,nw2);
                    int j2b = fconv_n(w1_s,nw2);
                    if(s%4==0 ){j2b = nw2-1-j2b;}
                    else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
                    if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;

                };//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'L'){

                    //K1:
                    site x(a,b,c);
                    int i1 = fconv_n(s,nw1);
                    if(i1<nw1/2 && i1>0){
                        x = site_switch(x.a,x.b,x.c);
                        i1 = nw1-i1;};
                    value +=  K1_vval(x.a,x.b,x.c,i1);

                    //K2
                     x.set(a,b,c);
                    int i2 = fconv_n(s,nw2);
                    int j2 = fconv_n(w1_s,nw2);
                    if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                    //K2b
                    x=site_switch(a,b,c);
                    int i2b = fconv_n(s,nw2);
                    int j2b = fconv_n(w1_s,nw2);
                    if(s%4==0 ){j2b = nw2-1-j2b;}
                    else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
                    if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                    if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2 = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                    value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                    //K3
                      x.set(a,b,c);
                    int i3 =fconv_n(s,nw3);
                    int j3 = fconv_n(w1_s,nw3);
                    int k3 = fconv_n(w2_s,nw3);
                    int i3_eff,j3_eff,k3_eff;
                    if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
                    if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                    if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                    else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                    value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
                };
            };
        };
#endif
        return value;}
    else{return 0;}
}
template<typename T0>

    double svert::vvalsmooth(int red_side,int map, int a, int b, int c, T0 q, T0 w1, T0 w2, char channel, int p, char f){


        //THIS FUNCTION IS NEEDED ONLY WHEN BUBBLE FUNCTIONS ARE USED WITH VERTEX OF TYPE "PARVERT" INSTEAD OF "FULLVERT". LEADS TO PREVIOUS FUNCTION.
        return vvalsmooth(a, b, c, q, w1,w2, channel, p, f);
    }
template<typename T0>
double svert::vvalsmooth(int a, int b, int c,  T0 q, T0 w1, T0 w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    if(distance(a,b,c) <= d_c){//cutoff distance

        T0  s,w1_s,w2_s;
        s = q;
        w1_s = w1;
        w2_s = w2;

#if temp==0
        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};

        double value=0;
        value = R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations
#elif temp==1
        double value=0;

        //K1:
        site x(a,b,c);
        int i1 = fconv_n(s,nw1);
        if(i1<nw1/2 && i1>0){
            x = site_switch(x.a,x.b,x.c);
            i1 = nw1-i1;};
        value +=  K1_vval(x.a,x.b,x.c,i1);

        //K2
         x.set(a,b,c);
        int i2 = fconv_n(s,nw2);
        int j2 = fconv_n(w1_s,nw2);
        if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
        if(s%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c);}
        else if(s%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
        value += K2_vval(x.a,x.b,x.c,i2,j2) ;

        //K2b
        x=site_switch(a,b,c);
        int i2b = fconv_n(s,nw2);
        int j2b = fconv_n(w1_s,nw2);
        if(s%4==0 ){j2b = nw2-1-j2b;}
        else if(s%4!=0 && j2b>0){j2b = nw2-j2b;};
        if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
        if(s%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c);}
        else if(s%4!=0 && j2b < nw2/2 && j2b>0){j2 = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
        value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


        //K3
          x.set(a,b,c);
        int i3 =fconv_n(s,nw3);
        int j3 = fconv_n(w1_s,nw3);
        int k3 = fconv_n(w2_s,nw3);
        int i3_eff,j3_eff,k3_eff;
        if(i3<nw3/2){x = site_switch(x.a,x.b,x.c);i3_eff = nw3-i3;}
        if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
        if(s%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
        else if(s%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
        value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
#endif
        return value;}
    else{return 0;}

}
class tvert{

    //    //K3:
    //    int K3_dim6=3;
    //    int K3_dim5=(nuc_eff+1)/2 * K3_dim6;
    //    int K3_dim4=nuc_eff  * K3_dim5;
    //    int K3_dim3=nw3_w2*K3_dim4;
    //    int K3_dim2=nw3_w1 * K3_dim3;
    //    int K3_dim1=nw3_q * K3_dim2;


    //    //K1:


    //    int K1_dim4=3 ;
    //    int K1_dim3=(nuc_eff+1)/2 * K1_dim4 ;
    //    int K1_dim2 = nuc_eff * K1_dim3;
    //    int K1_dim1=nw1_q*K1_dim2;


    //    //K2:

    //    int K2_dim5=3 ;
    //    int K2_dim4=(nuc_eff+1)/2 * K2_dim5 ;
    //    int K2_dim3 = nuc_eff * K2_dim4;
    //    int K2_dim2=nw2_w1*K2_dim3;
    //    int K2_dim1=nw2_q * K2_dim2;



    //K3:

    vector<double> K3 = vector<double>(K3_dim1);

    //K1:

    vector<double> K1 = vector<double>(K1_dim1);

    //K2:

    vector<double> K2 = vector<double>(K2_dim1);


#if sym==0
    //K2b:


    vector<double> K2b = vector<double>(K2_dim1);//same dimensionality as K2-clas
#endif



public:
    template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0, char);
        template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
        template<typename T0>
    double vvalsmooth(int,int,int, int, int, T0,T0,T0, char, int,char);//first two arguments: int red_side, int map
        template<typename T0>
    double vvalsmooth(int, int, int,T0,T0,T0);

    void K1_copy(vector<double>);
    vector<double> allout_K1(){return K1;}
    void K2_copy(vector<double>);
    vector<double> allout_K2(){return K2;}
#if sym==0
    void K2b_copy(vector<double>);
    vector<double> allout_K2b(){return K2b;};
#endif
    void K3_copy(vector<double>);
    vector<double> allout_K3(){return K3;}
    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int,int, double);
        void R_direct_set( int,  double);
        void K1_direct_set( int,  double);
        void K2_direct_set( int, double);
#if sym==0
    void K2b_setvert( int, int, int, int,int, double);
    void K2b_direct_set( int, double);
#endif
    double R_vval(int, int, int, int,int,int);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int,int);
    double K2b_vval(int, int, int, int,int);
    double R_vvalsmooth(int, int, int, double, double, double);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double,  double);
    double K2b_vvalsmooth(int, int, int, double,  double);
    friend tvert operator*(double alpha, const tvert& vertex);
    friend tvert operator*(const tvert& vertex, double alpha);
    friend tvert operator+(const tvert& vertex1, const tvert& vertex2);
    friend tvert abs_sum_tiny(const tvert& vertex1, const tvert& vertex2, double tiny);
    double K1_acc(int i){return K1[i];}
    double K2_acc(int i){return K2[i];}
#if sym ==0
    double K2b_acc(int i){return K2b[i];}
#endif
    double K3_acc(int i){return K3[i];}
};
template<typename T0>
            double tvert::vvalsmooth(int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45


            if(distance(a,b,c) <= d_c){//cutoff distance

            T0  t,w1_t,w2_t;
            if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2.;
            w2_t = (-w1-w2+q)/2.;
}
            else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
            else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2.;
            w2_t =(w1+w2+q)/2.;}
            else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2.;
            w2_t = w1+(q-w2)/2.;

};


            double value=0;
#if temp==0
            if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
            if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
            if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};



            value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t)+K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations

#elif temp==1
            //K1:
            site x(a,b,c);
            int i1 = fconv_n(t,nw1);
            if(i1<nw1/2 && i1>0){
                x = site_switch(x.a,x.b,x.c);
                i1 = nw1-i1;};
            value +=  K1_vval(x.a,x.b,x.c,i1);

            //K2
              x.set(a,b,c);
            int i2 = fconv_n(t,nw2);
            int j2 = fconv_n(w1_t,nw2);
            if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
            if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
            else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
            value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            //K2b
            x=site_switch(a,b,c);
            int i2b = fconv_n(t,nw2);
            int j2b = fconv_n(w1_t,nw2);
            i2b = nw2-i2b;

            if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
            if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
            else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
            value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            //K3
               x.set(a,b,c);
            int i3 =fconv_n(t,nw3);
            int j3 = fconv_n(w1_t,nw3);
            int k3 = fconv_n(w2_t,nw3);
            int i3_eff,j3_eff,k3_eff;
            if(i3<nw3/2){i3_eff = nw3-i3;}
            if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
            if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
            else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
            value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);



        #endif
            return value;}
            else{return 0;}


}
template<typename T0>
            double tvert::vvalsmooth(int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class


            if(distance(a,b,c) <= d_c){//cutoff distance

            T0 t,w1_t,w2_t;
            if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2.;
            w2_t = (-w1-w2+q)/2.;
}
            else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
            else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2.;
            w2_t =(w1+w2+q)/2.;}
            else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2.;
            w2_t = w1+(q-w2)/2.;

};


            double value=0;
#if temp==0
            if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
            if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
            if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};


            if(p==1){
            if(channel=='t'){
            if(f == 'R' || f == 'M'){value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2b_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'L'){value = K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='u'){
            if(f == 'R' || f== 'M'){value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);};
};
}

            else if(p==2){
            if(channel=='t'){
            if(f == 'R' || f == 'L'){value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'M'){value = K1_vvalsmooth(a,b,c,t) + K2b_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='u'){
            if(f == 'R' || f== 'L'){value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);

};
};
};
#elif temp==1
            if(p==1){
            if(channel=='t'){
            if(f == 'R' || f == 'M'){

                //K2b
                site x=site_switch(a,b,c);
                int i2b = fconv_n(t,nw2);
                int j2b = fconv_n(w1_t,nw2);
                i2b = nw2-i2b;

                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
                if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
                else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                  x.set(a,b,c);
                int i3 =fconv_n(t,nw3);
                int j3 = fconv_n(w1_t,nw3);
                int k3 = fconv_n(w2_t,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
                if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
                else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);


            }//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'L'){

                //K1:
                site x(a,b,c);
                int i1 = fconv_n(t,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(t,nw2);
                int j2 = fconv_n(w1_t,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
                if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
                else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            };//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='u'){
            if(f == 'R' || f== 'M'){

                //K1:
                site x(a,b,c);
                int i1 = fconv_n(t,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(t,nw2);
                int j2 = fconv_n(w1_t,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
                if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
                else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(t,nw2);
                int j2b = fconv_n(w1_t,nw2);
                i2b = nw2-i2b;

                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
                if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
                else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                 x.set(a,b,c);
                int i3 =fconv_n(t,nw3);
                int j3 = fconv_n(w1_t,nw3);
                int k3 = fconv_n(w2_t,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
                if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
                else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

            };
};
}

            else if(p==2){
            if(channel=='t'){
            if(f == 'R' || f == 'L'){

                //K2
                site x(a,b,c);
                int i2 = fconv_n(t,nw2);
                int j2 = fconv_n(w1_t,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
                if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
                else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                //K3
                  x.set(a,b,c);
                int i3 =fconv_n(t,nw3);
                int j3 = fconv_n(w1_t,nw3);
                int k3 = fconv_n(w2_t,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
                if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
                else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

            }//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'M'){
                //K1:
                site x(a,b,c);
                int i1 = fconv_n(t,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(t,nw2);
                int j2b = fconv_n(w1_t,nw2);
                i2b = nw2-i2b;

                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
                if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
                else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            };//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='u'){
            if(f == 'R' || f== 'L'){

                //K1:
                site x(a,b,c);
                int i1 = fconv_n(t,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(t,nw2);
                int j2 = fconv_n(w1_t,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
                if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
                else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(t,nw2);
                int j2b = fconv_n(w1_t,nw2);
                i2b = nw2-i2b;

                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
                if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
                else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                  x.set(a,b,c);
                int i3 =fconv_n(t,nw3);
                int j3 = fconv_n(w1_t,nw3);
                int k3 = fconv_n(w2_t,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
                if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
                else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

};
};
};
        #endif


            return value;}
            else{return 0;}


}
            //overload of previous function
template<typename T0>
            double tvert::vvalsmooth(int red_side, int map,int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel, int p, char f){


            return vvalsmooth( a, b, c, q, w1,  w2,  channel,p,  f);

}
template<typename T0>
            double tvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45



            if(distance(a,b,c) <= d_c){//cutoff distance

            T0  t,w1_t,w2_t;
            t = q;
            w1_t = w1;
            w2_t = w2;
            double value=0;

#if temp==0
            if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
            if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
            if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};

            value = R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations
        #elif temp==1
            //K1:
            site x(a,b,c);
            int i1 = fconv_n(t,nw1);
            if(i1<nw1/2 && i1>0){
                x = site_switch(x.a,x.b,x.c);
                i1 = nw1-i1;};
            value +=  K1_vval(x.a,x.b,x.c,i1);

            //K2
               x.set(a,b,c);
            int i2 = fconv_n(t,nw2);
            int j2 = fconv_n(w1_t,nw2);
            if(i2 < nw2/2 && i2>0){i2 = nw2-i2;};
            if(t%4==0 && j2 < nw2/2){j2 = nw2-1-j2; }
            else if(t%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2;};
            value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            //K2b
            x=site_switch(a,b,c);
            int i2b = fconv_n(t,nw2);
            int j2b = fconv_n(w1_t,nw2);
            i2b = nw2-i2b;

            if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b;};
            if(t%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; }
            else if(t%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b;};
            value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            //K3
             x.set(a,b,c);
            int i3 =fconv_n(t,nw3);
            int j3 = fconv_n(w1_t,nw3);
            int k3 = fconv_n(w2_t,nw3);
            int i3_eff,j3_eff,k3_eff;
            if(i3<nw3/2){i3_eff = nw3-i3;}
            if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;x = site_switch(x.a,x.b,x.c);};
            if(t%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;}
            else if(t%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;};
            value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

        #endif
            return value;}
            else{return 0;}


}
class uvert{
    //    int K3_dim6=3;
    //    int K3_dim5=(nuc_eff+1)/2 * K3_dim6;
    //    int K3_dim4=nuc_eff  * K3_dim5;
    //    int K3_dim3=nw3_w2*K3_dim4;
    //    int K3_dim2=nw3_w1 * K3_dim3;
    //    int K3_dim1=nw3_q * K3_dim2;


    //    //K1:


    //    int K1_dim4=3 ;
    //    int K1_dim3=(nuc_eff+1)/2 * K1_dim4 ;
    //    int K1_dim2 = nuc_eff * K1_dim3;
    //    int K1_dim1=nw1_q*K1_dim2;


    //    //K2:

    //    int K2_dim5=3 ;
    //    int K2_dim4=(nuc_eff+1)/2 * K2_dim5 ;
    //    int K2_dim3 = nuc_eff * K2_dim4;
    //    int K2_dim2=nw2_w1*K2_dim3;
    //    int K2_dim1=nw2_q * K2_dim2;



    //K3:

    vector<double> K3 = vector<double>(K3_dim1);

    //K1:

    vector<double> K1 = vector<double>(K1_dim1);

    //K2:

    vector<double> K2 = vector<double>(K2_dim1);


#if sym==0

    //K2b:


    vector<double> K2b = vector<double>(K2_dim1);//same dimensionality as K2-clas

#endif

public:
template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0, char);
    template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    template<typename T0>
    double vvalsmooth(int,int,int, int, int, T0,T0,T0, char, int,char);//first two arguments: int red_side, int map
    template<typename T0>
    double vvalsmooth(int, int, int, T0,T0,T0);


    void K1_copy(vector<double>);
    vector<double> allout_K1(){return K1;}
    void K2_copy(vector<double>);
    vector<double> allout_K2(){return K2;}
#if sym==0
    void K2b_copy(vector<double>);
    vector<double> allout_K2b(){return K2b;}
#endif
    void K3_copy(vector<double>);
    vector<double> allout_K3(){return K3;}
    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int,int, double);
        void R_direct_set( int,  double);
        void K1_direct_set( int,  double);
        void K2_direct_set( int, double);
#if sym==0
    void K2b_setvert( int, int, int, int,int, double);
    void K2b_direct_set( int, double);
#endif
    double R_vval(int, int, int, int,int,int);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int,int);
    double K2b_vval(int, int, int, int,int);
    double R_vvalsmooth(int, int, int, double, double, double);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double,  double);
    double K2b_vvalsmooth(int, int, int, double,  double);
    void bust(char);
    friend uvert operator*(double alpha, const uvert& vertex);
    friend uvert operator*(const uvert& vertex, double alpha);
    friend uvert operator+(const uvert& vertex1, const uvert& vertex2);
    friend uvert abs_sum_tiny(const uvert& vertex1, const uvert& vertex2, double tiny);
    double K1_acc(int i){return K1[i];}
    double K2_acc(int i){return K2[i];}
#if sym ==0
    double K2b_acc(int i){return K2b[i];}
#endif
    double K3_acc(int i){return K3[i];}

};
template<typename T0>
            double uvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45


            if(distance(a,b,c) <= d_c){//cutoff distance

            T0  u,w1_u,w2_u;
            if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2.;
            w2_u = (-w1+w2+q)/2.;}
            else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2.;
            w2_u = (w1+w2+q)/2.;}
            else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
            else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2.;
            w2_u = (w1+w2)/2.;};

            double value=0;
#if temp==0
            if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
            if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
            if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};


            value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u)  + K2b_vvalsmooth(a,b,c,u,w2_u)  ;//K2b is extracted from K2 by the symmetry relations
#elif temp==1
            //K1:
            site x(a,b,c);
            int i1 = fconv_n(u,nw1);
            if(i1<nw1/2 && i1>0){
                x = site_switch(x.a,x.b,x.c);
                i1 = nw1-i1;};
            value +=  K1_vval(x.a,x.b,x.c,i1);

            //K2
              x.set(a,b,c);
            int i2 = fconv_n(u,nw2);
            int j2 = fconv_n(w1_u,nw2);
            if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
            if(u%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c); }
            else if(u%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            //K2b
            x=site_switch(a,b,c);
            int i2b = fconv_n(u,nw2);
            int j2b = fconv_n(w1_u,nw2);
            i2b = nw2-i2b;
            if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
            if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
            else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            //K3
               x.set(a,b,c);
            int i3 =fconv_n(u,nw3);
            int j3 = fconv_n(w1_u,nw3);
            int k3 = fconv_n(w2_u,nw3);
            int i3_eff,j3_eff,k3_eff;
            if(i3<nw3/2){i3_eff = nw3-i3;x = site_switch(x.a,x.b,x.c);}
            if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
            if(u%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
            else if(u%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
            value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);

        #endif
            return value;  }
            else{return 0;}



}
template<typename T0>
            double uvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class


            if(distance(a,b,c) <= d_c){//cutoff distance



            T0  u,w1_u,w2_u;
            if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2.;
            w2_u = (-w1+w2+q)/2.;}
            else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2.;
            w2_u = (w1+w2+q)/2.;}
            else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
            else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2.;
            w2_u = (w1+w2)/2.;};

            double value=0;
#if temp==0
            if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
            if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
            if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

            if(p==1){
            if(channel=='u'){
            if(f == 'R' || f == 'M'){value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
            else if(f == 'K' || f == 'L'){value = K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='t'){
            if(f == 'R' || f== 'M'){value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);};
};
}

            else if(p==2){
            if(channel=='u'){
            if(f == 'R' || f == 'L'){value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'M'){value = K1_vvalsmooth(a,b,c,u) + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='t'){
            if(f == 'R' || f== 'L'){value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);
};
};
};

#elif temp==1
            if(p==1){
            if(channel=='u'){
            if(f == 'R' || f == 'M'){


                //K2b
                site x=site_switch(a,b,c);
                int i2b = fconv_n(u,nw2);
                int j2b = fconv_n(w1_u,nw2);
                i2b = nw2-i2b;
                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                  x.set(a,b,c);
                int i3 =fconv_n(u,nw3);
                int j3 = fconv_n(w1_u,nw3);
                int k3 = fconv_n(w2_u,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;x = site_switch(x.a,x.b,x.c);}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                if(u%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                else if(u%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
            }//if outer legs are conntected to different  vertex
            else if(f == 'K' || f == 'L'){

                //K1:
                site x(a,b,c);
                int i1 = fconv_n(u,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(u,nw2);
                int j2 = fconv_n(w1_u,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            };//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='t'){
            if(f == 'R' || f== 'M'){
                //K1:
                site x(a,b,c);
                int i1 = fconv_n(u,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(u,nw2);
                int j2 = fconv_n(w1_u,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(u,nw2);
                int j2b = fconv_n(w1_u,nw2);
                i2b = nw2-i2b;
                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                 x.set(a,b,c);
                int i3 =fconv_n(u,nw3);
                int j3 = fconv_n(w1_u,nw3);
                int k3 = fconv_n(w2_u,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;x = site_switch(x.a,x.b,x.c);}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                if(u%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                else if(u%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
            };
};
}

            else if(p==2){
            if(channel=='u'){
            if(f == 'R' || f == 'L'){value = R_vval(a,b,c,fconv_n(u,nw3),fconv_n(w1_u,nw3),fconv_n(w2_u,nw3)) + K2_vval(a,b,c,fconv_n(u,nw2),fconv_n(w2_u,nw2));}//if outer legs are conntected to different bare vertex
            else if(f == 'K' || f == 'M'){
                //K1:
                site x(a,b,c);
                int i1 = fconv_n(u,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);


                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(u,nw2);
                int j2b = fconv_n(w1_u,nw2);
                i2b = nw2-i2b;
                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            };//if outer legs are conntected to same bare vertex

}
            else if (channel=='s' || channel=='t'){
            if(f == 'R' || f== 'L'){
                //K1:
                site x(a,b,c);
                int i1 = fconv_n(u,nw1);
                if(i1<nw1/2 && i1>0){
                    x = site_switch(x.a,x.b,x.c);
                    i1 = nw1-i1;};
                value +=  K1_vval(x.a,x.b,x.c,i1);

                //K2
                  x.set(a,b,c);
                int i2 = fconv_n(u,nw2);
                int j2 = fconv_n(w1_u,nw2);
                if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2,j2) ;

                //K2b
                x=site_switch(a,b,c);
                int i2b = fconv_n(u,nw2);
                int j2b = fconv_n(w1_u,nw2);
                i2b = nw2-i2b;
                if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
                if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
                else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
                value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


                //K3
                  x.set(a,b,c);
                int i3 =fconv_n(u,nw3);
                int j3 = fconv_n(w1_u,nw3);
                int k3 = fconv_n(w2_u,nw3);
                int i3_eff,j3_eff,k3_eff;
                if(i3<nw3/2){i3_eff = nw3-i3;x = site_switch(x.a,x.b,x.c);}
                if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
                if(u%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
                else if(u%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
                value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
};
};
};

        #endif


            return value;  }
            else{return 0;}



}
            //overload of previous function
template<typename T0>
            double uvert::vvalsmooth(int red_side, int map,int a, int b, int c,  T0 q, T0 w1, T0 w2, char channel, int p, char f){


            return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
}
       template<typename T0>
            double uvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45


            if(distance(a,b,c) <= d_c){//cutoff distance


            T0  u,w1_u,w2_u;

            u = q;
            w1_u = w1;
            w2_u = w2;

            double value=0;
#if temp==0
            if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
            if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
            if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

            value = R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
           #elif temp==1
            //K1:
            site x(a,b,c);
            int i1 = fconv_n(u,nw1);
            if(i1<nw1/2 && i1>0){
                x = site_switch(x.a,x.b,x.c);
                i1 = nw1-i1;};
            value +=  K1_vval(x.a,x.b,x.c,i1);

            //K2
              x.set(a,b,c);
            int i2 = fconv_n(u,nw2);
            int j2 = fconv_n(w1_u,nw2);
            if(i2 < nw2/2 && i2>0){i2 = nw2-i2; x = site_switch(x.a,x.b,x.c);};
            if(u%4==0 && j2 < nw2/2){j2 = nw2-1-j2; x = site_switch(x.a,x.b,x.c); }
            else if(u%4!=0 && j2 < nw2/2 && j2>0){j2 = nw2-j2; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2,j2) ;

            //K2b
            x=site_switch(a,b,c);
            int i2b = fconv_n(u,nw2);
            int j2b = fconv_n(w1_u,nw2);
            i2b = nw2-i2b;
            if(i2b < nw2/2 && i2b>0){i2b = nw2-i2b; x = site_switch(x.a,x.b,x.c);};
            if(u%4==0 && j2b < nw2/2){j2b = nw2-1-j2b; x = site_switch(x.a,x.b,x.c); }
            else if(u%4!=0 && j2b < nw2/2 && j2b>0){j2b = nw2-j2b; x = site_switch(x.a,x.b,x.c);};
            value += K2_vval(x.a,x.b,x.c,i2b,j2b) ;


            //K3
               x.set(a,b,c);
            int i3 =fconv_n(u,nw3);
            int j3 = fconv_n(w1_u,nw3);
            int k3 = fconv_n(w2_u,nw3);
            int i3_eff,j3_eff,k3_eff;
            if(i3<nw3/2){i3_eff = nw3-i3;x = site_switch(x.a,x.b,x.c);}
            if( abs(j3-nw3/2) > abs(k3-nw3/2)){j3_eff = k3; k3_eff = j3;};
            if(u%4==0 && j3_eff<nw3/2){j3_eff = nw3-1-j3_eff; k3_eff = nw3-1-k3_eff;x = site_switch(x.a,x.b,x.c);}
            else if(u%4!=0 && j3_eff<nw3/2 && j3_eff>0 && k3_eff >0){j3_eff = nw3-j3_eff; k3_eff = nw3-k3_eff;x = site_switch(x.a,x.b,x.c);};
            value += R_vval(x.a,x.b,x.c,i3_eff,j3_eff,k3_eff);
        #endif

            return value;  }
            else{return 0;}



};
class irreducible{


    //    int irred_dim3=3;
    //    int irred_dim2=(nuc_eff+1)/2 * irred_dim3;
    //    int irred_dim1=nuc_eff*irred_dim2;

    vector<double> U_bare = vector<double>(irred_dim1);
    //the irreducible vertex is approximated by the bare interaction in the parquet approx
public:
    double vval(int, int, int);
    double vvalsmooth(int, int, int);
    double vvalsmooth(int, int, int,double,double,double,char,int,char);
    vector<double> allout(){return U_bare;}
    void copy(vector<double>);
    void setvert(int,int,int,double);
    void direct_set(int,double);

    friend irreducible operator*(double alpha, const irreducible & vertex);
    friend irreducible  operator*(const irreducible & vertex, double alpha);
    friend irreducible  operator+(const irreducible & vertex1, const irreducible & vertex2);
    friend irreducible  abs_sum_tiny(const irreducible & vertex1, const irreducible & vertex2, double tiny);
    double acc(int i){return U_bare[i];}



};
/***************************************************************************************************************************/


/*******************************define "fullvert"FULLVERT" as collection of all diagrams in all channels*********************************/
struct fullvert{//collection of all channels
    irreducible irred;
    svert svertex;
    tvert tvertex;
    uvert uvertex;
public:
template<typename T0>
    double vvalsmooth(int,int,int,T0,T0,T0,char);
    template<typename T0>
    double vvalsmooth(int,int,int,T0,T0,T0,char,int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    template<typename T0>
    double vvalsmooth(int,int,int,int,int,T0,T0,T0,char,int,char);//first two arguments: red_side, map

    friend fullvert operator*(double alpha, const fullvert& vertex);
    friend fullvert operator+(const fullvert& vertex1,const fullvert& vertex2);
    friend fullvert abs_sum_tiny(const fullvert& vertex1,const fullvert& vertex2, double tiny);

};


/******************CLASS FOR SELF ENERGY *************/
class self{
    vector<double > selfenergy =  vector<double >(nw1);
public:
    void setself(int, double);
    vector<double> allout(){return selfenergy;}
    void self_copy(vector<double>);
    double sval(int);
template<typename T0>
    double svalsmooth(T0);

    friend self operator+(const self& self1,const  self& self2);
    friend self operator+=(const self& self1,const  self& self2);
    friend self operator*(double alpha,const  self& self1);
    friend self operator*(const self& self1, double alpha);
    double acc(int i){return selfenergy[i];}
    friend self abs_sum_tiny(const self& self1, const self& self2, double tiny);

};
template<typename T0>
            double self::svalsmooth(T0 w){//smoothly interpolates for values between discrete frequ values of mesh


            double value=0;
            if(abs(w) < ffreqs[nw-1]){
            int W = fconv(w);
        #if temp==0
            value = ((selfenergy[W]*(ffreqs[W+1]-w)+selfenergy[W+1]*(-ffreqs[W]+w))/(ffreqs[W+1]-ffreqs[W]));
#elif temp==1
            value = (selfenergy[W]);
        #endif
};

            return value;
}

/*******PROPAGATOR FUNCTION***********/
            template<typename T0>
double propag(double Lambda, T0 w, self selfenergy,self diffselfenergy, char type);





/*******BUBBLE INTEGRATION FUNCTIONS***********/
//s-bubble:

#if temp==0
template<typename  T1,typename  T2>
struct sbubble_params{


    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;
    self& selfen;
    self& diffselfen;


    int a; int b; int c;

    int d; int e; int f;

    double s;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double sbubble_re_full(double w, void * p){
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);

    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  propag(Lambda,w+s/2,selfen, diffselfen,'g') * propag(Lambda,s/2-w,selfen,diffselfen,'g') ;
    return (1./(2*pi)*val);
}
template<typename T1,typename T2>
double sbubble_re_singsc(double w, void * p){
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);

    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  (propag(Lambda,w+s/2,selfen, diffselfen,'s') * propag(Lambda,s/2-w,selfen,diffselfen,'g')+propag(Lambda,w+s/2,selfen, diffselfen,'g') * propag(Lambda,s/2-w,selfen,diffselfen,'s')) ;
    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double sbubble_re_kat(double w, void * p){
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);

    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  (propag(Lambda,w+s/2,selfen, diffselfen,'k') * propag(Lambda,s/2-w,selfen,diffselfen,'g')+propag(Lambda,w+s/2,selfen, diffselfen,'g') * propag(Lambda,s/2-w,selfen,diffselfen,'k')) ;
    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double sbubble_re_ext(double w, void * p){
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);

    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  (propag(Lambda,w+s/2,selfen, diffselfen,'e') * propag(Lambda,s/2-w,selfen,diffselfen,'g')+propag(Lambda,w+s/2,selfen, diffselfen,'g') * propag(Lambda,s/2-w,selfen,diffselfen,'e')) ;
    return (1./(2*pi)*val);
}

#endif
template<typename T1,typename T2,typename T3>
double sbubble(int red_side,int map1, int map2, gsl_integration_workspace* w, double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, T3 s,T3 w1, T3 w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.


    //Note: factor 1/2 is included due to indistiguishibility of propagators in s-bubble

    double B = 0;


    double vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,w2,'s',1,h);
    double vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,wlimit,'s',2,h);





#if temp==0
    struct sbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};
    gsl_function F;
     double error_re=0;
     double upper_vert=0;
    double lower_vert=0;
    double fmax1=ffreqs[nw1-1];
    double fmax2 = ffreqs[(nw+nw2)/2-1];
    double fmax3= ffreqs[(nw+nw3)/2-1];
    double fmin1=ffreqs[0];
    double fmin2 = ffreqs[(nw-nw2)/2];
    double fmin3= ffreqs[(nw-nw3)/2];
    double bmax1=ffreqs[nw1-1];
    double bmax2 = ffreqs[(nw+nw2)/2-1];
    double bmax3= ffreqs[(nw+nw3)/2-1];
    double bmin1=ffreqs[0];
    double bmin2 = ffreqs[(nw-nw2)/2];
    double bmin3= ffreqs[(nw-nw3)/2];
#endif

#if temp==1
     int upper_vert=0;
    int lower_vert=0;
    int fmax1=nw1-1;
    int fmax2 = nw2-1;
    int fmax3= nw3-1;
    int fmin1=-(nw1-1);
    int fmin2 = -(nw2-1);
    int fmin3= -(nw3-1);
    int bmax1=nw1-2;
    int bmax2 = nw2-2;
    int bmax3= nw3-2;
    int bmin1=-(nw1-2);
    int bmin2 = -(nw2-2);
    int bmin3= -(nw3-2);

#endif


    if(p1=='g' && p2 =='g'){

        if(h=='R'){
            //upper bound:

            auto v1_K1_t_up = bmax1+w2;
            auto v1_K2_t_up = min({bmax2+w2,2*fmax2-w2-s},comp);//choose the minimum since beyond this value, this vertex cannot contribute any more
            auto v1_K2b_t_up = min({bmax2+w2,2*fmax2-w2+s},comp);
            auto v1_R_t_up = min({bmax3+w2,2*fmax3-w2-s,2*fmax3-w2+s},comp);
            auto v1_K1_u_up= bmax1-w2;
            auto v1_K2_u_up = min({bmax2-w2,2*fmax2+w2-s},comp);
            auto v1_K2b_u_up = min({bmax2-w2,2*fmax2+w2+s},comp);
            auto v1_R_u_up = min({bmax3-w2,2*fmax3+w2-s,2*fmax3+w2+s},comp);

            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_up = bmax1+w1;
            auto v2_K2_t_up = min({bmax2+w1,2*fmax2-w1-s},comp);
            auto v2_K2b_t_up = min({bmax2+w1,2*fmax2-w1+s},comp);
            auto v2_R_t_up = min({bmax3+w1,2*fmax3-w1+s,2*fmax3-w1-s},comp);
            auto v2_K1_u_up = bmax1-w1;
            auto v2_K2_u_up = min({bmax2-w1,2*fmax2+w1+s},comp);
            auto v2_K2b_u_up = min({bmax2-w1,2*fmax2+w1-s},comp);
            auto v2_R_u_up = min({bmax3-w1,2*fmax3+w1+s,2*fmax3+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_up = fmax2;
            auto v_R_s_up = fmax3;


            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_t_low = bmin1+w2;
            auto v1_K2_t_low = max({bmin2+w2,2*fmin2-w2-s},comp);//choose the minimum since bezond this value, this vertex cannot contribute any more
            auto v1_K2b_t_low = max({bmin2+w2,2*fmin2-w2+s},comp);
            auto v1_R_t_low = max({bmin3+w2,2*fmin3-w2-s,2*fmin3-w2+s},comp);
            auto v1_K1_u_low = bmin1-w2;
            auto v1_K2_u_low = max({bmin2-w2,2*fmin2+w2-s},comp);
            auto v1_K2b_u_low = max({bmin2-w2,2*fmin2+w2+s},comp);
            auto v1_R_u_low = max({bmin3-w2,2*fmin3+w2-s,2*fmin3+w2+s},comp);

            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_low =  bmin1+w1;
            auto v2_K2_t_low = max({bmin2+w1,2*fmin2-w1-s},comp);
            auto v2_K2b_t_low = max({bmin2+w1,2*fmin2-w1+s},comp);
            auto v2_R_t_low = max({bmin3+w1,2*fmin3-w1+s,2*fmin3-w1-s},comp);
            auto v2_K1_u_low = bmin1-w1;
            auto v2_K2_u_low =  max({bmin2-w1,2*fmin2+w1+s},comp);
            auto v2_K2b_u_low = max({bmin2-w1,2*fmin2+w1-s},comp);
            auto v2_R_u_low = max({bmin3-w1,2*fmin3+w1+s,2*fmin3+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_low = fmin2;
            auto v_R_s_low = fmin3;


            upper_vert = max({v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up},comp)  ;
            lower_vert = min({v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low},comp);

        }
        else if (h=='K'){
            //upper bound:

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_up = fmax2;

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_low = fmin2;

            upper_vert = v_K2_s_up;
            lower_vert = v_K2_s_low ;

        }
        else if(h=='L'){
            //upper bound:


            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_up = bmax1+w1;
            auto v2_K2_t_up = min({bmax2+w1,2*fmax2-w1-s},comp);
            auto v2_K2b_t_up = min({bmax2+w1,2*fmax2-w1+s},comp);
            auto v2_R_t_up = min({bmax3+w1,2*fmax3-w1+s,2*fmax3-w1-s},comp);
            auto v2_K1_u_up = bmax1-w1;
            auto v2_K2_u_up = min({bmax2-w1,2*fmax2+w1+s},comp);
            auto v2_K2b_u_up = min({bmax2-w1,2*fmax2+w1-s},comp);
            auto v2_R_u_up = min({bmax3-w1,2*fmax3+w1+s,2*fmax3+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_up = fmax2;
            auto v_R_s_up = fmax3;


            //lower bound:


            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_low =  bmin1+w1;
            auto v2_K2_t_low = max({bmin2+w1,2*fmin2-w1-s},comp);
            auto v2_K2b_t_low = max({bmin2+w1,2*fmin2-w1+s},comp);
            auto v2_R_t_low = max({bmin3+w1,2*fmin3-w1+s,2*fmin3-w1-s},comp);
            auto v2_K1_u_low = bmin1-w1;
            auto v2_K2_u_low =  max({bmin2-w1,2*fmin2+w1+s},comp);
            auto v2_K2b_u_low = max({bmin2-w1,2*fmin2+w1-s},comp);
            auto v2_R_u_low = max({bmin3-w1,2*fmin3+w1+s,2*fmin3+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_low = fmin2;
            auto v_R_s_low = fmin3;


            upper_vert = max({v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up},comp)  ;
            lower_vert = min({v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low},comp);

        }
        else if(h=='M'){
            //upper bound:

            auto v1_K1_t_up = bmax1+w2;
            auto v1_K2_t_up = min({bmax2+w2,2*fmax2-w2-s},comp);//choose the minimum since beyond this value, this vertex cannot contribute any more
            auto v1_K2b_t_up = min({bmax2+w2,2*fmax2-w2+s},comp);
            auto v1_R_t_up = min({bmax3+w2,2*fmax3-w2-s,2*fmax3-w2+s},comp);
            auto v1_K1_u_up= bmax1-w2;
            auto v1_K2_u_up = min({bmax2-w2,2*fmax2+w2-s},comp);
            auto v1_K2b_u_up = min({bmax2-w2,2*fmax2+w2+s},comp);
            auto v1_R_u_up = min({bmax3-w2,2*fmax3+w2-s,2*fmax3+w2+s},comp);


            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_up = fmax2;
            auto v_R_s_up =fmax3;


            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_t_low = bmin1+w2;
            auto v1_K2_t_low = max({bmin2+w2,2*fmin2-w2-s},comp);//choose the minimum since bezond this value, this vertex cannot contribute any more
            auto v1_K2b_t_low = max({bmin2+w2,2*fmin2-w2+s},comp);
            auto v1_R_t_low = max({bmin3+w2,2*fmin3-w2-s,2*fmin3-w2+s},comp);
            auto v1_K1_u_low = bmin1-w2;
            auto v1_K2_u_low = max({bmin2-w2,2*fmin2+w2-s},comp);
            auto v1_K2b_u_low = max({bmin2-w2,2*fmin2+w2+s},comp);
            auto v1_R_u_low = max({bmin3-w2,2*fmin3+w2-s,2*fmin3+w2+s},comp);


            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_s_low = fmin2;
            auto v_R_s_low = fmin3;


            upper_vert = max({v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_s_up,v_R_s_up},comp)  ;
            lower_vert = min({v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_s_low,v_R_s_low},comp);

        };

        #if temp==0
        F.function = &sbubble_re_full<T1,T2>;
        F.params = &params;


        //NUMERICAL PART:
        //boundaries at which the sharp cutoff has its step plus margin where regulators vanishes
        double upper_reg = abs(s/2) + 7*Lambda;
        double lower_reg = -abs(s/2) - 7*Lambda;

        double upper_self = ffreqs[nw-1]+abs(s)/2;
        double lower_self = ffreqs[0]-abs(s)/2;

        //boundaries which need to be computed numerically:

        double upper = max({upper_reg,upper_self,upper_vert},comp);
        double lower = min({lower_reg,lower_self,lower_vert},comp);


        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        // int status=1;
        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);
        //            if(status !=0){
        //                if( rel_error < 1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "s-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "sbubble-GG not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };

        gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);

        gsl_set_error_handler(old_handler);
        //ANALYTICAL PART:
        B += 1./(2*pi)*1./(2.*s)* ( -log((2*lower+s)/(2*lower-s)) + log((2*upper+s)/(2*upper-s))) * vert1_const * vert2_const;
#elif temp==1

int Lambda_int = static_cast<int>(abs(s/2)+Lambda/(pi*T));
if(Lambda_int%2 == (s/2)%2){Lambda_int -=1;};
if(abs(s/2.)+Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};

if(s%4==0 && upper_vert%2==0){
    upper_vert+=1;
}
else if(s%4!=0 && upper_vert%2!=0){
    upper_vert+=1;
};
if(s%4==0 && lower_vert%2==0){
   lower_vert-=1;
}
else if(s%4!=0 &&  lower_vert%2!=0){
    lower_vert-=1;
};

int selfmax_int = abs(s)/2+nw-1;

int upper = max({Lambda_int,selfmax_int,upper_vert},comp);
int lower = min({-Lambda_int,-selfmax_int,lower_vert},comp);



    for(int i=lower;i<upper+1;i+=2){
        B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,i,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,i,'s',2,h)  *  propag(Lambda,i+s/2,se, dse,'g') * propag(Lambda,s/2-i,se,dse,'g') ;

        B-=  -1./2 * vert1_const * vert2_const  *  1/((i+s/2)*pi*T) *  1/((s/2-i)*pi*T);//for s ==0, this is the summation to infinity on both sides
    };

    if (s==0){
        B+=1/pow(2*T,2)* vert1_const * vert2_const ;};





#endif

    }

    else if(p1=='s' ||p2 =='s'){
#if temp==0
#if reg==1 


        //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
        B = (-1)*1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,s+Lambda,se,dse,p2) ;
        B += (-1)*1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,s-Lambda, se,dse, p2);

        B +=   (-1)*1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2);
        B += (-1)*1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2);




#elif reg==2 || reg==3

        F.function = &sbubble_re_singsc<T1,T2>;
        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();


        double result;
        //  int status=1;
        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,-abs(s/2)-7*Lambda,7*Lambda+abs(s/2),abs_error,rel_error,1500,6,w, &result, &error_re);
        //            if(status ==1){
        //                if( rel_error < 1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "s-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "sbubble1 not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };

        gsl_integration_qag(&F,-abs(s/2)-7*Lambda,7*Lambda+abs(s/2),abs_error,rel_error,1500,6,w, &result, &error_re);


        B=result;
        gsl_set_error_handler(old_handler);

#endif

#elif temp==1
    int Lambda_int = static_cast<int>(Lambda/(pi*T));
    if(Lambda_int%2 == 0){Lambda_int -=1;};
    if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda_int-s/2,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda_int-s/2,'s',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,s+Lambda_int,se,dse,'g');
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda_int-s/2,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda_int-s/2,'s',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,s-Lambda_int,se,dse,'g') ;

    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda_int,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda_int,'s',2,h)  *  propag(Lambda,Lambda_int+s,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda_int,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda_int,'s',2,h)  *  propag(Lambda,-Lambda_int+s,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');
#endif
}

    else if(p1=='k'  || p2 =='k'){

#if temp==0

#if reg==1

        //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
        B =  (-1)*1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,s+Lambda,se,dse,p2) ;
        B += (-1)*1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,s-Lambda, se,dse, p2);

        B +=   (-1)*1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s') ;
        B += (-1)*1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s') ;



        //ADD KATANIN EXTENSION

        F.function = &sbubble_re_ext<T1,T2>;
        F.params = &params;
        double result_katanin=0;
        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        status=1;
        abs_error=1e-3; rel_error=1e-7;
        while(status)
        {
            status=gsl_integration_qag(&F,-abs(s/2)+ffreqs[0],abs(s/2)+ffreqs[nw-1],abs_error,rel_error,1500,6,w, &result_katanin, &error_re);
            if(status !=0){
                if( rel_error < 1e-1){
                    abs_error *= 1.2;
                    rel_error *= 1.2;
                    cout << "s-bubble-katanin: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}

                else {
                    cout << "sbubble-kat not computable with accepatbale precision" << endl;
                    exit (EXIT_FAILURE);
                };
            };
        };
        gsl_set_error_handler(old_handler);
        B += result_katanin;
#elif reg==2 || reg==3


        F.function = &sbubble_re_kat<T1,T2>;

        F.params = &params;
        double result;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        int status=1;


        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,-abs(s/2)-7*Lambda,abs(s/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        //            if(status !=0){
        //                if( rel_error < 1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "s-bubble1: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "sbubble1 not computable with accepatbale precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };
        gsl_integration_qag(&F,-abs(s/2)-7*Lambda,abs(s/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        gsl_set_error_handler(old_handler);

        B=result;



#endif
#elif temp==1
   //single scale contribution
    int Lambda_int = static_cast<int>(Lambda/(pi*T));
    if(Lambda_int%2 == 0){Lambda_int -=1;};
    if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda_int-s/2,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda_int-s/2,'s',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,s+Lambda_int,se,dse,'g');
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda_int-s/2,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda_int-s/2,'s',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,s-Lambda_int,se,dse,'g') ;

    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda_int,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda_int,'s',2,h)  *  propag(Lambda,Lambda_int+s,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
    B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda_int,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda_int,'s',2,h)  *  propag(Lambda,-Lambda_int+s,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');

//katanin extension




    for(int i=-s/2-(nw-1);i<-s/2+nw;i+=2){
        B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,i,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,i,'s',2,h)  *  propag(Lambda,i+s/2,se, dse,'e') * propag(Lambda,s/2-i,se,dse,'g') ;
        };



    for(int i=s/2-(nw-1);i<s/2+nw;i+=2){
        B+= (-1)* 1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,i,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,i,'s',2,h)  *  propag(Lambda,i+s/2,se, dse,'g') * propag(Lambda,s/2-i,se,dse,'e') ;
        };

#endif
    };

    if( abs(B) >1e6){ return 0;};


    return B;
}







//t-bubble:
template<typename  T1,typename  T2>
struct tbubble_params{
    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;

    self& selfen;
    self& diffselfen;

    int a; int b; int c;

    int d; int e; int f;

    double t;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double tbubble_re_full(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1.)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  propag(Lambda,w-t/2,selfen,diffselfen,'g') * propag(Lambda,w+t/2,selfen,diffselfen,'g') ;

    return (-1./(2*pi)*val);//minus sign is from definition of t-bubble
}

template<typename T1,typename T2>
double tbubble_re_singsc(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1.)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  (propag(Lambda,w-t/2,selfen,diffselfen,'s') * propag(Lambda,w+t/2,selfen,diffselfen,'g')+propag(Lambda,w-t/2,selfen,diffselfen,'g') * propag(Lambda,w+t/2,selfen,diffselfen,'s')) ;

    return (-1./(2*pi)*val);//minus sign is from definition of t-bubble
}
template<typename T1,typename T2>
double tbubble_re_kat(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1.)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) * (propag(Lambda,w-t/2,selfen,diffselfen,'k') * propag(Lambda,w+t/2,selfen,diffselfen,'g')+propag(Lambda,w-t/2,selfen,diffselfen,'g') * propag(Lambda,w+t/2,selfen,diffselfen,'k')) ;

    return (-1./(2*pi)*val);//minus sign is from definition of t-bubble
}
template<typename T1,typename T2>
double tbubble_re_ext(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1.)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  (propag(Lambda,w-t/2,selfen,diffselfen,'s') * propag(Lambda,w+t/2,selfen,diffselfen,'e')+propag(Lambda,w-t/2,selfen,diffselfen,'g') * propag(Lambda,w+t/2,selfen,diffselfen,'e')) ;
    return (-1./(2*pi)*val);//minus sign is from definition of t-bubble
}


template<typename T1,typename T2>
double tbubble_re_test(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1.)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  (propag(Lambda,w-t/2,selfen,diffselfen,'k') * propag(Lambda,w+t/2,selfen,diffselfen,'g')+propag(Lambda,w-t/2,selfen,diffselfen,'g') * propag(Lambda,w+t/2,selfen,diffselfen,'k')) ;
    string filename = "testfunc.txt";
    ofstream myfile;
    myfile.open (filename);

    myfile << w<< ' ' << (-1./(2*pi)*val) <<endl;
    myfile.close();

    return (-1./(2*pi)*val);//minus sign is from definition of t-bubble
}



template<typename T1,typename T2,typename T3>
double tbubble(int red_side, int map1, int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, T3 t,T3 w1, T3 w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.





    double B=0;



    double vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,w2,'t',1,h);
    double vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,wlimit,'t',2,h);





#if temp==0
    double error_re=0;
    double upper_vert=0;
    double lower_vert=0;
    struct tbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};
    gsl_function F;

#endif


    if(p1=='g' && p2 =='g'){
#if temp==0
        double fmax1=ffreqs[nw1-1];
        double fmax2 = ffreqs[(nw+nw2)/2-1];
        double fmax3= ffreqs[(nw+nw3)/2-1];
        double fmin1=ffreqs[0];
        double fmin2 = ffreqs[(nw-nw2)/2];
        double fmin3= ffreqs[(nw-nw3)/2];
        double bmax1=ffreqs[nw1-1];
        double bmax2 = ffreqs[(nw+nw2)/2-1];
        double bmax3= ffreqs[(nw+nw3)/2-1];
        double bmin1=ffreqs[0];
        double bmin2 = ffreqs[(nw-nw2)/2];
        double bmin3= ffreqs[(nw-nw3)/2];
#elif temp==1
       int upper_vert=0;
       int lower_vert=0;
       int fmax1=nw1-1;
       int fmax2 = nw2-1;
       int fmax3= nw3-1;
       int fmin1=-(nw1-1);
       int fmin2 = -(nw2-1);
       int fmin3= -(nw3-1);
       int bmax1=nw1-2;
       int bmax2 = nw2-2;
       int bmax3= nw3-2;
       int bmin1=-(nw1-2);
       int bmin2 = -(nw2-2);
       int bmin3= -(nw3-2);
#endif

        if(h=='R'){
            //upper bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_u_up = bmax1+w2;
            auto v1_K2_u_up = min({bmax2+w2,2*fmax2-w2+t},comp);
            auto v1_K2b_u_up = min({bmax2+w2,2*fmax2-w2-t},comp);
            auto v1_R_u_up = min({bmax3+w2,2*fmax3-w2+t,2*fmax3-w2-t},comp);
            auto v1_K1_s_up= bmax1-w2;
            auto v1_K2_s_up = min({bmax2-w2,2*fmax2+w2+t},comp);
            auto v1_K2b_s_up = min({bmax2-w2,2*fmax2+w2-t},comp);
            auto v1_R_s_up = min({bmax3-w2,2*fmax3+w2+t,2*fmax3+w2-t},comp);

            //conditions from vertex 2 (unprimed):
            auto v2_K1_u_up = bmax1+w1;
            auto v2_K2_u_up = min({bmax2+w1,2*fmax2-w1+t},comp);
            auto v2_K2b_u_up = min({bmax2+w1,2*fmax2-w1-t},comp);
            auto v2_R_u_up = min({bmax3+w1,2*fmax3-w1+t,2*fmax3-w1-t},comp);
            auto v2_K1_s_up = bmax1-w1;
            auto v2_K2_s_up = min({bmax2-w1,2*fmax2+w1-t},comp);
            auto v2_K2b_s_up = min({bmax2-w1,2*fmax2+w1+t},comp);
            auto v2_R_s_up = min({bmax3-w1,2*fmax3+w1+t,2*fmax3+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_up = fmax2;
            auto v_R_t_up = fmax3;



            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_u_low = bmin1+w2;
            auto v1_K2_u_low = max({bmin2+w2,2*fmin2-w2+t},comp);
            auto v1_K2b_u_low =  max({bmin2+w2,2*fmin2-w2-t},comp);
            auto v1_R_u_low =  max({bmin3+w2,2*fmin3-w2+t,2*fmin3-w2-t},comp);
            auto v1_K1_s_low = bmin1-w2;
            auto v1_K2_s_low =  max({bmin2-w2,2*fmin2+w2+t},comp);
            auto v1_K2b_s_low = max({bmin2-w2,2*fmin2+w2-t},comp);
            auto v1_R_s_low =  max({bmin3-w2,2*fmin3+w2+t,2*fmin3+w2-t},comp);

            //conditions from vertex 2 (unprimed):
            auto v2_K1_u_low = bmin1+w1;
            auto v2_K2_u_low = max({bmin2+w1,2*fmin2-w1+t},comp);
            auto v2_K2b_u_low = max({bmin2+w1,2*fmin2-w1-t},comp);
            auto v2_R_u_low = max({bmin3+w1,2*fmin3-w1+t,2*fmin3-w1-t},comp);
            auto v2_K1_s_low =  bmin1-w1;
            auto v2_K2_s_low = max({bmin2-w1,2*fmin2+w1-t},comp);
            auto v2_K2b_s_low = max({bmin2-w1,2*fmin2+w1+t},comp);
            auto v2_R_s_low =max({bmin3-w1,2*fmin3+w1+t,2*fmin3+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_low = fmin2;
            auto v_R_t_low = fmin3;

            upper_vert = max({v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up},comp)  ;
            lower_vert = min({v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low},comp);
        }

        else if(h=='K'){
            //upper bound:

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_up = fmax2;



            //lower bound:

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_low = fmin2;


            upper_vert = v_K2_t_up;
            lower_vert = v_K2_t_low;
        }
        else if(h=='L'){
            //upper bound:


            //conditions from vertex 2 (unprimed):
            auto v2_K1_u_up = bmax1+w1;
            auto v2_K2_u_up = min({bmax2+w1,2*fmax2-w1+t},comp);
            auto v2_K2b_u_up = min({bmax2+w1,2*fmax2-w1-t},comp);
            auto v2_R_u_up = min({bmax3+w1,2*fmax3-w1+t,2*fmax3-w1-t},comp);
            auto v2_K1_s_up = bmax1-w1;
            auto v2_K2_s_up = min({bmax2-w1,2*fmax2+w1-t},comp);
            auto v2_K2b_s_up = min({bmax2-w1,2*fmax2+w1+t},comp);
            auto v2_R_s_up = min({bmax3-w1,2*fmax3+w1+t,2*fmax3+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_up = fmax2;
            auto v_R_t_up = fmax3;

            //lower bound:

            //conditions from vertex 2 (unprimed):
            auto v2_K1_u_low = bmin1+w1;
            auto v2_K2_u_low = max({bmin2+w1,2*fmin2-w1+t},comp);
            auto v2_K2b_u_low = max({bmin2+w1,2*fmin2-w1-t},comp);
            auto v2_R_u_low = max({bmin3+w1,2*fmin3-w1+t,2*fmin3-w1-t},comp);
            auto v2_K1_s_low =  bmin1-w1;
            auto v2_K2_s_low = max({bmin2-w1,2*fmin2+w1-t},comp);
            auto v2_K2b_s_low = max({bmin2-w1,2*fmin2+w1+t},comp);
            auto v2_R_s_low =max({bmin3-w1,2*fmin3+w1+t,2*fmin3+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_low = fmin2;
            auto v_R_t_low = fmin3;

            upper_vert = max({v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up},comp)  ;
            lower_vert = min({v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low},comp);
        }

        else if(h=='M'){
            //upper bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_u_up = bmax1+w2;
            auto v1_K2_u_up = min({bmax2+w2,2*fmax2-w2+t},comp);
            auto v1_K2b_u_up = min({bmax2+w2,2*fmax2-w2-t},comp);
            auto v1_R_u_up = min({bmax3+w2,2*fmax3-w2+t,2*fmax3-w2-t},comp);
            auto v1_K1_s_up= bmax1-w2;
            auto v1_K2_s_up = min({bmax2-w2,2*fmax2+w2+t},comp);
            auto v1_K2b_s_up = min({bmax2-w2,2*fmax2+w2-t},comp);
            auto v1_R_s_up = min({bmax3-w2,2*fmax3+w2+t,2*fmax3+w2-t},comp);


            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_up = fmax2;
            auto v_R_t_up = fmax3;



            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_u_low = bmin1+w2;
            auto v1_K2_u_low = max({bmin2+w2,2*fmin2-w2+t},comp);
            auto v1_K2b_u_low =  max({bmin2+w2,2*fmin2-w2-t},comp);
            auto v1_R_u_low =  max({bmin3+w2,2*fmin3-w2+t,2*fmin3-w2-t},comp);
            auto v1_K1_s_low = bmin1-w2;
            auto v1_K2_s_low =  max({bmin2-w2,2*fmin2+w2+t},comp);
            auto v1_K2b_s_low = max({bmin2-w2,2*fmin2+w2-t},comp);
            auto v1_R_s_low =  max({bmin3-w2,2*fmin3+w2+t,2*fmin3+w2-t},comp);


            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            auto v_K2_t_low = fmin2;
            auto v_R_t_low = fmin3;

            upper_vert = max({v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_t_up,v_R_t_up},comp)  ;
            lower_vert = min({v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_t_low,v_R_t_low},comp);
        };
#if temp==0
        F.function = &tbubble_re_full<T1,T2>;
        F.params = &params;


        //NUMERICAL PART:
        //boundaries at which the sharp cutoff has its step plus margin where regulators vanishes
        double upper_reg = abs(t/2) + 7*Lambda;
        double lower_reg = -abs(t/2) - 7*Lambda;

        double upper_self = ffreqs[nw-1]+abs(t)/2;
        double lower_self = ffreqs[0]-abs(t)/2;

        //boundaries which need to be computed numerically:

        double upper = max({upper_reg,upper_self,upper_vert},comp);
        double lower = min({lower_reg,lower_self,lower_vert},comp);


        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //   int status=1;
        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);
        //            if(status !=0){
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "t-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "tbubble-GG not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };

        gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);
        gsl_set_error_handler(old_handler);

        //ANALYTICAL PART:
        B += 1./(2*pi)*1./t *( -log((2*lower+t)/(2*lower-t)) + log((2*upper+t)/(2*upper-t))) * vert1_const * vert2_const;

#elif temp==1
        int Lambda_int = static_cast<int>(abs(t/2)+Lambda/(pi*T));
        if(Lambda_int%2 == (t/2)%2){Lambda_int -=1;};
        if(abs(t/2.)+Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};


        if(t%4==0 && upper_vert%2==0){
            upper_vert+=1;
        }
        else if(t%4!=0 && upper_vert%2!=0){
            upper_vert+=1;
        };
        if(t%4==0 && lower_vert%2==0){
           lower_vert-=1;
        }
        else if(t%4!=0 &&  lower_vert%2!=0){
            lower_vert-=1;
        };

       int selfmax_int = abs(t)/2+nw-1;

        int upper = max({Lambda_int,selfmax_int,upper_vert},comp);
        int lower = min({-Lambda_int,-selfmax_int,lower_vert},comp);



            for(int i=lower;i<upper+1;i+=2){
                B+=(-1)* (-1)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,i,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,i,'t',2,h)  *  propag(Lambda,i-t/2,se, dse,'g') * propag(Lambda,t/2+i,se,dse,'g') ;
  
                B-=  (-1)* (-1)* vert1_const * vert2_const  *  1/((i-t/2)*pi*T) *  1/((t/2+i)*pi*T);
            };
            if (t==0){
                B+=1/pow(2*T,2)* vert1_const * vert2_const ;};
        
        
#endif
    }

    else if(p1=='s'|| p2 =='s'){
#if temp==0
#if reg==1

        //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
        B = (-1)* (-1.) *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) * propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,'g') ;
        B += (-1)*(-1.) *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) * propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,'g') ;

        B += (-1)*(-1.) * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) * propag(Lambda,-Lambda-t,se,dse,'g') * propag(Lambda,-Lambda,se,dse,'s') ;
        B += (-1)*(-1.) *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) * propag(Lambda,Lambda-t,se,dse,'g') * propag(Lambda,Lambda,se,dse,'s') ;


#elif reg==2 || reg==3

        F.function = &tbubble_re_singsc<T1,T2>;
        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //     int status=1;
        double result;
        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {

        //            status=gsl_integration_qag(&F,-abs(t/2)-7*Lambda,abs(t/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        //            if (status !=0){
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "t-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "tbubble1 not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };
        gsl_integration_qag(&F,-abs(t/2)-7*Lambda,abs(t/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        gsl_set_error_handler(old_handler);
        B=result;




#endif
#elif temp==1

int Lambda_int = static_cast<int>(Lambda/(pi*T));
if(Lambda_int%2 == 0){Lambda_int -=1;};
if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};
B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda_int+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda_int+t/2,'t',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,t-Lambda_int,se,dse,'g');
B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda_int+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda_int+t/2,'t',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,t+Lambda_int,se,dse,'g') ;

B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-t/2-Lambda_int,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-t/2-Lambda_int,'t',2,h)  *  propag(Lambda,Lambda_int-t,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-t/2+Lambda_int,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-t/2+Lambda_int,'t',2,h)  *  propag(Lambda,-Lambda_int-t,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');


#endif
    }

    else if(p1=='k' || p2 =='k'){

#if temp==0

#if reg==1

        //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
        B =   (-1)* (-1.) *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) * propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2) ;
        B += (-1)*(-1.) *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) * propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2) ;

        B +=  (-1)*(-1.) * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) * propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s') ;
        B += (-1)*(-1.) *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) * propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s') ;

        //ADD KATANIN EXTENSION

        F.function = &tbubble_re_ext<T1,T2>;
        F.params = &params;
        double result_katanin=0;
        abs_error=1e-3; rel_error=1e-7;
        status=1;
        while(status)
        {
            status=gsl_integration_qag(&F,-abs(t/2)+ffreqs[0],abs(t/2)+ffreqs[nw-1],abs_error,rel_error,1500,6,w, &result_katanin, &error_re);
            if(status !=0){
                if(rel_error<1e-1){
                    abs_error *= 1.2;
                    rel_error *= 1.2;
                    cout << "t-bubble-katanin: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
                else {
                    cout << "tbubble-kat not computable with acceptable precision" << endl;
                    exit (EXIT_FAILURE);
                };
            };
        };
        gsl_set_error_handler(old_handler);
        B += result_katanin;
#elif reg==2 || reg==3


        F.function = &tbubble_re_kat<T1,T2>;

        F.params = &params;
        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //    int status=1;

        double abs_error=1e-3, rel_error=1e-7;
        double result;

        //        while(status)
        //        {

        //            status=gsl_integration_qag(&F,-abs(t/2)-7*Lambda,abs(t/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        //         if(status != 0){
        //               F.function = &tbubble_re_test<T1,T2>;
        //               status=gsl_integration_qag(&F,-abs(t/2)-7*Lambda,abs(t/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "t-bubble1: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "tbubble1 not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //      };
        gsl_integration_qag(&F,-abs(t/2)-7*Lambda,abs(t/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        gsl_set_error_handler(old_handler);
        B=result;


#endif

#elif temp==1
    //single scale contribution
              int Lambda_int = static_cast<int>(Lambda/(pi*T));
              if(Lambda_int%2 == 0){Lambda_int -=1;};
              if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};
              B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda_int+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda_int+t/2,'t',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,t-Lambda_int,se,dse,'g');
              B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda_int+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda_int+t/2,'t',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,t+Lambda_int,se,dse,'g') ;

              B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-t/2-Lambda_int,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-t/2-Lambda_int,'t',2,h)  *  propag(Lambda,Lambda_int-t,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
              B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,-t/2+Lambda_int,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-t/2+Lambda_int,'t',2,h)  *  propag(Lambda,-Lambda_int-t,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');


//katanin extension

    for(int i=t/2-(nw-1);i<t/2+nw;i+=2){
        B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,i,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,i,'t',2,h)  *  propag(Lambda,i-t/2,se, dse,'e') * propag(Lambda,t/2+i,se,dse,'g') ;
        };



    for(int i=-t/2-(nw-1);i<-t/2+nw;i+=2){
        B+= (-1)* (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,t,i,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,i,'t',2,h)  *  propag(Lambda,i-t/2,se, dse,'g') * propag(Lambda,t/2+i,se,dse,'e') ;
        };

#endif
    };

    if( abs(B) >1e6){ return 0;};
    return B;
}


//u-bubble:
template<typename  T1,typename  T2>
struct ububble_params{
    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;

    self& selfen;
    self& diffselfen;

    int a; int b; int c;

    int d; int e; int f;

    double u;double w1; double w2;
    char h;//class of bubble (K = K1, L=K2, M = K2b, R = R)

};

template<typename T1,typename T2>
double ububble_re_full(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double u = (params->u);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)*vert1.vvalsmooth(red_side,map1,a,b,c,u,w,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,w,'u',2,h) * propag(Lambda,w-u/2,selfen,diffselfen,'g') * propag(Lambda,w+u/2,selfen,diffselfen,'g') ;
    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double ububble_re_singsc(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double u = (params->u);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)*vert1.vvalsmooth(red_side,map1,a,b,c,u,w,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,w,'u',2,h) * ( propag(Lambda,w-u/2,selfen,diffselfen,'g') * propag(Lambda,w+u/2,selfen,diffselfen,'s') + propag(Lambda,w-u/2,selfen,diffselfen,'s') * propag(Lambda,w+u/2,selfen,diffselfen,'g') );
    return (1./(2*pi)*val);
}
template<typename T1,typename T2>
double ububble_re_kat(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double u = (params->u);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)*vert1.vvalsmooth(red_side,map1,a,b,c,u,w,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,w,'u',2,h) * ( propag(Lambda,w-u/2,selfen,diffselfen,'g') * propag(Lambda,w+u/2,selfen,diffselfen,'k') + propag(Lambda,w-u/2,selfen,diffselfen,'k') * propag(Lambda,w+u/2,selfen,diffselfen,'g') );
    return (1./(2*pi)*val);
}
template<typename T1,typename T2>
double ububble_re_ext(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);

    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double u = (params->u);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
    double val = (-1)*vert1.vvalsmooth(red_side,map1,a,b,c,u,w,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,w,'u',2,h)   * ( propag(Lambda,w-u/2,selfen,diffselfen,'g') * propag(Lambda,w+u/2,selfen,diffselfen,'e') + propag(Lambda,w-u/2,selfen,diffselfen,'e') * propag(Lambda,w+u/2,selfen,diffselfen,'g') );
    return (1./(2*pi)*val);
}




template<typename T1,typename T2,typename T3>
double ububble(int red_side,int map1,int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, T3 u,T3 w1, T3 w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.



    double B=0;






    double vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,w2,'u',1,h);
    double vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,wlimit,'u',2,h);
#if temp==0
    struct ububble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};
    gsl_function F;
    double error_re=0;

#endif

    if(p1=='g' && p2 =='g'){

 
#if temp==0
  double upper_vert=0;
    double lower_vert=0;
        double fmax1=ffreqs[nw1-1];
        double fmax2 = ffreqs[(nw+nw2)/2-1];
        double fmax3= ffreqs[(nw+nw3)/2-1];
        double fmin1=ffreqs[0];
        double fmin2 = ffreqs[(nw-nw2)/2];
        double fmin3= ffreqs[(nw-nw3)/2];
        double bmax1=ffreqs[nw1-1];
        double bmax2 = ffreqs[(nw+nw2)/2-1];
        double bmax3= ffreqs[(nw+nw3)/2-1];
        double bmin1=ffreqs[0];
        double bmin2 = ffreqs[(nw-nw2)/2];
        double bmin3= ffreqs[(nw-nw3)/2];
#elif temp==1
       int upper_vert=0;
       int lower_vert=0;
       int fmax1=nw1-1;
       int fmax2 = nw2-1;
       int fmax3= nw3-1;
       int fmin1=-(nw1-1);
       int fmin2 = -(nw2-1);
       int fmin3= -(nw3-1);
       int bmax1=nw1-2;
       int bmax2 = nw2-2;
       int bmax3= nw3-2;
       int bmin1=-(nw1-2);
       int bmin2 = -(nw2-2);
       int bmin3= -(nw3-2);
#endif
        if(h=='R'){
            //upper bound:

            //conditions from vertex 1 (primed):
            auto v1_K1_t_up = bmax1+w2;
            auto v1_K2_t_up = min({bmax2+w2,2*fmax2-w2+u},comp);
            auto v1_K2b_t_up = min({bmax2+w2,2*fmax2-w2-u},comp);
            auto v1_R_t_up = min({bmax3+w2,2*fmax3-w2+u,2*fmax3-w2-u},comp);;
            auto v1_K1_s_up= bmax1-w2;
            auto v1_K2_s_up = min({bmax2-w2,2*fmax2+w2+u},comp);
            auto v1_K2b_s_up =  min({bmax2-w2,2*fmax2+w2-u},comp);
            auto v1_R_s_up = min({bmax3-w2,2*fmax3+w2-u,2*fmax3+w2+u},comp);

            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_up = bmax1+w1;
            auto v2_K2_t_up = min({bmax2+w1,2*fmax2-w1+u},comp);
            auto v2_K2b_t_up = min({bmax2+w1,2*fmax2-w1-u},comp);
            auto v2_R_t_up = min({bmax3+w1,2*fmax3-w1-u,2*fmax3-w1+u},comp);
            auto v2_K1_s_up = bmax1-w1;
            auto v2_K2_s_up = min({bmax2-w1,2*fmax2+w1-u},comp);
            auto v2_K2b_s_up = min({bmax2-w1,2*fmax2+w1+u},comp);
            auto v2_R_s_up = min({bmax3-w1,2*fmax3+w1+u,2*fmax3+w1-u},comp);;
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            auto v_K2_u_up = fmax2;
            auto v_R_u_up = fmax3;



            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_t_low = bmin1+w2;
            auto v1_K2_t_low = max({bmin2+w2,2*fmin2-w2+u},comp);
            auto v1_K2b_t_low =max({bmin2+w2,2*fmin2-w2-u},comp);
            auto v1_R_t_low = max({bmin3+w2,2*fmin3-w2+u,2*fmin3-w2-u},comp);;
            auto v1_K1_s_low = bmin1-w2;
            auto v1_K2_s_low = max({bmin2-w2,2*fmin2+w2+u},comp);
            auto v1_K2b_s_low = max({bmin2-w2,2*fmin2+w2-u},comp);
            auto v1_R_s_low = max({bmin3-w2,2*fmin3+w2-u,2*fmin3+w2+u},comp);


            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_low = bmin1+w1;
            auto v2_K2_t_low = max({bmin2+w1,2*fmin2-w1+u},comp);
            auto v2_K2b_t_low = max({bmin2+w1,2*fmin2-w1-u},comp);
            auto v2_R_t_low =  max({bmin3+w1,2*fmin3-w1-u,2*fmin3-w1+u},comp);
            auto v2_K1_s_low = bmin1-w1;
            auto v2_K2_s_low =  max({bmin2-w1,2*fmin2+w1-u},comp);
            auto v2_K2b_s_low = max({bmin2-w1,2*fmin2+w1+u},comp);
            auto v2_R_s_low =  max({bmin3-w1,2*fmin3+w1+u,2*fmin3+w1-u},comp);;
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_u_low = fmin2;
            auto v_R_u_low = fmin3;

            upper_vert = max({v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v_K2_u_up,v_R_u_up},comp)  ;
            lower_vert = min({v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v_K2_u_low,v_R_u_low},comp);
        }
        else if(h=='K'){
            //upper bound:

            //conditions from vertex 1 (primed):

            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            auto v_K2_u_up = fmax2;

            //lower bound:

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_u_low = fmin2;


            upper_vert =v_K2_u_up;
            lower_vert =v_K2_u_low;
        }
        else if(h=='L'){
            //upper bound:



            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_up = bmax1+w1;
            auto v2_K2_t_up = min({bmax2+w1,2*fmax2-w1+u},comp);
            auto v2_K2b_t_up = min({bmax2+w1,2*fmax2-w1-u},comp);
            auto v2_R_t_up = min({bmax3+w1,2*fmax3-w1-u,2*fmax3-w1+u},comp);
            auto v2_K1_s_up = bmax1-w1;
            auto v2_K2_s_up = min({bmax2-w1,2*fmax2+w1-u},comp);
            auto v2_K2b_s_up = min({bmax2-w1,2*fmax2+w1+u},comp);
            auto v2_R_s_up = min({bmax3-w1,2*fmax3+w1+u,2*fmax3+w1-u},comp);;
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            auto v_K2_u_up = fmax2;
            auto v_R_u_up = fmax3;



            //lower bound:

            //conditions from vertex 2 (unprimed):
            auto v2_K1_t_low = bmin1+w1;
            auto v2_K2_t_low = max({bmin2+w1,2*fmin2-w1+u},comp);
            auto v2_K2b_t_low = max({bmin2+w1,2*fmin2-w1-u},comp);
            auto v2_R_t_low =  max({bmin3+w1,2*fmin3-w1-u,2*fmin3-w1+u},comp);
            auto v2_K1_s_low = bmin1-w1;
            auto v2_K2_s_low =  max({bmin2-w1,2*fmin2+w1-u},comp);
            auto v2_K2b_s_low = max({bmin2-w1,2*fmin2+w1+u},comp);
            auto v2_R_s_low =  max({bmin3-w1,2*fmin3+w1+u,2*fmin3+w1-u},comp);;
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_u_low = fmin2;
            auto v_R_u_low = fmin3;

            upper_vert = max({v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v_K2_u_up,v_R_u_up},comp)  ;
            lower_vert = min({v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v_K2_u_low,v_R_u_low},comp);
        }
        else if(h=='M'){
            //upper bound:

            //conditions from vertex 1 (primed):
            auto v1_K1_t_up = bmax1+w2;
            auto v1_K2_t_up = min({bmax2+w2,2*fmax2-w2+u},comp);
            auto v1_K2b_t_up = min({bmax2+w2,2*fmax2-w2-u},comp);
            auto v1_R_t_up = min({bmax3+w2,2*fmax3-w2+u,2*fmax3-w2-u},comp);;
            auto v1_K1_s_up= bmax1-w2;
            auto v1_K2_s_up = min({bmax2-w2,2*fmax2+w2+u},comp);
            auto v1_K2b_s_up =  min({bmax2-w2,2*fmax2+w2-u},comp);
            auto v1_R_s_up = min({bmax3-w2,2*fmax3+w2-u,2*fmax3+w2+u},comp);


            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            auto v_K2_u_up = fmax2;
            auto v_R_u_up = fmax3;



            //lower bound:
            //conditions from vertex 1 (primed):
            auto v1_K1_t_low = bmin1+w2;
            auto v1_K2_t_low = max({bmin2+w2,2*fmin2-w2+u},comp);
            auto v1_K2b_t_low =max({bmin2+w2,2*fmin2-w2-u},comp);
            auto v1_R_t_low = max({bmin3+w2,2*fmin3-w2+u,2*fmin3-w2-u},comp);;
            auto v1_K1_s_low = bmin1-w2;
            auto v1_K2_s_low = max({bmin2-w2,2*fmin2+w2+u},comp);
            auto v1_K2b_s_low = max({bmin2-w2,2*fmin2+w2-u},comp);
            auto v1_R_s_low = max({bmin3-w2,2*fmin3+w2-u,2*fmin3+w2+u},comp);



            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            auto v_K2_u_low = fmin2;
            auto v_R_u_low = fmin3;

            upper_vert = max({v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v_K2_u_up,v_R_u_up},comp)  ;
            lower_vert = min({v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v_K2_u_low,v_R_u_low},comp);
        }
        #if temp==0
        F.function = &ububble_re_full<T1,T2>;
        F.params = &params;


        //NUMERICAL PART:
        //boundaries at which the sharp cutoff has its step plus margin where regulators vanishes
        double upper_reg = abs(u/2) + 7*Lambda;
        double lower_reg = -abs(u/2) - 7*Lambda;

        double upper_self = ffreqs[nw-1]+abs(u)/2;
        double lower_self = ffreqs[0]-abs(u)/2;

        //boundaries which need to be computed numerically:

        double upper = max({upper_reg,upper_self,upper_vert},comp);
        double lower = min({lower_reg,lower_self,lower_vert},comp);

        double abs_error=1e-3, rel_error=1e-7;
        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        // int status=1;
        //        while(status)
        //        {

        //            status=gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);

        //            if(status !=0){
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "u-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "ububble_GG not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };
        //            };
        //        };

        gsl_integration_qag(&F,lower,upper,abs_error,rel_error,1500,6,w, &B, &error_re);
        gsl_set_error_handler(old_handler);

        //ANALYTICAL PART:
        B += (-1.)/(2*pi)*1./u *( -log((2*lower+u)/(2*lower-u)) + log((2*upper+u)/(2*upper-u))) * vert1_const * vert2_const;
#elif temp==1
    int Lambda_int = static_cast<int>(abs(u/2)+Lambda/(pi*T));
    if(Lambda_int%2 == (u/2)%2){Lambda_int -=1;};
    if(abs(u/2.)+Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};


    int selfmax_int = abs(u)/2+nw-1;


   if(u%4==0 && upper_vert%2==0){
            upper_vert+=1;
        }
        else if(u%4!=0 && upper_vert%2!=0){
            upper_vert+=1;
        };
        if(u%4==0 && lower_vert%2==0){
           lower_vert-=1;
        }
        else if(u%4!=0 &&  lower_vert%2!=0){
            lower_vert-=1;
        };

     

        int upper = max({Lambda_int,selfmax_int,upper_vert},comp);
        int lower = min({-Lambda_int,-selfmax_int,lower_vert},comp);



        for(int i=lower;i<upper+1;i+=2){

            B+= (-1)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,i,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,i,'u',2,h)  *  propag(Lambda,i-u/2,se, dse,'g') * propag(Lambda,u/2+i,se,dse,'g') ;
  
         B-=   (-1)* vert1_const * vert2_const  *  1/((i-u/2)*pi*T) *  1/((u/2+i)*pi*T);
        };
        if (u==0){
                B+=-1/pow(2*T,2)* vert1_const * vert2_const ;};


#endif

    }

    else if(p1=='s' || p2 =='s'){
#if temp==0
#if reg==1

        B = -1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h)* propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,-Lambda+u,se,dse,p2) ;
        B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,Lambda+u,se,dse,p2) ;


        B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2) ;
        B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2) ;


#elif reg==2 || reg==3

        F.function = &ububble_re_singsc<T1,T2>;
        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //   int status=1;
        double abs_error=1e-3, rel_error=1e-7;
        double result;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,-abs(u/2)-7*Lambda,abs(u/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);

        //            if(status !=0){
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "u-bubble: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "ububble1 not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };

        //            };
        //        };

        gsl_integration_qag(&F,-abs(u/2)-7*Lambda,abs(u/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);

        gsl_set_error_handler(old_handler);
        B= result;



#endif
#elif temp==1

    int Lambda_int = static_cast<int>(Lambda/(pi*T));
    if(Lambda_int%2 == 0){Lambda_int -=1;};
    if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};
    B+= (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda_int+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda_int+u/2,'u',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,u-Lambda_int,se,dse,'g');
    B+=  (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda_int+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda_int+u/2,'u',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,u+Lambda_int,se,dse,'g') ;

    B+=  (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-u/2-Lambda_int,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-u/2-Lambda_int,'u',2,h)  *  propag(Lambda,Lambda_int-u,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
    B+=  (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-u/2+Lambda_int,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-u/2+Lambda_int,'u',2,h)  *  propag(Lambda,-Lambda_int-u,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');

#endif
    }

    else if(p1=='k' || p2 =='k'){
#if temp==0

#if reg==1
        if(p1=='k'){
            //factor (-1) since two (physically) purely imaginary propagators are multiplied which are implemented as (real) double
            B = -1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h)* propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,-Lambda+u,se,dse,p2) ;
            B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,Lambda+u,se,dse,p2) ;

            B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2) ;
            B += -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2) ;

            //ADD KATANIN EXTENSION

            F.function = &ububble_re_ext<T1,T2>;
            F.params = &params;
            double result_katanin;
            gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
            status=1;
            abs_error=1e-3; rel_error=1e-7;
            while(status)
            {
                status=gsl_integration_qag(&F,-abs(u/2)+ffreqs[0],abs(u/2)+ffreqs[nw-1],abs_error,rel_error,1500,6,w, &result_katanin, &error_re);
                if(status !=0){
                    if(rel_error<1e-1){
                        abs_error *= 1.2;
                        rel_error *= 1.2;
                        cout << "u-bubble_katanin: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
                    else {
                        cout << "ububble_kat not computable with acceptable precision" << endl;
                        exit (EXIT_FAILURE);
                    };
                };

            };
            B += result_katanin;
            gsl_set_error_handler(old_handler);
#elif reg==2 || reg==3

        F.function = &ububble_re_kat<T1,T2>;


        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        double result;
        // int status=1;
        double abs_error=1e-3, rel_error=1e-7;
        //        while(status)
        //        {
        //            status=gsl_integration_qag(&F,-abs(u/2)-7*Lambda,abs(u/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);


        //            if(status !=0){
        //                if(rel_error<1e-1){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "u-bubble1: errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;}
        //                else {
        //                    cout << "ububble1 not computable with acceptable precision" << endl;
        //                    exit (EXIT_FAILURE);
        //                };

        //            };
        //        };

        gsl_integration_qag(&F,-abs(u/2)-7*Lambda,abs(u/2)+7*Lambda,abs_error,rel_error,1500,6,w, &result, &error_re);
        gsl_set_error_handler(old_handler);

        B=result;



#endif
#elif temp==1
    //single scale contribution
              int Lambda_int = static_cast<int>(Lambda/(pi*T));
              if(Lambda_int%2 == 0){Lambda_int -=1;};
              if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};

              B+= (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda_int+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda_int+u/2,'u',2,h)  *  propag(Lambda,-Lambda_int,se, dse,'s') * propag(Lambda,u-Lambda_int,se,dse,'g');
              B+= (-1)* vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda_int+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda_int+u/2,'u',2,h)  *  propag(Lambda,Lambda_int,se, dse,'s') * propag(Lambda,u+Lambda_int,se,dse,'g') ;

              B+=  (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-u/2-Lambda_int,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-u/2-Lambda_int,'u',2,h)  *  propag(Lambda,Lambda_int-u,se, dse,'g') * propag(Lambda,-Lambda_int,se,dse,'s');
              B+= (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,-u/2+Lambda_int,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-u/2+Lambda_int,'u',2,h)  *  propag(Lambda,-Lambda_int-u,se, dse,'g') * propag(Lambda,Lambda_int,se,dse,'s');


//katanin extension




    for(int i=u/2-(nw-1);i<u/2+nw;i+=2){

        B+=  (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,i,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,i,'u',2,h)  *  propag(Lambda,i-u/2,se, dse,'e') * propag(Lambda,u/2+i,se,dse,'g') ;
        };



    for(int i=-u/2-(nw-1);i<-u/2+nw;i+=2){
        B+= (-1) * vert1.vvalsmooth(red_side,map1,a,b,c,u,i,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,i,'u',2,h)  *  propag(Lambda,i-u/2,se, dse,'g') * propag(Lambda,u/2+i,se,dse,'e') ;
        };
#endif

    };

    if( abs(B) >1e6){ return 0;};


    return B;
}



//define parvert as touple of spin and density vertex
template<class T>
struct parvert{//define a touple for parametrized vertices that contains one spin vertex and one density vertex
    T spinvertex;
    T densvertex;
};




//define operators for parvert
parvert<svert> operator+(parvert<svert> ,parvert<svert>);
parvert<tvert> operator+(parvert<tvert> ,parvert<tvert>);
parvert<uvert> operator+(parvert<uvert> ,parvert<uvert>);
parvert<svert> abs_sum_tiny(parvert<svert> ,parvert<svert>, double tiny);
parvert<tvert> abs_sum_tiny(parvert<tvert> ,parvert<tvert>, double tiny);
parvert<uvert> abs_sum_tiny(parvert<uvert> ,parvert<uvert>, double tiny);
parvert<irreducible> operator+(parvert<irreducible> ,parvert<irreducible>);
parvert<irreducible> abs_sum_tiny(parvert<irreducible> ,parvert<irreducible>, double tiny);
parvert<svert> operator+=(parvert<svert> ,parvert<svert>);
parvert<tvert> operator+=(parvert<tvert> ,parvert<tvert>);
parvert<uvert> operator+=(parvert<uvert> ,parvert<uvert>);
parvert<fullvert> operator+(parvert<fullvert> ,parvert<fullvert>);
parvert<fullvert> abs_sum_tiny(parvert<fullvert> ,parvert<fullvert>, double tiny);
parvert<irreducible> operator+=(parvert<irreducible> ,parvert<irreducible>);
parvert<svert> operator*(double  ,parvert<svert>&);
parvert<svert> operator*(parvert<svert>& ,double );
parvert<tvert> operator*(double  ,parvert<tvert>&);
parvert<tvert> operator*(parvert<tvert>& ,double );
parvert<uvert> operator*(double  ,parvert<uvert>&);
parvert<uvert> operator*(parvert<uvert>& ,double );
parvert<irreducible> operator*(double  ,parvert<irreducible>&);
parvert<irreducible> operator*(parvert<irreducible>& ,double );


#if temp==0
////loop function:
template<typename  T>
struct loop_params{
    double Lambda;
    T& vert;
    char ptype;
    self& selfen;
    self& diffselfen;

    int a; int b; int c;


    double w1;
    int red;

};


//}

template<typename T0, typename T1>
double loop_im(T1 w, void * p){
    struct loop_params<T0> * params
            = static_cast< struct loop_params<T0> *>(p);
    char ptype = (params ->ptype);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T0& vert = (params-> vert);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);
    double w1 = (params->w1);
    double red = (params->red);
    //ATTENTION: The result is a double but it represents the IMAGINARY PART of the loop function (which has no real part due to the symmetries of the system)
    return(1./(2*pi)* (-1.) * (vert.vvalsmooth(red,0,a,b,c,0.,w1,w,'t',1,'M') )*propag(Lambda,w,selfen,diffselfen,ptype));
}

/**same loop function but with changed arguments in vertex (needed to parvert-loops), see SB1, p.76*/


template<typename T0, typename T1>
double loop_im_sw(T1 w, void * p){
    struct loop_params<T0> * params
            = static_cast< struct loop_params<T0> *>(p);
    char ptype = (params ->ptype);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T0& vert = (params-> vert);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);
    double w1 = (params->w1);
    double red = (params->red);
    //ATTENTION: The result is a double but it represents the IMAGINARY PART of the loop function (which has no real part due to the symmetries of the system)
    return(1./(2*pi)*  (-1.) * (vert.vvalsmooth(red,0,a,b,c,0.,w1,w,'u',1,'M'))*propag(Lambda,w,selfen,diffselfen,ptype));

}



#endif


template<class T0, typename T1>
//loop Attention: Only use for freqeuency-DEPENDENT vertices. Otherwise, the loop yields simply the fermionic density.
//ATTENTION: The result is a double but it represents the IMAGINARY PART of the loop function (which has no real part due to the symmetries of the system)
double loop_std(int red,double Lambda, T0& vertex, char p, self& se, self& dse,int a, int b, int c,T1 w1){

    double result=0;



#if temp==0
double  error_im;
 double abs_error = 1e-2,rel_error=1e-7;
    struct loop_params<T0> params = {Lambda,vertex, p,se, dse, a,  b,  c,  w1,red};
    gsl_function F;
#endif


    if(p =='s'){
#if temp==0
#if reg==1


        result = -1. *1./(2*pi)* vertex.vvalsmooth(red,0,a,b,c,0,w1,-Lambda,'t',1,'M')*propag(Lambda,-Lambda,se,dse,p);
        result += -1. *1./(2*pi)* vertex.vvalsmooth(red,0,a,b,c,0,w1,Lambda,'t',1,'M')*propag(Lambda,Lambda,se,dse,p);
        return result;
    //the contribution from the constant part of the vertex cancels


#elif reg==2 || reg==3

        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        double upper_bound = 7*Lambda;
        double lower_bound = -7*Lambda;


        F.params = &params;
        F.function = &loop_im<T0>;
        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //            int status=1;
        //            while(status)
        //            {
        //                status= gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        //                if(status !=0){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "loop errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;};
        //            };

        gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        gsl_set_error_handler(old_handler);




        gsl_integration_workspace_free(w);

      #endif
#elif temp==1
        int Lambda_int = static_cast<int>(Lambda/(pi*T));
        if(Lambda_int%2 == 0){Lambda_int -=1;};
        if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};

        result = -1. * vertex.vvalsmooth(red,0,a,b,c,0,w1,-Lambda_int,'t',1,'M')*propag(Lambda,-Lambda_int,se,dse,p);
        result += -1. * vertex.vvalsmooth(red,0,a,b,c,0,w1,Lambda_int,'t',1,'M')*propag(Lambda,Lambda_int,se,dse,p);
#endif
    }
    else if(p =='e'){
#if temp==0
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        double upper_bound = ffreqs[nw-1];
        double lower_bound = ffreqs[0];

        F.params = &params;
        F.function = &loop_im<T0>;
        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //            int status=1;
        //            while(status)
        //            {
        //                status= gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        //                if(status !=0){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "loop errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;};
        //            };

        gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        gsl_set_error_handler(old_handler);

        gsl_integration_workspace_free(w);
#elif temp==1
         for(int i=-(nw-1); i<nw+1; i+=2){
 result += -1. * vertex.vvalsmooth(red,0,a,b,c,0,w1,i,'t',1,'M')*propag(Lambda,i,se,dse,p);};
#endif
    };


    return result;

}




template<class T0, typename T1>
//loop Attention: Only use for freqeuency-DEPENDENT vertices. Otherwise, the loop yields simply the fermionic density.
//ATTENTION: The result is a double but it represents the IMAGINARY PART of the loop function (which has no real part due to the symmetries of the system)
double loop_sw(int red, double Lambda, T0& vertex, char p, self& se, self& dse,int a, int b, int c, T1 w1){

    double result=0;


#if temp==0
      double abs_error = 1e-2,rel_error=1e-7;
       double  error_im;
    struct loop_params<T0> params = {Lambda,vertex, p,se, dse, a,  b,  c,  w1, red};
    gsl_function F;
#endif



    if(p =='s'){
#if temp==0

#if reg==1
    result += -1. *1./(2*pi)* vertex.vvalsmooth(red,0,a,b,c,0,w1,-Lambda,'u',1,'M')*propag(Lambda,-Lambda,se,dse,p);
    result += -1. *1./(2*pi)* vertex.vvalsmooth(red,0,a,b,c,0,w1,Lambda,'u',1,'M')*propag(Lambda,Lambda,se,dse,p);
    return result;
#elif reg==2 || reg==3


        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);
        double upper_bound = 7*Lambda;
        double lower_bound = -7*Lambda;

        F.function = &loop_im<T0>;
        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //            int status=1;
        //            while(status)
        //            {
        //                status=  gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        //                if(status !=0){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "loop errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;};
        //            };
        gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        gsl_set_error_handler(old_handler);

        gsl_integration_workspace_free(w);
#endif
#elif temp==1

        int Lambda_int = static_cast<int>(Lambda/(pi*T));
        if(Lambda_int%2 == 0){Lambda_int -=1;};
        if(Lambda/(pi*T)-Lambda_int>1.){Lambda_int +=2;};

        result = (-1.) * vertex.vvalsmooth(red,0,a,b,c,0,w1,Lambda_int,'u',1,'M')*propag(Lambda,Lambda_int,se,dse,p);
        result += (-1.) * vertex.vvalsmooth(red,0,a,b,c,0,w1,-Lambda_int,'u',1,'M')*propag(Lambda,-Lambda_int,se,dse,p);
#endif
    }

    else if(p =='e'){
#if temp==0
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);
        double upper_bound = ffreqs[nw-1];
        double lower_bound = ffreqs[0];


        F.function = &loop_im<T0>;
        F.params = &params;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        //            int status=1;
        //            while(status)
        //            {
        //                   status=  gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);

        //                if(status !=0){
        //                    abs_error *= 1.2;
        //                    rel_error *= 1.2;
        //                    cout << "loop errors corrected: rel_error = " << rel_error << " , abs_error = " << abs_error << endl;};
        //            };

        gsl_integration_qag(&F,lower_bound,upper_bound,abs_error, rel_error,1500,6,w, &result, &error_im);
        gsl_set_error_handler(old_handler);

        gsl_integration_workspace_free(w);

    #elif temp==1
        for(int i=-(nw-1); i<nw+1; i+=2){
         result += (-1.) * vertex.vvalsmooth(red,0,a,b,c,0,w1,i,'u',1,'M')*propag(Lambda,i,se,dse,p);
        };
#endif
    };
    return result;

}





template<class T0>
self loop(int red,double Lambda, parvert<T0>& vertex, char p, self& se, self& dse){//see SB1, p. 76 for details. The last argument specifies if only tbar contributions are returned from the vetex as necessary in higher order corrections

    self self1, self2, self3;
    

    //**************first contr.********//
#pragma omp parallel for
    for(int i=0; i<nw; i++){//for all different external freqs..
        double k= 0,l=0,m=0;
#if temp==0
                        double w = ffreqs[i];
#elif temp==1
                        int w = -nw1+1+i*2;
#endif

        for(int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1;a++){
            for(int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1;b++){
                for(int c= 1;c<4;c++){//..iterate through all sites..
                    if(distance(a,b,c) <= d_c){//cutoff distance


                        k += loop_std(red,Lambda,vertex.densvertex,p,se,dse,a,b,c,w);


                    };};};};

        self1.setself(i,2.*k);//factor 2 is a combinatorial factor

        //**************second contr.********//

        l = loop_sw(red,Lambda,vertex.spinvertex,p,se,dse,0,0,1,w);
        //cout << "2ok" << endl;
        self2.setself(i,-3.*l);


        //**************third contr.********//


        m = loop_sw(red,Lambda,vertex.densvertex,p,se,dse,0,0,1,w);

        self3.setself(i,-1.*m);
    };//closes loop over freqs i.



    return (self1 + self2 + self3);
}





class Susc{
    vector<double>Sus =vector<double>(irred_dim1);

//    vector<vector<vector<double> >  >   Sus =
//            vector<vector<vector<double>  > >
//            (nuc_eff,vector<vector<double>  >
//             ((nuc_eff+1)/2, vector<double>(3)
//              ));//three atoms per unit cell

public:
    double vval(int, int, int);//arguments 1-3: site, argument 4, bosonic frequency
       double acc(int);//arguments 1-3: site, argument 4, bosonic frequency
    void write(int,int,int, double );//first argument: Lambda
    void add_write(int,int,int, double );//first argument: Lambda
};


Susc suscept(double , parvert<fullvert>& ,  self );



//define a struct object which includes the self energy and the vertex which are needed to evaluate the RHS of the flow equations.

struct state{
    // double Lambda;
    self selfenergy;
    //Susc sus;
    parvert<fullvert>  vertex;};


state operator+(state, state);
state abs_sum_tiny(state, state, double);
state operator*(double,state);
state operator*(state, double);
double max_err(state&,state&);






void write_hdf(const H5std_string FILE_NAME,double Lambdas_i, long Lambda_size,int total_iterations,state& vertex);
void add_hdf(const H5std_string FILE_NAME,int Lambdas_it, long Lambda_size,state& vertex);
void add_hdf_sus(const H5std_string FILE_NAME,int iteration, vector<double>& Lambdas, long total_iterations,Susc& susceptibility);
state read_hdf(const H5std_string FILE_NAME,int Lambda,long Lambda_size,long total_iterations, vector<double> &Lambdas);

void state_bcast(state & data);

#endif
