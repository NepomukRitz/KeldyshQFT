#ifndef _kagome_hpp_
#define _kagome_hpp_

#include <complex>
#include <vector>
#include "H5Cpp.h"
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <omp.h>

using namespace std;



const H5std_string	DATASET_R("R");
const H5std_string	DATASET_K1("K1");
const H5std_string	DATASET_K2("K2");
const H5std_string	DATASET_K2b("K2b");
const H5std_string	DATASET_irred("irred");
const H5std_string	DATASET_sus("sus");
const H5std_string	FERM_FREQS_LIST("ferm_freqslist");
const H5std_string	BOS_FREQS_LIST("bos_freqslist");
const H5std_string	SELF_LIST("selflist");
const H5std_string	LAMBDA_LIST("lambdas");
const H5std_string      MEMBER1( "spin_re" );
//const H5std_string      MEMBER2( "spin_im" );
const H5std_string      MEMBER2( "dens_re" );
//const H5std_string      MEMBER4( "dens_im" );
const H5std_string      MEMBER5( "re" );
const H5std_string      MEMBER6( "im" );



//IMPORTANT: COMPILE WITH C++11 !!
//the first two of the three spacial lattice indices on the kagome lattice specify the unit cell and the third one denotes the site within this UC (nuc,nuc,3) components. the first two can take the values [-(nuc-1)/2,(nuc-1)/2] and the third can take the values 1,2,3
// For different lattice structures that can be described by three independent indices, one only needs to change to conversion from input to lattice values in definition of vertex classes and the site-conversion in the t-bubble.
extern const double pi;
extern const complex<double> ci; //complex i
extern const int reg;
extern const int sym;//with all symmetries if sym = 2, with only diagrammatic symmetries if sym=1 and without if sym =0.
extern const double sharp; // governs the sharpness of the second regularizer
extern const double wt ;//transition frequency between log and lin grid.
extern const double w0;
extern const double wmax;
extern const int nlog;
extern const double k;
extern const double delw;
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
extern const complex<double> zero;
const double wlimit = 1e10;//frequency used to perform the numerical limits for the asymptotic functions
int fconv(double w); //function that converts physical freqs to the lattice site on the freq mesh that corresponds to the closest grid frequency below
int fconv_n(double w, int n); //function that converts physical freqs to the lattice site on a mesh of size n that corresponds to the closest grid frequency below

//int bconv(double q); //function that converts physical freqs to the lattice site on the freq mesh that corresponds to the closest grid frequency below
//int bconv_n(double q, int n); //function that converts physical freqs to the lattice site on a mesh of size n that corresponds to the closest grid frequency below

bool comp(int , int);

struct site{
    int a;
    int b;
    int c;
public:
    site(int x,int y,int z){a = x; b = y; c = z;}
    site(){}
};
site site_switch(int a, int b, int c);
site site_project(int a, int b, int c);//projects any site to its equivalent site in the upper half of the lattice

double distance(int a, int b, int c);//yields the distance of an arbitrary site to the reference site at the origin.



/*******************************CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX*********************************/

class svert{


    //rest function:

    vector<vector<vector<vector<vector<vector<double > > > > > > R =
            vector<vector<vector<vector<vector<vector<double > > > > > >
            (nuc_eff,vector<vector<vector<vector<vector<double > > > > >
             ((nuc_eff+1)/2, vector<vector<vector<vector<double > > > >
              (3,vector<vector<vector<double > > >//three atoms per unit cell
               (nw3_q, vector<vector<double > >//nw frequency entries for each of the three frequencies
                (nw3_w1, vector<double >(nw3_w2))))));

    //K1:
    vector<vector<vector<vector<double > > > >  K1 =
            vector<vector<vector<vector<double > > > >
            (nuc_eff,vector<vector<vector<double > > >
             ((nuc_eff+1)/2, vector<vector<double > >
              (3,vector<double >(nw1_q))));//three atoms per unit cell

    //K2:
    vector<vector<vector<vector<vector<double > > > > > K2 =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));


    //K2b:
    vector<vector<vector<vector<vector<double > > > > > K2b =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));





public:

    double vvalsmooth(int, int, int, double, double, double, char);
    double vvalsmooth(int, int, int, double, double, double, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int, int, int, double, double, double, char, int,char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int, int, double);
    void K2b_setvert( int, int, int, int,int, double);
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

};
class tvert{



    vector<vector<vector<vector<vector<vector<double > > > > > > R =
            vector<vector<vector<vector<vector<vector<double > > > > > >
            (nuc_eff,vector<vector<vector<vector<vector<double > > > > >
             ((nuc_eff+1)/2, vector<vector<vector<vector<double > > > >
              (3,vector<vector<vector<double > > >//three atoms per unit cell
               (nw3_q, vector<vector<double > >//nw frequency entries for each of the three frequencies
                (nw3_w1, vector<double >(nw3_w2))))));

    //K1:
    vector<vector<vector<vector<double > > > >  K1 =
            vector<vector<vector<vector<double > > > >
            (nuc_eff,vector<vector<vector<double > > >
             ((nuc_eff+1)/2, vector<vector<double > >
              (3,vector<double >(nw1_q))));//three atoms per unit cell

    //K2:
    vector<vector<vector<vector<vector<double > > > > > K2 =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));


    //K2b:
    vector<vector<vector<vector<vector<double > > > > > K2b =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));





public:
    double vvalsmooth(int, int, int, double, double, double,char);
    double vvalsmooth(int, int, int, double, double, double, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int, int, int, double, double, double, char, int,char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int,int, double);
    void K2b_setvert( int, int, int, int,int, double);
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

};
class uvert{


    //rest function:
    vector<vector<vector<vector<vector<vector<double > > > > > > R =
            vector<vector<vector<vector<vector<vector<double > > > > > >
            (nuc_eff,vector<vector<vector<vector<vector<double > > > > >
             ((nuc_eff+1)/2, vector<vector<vector<vector<double > > > >
              (3,vector<vector<vector<double > > >//three atoms per unit cell
               (nw3_q, vector<vector<double > >//nw frequency entries for each of the three frequencies
                (nw3_w1, vector<double >(nw3_w2))))));

    //K1:
    vector<vector<vector<vector<double > > > >  K1 =
            vector<vector<vector<vector<double > > > >
            (nuc_eff,vector<vector<vector<double > > >
             ((nuc_eff+1)/2, vector<vector<double > >
              (3,vector<double >(nw1_q))));//three atoms per unit cell

    //K2:
    vector<vector<vector<vector<vector<double > > > > > K2 =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));


    //K2b:
    vector<vector<vector<vector<vector<double > > > > > K2b =
            vector<vector<vector<vector<vector<double > > > > >
            (nuc_eff,vector<vector<vector<vector<double > > > >
             ((nuc_eff+1)/2, vector<vector<vector<double > > >
              (3,vector<vector<double > >//three atoms per unit cell
               (nw2_q, vector<double >(nw2_w1)))));





public:
    double vvalsmooth(int, int, int, double, double, double, char);
    double vvalsmooth(int, int, int, double, double, double, char, int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int, int, int, double, double, double, char, int,char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void R_setvert( int, int, int, int,int,int, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int,int, double);
    void K2b_setvert( int, int, int, int,int, double);
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

};
class irreducible{
    vector<vector<vector<double > > >  U_bare =
            vector<vector<vector<double > > >
            (nuc_eff,vector<vector<double > >
             ((nuc_eff+1)/2, vector<double >(3)));//three atoms per unit cell
    //the irreducible vertex is approximated by the bare interaction in the parquet approx
public:
    double vval(int, int, int);
    double vvalsmooth(int, int, int);
    double vvalsmooth(int, int, int,double,double,double,char,int,char);

    void setvert(int,int,int,double);

    friend irreducible operator*(double alpha, const irreducible & vertex);
    friend irreducible  operator*(const irreducible & vertex, double alpha);
    friend irreducible  operator+(const irreducible & vertex1, const irreducible & vertex2);


};
/***************************************************************************************************************************/


/*******************************define "fullvert"FULLVERT" as collection of all diagrams in all channels*********************************/
struct fullvert{//collection of all channels
    irreducible irred;
    svert svertex;
    tvert tvertex;
    uvert uvertex;
public:

    double vvalsmooth(int,int,int,double,double,double,char);
    double vvalsmooth(int,int,int,double,double,double,char,int,char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map
    friend fullvert operator*(double alpha, const fullvert& vertex);
    friend fullvert operator+(const fullvert& vertex1,const fullvert& vertex2);

};


/******************CLASS FOR SELF ENERGY *************/
class self{
    vector<complex<double> > selfenergy =  vector<complex<double > >(nw);
public:
    void setself(int, complex<double>);
    complex<double> sval(int);
    complex<double> svalsmooth(double);
    friend self operator+(const self& self1,const  self& self2);
    friend self operator+=(const self& self1,const  self& self2);
    friend self operator*(double alpha,const  self& self1);
    friend self operator*(const self& self1, double alpha);
};


/*******PROPAGATOR FUNCTION***********/
complex<double> propag(double Lambda, double w, self selfenergy,self diffselfenergy, char type);





/*******BUBBLE INTEGRATION FUNCTIONS***********/
//s-bubble:
template<typename  T1,typename  T2>

struct sbubble_params{


    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;
    char ptype1;
    char ptype2;
    self& selfen;
    self& diffselfen;


    int a; int b; int c;

    int d; int e; int f;

    double s;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double sbubble_re(double w, void * p){
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
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
    double val = real(1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  propag(Lambda,w+s/2,selfen, diffselfen,ptype1) * propag(Lambda,s/2-w,selfen,diffselfen,ptype2) );
    return (1./(2*pi)*val);
}


template<typename T1,typename T2>
double sgreensfunc_re(double w, void * p){//intergrates only propagators - used if vertices are constant in respective integration interval
    struct sbubble_params<T1,T2> * params
            = static_cast< struct sbubble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);

    double s = (params->s);

    double val = real(1./2 * propag(Lambda,w+s/2,selfen, diffselfen,ptype1) * propag(Lambda,s/2-w,selfen,diffselfen,ptype2) );

    return (1./(2*pi)*val);
}


template<typename T1,typename T2>
double sbubble(int red_side,int map1, int map2, gsl_integration_workspace* w, double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, double s,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.
    double abs_error=1e-2, rel_error=1e-4;
   double abs_error_bare=1e-4,rel_error_bare=1e-4;
    //Note: factor 1/2 is included due to indistiguishibility of propagators in s-bubble

    double B =0;

    if((reg ==1 && p1 =='s') ){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,s+Lambda,se,dse,p2)) ;
        B += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,s-Lambda, se,dse, p2));
    }
    else if((reg ==1 && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration
        B =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2)) ;
        B += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2)) ;

    }


    else{

        double limit_low,limit_low1,limit_low2;
        double limit_up,limit_up1,limit_up2;



        double vert1_const,vert2_const;



        int mode =1;

        if(h=='R'){
            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-s},comp);//choose the minimum since beyond this value, this vertex cannot contribute any more
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+s},comp);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2-s,2*ffreqs[(nw+nw3)/2-1]-w2+s},comp);
            double v1_K1_u_up= bfreqs[nw1-1]-w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-s},comp);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+s},comp);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-s,2*ffreqs[(nw+nw3)/2-1]+w2+s},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1-s},comp);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1+s},comp);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw+nw3)/2-1]-w1+s,2*ffreqs[(nw+nw3)/2-1]-w1-s},comp);
            double v2_K1_u_up = bfreqs[nw1-1]-w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+s},comp);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-s},comp);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+s,2*ffreqs[(nw+nw3)/2-1]+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-s},comp);//choose the minimum since bezond this value, this vertex cannot contribute any more
            double v1_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+s},comp);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2-s,2*ffreqs[(nw-nw3)/2]-w2+s},comp);
            double v1_K1_u_low = bfreqs[0]-w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-s},comp);
            double v1_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+s},comp);
            double v1_R_u_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-s,2*ffreqs[(nw-nw3)/2]+w2+s},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low =  bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw-+nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1-s},comp);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1+s},comp);
            double v2_R_t_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw-nw3)/2]-w1+s,2*ffreqs[(nw-nw3)/2]-w1-s},comp);
            double v2_K1_u_low = bfreqs[0]-w1;
            double v2_K2_u_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+s},comp);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-s},comp);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+s,2*ffreqs[(nw-nw3)/2]+w1-s},comp);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];

            vector<double> upperR_v1{v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperR_v1.begin(),upperR_v1.end());
            vector<double> lowerR_v1{v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerR_v1.begin(),lowerR_v1.end());

            vector<double> upperR_v2{v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperR_v2.begin(),upperR_v2.end());
            vector<double> lowerR_v2{v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerR_v2.begin(),lowerR_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,w2,'s',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,wlimit,'s',2,h);

            for(int i=upperR_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,upperR_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperR_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperR_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperR_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperR_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperR_v1[i],'s',2,h))>0 ){limit_up1 = upperR_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerR_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerR_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerR_v1[i],'s',2,h))>0){limit_low1 = lowerR_v1[i];break;};
                if(i==lowerR_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerR_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerR_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerR_v2[i],w2,'s',1,h))>0){limit_low2 = lowerR_v2[i];break;};
                if(i==lowerR_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h=='L'){//in this case, only Gamma (vert2) sets the limits

            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1-s},comp);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1+s},comp);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw+nw3)/2-1]-w1+s,2*ffreqs[(nw+nw3)/2-1]-w1-s},comp);
            double v2_K1_u_up = bfreqs[nw1-1]-w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+s},comp);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-s},comp);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+s,2*ffreqs[(nw+nw3)/2-1]+w1-s},comp);


            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low =  bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw-+nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1-s},comp);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1+s},comp);
            double v2_R_t_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw-nw3)/2]-w1+s,2*ffreqs[(nw-nw3)/2]-w1-s},comp);
            double v2_K1_u_low = bfreqs[0]-w1;
            double v2_K2_u_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+s},comp);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-s},comp);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+s,2*ffreqs[(nw-nw3)/2]+w1-s},comp);

            vector<double> upperK2_v1{v_K2_s_up};
            sort(upperK2_v1.begin(),upperK2_v1.end());
            vector<double> lowerK2_v1{v_K2_s_low};
            sort(lowerK2_v1.begin(),lowerK2_v1.end());

            vector<double> upperK2_v2{v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,wlimit,'s',1,h);
            vert2_const = vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,wlimit,'s',2,h);


            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                        if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK2_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2_v1.size()-1; i>-1 ; i--){

                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2_v1[i],'s',2,h))>0 ){limit_up1 = upperK2_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};

            for(int i=0; i< lowerK2_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2_v1[i],'s',2,h))>0){limit_low1 = lowerK2_v1[i];break;};
                if(i==lowerK2_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK2_v2[i];break;};
                if(i==lowerK2_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};


            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};

        }

        else if(h=='M'){//in this case, only Gamma' (vert1) sets the limits

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-s},comp);//choose the minimum since beyond this value, this vertex cannot contribute any more
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+s},comp);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2-s,2*ffreqs[(nw+nw3)/2-1]-w2+s},comp);
            double v1_K1_u_up= bfreqs[nw1-1]-w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-s},comp);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+s},comp);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-s,2*ffreqs[(nw+nw3)/2-1]+w2+s},comp);

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-s},comp);//choose the minimum since bezond this value, this vertex cannot contribute any more
            double v1_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+s},comp);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2-s,2*ffreqs[(nw-nw3)/2]-w2+s},comp);
            double v1_K1_u_low = bfreqs[0]-w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-s},comp);
            double v1_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+s},comp);
            double v1_R_u_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-s,2*ffreqs[(nw-nw3)/2]+w2+s},comp);

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];

            vector<double> upperK2b_v1{v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_s_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_s_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,w2,'s',1,h);
            vert2_const =vert1.vvalsmooth(red_side,map2,d,e,f,s,wlimit,wlimit,'s',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2b_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2b_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2b_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2b_v1[i],'s',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;};
                if(i==0){mode1_up=0;};
            };


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2b_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2b_v1[i],'s',2,h))>0){limit_low1 = lowerK2b_v1[i];break;};
                if(i==lowerK2b_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2b_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2b_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK2b_v2[i];break;};
                if(i==lowerK2b_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};

        }

        else if(h=='K'){ //in this case the only contribution to the non-constat part comes from K2,t

            double v_K2_s_low = ffreqs[(nw1-nw2)/2];

            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];


            vector<double> upperK1_v1{v_K2_s_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_s_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_s_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_s_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,wlimit,'s',1,h);
            vert2_const = vert1.vvalsmooth(red_side,map2,d,e,f,s,wlimit,wlimit,'s',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK1_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK1_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK1_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK1_v1[i],'s',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK1_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK1_v1[i],'s',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK1_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK1_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };
            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};


        };

        double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct sbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};
        gsl_function F;

        if(mode==0 && abs(vert1_const * vert2_const)>0){

            if(reg==1 && p1=='g' && p2 =='g'){
                F.function = &sgreensfunc_re<T1,T2>;
                F.params = &params;

              //even integrand --> sufficient to intergrate one side and multiply result by two
                double upper = abs(s/2) + Lambda;


                F.function = &sgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);             
                double compl_real = 2. * resultgfu_re *vert1_const * vert2_const;

                B = compl_real;

            }


            else if(reg==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(s/2)-Lambda;
                double upper = abs(s/2) + Lambda;
                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;




                if(p1=='k'){
                    singsc =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,s+Lambda,se,dse,p2)) ;
                    singsc += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,s-Lambda, se,dse, p2));
                         p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]-s/2;dse_cutoff_low= ffreqs[0]-s/2;}

                else if(p2=='k'){
                    singsc =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                                 p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]+s/2;dse_cutoff_low= ffreqs[0]+s/2;};

                struct sbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;};

                    F.function = &sgreensfunc_re<T1,T2>;
                    F.params = &params_mod;

                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = resultgfl_re *vert1_const * vert2_const;                   
                    B= singsc+ compl_real ;

                }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};

                    F.function = &sgreensfunc_re<T1,T2>;
                    F.params = &params_mod;

                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                            double compl_real = resultgfl_re *vert1_const * vert2_const;
                        B= singsc + compl_real;
                   }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){
                    F.function = &sgreensfunc_re<T1,T2>;
                    F.params = &params_mod;


                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                           F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);
                           double compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    B= singsc+ compl_real ;
                   };}


            else if(reg==2 && (p1=='g' &&p2=='g')){


                F.params = &params;
                F.function = &sgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                      w, &result_re, &errorgfu_re);

                double compl_real = 2*result_re *vert1_const * vert2_const;           
                B=compl_real;

            }
            else if(reg==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0;

                if(p1=='s'){bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;}
                else if(p2=='s'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;}
                F.params = &params;
                F.function = &sgreensfunc_re<T1,T2>;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                    w, &result_re, &error_re);
                     double compl_real = result_re *vert1_const * vert2_const;
                  B=compl_real;
             }

            else if(reg==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0, dse_cutoff_low=0, dse_cutoff_up=0;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;
                    dse_cutoff_low = ffreqs[0]-s/2; dse_cutoff_up = ffreqs[nw-1]-s/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;
                    dse_cutoff_low = ffreqs[0]+s/2; dse_cutoff_up = ffreqs[nw-1] +s/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct sbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};
                F.params = &params_singsc;
                F.function = &sgreensfunc_re<T1,T2>;//integration of single scale
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct sbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of katanin extension
                F.function = &sgreensfunc_re<T1,T2>;
                double result_re2=0,error_re2=0;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error, rel_error,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real = (result_re + result_re2 ) *vert1_const * vert2_const;

                B=compl_real;

            };
        }

        else if(mode==1){


                if(reg ==1 && p1=='g' && p2=='g'){
                    double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                    double lower = -abs(s/2)-Lambda;
                    double upper = abs(s/2) + Lambda;

                    if(lower <= limit_low){//left side of sharp cutoff

                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfl_re, &errorgfl_re);

                            double compl_real = resultgfl_re  *vert1_const * vert2_const;

                            lhs=compl_real;

                        };}
                    else if(lower <= limit_up){
                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfl_re, &errorgfl_re);
                                   compl_real = resultgfl_re *vert1_const * vert2_const;

                        };

                        F.params = &params;
                        F.function = &sbubble_re<T1,T2>;

                        gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                               lhs=result_re + compl_real;

                    }

                    else if(lower > limit_up){

                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfl_re, &errorgfl_re);

                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                                          compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                         };
                        F.params = &params;
                        F.function = &sbubble_re<T1,T2>;

                        gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                                lhs=result_re + compl_real;

                    };

                    if(upper >= limit_up){//right side of sharp cutoff

                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfu_re, &errorgfu_re);
                                  double compl_real = resultgfu_re *vert1_const * vert2_const;
                                 rhs=compl_real;
                            };}
                    else if(upper >= limit_low){
                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfu_re, &errorgfu_re);
                                    compl_real =resultgfu_re  *vert1_const * vert2_const;
                         };
                        F.params = &params;
                        F.function = &sbubble_re<T1,T2>;

                        gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                                 rhs=result_re + compl_real;

                    }

                    else if (upper < limit_low){
                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0){
                            F.params = &params;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfl_re, &errorgfl_re);
                                    F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                                  w, &resultgfu_re, &errorgfu_re);
                                       compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                          };
                        F.params = &params;
                        F.function = &sbubble_re<T1,T2>;

                        gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                                rhs = result_re + compl_real;
                        };

                    B =rhs + lhs;
                   }

            else if(reg==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){

                double lower = -abs(s/2)-Lambda;
                double upper = abs(s/2) + Lambda;
                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;


                if(p1=='k'){
                    singsc =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,s+Lambda,se,dse,p2)) ;
                    singsc += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,s-Lambda, se,dse, p2));
                         p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]-s/2;dse_cutoff_low= ffreqs[0]-s/2;}

                else if(p2=='k'){
                    singsc =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                                 p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]+s/2;dse_cutoff_low= ffreqs[0]+s/2;};


                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct sbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                              double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    lhs=compl_real;

                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                             compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &sbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                              lhs=result_re + compl_real;
                      }

                    else if(lower > limit_up){

                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                            double up_eff=0;
                            if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                            F.params = &params_mod;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfl_re, &errorgfl_re);
                                    compl_real = resultgfl_re  *vert1_const * vert2_const;
                         };


                        if(dse_cutoff_low <limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &sbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);
                               };


                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                            double up_eff=0;
                            if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                            F.params = &params_mod;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfl_re, &errorgfl_re);
                                 compl_real += resultgfl_re  *vert1_const * vert2_const;

                        };

                        lhs=result_re + compl_real;
                      };

                    if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                               double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                        rhs=compl_real;

                    }
                    else if(upper >= limit_low ){
                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                               compl_real =(resultgfu_re  )*vert1_const * vert2_const;
                         };
                        if(dse_cutoff_up > upper){
                            double low_eff=0;
                            if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &sbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);

                        };

                        rhs=result_re + compl_real;

                    }

                    else if(upper < limit_low){

                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                            double low_eff=0;
                            if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                            F.params = &params_mod;

                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                                   compl_real = (resultgfu_re  )*vert1_const * vert2_const;
                         };


                        if(dse_cutoff_up >limit_low){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &sbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);
                             };


                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                            F.params = &params_mod;
                            F.function = &sgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                                 compl_real += (resultgfu_re  )*vert1_const * vert2_const;

                        };

                        rhs=result_re + compl_real;

                    };


                    B = rhs + lhs + singsc;
                        };}


            else if(reg==2 && p1=='g' && p2=='g'){

                F.function = &sbubble_re<T1,T2>;
                F.params = &params;

                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                double compl_real=0;

                if(abs(vert1_const * vert2_const)>0){
                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &sgreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=result_re + compl_real ;
              }

            else if(reg==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){

                   double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                  if(p1=='s'){bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;}
                else if(p2=='s'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                   if(m1==1 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                             //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                       if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };




                    B=result_re + compl_real ;
                  }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                           //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                              compl_real = (resultgfl_re )*vert1_const * vert2_const;
                      };


                    B = result_re + compl_real ;
                 }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;
                    };



                    B=result_re + compl_real ;
                  }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                               compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=compl_real ;
                 }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    B=result_re  ;
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real =(resultgfu_re )*vert1_const * vert2_const;
                      };

                    B=compl_real ;
                 }

            }

            else if(reg==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct sbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                             compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                      };




                    singsc=result_re + compl_real;
                 }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re)*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                         compl_real = (resultgfu_re )*vert1_const * vert2_const;
                     };



                    singsc=result_re + compl_real ;
                  }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                         compl_real = (resultgfl_re )*vert1_const * vert2_const;
                     };


                    singsc= compl_real ;
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    singsc=result_re  ;
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real = (resultgfu_re)*vert1_const * vert2_const;

                    };

                    singsc= compl_real ;
                 }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]-s/2; bound_up = ffreqs[nw-1]-s/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]+s/2; bound_up = ffreqs[nw-1]+s/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct sbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                         //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                         compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    kat=result_re + compl_real ;
                 }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                           //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                         compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=(result_re + compl_real );
                 }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){


                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=(result_re + compl_real );
                 }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    //    double compl_imag;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                         compl_real =(resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=( compl_real );
                 }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &sbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &sgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                       compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat= compl_real ;
                 }
                B=kat + singsc;

            };


        };

    };
    if(abs(B)<1e-20){B=0;};
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
    char ptype1;
    char ptype2;
    self& selfen;
    self& diffselfen;

    int a; int b; int c;

    int d; int e; int f;

    double t;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double tbubble_re(double w, void * p){
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
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

    double val = real(-1.*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  propag(Lambda,w-t/2,selfen,diffselfen,ptype1) * propag(Lambda,w+t/2,selfen,diffselfen,ptype2) );

    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double tgreensfunc_re(double w, void * p){//returns the value of all constituents of a bubble function at a specific value of the integration parameter WITHOUT contribution of bare vertices!
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);

    double t = (params->t);


    double val = real( propag(Lambda,w-t/2,selfen,diffselfen,ptype1) * propag(Lambda,w+t/2,selfen,diffselfen,ptype2) );


    return (-1./(2.*pi)*val);//minus sign is from definition of t-bubble
}





template<typename T1,typename T2>
double tbubble(int red_side, int map1, int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, double t,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.

    double abs_error=1e-2, rel_error=1e-4;

     double abs_error_bare=1e-2,rel_error_bare=1e-4;

    double B=0;


    if((reg ==1 && p1 =='s' && p2 =='g' ) ){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B = real( -1. *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) * propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
        B += real(-1. *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) * propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
    }
    else if((reg ==1 && p1=='g' && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration

        B = real(-1. * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) * propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
        B += real(-1. *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) * propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
    }

    else{

        double limit_low=0,limit_low1=0,limit_low2=0;
        double limit_up=0,limit_up1=0,limit_up2=0;

        double vert1_const,vert2_const;





        int mode =1;

        if(h == 'R'){

            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_u_up = bfreqs[nw1-1]+w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+t},comp);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-t},comp);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+t,2*ffreqs[(nw+nw3)/2-1]-w2-t},comp);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+t},comp);
            double v1_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-t},comp);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2+t,2*ffreqs[(nw+nw3)/2-1]+w2-t},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_up = bfreqs[nw1-1]+w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+t},comp);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-t},comp);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1+t,2*ffreqs[(nw1+nw3)/2-1]-w1-t},comp);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1-t},comp);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1+t},comp);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw1+nw3)/2-1]+w1+t,2*ffreqs[(nw1+nw3)/2-1]+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_u_low = bfreqs[0]+w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+t},comp);
            double v1_K2b_u_low =  max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-t},comp);
            double v1_R_u_low =  max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+t,2*ffreqs[(nw-nw3)/2]-w2-t},comp);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+t},comp);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-t},comp);
            double v1_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2+t,2*ffreqs[(nw-nw3)/2]+w2-t},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_low = bfreqs[0]+w1;
            double v2_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+t},comp);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-t},comp);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1+t,2*ffreqs[(nw1-nw3)/2]-w1-t},comp);
            double v2_K1_s_low =  bfreqs[0]-w1;
            double v2_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1-t},comp);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1+t},comp);
            double v2_R_s_low =max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw1-nw3)/2]+w1+t,2*ffreqs[(nw1-nw3)/2]+w1-t},comp);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];




            vector<double> upperR_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperR_v1.begin(),upperR_v1.end());
            vector<double> lowerR_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerR_v1.begin(),lowerR_v1.end());

            vector<double> upperR_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperR_v2.begin(),upperR_v2.end());
            vector<double> lowerR_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerR_v2.begin(),lowerR_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,w2,'t',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,wlimit,'t',2,h);

            for(int i=upperR_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperR_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperR_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperR_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperR_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperR_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperR_v1[i],'t',2,h))>0 ){limit_up1 = upperR_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};



            for(int i=0; i< lowerR_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerR_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerR_v1[i],'t',2,h))>0){limit_low1 = lowerR_v1[i];break;};
                if(i==lowerR_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerR_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerR_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerR_v2[i],w2,'t',1,h))>0){limit_low2 = lowerR_v2[i];break;};
                if(i==lowerR_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
              if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h == 'L'){ //constant vertices


            //conditions from vertex 2 (unprimed):
            double v2_K1_u_up = bfreqs[nw1-1]+w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+t},comp);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-t},comp);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1+t,2*ffreqs[(nw1+nw3)/2-1]-w1-t},comp);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1-t},comp);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1+t},comp);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw1+nw3)/2-1]+w1+t,2*ffreqs[(nw1+nw3)/2-1]+w1-t},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_low = bfreqs[0]+w1;
            double v2_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+t},comp);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-t},comp);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1+t,2*ffreqs[(nw1-nw3)/2]-w1-t},comp);
            double v2_K1_s_low =  bfreqs[0]-w1;
            double v2_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1-t},comp);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1+t},comp);
            double v2_R_s_low =max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw1-nw3)/2]+w1+t,2*ffreqs[(nw1-nw3)/2]+w1-t},comp);



            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];




            vector<double> upperK2_v1{v_K2_t_up};

            vector<double> lowerK2_v1{v_K2_t_low};


            vector<double> upperK2_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,wlimit,'t',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,wlimit,'t',2,h);

            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                 if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK2_v2[i];break;};
                if(i==0){mode2_up=0;};


            };



            if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2_v1[0],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2_v1[0],'t',2,h))>0 ){limit_up1 = upperK2_v1[0];}
            else{mode1_up=0;};

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2_v1[0],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2_v1[0],'t',2,h))>0){limit_low1 = lowerK2_v1[0];}
                else{mode1_low=0;};


            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK2_v2[i];break;};
                if(i==lowerK2_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
             if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};


        }

        else if(h == 'M'){//in this case, only Gamma' (vert1) sets the limits

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_u_low = bfreqs[0]+w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+t},comp);
            double v1_K2b_u_low =  max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-t},comp);
            double v1_R_u_low =  max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+t,2*ffreqs[(nw-nw3)/2]-w2-t},comp);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+t},comp);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-t},comp);
            double v1_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2+t,2*ffreqs[(nw-nw3)/2]+w2-t},comp);

            //conditions from vertex 1 (primed):
            double v1_K1_u_up = bfreqs[nw1-1]+w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+t},comp);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-t},comp);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+t,2*ffreqs[(nw+nw3)/2-1]-w2-t},comp);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+t},comp);
            double v1_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-t},comp);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2+t,2*ffreqs[(nw+nw3)/2-1]+w2-t},comp);


            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];


            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];

            vector<double> upperK2b_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_t_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_t_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;

            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,w2,'t',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,t,wlimit,wlimit,'t',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2b_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2b_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2b_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2b_v1[i],'t',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;};
                if(i==0){mode1_up=0;};
            };


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};

            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2b_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2b_v1[i],'t',2,h))>0){limit_low1 = lowerK2b_v1[i];break;};
                if(i==lowerK2b_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2b_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2b_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK2b_v2[i];break;};
                if(i==lowerK2b_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
             if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h == 'K'){ //in this case the only contribution to the non-constat part comes from K2,t

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];



            vector<double> upperK1_v1{v_K2_t_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_t_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_t_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_t_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,wlimit,'t',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,t,wlimit,wlimit,'t',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK1_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK1_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK1_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK1_v1[i],'t',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK1_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK1_v1[i],'t',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK1_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK1_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};

            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};



        };


      double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct tbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};
        gsl_function F;





        if(mode==0 && abs(vert1_const * vert2_const)>0){

            if(reg==1 && p1=='g' && p2 =='g'){

                F.function = &tgreensfunc_re<T1,T2>;
                F.params = &params;


                double upper = abs(t/2) + Lambda;

//even integrand -> integrate only upper half and multiply by two


                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);

                double compl_real = 2*resultgfu_re *vert1_const * vert2_const;

                B= compl_real ;

            }
            else if(reg==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;



                char p1_new,p2_new;
                double dse_cutoff_up, dse_cutoff_low;//cutoff due to finite number of saved value in differentiated self energy
                double singsc;

                if(p1=='k'){
                    singsc = real( -1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+t/2;dse_cutoff_low= ffreqs[0]+t/2;}

                else if(p2=='k'){
                    singsc =  real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) *  propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) *  propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-t/2;dse_cutoff_low= ffreqs[0]-t/2;};

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;};


                    F.params = &params_mod;

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    B=singsc+ compl_real ;
                }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};


                    F.params = &params_mod;

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    B= singsc + compl_real ;
                 }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){

                    F.params = &params_mod;


                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);


                    double compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                      B= singsc ;//+ compl_real ;
                  };}


            else if(reg==2 && (p1=='g' &&p2=='g')){
                   F.params = &params;


                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                      w, &result_re, &errorgfu_re);


                double compl_real = 2*result_re *vert1_const * vert2_const;

                B=compl_real;

            }
            else if(reg==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low, bound_up;
                if(p1=='s'){bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;}
                else if(p2=='s'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;}
                F.params = &params;
                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                    w, &result_re, &error_re);
                   double compl_real = (result_re) *vert1_const * vert2_const;

                B=compl_real;

            }

            else if(reg==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low, bound_up, dse_cutoff_low, dse_cutoff_up;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;
                    dse_cutoff_low = ffreqs[0]+t/2; dse_cutoff_up = ffreqs[nw-1]+t/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;
                    dse_cutoff_low = ffreqs[0]-t/2; dse_cutoff_up = ffreqs[nw-1] -t/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct tbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};
                F.params = &params_singsc;
                F.function = &tgreensfunc_re<T1,T2>;//integration of single scale
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of katanin extension
                F.function = &tgreensfunc_re<T1,T2>;
                double result_re2,error_re2,result_im2,error_im2;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare, rel_error_bare,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real =(result_re + result_re2 ) *vert1_const * vert2_const;

                B=compl_real;

            };
        }
        else if(mode==1){


            if(reg==1 && p1=='g' && p2=='g'){
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;

                  if(lower <= limit_low){//left side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                         double compl_real = (resultgfl_re )*vert1_const * vert2_const;
                         lhs=compl_real;
                     };}
                else if(lower <= limit_up){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                           compl_real = (resultgfl_re  )*vert1_const * vert2_const;
                     };

                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      lhs= result_re + compl_real;
                 }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     lhs= result_re + compl_real;


                };

                if(upper >= limit_up){//right side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                          double compl_real =(resultgfu_re  )*vert1_const * vert2_const;

                        rhs=(compl_real);

                    };}
                else if(upper >= limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                          compl_real = (resultgfu_re  )*vert1_const * vert2_const;
                     };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;

                }

                else if (upper < limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                          compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;

                };
                B= rhs +lhs;
            }
            else if((reg==1 && p1=='k' && p2=='g') || (reg==1 && p1=='g' && p2=='k') ){


                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;


                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;

                if(p1=='k'){

                    singsc = real(- 1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+t/2;dse_cutoff_low= ffreqs[0]+t/2;}

                else if(p2=='k'){

                    singsc = real( -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) *  propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) *  propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-t/2;dse_cutoff_low= ffreqs[0]-t/2;};
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    lhs=compl_real;

                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                          compl_real = (resultgfl_re )*vert1_const * vert2_const;
                     };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                     };

                    lhs=result_re + compl_real;
                 }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           compl_real = resultgfl_re*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_low <limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                     };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                          compl_real += (resultgfl_re  )*vert1_const * vert2_const;

                    };

                    lhs=result_re + compl_real;
                 };

                if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                    double low_eff=0;
                    if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                    F.params = &params_mod;
                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);

                    double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    rhs=compl_real;
                }
                else if(upper >= limit_low ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                     };

                    rhs=result_re + compl_real;

                }

                else if(upper < limit_low){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                           compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_up >limit_low){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                      };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real +=(resultgfu_re  )*vert1_const * vert2_const;
                      };

                    rhs=result_re + compl_real;
                  };

                B= rhs + lhs + singsc;
              }

            else if(reg==2 && p1=='g' && p2=='g'){

                F.function = &tbubble_re<T1,T2>;
                F.params = &params;



                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);
                  //integration of greens function outside of w-dependent interval

                double compl_real=0;
                 if(abs(vert1_const * vert2_const)>0){
                    F.function = tgreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=result_re + compl_real ;
              }
            else if(reg==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){

                double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                if(p1=='s'){bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;}
                else if(p2=='s'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                         //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                             compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    B=result_re + compl_real ;
                 }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    B=(result_re + compl_real );
                  }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                          compl_real =(resultgfu_re )*vert1_const * vert2_const;
                    };



                    B=result_re + compl_real ;
                 }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real =(resultgfl_re )*vert1_const * vert2_const;
                       };


                    B=( compl_real );
                  }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    B=result_re  ;
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };

                    B=( compl_real );
                  }

            }

            else if(reg==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct tbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                            //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                                 compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                     };




                    singsc=result_re + compl_real ;
                 }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                          compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                  }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                           compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    singsc=result_re + compl_real ;
                  }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                          compl_real =(resultgfl_re )*vert1_const * vert2_const;
                       };


                    singsc=compl_real ;
                 }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                        //integration of greens function outside of w-dependent interval

                    singsc=result_re ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                          compl_real =(resultgfu_re )*vert1_const * vert2_const;
                       };

                    singsc = compl_real ;
                  }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]+t/2; bound_up = ffreqs[nw-1]+t/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]-t/2; bound_up = ffreqs[nw-1]-t/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                        //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re)*vert1_const * vert2_const;

                    };




                    kat=result_re + compl_real ;
                  }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                         //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    kat=result_re + compl_real ;
                 }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){


                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re)*vert1_const * vert2_const;

                    };


                    kat= compl_real ;
                  }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                            compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat= compl_real ;
                 }
                B=kat + singsc;

            };


        };


     };



    if(abs(B)<1e-20){B=0;};//prevents errors from unprecise cancellation of diagrams
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
    char ptype1;
    char ptype2;
    self& selfen;
    self& diffselfen;

    int a; int b; int c;

    int d; int e; int f;

    double u;double w1; double w2;
    char h;//class of bubble (K = K1, L=K2, M = K2b, R = R)

};

template<typename T1,typename T2>
double ububble_re(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
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

    double val = real(vert1.vvalsmooth(red_side,map1,a,b,c,u,w,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,w,'u',2,h) *  propag(Lambda,w-u/2,selfen,diffselfen,ptype1) * propag(Lambda,w+u/2,selfen,diffselfen,ptype2) );
    return (1./(2*pi)*val);
}



template<typename T1,typename T2>
double ugreensfunc_re(double w, void * p){
    struct ububble_params<T1,T2> * params
            = static_cast< struct ububble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);

    double u = (params->u);


    double val = real( propag(Lambda,w-u/2,selfen,diffselfen,ptype1) * propag(Lambda,w+u/2,selfen,diffselfen,ptype2) );
    return (1./(2*pi)*val);
}





template<typename T1,typename T2>
double ububble(int red_side,int map1,int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, self& se, self& dse, double u,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.

    double abs_error=1e-2, rel_error=1e-4;

     double abs_error_bare=1e-4,rel_error_bare=1e-4;
    double B=0;

    if((reg ==1 && p1 =='s')){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B = real(1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h)* propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,-Lambda+u,se,dse,p2)) ;
        B += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,Lambda+u,se,dse,p2)) ;
    }
    else if((reg ==1 && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration

        B =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2)) ;
        B += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2)) ;
    }





    else{
        double limit_low,limit_low1,limit_low2;
        double limit_up,limit_up1,limit_up2;



        double vert1_const,vert2_const;




        int mode =1;

        if(h=='R'){


            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+u},comp);
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-u},comp);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+u,2*ffreqs[(nw+nw3)/2-1]-w2-u},comp);;
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+u},comp);
            double v1_K2b_s_up =  min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-u},comp);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-u,2*ffreqs[(nw+nw3)/2-1]+w2+u},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+u},comp);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-u},comp);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1-u,2*ffreqs[(nw1+nw3)/2-1]-w1+u},comp);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-u},comp);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+u},comp);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+u,2*ffreqs[(nw+nw3)/2-1]+w1-u},comp);;
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_u_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+u},comp);
            double v1_K2b_t_low =max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-u},comp);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+u,2*ffreqs[(nw-nw3)/2]-w2-u},comp);;
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+u},comp);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-u},comp);
            double v1_R_s_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-u,2*ffreqs[(nw-nw3)/2]+w2+u},comp);


            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low = bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+u},comp);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-u},comp);
            double v2_R_t_low =  max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1-u,2*ffreqs[(nw1-nw3)/2]-w1+u},comp);
            double v2_K1_s_low = bfreqs[0]-w1;
            double v2_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-u},comp);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+u},comp);
            double v2_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+u,2*ffreqs[(nw-nw3)/2]+w1-u},comp);;
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_R_u_low = ffreqs[(nw1-nw3)/2];


            vector<double> upperR_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v_K2_u_up,v_R_u_up};
            sort(upperR_v1.begin(),upperR_v1.end());
            vector<double> lowerR_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v_K2_u_low,v_R_u_low};
            sort(lowerR_v1.begin(),lowerR_v1.end());

            vector<double> upperR_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v_K2_u_up,v_R_u_up};
            sort(upperR_v2.begin(),upperR_v2.end());
            vector<double> lowerR_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v_K2_u_low,v_R_u_low};
            sort(lowerR_v2.begin(),lowerR_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,w2,'u',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,wlimit,'u',2,h);

            for(int i=upperR_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperR_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperR_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperR_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperR_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperR_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperR_v1[i],'u',2,h))>0 ){limit_up1 = upperR_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};



            for(int i=0; i< lowerR_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerR_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerR_v1[i],'u',2,h))>0){limit_low1 = lowerR_v1[i];break;};
                if(i==lowerR_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerR_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerR_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerR_v2[i],w2,'u',1,h))>0){limit_low2 = lowerR_v2[i];break;};
                if(i==lowerR_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
               if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}


        else if(h=='L'){ //constant vertices



            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_R_u_low = ffreqs[(nw1-nw3)/2];
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_u_up = ffreqs[(nw1+nw3)/2-1];

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low = bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+u},comp);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-u},comp);
            double v2_R_t_low =  max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1-u,2*ffreqs[(nw1-nw3)/2]-w1+u},comp);
            double v2_K1_s_low = bfreqs[0]-w1;
            double v2_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-u},comp);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+u},comp);
            double v2_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+u,2*ffreqs[(nw-nw3)/2]+w1-u},comp);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+u},comp);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-u},comp);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1-u,2*ffreqs[(nw1+nw3)/2-1]-w1+u},comp);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-u},comp);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+u},comp);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+u,2*ffreqs[(nw+nw3)/2-1]+w1-u},comp);

            vector<double> upperK2_v1{v_K2_u_up};
            sort(upperK2_v1.begin(),upperK2_v1.end());
            vector<double> lowerK2_v1{v_K2_u_low};
            sort(lowerK2_v1.begin(),lowerK2_v1.end());

            vector<double> upperK2_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v_K2_u_up,v_R_u_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v_K2_u_low,v_R_u_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,wlimit,'u',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,wlimit,'u',2,h);


            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK2_v2[i];break;};
                if(i==0){mode2_up=0;};

            };

            for(int i=upperK2_v1.size()-1; i>-1 ; i--){

                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2_v1[i],'u',2,h))>0 ){limit_up1 = upperK2_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};

            for(int i=0; i< lowerK2_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2_v1[i],'u',2,h))>0){limit_low1 = lowerK2_v1[i];break;};
                if(i==lowerK2_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK2_v2[i];break;};
                if(i==lowerK2_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};

            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};
          }

        else if(h=='M'){//in this case, only Gamma' (vert1) sets the limits

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_R_u_low = ffreqs[(nw1-nw3)/2];

            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_u_up = ffreqs[(nw1+nw3)/2-1];


            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+u},comp);
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-u},comp);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+u,2*ffreqs[(nw+nw3)/2-1]-w2-u},comp);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+u},comp);
            double v1_K2b_s_up =  min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-u},comp);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-u,2*ffreqs[(nw+nw3)/2-1]+w2+u},comp);

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+u},comp);
            double v1_K2b_t_low =max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-u},comp);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+u,2*ffreqs[(nw-nw3)/2]-w2-u},comp);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+u},comp);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-u},comp);
            double v1_R_s_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-u,2*ffreqs[(nw-nw3)/2]+w2+u},comp);


            vector<double> upperK2b_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v_K2_u_up,v_R_u_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v_K2_u_low,v_R_u_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_u_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_u_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;

            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,w2,'u',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,u,wlimit,wlimit,'u',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2b_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2b_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2b_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2b_v1[i],'u',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;};
                if(i==0){mode1_up=0;};
            };


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2b_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2b_v1[i],'u',2,h))>0){limit_low1 = lowerK2b_v1[i];break;};
                if(i==lowerK2b_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2b_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2b_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK2b_v2[i];break;};
                if(i==lowerK2b_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
              if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h=='K'){ //in this case the only contribution to the non-constat part comes from K2,u


            vector<double> upperK1_v1{v_K2_u_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_u_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_u_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_u_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,wlimit,'u',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,u,wlimit,wlimit,'u',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK1_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK1_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK1_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK1_v1[i],'u',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK1_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK1_v1[i],'u',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK1_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK1_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
             if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};};

        double v_K2_u_low = ffreqs[(nw1-nw2)/2];
         double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];

        double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct ububble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};
        gsl_function F;

        if(mode==0 && abs(vert1_const * vert2_const)>0){


            if(reg==1 && p1=='g' && p2 =='g'){
                F.function = &ugreensfunc_re<T1,T2>;
                F.params = &params;


                double upper = abs(u/2) + Lambda;

                //integrand is even: sufficient to integrate positive part and multiply result by 2



                F.function = &ugreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);


                double compl_real = 2*resultgfu_re *vert1_const * vert2_const;

                B= compl_real ;
              }


            else if(reg==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;

                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;




                if(p1=='k'){
                    singsc = real(1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+u,se,dse,p2)) ;
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+u,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+u/2;dse_cutoff_low= ffreqs[0]+u/2;}

                else if(p2=='k'){
                    singsc =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s') );
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-u/2;dse_cutoff_low= ffreqs[0]-u/2;};

                struct ububble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;};


                    F.params = &params_mod;

                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                         double compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    B= singsc+ compl_real ;
                 }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};


                    F.params = &params_mod;

                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                          double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    B=singsc + compl_real ;
                 }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){
                    F.function = &ugreensfunc_re<T1,T2>;
                    F.params = &params_mod;


                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);

                    double compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                      B=singsc + compl_real ;
                 };}


            else if(reg==2 && (p1=='g' &&p2=='g')){

                F.function = &ugreensfunc_re<T1,T2>;
                  F.params = &params;


                  gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                        w, &result_re, &errorgfu_re);


                  double compl_real = 2*result_re *vert1_const * vert2_const;

                B=compl_real;


            }
            else if(reg==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0;
                if(p1=='s'){bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;}
                else if(p2=='s'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;}
                F.params = &params;
                F.function = &ugreensfunc_re<T1,T2>;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                double compl_real =(result_re) *vert1_const * vert2_const;

                B=compl_real;

            }

            else if(reg==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially

                double bound_low=0, bound_up=0, dse_cutoff_low=0, dse_cutoff_up=0;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;
                    dse_cutoff_low = ffreqs[0]+u/2; dse_cutoff_up = ffreqs[nw-1]+u/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;
                    dse_cutoff_low = ffreqs[0]-u/2; dse_cutoff_up = ffreqs[nw-1] -u/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct ububble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};
                F.params = &params_singsc;
                F.function = &ugreensfunc_re<T1,T2>;//integration of single scale
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct ububble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of katanin extension
                F.function = &ugreensfunc_re<T1,T2>;
                double result_re2=0,error_re2=0;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error, rel_error,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real =(result_re + result_re2 ) *vert1_const * vert2_const;
                  B=compl_real;
             };
        }

        else if(mode==1){


            if(reg==1 && p1=='g' && p2=='g'){
                 double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;
                if(lower <= limit_low){//left side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        double compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                        lhs=compl_real;

                    };}
                else if(lower <= limit_up){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                           compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    };

                    F.params = &params;
                    F.function = &ububble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     lhs=result_re + compl_real;

                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                             F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };
                    F.params = &params;
                    F.function = &ububble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     lhs=result_re + compl_real;


                };

                if(upper >= limit_up){//right side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                         double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                        rhs=compl_real;

                    };}
                else if(upper >= limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                         compl_real =(resultgfu_re  )*vert1_const * vert2_const;
                     };
                    F.params = &params;
                    F.function = &ububble_re<T1,T2>;

                    gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                     rhs=result_re + compl_real;

                }

                else if (upper < limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                     };
                    F.params = &params;
                    F.function = &ububble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;


                };

                B= rhs + lhs;


            }
            else if(reg==1 && ((p1=='k' &&p2=='g') || (p1=='g'&& p2 =='k'))){

                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;

                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;

                if(p1=='k'){
                    singsc = real( 1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+u,se,dse,p2) );
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+u,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+u/2;dse_cutoff_low= ffreqs[0]+u/2;}

                else if(p2=='k'){
                    singsc =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-u/2;dse_cutoff_low= ffreqs[0]-u/2;};
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct ububble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                         double compl_real = (resultgfl_re )*vert1_const * vert2_const;
                      lhs=compl_real;
                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                         compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &ububble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                       };

                    lhs=(result_re + compl_real);

                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_low <limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &ububble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                        };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                          compl_real += (resultgfl_re )*vert1_const * vert2_const;
                     };

                    lhs=(result_re + compl_real);


                };

                if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                    double low_eff=0;
                    if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                    F.params = &params_mod;
                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);
                       double compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    rhs=compl_real;
                }
                else if(upper >= limit_low ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                          compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &ububble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };

                    rhs=(result_re + compl_real);

                }

                else if(upper < limit_low){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;

                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                         compl_real =(resultgfu_re  )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_up >limit_low){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &ububble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                        };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                             compl_real +=(resultgfu_re )*vert1_const * vert2_const;

                    };

                    rhs=(result_re + compl_real);


                };

                B=rhs +lhs+ singsc;
             }

            else if(reg==2 && p1=='g' && p2=='g'){

                F.function = &ububble_re<T1,T2>;
                F.params = &params;

                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //integration of greens function outside of w-dependent interval

                double compl_real=0;
                  if(abs(vert1_const * vert2_const)>0){
                    F.function = ugreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &ugreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=(result_re + compl_real );
              }

            else if(reg==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){


                double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                if(p1=='s'){bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;}
                else if(p2=='s'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    B=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real =(resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=(result_re + compl_real );
                 }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &result_re, &error_re);
                     //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);


                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    B=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    B=(result_re  );
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };

                    B=( compl_real );
                }

            }

            else if(reg==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct ububble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    singsc=(result_re + compl_real );
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    singsc=(result_re + compl_real);
                 }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                           compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    singsc=( compl_real );
                 }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    singsc=(result_re  );
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    singsc=( compl_real );
                }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]+u/2; bound_up = ffreqs[nw-1]+u/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]-u/2; bound_up = ffreqs[nw-1]-u/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct ububble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                      if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re)*vert1_const * vert2_const;

                    };




                    kat=(result_re + compl_real );
                 }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=(result_re + compl_real );
                 }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                        //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){


                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=(result_re + compl_real );
                  }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                     if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);


                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=( compl_real );
                 }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &ububble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                      //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                 }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    //   double compl_imag;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &ugreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat=( compl_real );
                 }
                B=(kat + singsc);

            };


        };



          //cout << " ok" << endl;



    };
 if(abs(B)<1e-20){B=0;};

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
parvert<irreducible> operator+(parvert<irreducible> ,parvert<irreducible>);
parvert<svert> operator+=(parvert<svert> ,parvert<svert>);
parvert<tvert> operator+=(parvert<tvert> ,parvert<tvert>);
parvert<uvert> operator+=(parvert<uvert> ,parvert<uvert>);
parvert<irreducible> operator+=(parvert<irreducible> ,parvert<irreducible>);
parvert<svert> operator*(double  ,parvert<svert>&);
parvert<svert> operator*(parvert<svert>& ,double );
parvert<tvert> operator*(double  ,parvert<tvert>&);
parvert<tvert> operator*(parvert<tvert>& ,double );
parvert<uvert> operator*(double  ,parvert<uvert>&);
parvert<uvert> operator*(parvert<uvert>& ,double );
parvert<irreducible> operator*(double  ,parvert<irreducible>&);
parvert<irreducible> operator*(parvert<irreducible>& ,double );



//loop function:
template<typename  T>
struct loop_params{
    double Lambda;
    T& vert;
    char ptype;
    self& selfen;
    self& diffselfen;

    int a; int b; int c;


    double w1;

};


//}

template<typename T>
double loop_im(double w, void * p){
    struct loop_params<T> * params
            = static_cast< struct loop_params<T> *>(p);
    char ptype = (params ->ptype);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T& vert = (params-> vert);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);
    double w1 = (params->w1);
    return(1./(2*pi)* imag( -1. * (vert.vvalsmooth(a,b,c,w1,w,w1,'v')-vert.vvalsmooth(a,b,c,w1,wlimit,w1,'v') )*propag(Lambda,w,selfen,diffselfen,ptype)));//*exp(ci*w*1e-16)));
}

/**same loop function but with changed arguments in vertex (needed to parvert-loops), see SB1, p.76*/


template<typename T>
double loop_im_sw(double w, void * p){
    struct loop_params<T> * params
            = static_cast< struct loop_params<T> *>(p);
    char ptype = (params ->ptype);
    self& selfen = (params->selfen);
    self& diffselfen = (params->diffselfen);
    double Lambda = (params->Lambda);
    T& vert = (params-> vert);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);
    double w1 = (params->w1);
    return(1./(2*pi)* imag( -1. * (vert.vvalsmooth(a,b,c,w1,w,w,'v')-vert.vvalsmooth(a,b,c,w1,wlimit,wlimit,'v'))*propag(Lambda,w,selfen,diffselfen,ptype)));//*exp(ci*w*1e-16)));

}






template<class T>
//loop Attention: Only use for freqeuency-DEPENDENT vertices. Otherwise, the loop yields simply the fermionic density.
complex<double> loop_std(double Lambda, T& vertex, char p, self& se, self& dse,int a, int b, int c, double w1){
    double abs_error = 1e-4,rel_error=1e-4;
    complex<double> result(0.0,0.0);
    if((reg ==1 && p =='s') ){

        result = -1. *1./(2*pi)* vertex.vvalsmooth(a,b,c,w1,-Lambda,w1,'v')*propag(Lambda,-Lambda,se,dse,p);
        result += -1. *1./(2*pi)* vertex.vvalsmooth(a,b,c,w1,Lambda,w1,'v')*propag(Lambda,Lambda,se,dse,p);
        //the contribution from the constant part of the vertex cancels
    }

    else{
        int mode = 1;
        double limit_low;
        double limit_up;

        //upper bound:

        //conditions from vertex :

        double K1_s_up= bfreqs[nw1-1]-w1;
        double K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,w1-2*ffreqs[(nw1-nw2)/2]},comp);
        double R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,w1-2*ffreqs[(nw1-nw3)/2]},comp);

        double K1_u_up = bfreqs[nw1-1]+w1;
        double K2_u_up =  min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1},comp);
        double R_u_up =  min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1},comp);

        double K2_t_up = ffreqs[(nw1+nw2)/2-1];
        double R_t_up = ffreqs[(nw1+nw3)/2-1];

        //lower bound:
        //conditions from vertex :
        double K1_s_low= bfreqs[0]-w1;
        double K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,w1-2*ffreqs[(nw1+nw2)/2-1]},comp);
        double R_s_low = max({bfreqs[(nw1-nw3)/2]-w1,w1-2*ffreqs[(nw1+nw3)/2-1]},comp);


        double K1_u_low = bfreqs[0]+w1;
        double K2_u_low =  max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1},comp);
        double R_u_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1},comp);

        double K2_t_low = ffreqs[(nw1-nw2)/2];
        double R_t_low = ffreqs[(nw1-nw3)/2];

        vector<double> upper{K1_s_up,K2_s_up,R_s_up,K1_u_up,K2_u_up,R_u_up,K2_t_up,R_t_up};
        sort(upper.begin(),upper.end());
        vector<double> lower{K1_s_low,K2_s_low,R_s_low,K1_u_low,K2_u_low,R_u_low,K2_t_low,R_t_low};
        sort(lower.begin(),lower.end());

        complex<double> vert_const = vertex.vvalsmooth(a,b,c,w1,wlimit,w1,'v');

        int mode_up =1;
        int mode_low =1;

        for(int i=upper.size()-1; i>-1 ; i--){
            if(abs(vertex.vvalsmooth(a,b,c,w1,upper[i],w1,'v')-vert_const)>0 ){limit_up = upper[i];break;};
            if(i==0){mode_up=0;};
        };

        for(int i=0; i< lower.size() ; i++){
            if(abs(vertex.vvalsmooth(a,b,c,w1,lower[i],w1,'v')-vert_const)>0){limit_low = lower[i];break;};
            if(i==lower.size()-1){mode_low=0;};
        };


        if(mode_low ==0 && mode_up ==0){mode =0;}; // in this case, the differentiated self energy has no contribution (contant contributions are not computed here)

         gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        double result_re, error_re, result_im, error_im;


        struct loop_params<T> params = {Lambda,vertex, p,se, dse, a,  b,  c,  w1};
        gsl_function F;



        if(mode==1){

            if(p =='s' && reg ==2){

                double upper_bound = 7*Lambda;
                double lower_bound = -7*Lambda;

                if(lower_bound < limit_up && upper_bound > limit_low){
                    double limit_low_eff = max({limit_low,lower_bound},comp);
                    double limit_up_eff = min({limit_up,upper_bound},comp);

                    F.params = &params;
                    F.function = &loop_im<T>;
                    gsl_integration_qag(&F,limit_low_eff,limit_up_eff,abs_error, rel_error,1500,2,
                                        w, &result_im, &error_im);



                    result.imag( result_im );};//loop has only imaginary part



            };};
        gsl_integration_workspace_free(w);};


    return result;

}




template<class T>
//loop Attention: Only use for freqeuency-DEPENDENT vertices. Otherwise, the loop yields simply the fermionic density.
complex<double> loop_sw(double Lambda, T& vertex, char p, self& se, self& dse,int a, int b, int c, double w1){
    double abs_error = 1e-4,rel_error=1e-4;
    complex<double> result(0.0,0.0);
    if((reg ==1 && p =='s') ){

        result += -1. *1./(2*pi)* vertex.vvalsmooth(a,b,c,w1,-Lambda,-Lambda,'v')*propag(Lambda,-Lambda,se,dse,p);
        result += -1. *1./(2*pi)* vertex.vvalsmooth(a,b,c,w1,Lambda,Lambda,'v')*propag(Lambda,Lambda,se,dse,p);
    }

    else{
        int mode =1;
        double limit_low;
        double limit_up;

        //upper bound:



        double K1_s_up= bfreqs[nw1-1]-w1;
        double K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,w1-2*ffreqs[(nw1-nw2)/1]},comp);
        double R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,w1-2*ffreqs[(nw1-nw3)/1]},comp);

        double K1_t_up = bfreqs[nw1-1]+w1;
        double K2_t_up =  min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1},comp);
        double R_t_up =  min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1},comp);

        double K2_u_up = ffreqs[(nw1+nw2)/2-1];
        double R_u_up = ffreqs[(nw1+nw3)/2-1];

        //lower bound:
        //conditions from vertex :
        double K1_s_low= bfreqs[0]-w1;
        double K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,w1-2*ffreqs[(nw1+nw2)/2-1]},comp);
        double R_s_low = max({bfreqs[(nw1-nw3)/2]-w1,w1-2*ffreqs[(nw1+nw3)/2-1]},comp);


        double K1_t_low = bfreqs[0]+w1;
        double K2_t_low =  max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1},comp);
        double R_t_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1},comp);

        double K2_u_low = ffreqs[(nw1-nw2)/2];
        double R_u_low = ffreqs[(nw1-nw3)/2];









        vector<double> upper{K1_s_up,K2_s_up,R_s_up,K1_t_up,K2_t_up,R_t_up,K2_u_up,R_u_up};
        sort(upper.begin(),upper.end());
        vector<double> lower{K1_s_low,K2_s_low,R_s_low,K1_t_low,K2_t_low,R_t_low,K2_u_low,R_u_low};
        sort(lower.begin(),lower.end());

        complex<double> vert_const = vertex.vvalsmooth(a,b,c,w1,wlimit,wlimit,'v');


        int mode_up;
        int mode_low;
        for(int i=upper.size()-1; i>-1 ; i--){
            if(abs(vertex.vvalsmooth(a,b,c,w1,upper[i],upper[i],'v')-vert_const)>0 ){limit_up = upper[i];break;};
            if(i==0){mode_up=0;};
        };

        for(int i=0; i< lower.size() ; i++){
            if(abs(vertex.vvalsmooth(a,b,c,w1,lower[i],lower[i],'v')-vert_const)>0){limit_low = lower[i];break;};
            if(i==lower.size()-1){mode_low=0;};
        };


        if(mode_low ==0 && mode_up ==0){mode =0;};


        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        double result_re, error_re, result_im, error_im;


        struct loop_params<T> params = {Lambda,vertex, p,se, dse, a,  b,  c,  w1};
        gsl_function F;

        

        if(mode==1){

            if(p =='s' && reg ==2){


                double upper_bound = 7*Lambda;
                double lower_bound = -7*Lambda;

                if(lower_bound < limit_up && upper_bound > limit_low){
                    double limit_low_eff = max({limit_low,lower_bound},comp);
                    double limit_up_eff = min({limit_up,upper_bound},comp);


                    F.function = &loop_im<T>;
                    F.params = &params;
                    gsl_integration_qag(&F,limit_low_eff,limit_up_eff,abs_error, rel_error,1500,2,
                                        w, &result_im, &error_im);


                    result.imag( result_im );};

            };};
        gsl_integration_workspace_free(w);};


    return result;

}





template<class T>
self loop(double Lambda, parvert<T>& vertex, char p, self& se, self& dse){//see SB1, p. 76 for details

    self self1, self2, self3;
    

    //**************first contr.********//
#pragma omp parallel
    for(int i=0; i<nw; i++){//for all different external freqs..
        complex<double> k= zero,l=zero,m=zero;


        for(int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1;a++){
            for(int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1;b++){
                for(int c= 1;c<4;c++){//..iterate through all sites..
                    if(distance(a,b,c) <= d_c){//cutoff distance


                        k += loop_std(Lambda,vertex.densvertex,p,se,dse,a,b,c,ffreqs[i]);


                    };};};};

        self1.setself(i,2.*k);//factor 2 is a combinatorial factor

        //**************second contr.********//

        l = loop_sw(Lambda,vertex.spinvertex,p,se,dse,0,0,1,ffreqs[i]);
        //cout << "2ok" << endl;
        self2.setself(i,-3.*l);


        //**************third contr.********//


        m = loop_sw(Lambda,vertex.densvertex,p,se,dse,0,0,1,ffreqs[i]);

        self3.setself(i,-1.*m);
    };//closes loop over freqs i.



    return (self1 + self2 + self3);
}



//define bubble functions with parvert as argument and output:


//for site-spin-freq parametrization from Reuther, define bubble function of vertices of the form (Gamma_spin *(Pauli)^2 + Gamma_dens * del^2). See Sb, p. 71
template<class T1, class T2> //last arguments defindes the object of type parvert<T> to which the result is written
parvert<svert> sbubble(int red_side, double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:

    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;//concerns bos. freq.
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if( sym== 2){sum_K2_i = nw/2;};

    int sum_K2_j;//concerns ferm. freq
    if(sym==0 ){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);




    parvert<svert> result;
    vector<site> upper_sites;

    int number=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};

#pragma omp parallel for
    for(int m=0; m<number; m++){

        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        int a = upper_sites[m].a;
        int b = upper_sites[m].b;
        int c = upper_sites[m].c;
        for(int i=sum_K1 ; i< (nw+nw1)/2 ; i++){


            double spinspin;
            double spindens;
            double densspin;
            double densdens;


            //note that the asymptotic behaviour is centered around s/2 (half the bosonic transfer freq in the s-channel)
            double s = bfreqs[i];

            spinspin = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
            spindens = sbubble(red_side,0,0,w, Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');
            densspin = sbubble(red_side,0,0,w, Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');
            densdens = sbubble(red_side,0,0,w, Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');
                     //set spin vertex. Combinatorial factors, see SBI
            result.spinvertex.K1_setvert(a,b,c,i,-4. * spinspin + 2. * densspin + 2. * spindens);//combine vertices with and without tilde since u- and t-channel always appear as a sum

            //set density vertex
             result.densvertex.K1_setvert(a,b,c,i,6. * spinspin + 2. * densdens);

        };

        // K2 and K2b



        for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){


                double spinspin;
                double spindens;
                double densspin;
                double densdens;


                double s = bfreqs[i];
                double  w1 = ffreqs[j];



                spinspin = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                spindens = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');
                densspin = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');
                densdens = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');                       //set spin vertex. Combinatorial factors, see SBI
                      result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-4. * spinspin + 2. * densspin + 2. * spindens);
                result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,6. * spinspin + 2. * densdens);


                if(sym==0){//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                    double bspinspin;
                    double bspindens;
                    double bdensspin;
                    double bdensdens;

                    bspinspin = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                    bspindens = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                    bdensspin = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                    bdensdens = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                      result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-4. * bspinspin + 2. * bdensspin + 2. * bspindens );
                     result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,6. * bspinspin + 2. * bdensdens );
                };
            };
        };

        //R (rest function):


//                for(int i=sum_R ; i<(nw+nw3)/2; i++){
//                    for(int j=sum_R ; j<(nw+nw3)/2; j++){
//                        for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble

//                            if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){

//                                double spinspin;
//                                double spindens;
//                                double densspin;
//                                double densdens;

//                                double s = bfreqs[i];
//                                double w1 = ffreqs[j];
//                                double w2 = ffreqs[k];

//                                spinspin = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
//                                spindens = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
//                                densspin = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
//                                densdens = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');

//                                //set spin vertex. Combinatorial factors, see SBI

//                                result.spinvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2  ,-4. * spinspin + 2. * densspin + 2. * spindens);
//                                //set density vertex



//                                result.densvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2 ,6. * spinspin + 2. * densdens  );

//                            };};};};

        gsl_integration_workspace_free(w);
    };

    return result;
}

//****************three types of t-bubbles in full parametrization ***********************
//NOTE: the last argument is a parvert<T> vertex to which is result is written in the ZEROTH component
template<class T1, class T2>//this is the RPA term, compare to SB1, page 73.
parvert<tvert> tbubble1(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;
    //mixed combinations like spindens or densspin do not contribute since they cancel out (see Sb1, p. 66)

#pragma omp parallel for
    //  K1:
    for(int i=sum_K1 ; i< (nw+nw1)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){

                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspinval, densdensval;
                        double t = bfreqs[i] ;
                        for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.


                            if(f==1){
                                for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                spinspinval += tbubble(red_side, 0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                         };};};};}

                            else if(f==2){
                                for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        int ap = a-d;
                                        int bp = b-e;//translate such that j lies in first Wigner Seitz cell

                                        int app = -ap-bp;
                                        int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                            };};};};}
                            else if(f==3){
                                for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        int ap = a - d;
                                        int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                        int app = bp;
                                        int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                            };

                                        };};};};
                        };
                        result.spinvertex.K1_setvert(a,b,c,i,2.*spinspinval);
                        result.densvertex.K1_setvert(a,b,c,i-(nw-nw1)/2,2.*densdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                        spinspinval = 0;//reset values for next iteration step with different frequencies
                        densdensval = 0;
                    };

                };};};
        gsl_integration_workspace_free(w);
    };


    //   K2 and K2b:


#pragma omp parallel for

    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){

                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspinval, densdensval;
                            double t = bfreqs[i];
                            double  w1 = ffreqs[j];
                            for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                if(f==1){
                                    for(int d =-(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};};};}

                                else if(f==2){
                                    for(int d =-(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a- d;
                                            int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                            int app = -ap-bp;
                                            int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                       spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};};};}
                                else if(f==3){
                                    for(int d = -(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a - d;
                                            int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                            int app = bp;
                                            int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};

                                        };};};
                            };

                            result.spinvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*spinspinval );
                            result.densvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*densdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                            spinspinval = 0;//reset values for next iteration step with different frequencies
                            densdensval =0;


                            //K2b (interchange k and j in bubble calculations):


                            if(sym==0){//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                                double bspinspinval, bdensdensval;
                                for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                    if(f==1){
                                        for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                            for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                if(distance(d,e,f) <= d_c){//cutoff distance
                                                    if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                        bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                        bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    };};};};}

                                    else if(f==2){
                                        for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                            for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                int ap = a- d;
                                                int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                                int app = -ap-bp;
                                                int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                                if(distance(d,e,f) <= d_c){//cutoff distance
                                                    if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                        bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                        bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    };};};};}
                                    else if(f==3){
                                        for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                            for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                int ap = a - d;
                                                int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                                int app = bp;
                                                int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                                if(distance(d,e,f) <= d_c){//cutoff distance
                                                    if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                        bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                        bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    };};

                                            };};};
                                };
                                result.spinvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bspinspinval);//
                                result.densvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bdensdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                                bspinspinval = 0;//reset values for next iteration step with different frequencies
                                bdensdensval = 0;


                            };};};

                };};};
        gsl_integration_workspace_free(w);
    };

    // R (rest function):

#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){

                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspinval, densdensval;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];

                                    double w2 = ffreqs[k];

                                    for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                        if(f==1){
                                            for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};};};}

                                        else if(f==2){
                                            for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    int ap = a- d;
                                                    int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                                    int app = -ap-bp;
                                                    int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                    int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};};};}
                                        else if(f==3){
                                            for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    int ap = a - d;
                                                    int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                                    int app = bp;
                                                    int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                    int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};

                                                };};};
                                    };
                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspinval );
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*densdensval  );

                                    spinspinval = 0;//reset values for next iteration step with different frequencies
                                    densdensval = 0;


                                };};};};};};};
        gsl_integration_workspace_free(w);
    };





    return result;
}

template<class T1, class T2>
parvert<tvert> tbubble2(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2,char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;



#pragma omp parallel for
    //K1:
    for(int i=sum_K1; i< (nw+nw1)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){


                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspin, spindens, densspin, densdens;
                        double t = bfreqs[i];
                        spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                        spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                        result.spinvertex.K1_setvert(a,b,c,i,spinspin - densspin);
                        result.densvertex.K1_setvert(a,b,c,i-(nw-nw1)/2,-3.*spindens -1.*densdens);
                                };
                };};};
        gsl_integration_workspace_free(w);
    };

#pragma omp parallel for
    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_K2_j; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){


                        //K2 and K2b:

                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspin, spindens, densspin, densdens;
                            double t = bfreqs[i];
                            double w1 = ffreqs[j] ;
                            spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                            spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                            result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*spinspin -1.*densspin );
                            result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*spindens -1.*densdens);

                            if(sym==0){//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                                double bspinspin, bspindens, bdensspin, bdensdens;
                                bspinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                                bspindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');

                                result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*bspinspin -1.*bdensspin );
                                result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*bspindens -1.*bdensdens);
                              };};};
                };};};
        gsl_integration_workspace_free(w);
    };

#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in t-bubble
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){
                                //rest function:

                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspin, spindens, densspin, densdens;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
                                    spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,1.*spinspin -1.*densspin);
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,-3.*spindens -1.*densdens);
                                   };};};

                    };};};};
        gsl_integration_workspace_free(w);
    };

    return result;
}

template<class T1, class T2>
parvert<tvert> tbubble3(int red_side, double Lambda, parvert<T1>& vert1, parvert<T2>& vert2,char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;



#pragma omp parallel for

    //K1:
    for(int i=sum_K1; i< (nw+nw1)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){


                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspin, spindens, densspin, densdens;
                        double t = bfreqs[i];
                        spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                        spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                        result.spinvertex.K1_setvert(a,b,c,i,spinspin - spindens);
                        result.densvertex.K1_setvert(a,b,c,i,-3.*densspin -1.*densdens);
                       };
                };};};
        gsl_integration_workspace_free(w);
    };
      //K2 and K2b:

#pragma omp parallel for
    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){



                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspin, spindens, densspin, densdens;
                            double t = bfreqs[i];
                            double w1 = ffreqs[j] ;
                            spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                            spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                            result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*spinspin -1.*spindens );
                            result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*densspin -1.*densdens);

                            if(sym==0){//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                                double bspinspin, bspindens, bdensspin, bdensdens;
                                bspinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                                bspindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');

                                result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*bspinspin -1.*bspindens);
                                result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*bdensspin -1.*bdensdens );
                                    };};};
                };};};
        gsl_integration_workspace_free(w);
    };

    //R:


#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){



                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspin, spindens, densspin, densdens;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
                                    spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,1.*spinspin -1.*spindens );
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,-3.*densspin -1.*densdens );

                                };};};};};};};
        gsl_integration_workspace_free(w);
    };

    return result;
}

//template<class T1, class T2>
//parvert<tvert> tbubble(double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
//    return(tbubble1(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + tbubble2(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + tbubble3(Lambda, vert1, vert2,  p1,p2,selfenergy, diffselfenergy) );
//}





template<class T1, class T2>//this is the RPA term, compare to SB1, page 73.
parvert<tvert> tbubble(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:

    vector<site> all_sites,upper_sites;

    int number_all=0, number_upper=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= -(nuc_eff-1)/2; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){

                if( distance(a,b,c) <= d_c){
                    number_all +=1;
                    site site_here(a,b,c);
                    all_sites.push_back(site_here);};};};};


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){

                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number_upper +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};



    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;
    //mixed combinations like spindens or densspin do not contribute since they cancel out (see Sb1, p. 66)

#pragma omp parallel for
    for(int m=0; m<number_upper;m++){
        int a = upper_sites[m].a;
        int b = upper_sites[m].b;
        int c = upper_sites[m].c;

        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);

        //K1 class
        for(int i=sum_K1 ; i< (nw+nw1)/2; i++){
            double t = bfreqs[i] ;


            //first diagrammatic type (RPA)
            double spinspinval=0, densdensval=0;
            double spinspin_2=0,spindens_2=0,densspin_2=0,densdens_2=0;double spinspin_3=0,spindens_3=0,densspin_3=0,densdens_3=0;

            for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                int d = all_sites[s].a;
                int e = all_sites[s].b;
                int f = all_sites[s].c;

                if(f==1 && distance(a-d,b-e,c) <= d_c){//cutoff distance

                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                }

                else if(f==2){

                    int ap = a-d;
                    int bp = b-e;//translate such that j lies in first Wigner Seitz cell

                    int app = -ap-bp;
                    int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                    int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)


                    if(distance(app,bpp,cp) <= d_c){//cutoff distance
                        spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    };}
                else if(f==3){

                    int ap = a - d;
                    int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                    int app = bp;
                    int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                    int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                    if(distance(app,bpp,cp) <= d_c){//cutoff distance

                        spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    };

                };
            };



            //second diagrammatic type

            spinspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the fisrt vertex is taken with the distance vector 0.(interaction of a site with itself)
            spindens_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
            densspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
            densdens_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

            //third diagrammatic type

            spinspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
            spindens_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
            densspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
            densdens_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');


            result.spinvertex.K1_setvert(a,b,c,i,2.*spinspinval+spinspin_2-densspin_2+spinspin_3-spindens_3);
            result.densvertex.K1_setvert(a,b,c,i-(nw-nw1)/2,2.*densdensval-3*spindens_2-densdens_2-3*densspin_3-densdens_3);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed


        };


        //K2-class:
        for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){

                double spinspinval=0, densdensval=0;
                double bspinspinval, bdensdensval;
                double spinspin_2=0,spindens_2=0,densspin_2=0,densdens_2=0;double spinspin_3=0,spindens_3=0,densspin_3=0,densdens_3=0;
                double bspinspin_2=0,bspindens_2=0,bdensspin_2=0,bdensdens_2=0;double bspinspin_3=0,bspindens_3=0,bdensspin_3=0,bdensdens_3=0;


                double t = bfreqs[i];
                double  w1 = ffreqs[j];
                //first diagrammatic contribution (RPA):
                for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                    int d = upper_sites[s].a;
                    int e = upper_sites[s].b;
                    int f = upper_sites[s].c;
                    if(f==1 && distance(a-d,b-e,c) <= d_c){//cutoff distance
                        spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                        densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                        if(sym==0){
                            bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');

                        };}

                    else if(f==2){

                        int ap = a- d;
                        int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                        int app = -ap-bp;
                        int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                        int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            if(sym==0){
                                bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            };
                        };}
                    else if(f==3){

                        int ap = a - d;
                        int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                        int app = bp;
                        int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                        int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            if(sym==0){
                                bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            };
                        };};



                };

                //second diagrammatic contribution

                spinspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                  spindens_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                densspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                densdens_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                if(sym==0){
                    bspinspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                    bspindens_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                    bdensspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                    bdensdens_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');};

                //third diagrammatic contribution
                spinspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                spindens_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                densspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                densdens_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                if(sym==0){
                    bspinspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                    bspindens_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                    bdensspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                    bdensdens_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');};


                if(sym==0){
                    result.spinvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bspinspinval+bspinspin_2-bdensspin_2+1.*bspinspin_3 -1.*bspindens_3);//
                    result.densvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bdensdensval-3*bspindens_2 - bdensdens_2- 3.*bdensspin_3 -1.*bdensdens_3);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed
                };

                result.spinvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*spinspinval + spinspin_2 - densspin_2 + 1.*spinspin_3 - 1.*spindens_3 );
                result.densvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*densdensval - 3.*spindens_2 - densdens_2 - 3.*densspin_3 - 1.*densdens_3);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed


            };};





        // R (rest function):

//                for(int i=sum_R ; i<(nw+nw3)/2; i++){
//                    for(int j=sum_R ; j<(nw+nw3)/2; j++){
//                        for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){
//                            if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
//                                double t = bfreqs[i];
//                                double w1 = ffreqs[j];
//                                double w2 = ffreqs[k];


//                                double spinspinval=0, densdensval=0;
//                                double spinspin_2=0,spindens_2=0,densspin_2=0,densdens_2=0;double spinspin_3=0,spindens_3=0,densspin_3=0,densdens_3=0;
//                                //first diagrammatic contribution (RPA):

//                                for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
//                                    int d = upper_sites[s].a;
//                                    int e = upper_sites[s].b;
//                                    int f = upper_sites[s].c;
//                                    if(f==1){
//                                        if(distance(a-d,b-e,c) <= d_c){//cutoff distance
//                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                        };}

//                                    else if(f==2){

//                                        int ap = a- d;
//                                        int bp = b- e;//translate such that j lies in first Wigner Seitz cell

//                                        int app = -ap-bp;
//                                        int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
//                                        int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

//                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

//                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                        };}
//                                    else if(f==3){

//                                        int ap = a - d;
//                                        int bp = b - e;//translate such that j lies in first Wigner Seitz cell
//                                        int app = bp;
//                                        int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
//                                        int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)


//                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance
//                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


//                                        };};
//                                };

//                                //second diagrammatic conmtribution
//                                spinspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
//                                spindens_2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                densspin_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                densdens_2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

//                                //third diagrammatic contribution
//                                spinspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
//                                spindens_3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                densspin_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
//                                densdens_3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


//                                result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspinval + spinspin_2 - densspin_2 + 1.*spinspin_3 - 1.*spindens_3 );
//                                result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*densdensval - 3*spindens_2 - densdens_2 - 3.*densspin_3 - 1.*densdens_3 );

//                                spinspinval = 0;//reset values for next iteration step with different frequencies
//                                densdensval = 0;
//                                spinspin_2=0;spindens_2=0;densspin_2=0;densdens_2=0; spinspin_3=0;spindens_3=0;densspin_3=0;densdens_3=0;

//                            };};};};
        gsl_integration_workspace_free(w);
    };
    return result;
}








//*******************************************************************************************
template<class T1, class T2>
parvert<uvert> ububble(int reg_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){


    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?0:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<uvert> result;
    vector<site> upper_sites;

    int number=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};



#pragma omp parallel for
    for(int m=0; m<number ; m++){
        gsl_integration_workspace * w
                = gsl_integration_workspace_alloc (1500);
        int a = upper_sites[m].a;
        int b = upper_sites[m].b;
        int c = upper_sites[m].c;
        for(int i=sum_K1; i< (nw+nw1)/2; i++){

            //K1:



            double spinspin, spindens, densspin, densdens;
            double u = bfreqs[i] ;
            //compute all four possible vertex combinations in s-bubble
            spinspin = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
            spindens = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
            densspin = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
            densdens = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
            result.spinvertex.K1_setvert(a,b,c,i,2.*spinspin + densspin + spindens);
            result.densvertex.K1_setvert(a,b,c,i, 3.* spinspin + densdens);
        };


        //K2 and K2b:



        for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){

                double spinspin, spindens, densspin, densdens;
                double u = bfreqs[i];
                double w1 = ffreqs[j] ;
                spinspin = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                spindens = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                densspin = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                densdens = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,2.*spinspin  +  densspin + spindens);
                result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2, 3.* spinspin + densdens );

                if(sym==0){//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                    double bspinspin, bspindens, bdensspin, bdensdens;

                    bspinspin = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                    bspindens = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                    bdensspin = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                    bdensdens = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                    result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,2.*bspinspin +  bdensspin + bspindens );
                    result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2, 3.* bspinspin + bdensdens );


                };};};





//                for(int i=sum_R ; i<(nw+nw3)/2; i++){
//                    for(int j=sum_R ; j<(nw+nw3)/2; j++){
//                        for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble
//                            if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){


//                                double spinspin, spindens, densspin, densdens;
//                                double u = bfreqs[i];
//                                double w1 = ffreqs[j];
//                                double w2 = ffreqs[k];
//                                spinspin = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
//                                spindens = ububble(reg_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
//                                densspin = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
//                                densdens = ububble(reg_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
//                                result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspin +  densspin + spindens );
//                                result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2, 3.* spinspin + densdens );
//                                //Note: The result i written into th Lambda=0 component since the bubble functions are called within the function evaluate where the vertex "result" only contains one Lambda-layer
//                            };};};};

        gsl_integration_workspace_free(w);
    };





    return result;
}



class Susc{
    vector<vector<vector<double> >  >   Sus =
            vector<vector<vector<double>  > >
            (nuc_eff,vector<vector<double>  >
             ((nuc_eff+1)/2, vector<double>(3)
              ));//three atoms per unit cell

public:
    double vval(int, int, int);//arguments 1-3: site, argument 4, bosonic frequency
    void write(int,int,int, double );//first argument: Lambda
    void add_write(int,int,int, double );//first argument: Lambda
};


Susc suscept(double , parvert<fullvert>& ,  self );



//define a struct object which includes the self energy and the vertex which are needed to evaluate the RHS of the flow equations.

struct state{
    double Lambda;
    self selfenergy;
    Susc sus;
    parvert<fullvert>  vertex;};


state operator+(state, state);
state operator*(double,state);
state operator*(state, double);






void write_hdf(const H5std_string FILE_NAME,int Lambdas_it, long Lambda_size,state& vertex);
void add_hdf(const H5std_string FILE_NAME,int Lambdas_it, long Lambda_size,state& vertex);
state read_hdf(const H5std_string FILE_NAME,int Lambda,long Lambda_size);


#endif
