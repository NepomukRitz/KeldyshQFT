//
// Created by Sa.Aguirre on 7/18/19.
//

// TODO: documentation!!

#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <vector>

#include "parameters.h"
#include "data_structures.h"

using namespace std;

/*******************************CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX*********************************/
template <typename Q>
class avert{
public:

    //K1:
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wa * n_in);

    //K2:
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wa * nw2_nua * n_in);

    //K3:
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in);



//    avert() {
//        K1 = vector<Q> (nw1_wa);
//        K2 = vector<vector<Q> > (nw2_wa, vector<Q> (nw2_nua));
//        K3 = vector<vector<vector<Q> > > (nw3_wa, vector<vector<Q> > (nw3_nua, vector<Q> (nw3_nuap)));
//    }
//
//    avert(int input[]) {
//        K1 = vector<Q> (nw1_wa, Q(input));
//        K2 = vector<vector<Q> > (nw2_wa, vector<Q> (nw2_nua, Q(input)));
//        K3 = vector<vector<vector<Q> > > (nw3_wa, vector<vector<Q> > (nw3_nua, vector<Q> (nw3_nuap, Q(input))));
//    }


/* TODO: check all member functions
public:
    double vvalsmooth(int, int, int, double, double, double, char);
    double vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int, int, double);
    void K3_setvert( int, int, int, int, int, int, double);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int, int);
    double K3_vval(int, int, int, int, int, int);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double, double);
    double K3_vvalsmooth(int, int, int, double, double, double);
    void bust(char);
    friend avert operator*(double alpha, const avert& vertex);
    friend avert operator*(const avert& vertex, double alpha);
    friend avert operator+(const avert& vertex1, const avert& vertex2);

*/
};

template <class Q>
class pvert{
public:

    //K1:
    //vector<Q> K1;
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wp * n_in);
//    vector<vector<vector<vector<double > > > >  K1 =
//    vector<vector<vector<vector<double > > > >
//    (nuc_eff,vector<vector<vector<double > > >
//             ((nuc_eff+1)/2, vector<vector<double > >
//                             (3,vector<double >(nw1_q))));//three atoms per unit cell

    //K2:
    //vector<vector<Q> > K2;
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);
//    vector<vector<vector<vector<vector<double > > > > > K2 =
//    vector<vector<vector<vector<vector<double > > > > >
//    (nuc_eff,vector<vector<vector<vector<double > > > >
//             ((nuc_eff+1)/2, vector<vector<vector<double > > >
//                             (3,vector<vector<double > >//three atoms per unit cell
//                                (nw2_q, vector<double >(nw2_w1)))));


    //K3:
    //vector<vector<vector<Q> > > K3;
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);
//    vector<vector<vector<vector<vector<vector<double > > > > > > K3 =
//    vector<vector<vector<vector<vector<vector<double > > > > > >
//    (nuc_eff,vector<vector<vector<vector<vector<double > > > > >
//             ((nuc_eff+1)/2, vector<vector<vector<vector<double > > > >
//                             (3,vector<vector<vector<double > > >//three atoms per unit cell
//                                (nw3_q, vector<vector<double > >//nw frequency entries for each of the three frequencies
//                                        (nw3_w1, vector<double >(nw3_w2))))));


//    pvert() {
//        K1 = vector<Q> (nw1_wp);
//        K2 = vector<vector<Q> > (nw2_wp, vector<Q> (nw2_nup));
//        K3 = vector<vector<vector<Q> > > (nw3_wp, vector<vector<Q> > (nw3_nup, vector<Q> (nw3_nupp)));
//    }
//
//    pvert(int input[]) {
//        K1 = vector<Q> (nw1_wp, Q(input));
//        K2 = vector<vector<Q> > (nw2_wp, vector<Q> (nw2_nup, Q(input)));
//        K3 = vector<vector<vector<Q> > > (nw3_wp, vector<vector<Q> > (nw3_nup, vector<Q> (nw3_nupp, Q(input))));
//    }



/* TODO: check all member functions
public:

    double vvalsmooth(int, int, int, double, double, double, char);
    double vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int, int, double);
    void K3_setvert( int, int, int, int, int, int, double);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int, int);
    double K3_vval(int, int, int, int, int, int);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double, double);
    double K3_vvalsmooth(int, int, int, double, double, double);

    friend pvert operator*(double alpha, const pvert& vertex);
    friend pvert operator*(const pvert& vertex, double alpha);
    friend pvert operator+(const pvert& vertex1, const pvert& vertex2);

*/
};

template <class Q>
class tvert{

    //K1:
    //vector<Q> K1;
    vec<Q> vec_K1 = vec<Q> (nK_K1 * nw1_wt * n_in);
//    vector<vector<vector<vector<double > > > >  K1 =
//    vector<vector<vector<vector<double > > > >
//    (nuc_eff,vector<vector<vector<double > > >
//             ((nuc_eff+1)/2, vector<vector<double > >
//                             (3,vector<double >(nw1_q))));//three atoms per unit cell

    //K2:
    //vector<vector<Q> > K2;
    vec<Q> vec_K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);
//    vector<vector<vector<vector<vector<double > > > > > K2 =
//    vector<vector<vector<vector<vector<double > > > > >
//    (nuc_eff,vector<vector<vector<vector<double > > > >
//             ((nuc_eff+1)/2, vector<vector<vector<double > > >
//                             (3,vector<vector<double > >//three atoms per unit cell
//                                (nw2_q, vector<double >(nw2_w1)))));


    //K3
    //vector<vector<vector<Q> > > K3;
    vec<Q> vec_K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);
//    vector<vector<vector<vector<vector<vector<double > > > > > > K3 =
//    vector<vector<vector<vector<vector<vector<double > > > > > >
//    (nuc_eff,vector<vector<vector<vector<vector<double > > > > >
//             ((nuc_eff+1)/2, vector<vector<vector<vector<double > > > >
//                             (3,vector<vector<vector<double > > >//three atoms per unit cell
//                                (nw3_q, vector<vector<double > >//nw frequency entries for each of the three frequencies
//                                        (nw3_w1, vector<double >(nw3_w2))))));

public:
    Q K1(int i_K, int i_w1) {
        return vec_K1[i_K*nw1_wt + i_w1];
    }
    Q K2(int i_K, int i_w1, int i_w2) {
        return vec_K2[i_K*nw2_wt*nw2_nut + i_w1*nw2_nut + i_w2];
    }
    Q K3(int i_K, int i_w1, int i_w2, int i_w3) {
        return vec_K3[i_K*nw3_wt*nw3_nut*nw3_nutp + i_w1*nw3_nut*nw3_nutp + i_w2*nw3_nutp + i_w3];
    }
    // TODO: repeat for channels a, p

//    Q K1(int i_K, int i_w1, int i_in[]) {
//        int index = 0;
//        for (int i=0; i<d_in; ++i) {
//            int mult = 1;
//            for (int j=i+1; j<d_in; ++j) {
//                mult *= j;
//            }
//            index += mult * i_in[i];
//        }
//        return vec_K1[i_K*nw1_wt*i_in + i_w1*n_in + i_in];
//    }


//    tvert() {
//        K1 = vector<Q> (nw1_wt);
//        K2 = vector<vector<Q> > (nw2_wt, vector<Q> (nw2_nut));
//        K3 = vector<vector<vector<Q> > > (nw3_wt, vector<vector<Q> > (nw3_nut, vector<Q> (nw3_nutp)));
//    }
//
//    tvert(int input[]) {
//        K1 = vector<Q> (nw1_wt, Q(input));
//        K2 = vector<vector<Q> > (nw2_wt, vector<Q> (nw2_nut, Q(input)));
//        K3 = vector<vector<vector<Q> > > (nw3_wt, vector<vector<Q> > (nw3_nut, vector<Q> (nw3_nutp, Q(input))));
//    }



/* TODO: check all member functions
public:
    double vvalsmooth(int, int, int, double, double, double, char);
    double vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    double vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, double);
    void K2_setvert( int, int, int, int, int, double);
    void K3_setvert( int, int, int, int, int, int, double);
    double K1_vval(int, int, int, int);
    double K2_vval(int, int, int, int, int);
    double K3_vval(int, int, int, int, int, int);
    double K1_vvalsmooth(int, int, int, double);
    double K2_vvalsmooth(int, int, int, double, double);
    double K3_vvalsmooth(int, int, int, double, double, double);
    friend tvert operator*(double alpha, const tvert& vertex);
    friend tvert operator*(const tvert& vertex, double alpha);
    friend tvert operator+(const tvert& vertex1, const tvert& vertex2);

*/
};


template <class Q>
class irreducible{
public:
    Q U_bare;
//    vector<vector<vector<double > > >  U_bare =
//    vector<vector<vector<double > > >
//    (nuc_eff,vector<vector<double > >
//             ((nuc_eff+1)/2, vector<double >(3)));//three atoms per unit cell
    //the irreducible vertex is approximated by the bare interaction in the parquet approx

//    irreducible() {};
//
//    irreducible(int input[]) {
//        U_bare = Q (input);
//    }

/* TODO: check all member functions
public:
    double vval(int, int, int);
    double vvalsmooth(int, int, int);
    double vvalsmooth(int, int, int,double,double,double,char,int,char);

    void setvert(int,int,int,double);

    friend irreducible operator*(double alpha, const irreducible & vertex);
    friend irreducible  operator*(const irreducible & vertex, double alpha);
    friend irreducible  operator+(const irreducible & vertex1, const irreducible & vertex2);

*/
};
/***************************************************************************************************************************/

/*******************************define "fullvert"FULLVERT" as collection of all diagrams in all channels*********************************/

template <class Q>
class fullvert{//collection of all channels
public:
    irreducible<Q> irred;
    avert<Q> avertex;
    pvert<Q> pvertex;
    tvert<Q> tvertex;

//    fullvert() {
//        irred = irreducible<Q> ();
//        avertex = avert<Q> ();
//        pvertex = pvert<Q> ();
//        tvertex = tvert<Q> ();
//    }
//    fullvert(int input[]) {
//        irred = irreducible<Q> (input);
//        avertex = avert<Q> (input);
//        pvertex = pvert<Q> (input);
//        tvertex = tvert<Q> (input);
//    }

/* TODO: check all member functions
public:

    double vvalsmooth(int,int,int,double,double,double,char);
    double vvalsmooth(int,int,int,double,double,double,char,int,char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map
    friend fullvert operator*(double alpha, const fullvert& vertex);
    friend fullvert operator+(const fullvert& vertex1,const fullvert& vertex2);
*/
};

//template <class T>
//class Keldyshcomp {
//public:
//    T PsiA;
//    T PsiB;
//    T PhiA;
//    T PhiB;
//    T PhiC;
//    T PhiD;
//
//    Keldyshcomp() {
//        PsiA = T ();
//        PsiB = T ();
//        PhiA = T ();
//        PhiB = T ();
//        PhiC = T ();
//        PhiD = T ();
//    }
//    Keldyshcomp(int input[]) {
//        PsiA = T (input);
//        PsiB = T (input);
//        PhiA = T (input);
//        PhiB = T (input);
//        PhiC = T (input);
//        PhiD = T (input);
//    }
//};


//define parvert as tuple of spin and density vertex
template <class T>
class parvert{//define a tuple for parametrized vertices that contains one spin vertex and one density vertex
public:
    T spinvertex;
    T densvertex;

//    Keldyshcomp<T> spinvertex;
//    Keldyshcomp<T> densvertex;
//
//    parvert() {
//        spinvertex = Keldyshcomp<T> () ;
//        densvertex = Keldyshcomp<T> () ;
//    }
//    parvert(int input[]) {
//        spinvertex = Keldyshcomp<T> (input) ;
//        densvertex = Keldyshcomp<T> (input) ;
//    }
};

/* TODO: check all those functions
//define operators for parvert
parvert<avert> operator+(parvert<avert> ,parvert<avert>);
parvert<pvert> operator+(parvert<pvert> ,parvert<pvert>);
parvert<tvert> operator+(parvert<tvert> ,parvert<tvert>);
parvert<irreducible> operator+(parvert<irreducible> ,parvert<irreducible>);
parvert<avert> operator+=(parvert<avert> ,parvert<avert>);
parvert<pvert> operator+=(parvert<pvert> ,parvert<pvert>);
parvert<tvert> operator+=(parvert<tvert> ,parvert<tvert>);
parvert<irreducible> operator+=(parvert<irreducible> ,parvert<irreducible>);
parvert<avert> operator*(double  ,parvert<avert>&);
parvert<avert> operator*(parvert<avert>& ,double );
parvert<pvert> operator*(double  ,parvert<pvert>&);
parvert<pvert> operator*(parvert<pvert>& ,double );
parvert<tvert> operator*(double  ,parvert<tvert>&);
parvert<tvert> operator*(parvert<tvert>& ,double );
parvert<irreducible> operator*(double  ,parvert<irreducible>&);
parvert<irreducible> operator*(parvert<irreducible>& ,double );
*/


#endif //KELDYSH_MFRG_VERTEX_H
