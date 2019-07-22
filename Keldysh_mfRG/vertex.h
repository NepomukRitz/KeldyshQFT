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

    //K1:
    vec<Q> vec_K1 = vec<Q> (nK_K1 * nw1_wa * n_in);

    //K2:
    vec<Q> vec_K2 = vec<Q> (nK_K2 * nw2_wa * nw2_nua * n_in);

    //K3:
    vec<Q> vec_K3 = vec<Q> (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in);
public:
    Q K1(int i_K, int i_w1, int i_in) {
        return vec_K1[i_K*nw1_wa*n_in + i_w1*n_in + i_in];
    }
    Q K2(int i_K, int i_w1, int i_w2, int i_in) {
        return vec_K2[i_K*nw2_wa*nw2_nua*n_in + i_w1*nw2_nua*n_in + i_w2*n_in + i_in];
    }
    Q K3(int i_K, int i_w1, int i_w2, int i_w3, int i_in) {
        return vec_K3[i_K*nw3_wa*nw3_nua*nw3_nuap*n_in + i_w1*nw3_nua*nw3_nuap*n_in + i_w2*nw3_nuap*n_in + i_w3*n_in + i_in];
    }


// TODO: check all member functions
public:
    Q vvalsmooth(int, double, double, double, int, char);
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, double, double, double, int);
    void K1_setvert(int, int, int, Q);
    void K2_setvert(int, int, int, int, Q);
    void K3_setvert(int, int, int, int, int, Q);
    Q K1_vval(int, int, int);
    Q K2_vval(int, int, int, int);
    Q K3_vval(int, int, int, int, int);
    Q K1_vvalsmooth(int, double, int);
    Q K2_vvalsmooth(int, double, double, int);
    Q K3_vvalsmooth(int, double, double, double, int);

    friend avert operator*(double alpha, const avert& vertex);
    friend avert operator*(const avert& vertex, double alpha);
    friend avert operator+(const avert& vertex1, const avert& vertex2);

    /*
    Q vvalsmooth(int, int, int, double, double, double, char);
    Q vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    Q vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, Q);
    void K2_setvert( int, int, int, int, int, Q);
    void K3_setvert( int, int, int, int, int, int, Q);
    Q K1_vval(int, int, int, int);
    Q K2_vval(int, int, int, int, int);
    Q K3_vval(int, int, int, int, int, int);
    Q K1_vvalsmooth(int, int, int, double);
    Q K2_vvalsmooth(int, int, int, double, double);
    Q K3_vvalsmooth(int, int, int, double, double, double);
    */
};

template <class Q>
class pvert{

    //K1:
    vec<Q> vec_K1 = vec<Q> (nK_K1 * nw1_wp * n_in);

    //K2:
    vec<Q> vec_K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);

    //K3:
    vec<Q> vec_K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);

public:
    Q K1(int i_K, int i_w1, int i_in) {
        return vec_K1[i_K*nw1_wp*n_in + i_w1*n_in + i_in];
    }
    Q K2(int i_K, int i_w1, int i_w2, int i_in) {
        return vec_K2[i_K*nw2_wp*nw2_nup*n_in + i_w1*nw2_nup*n_in + i_w2*n_in + i_in];
    }
    Q K3(int i_K, int i_w1, int i_w2, int i_w3, int i_in) {
        return vec_K3[i_K*nw3_wp*nw3_nup*nw3_nupp*n_in + i_w1*nw3_nup*nw3_nupp*n_in + i_w2*nw3_nupp*n_in + i_w3*n_in + i_in];
    }


// TODO: check all member functions
public:
    Q vvalsmooth(int, double, double, double, int, char);
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, double, double, double, int);
    void K1_setvert(int, int, int, Q);
    void K2_setvert(int, int, int, int, Q);
    void K3_setvert(int, int, int, int, int, Q);
    Q K1_vval(int, int, int);
    Q K2_vval(int, int, int, int);
    Q K3_vval(int, int, int, int, int);
    Q K1_vvalsmooth(int, double, int);
    Q K2_vvalsmooth(int, double, double, int);
    Q K3_vvalsmooth(int, double, double, double, int);

    friend pvert operator*(double alpha, const pvert& vertex);
    friend pvert operator*(const pvert& vertex, double alpha);
    friend pvert operator+(const pvert& vertex1, const pvert& vertex2);
    /*
    Q vvalsmooth(int, int, int, double, double, double, char);
    Q vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    Q vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, Q);
    void K2_setvert( int, int, int, int, int, Q);
    void K3_setvert( int, int, int, int, int, int, Q);
    Q K1_vval(int, int, int, int);
    Q K2_vval(int, int, int, int, int);
    Q K3_vval(int, int, int, int, int, int);
    Q K1_vvalsmooth(int, int, int, double);
    Q K2_vvalsmooth(int, int, int, double, double);
    Q K3_vvalsmooth(int, int, int, double, double, double);
    */
};

template <class Q>
class tvert{

    //K1:
    vec<Q> vec_K1 = vec<Q> (nK_K1 * nw1_wt * n_in);

    //K2:
    vec<Q> vec_K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);

    //K3
    vec<Q> vec_K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);

public:
    Q K1(int i_K, int i_w1, int i_in) {
        return vec_K1[i_K*nw1_wt*n_in + i_w1*n_in + i_in];
    }
    Q K2(int i_K, int i_w1, int i_w2, int i_in) {
        return vec_K2[i_K*nw2_wt*nw2_nut*n_in + i_w1*nw2_nut*n_in + i_w2*n_in + i_in];
    }
    Q K3(int i_K, int i_w1, int i_w2, int i_w3, int i_in) {
        return vec_K3[i_K*nw3_wt*nw3_nut*nw3_nutp*n_in + i_w1*nw3_nut*nw3_nutp*n_in + i_w2*nw3_nutp*n_in + i_w3*n_in + i_in];
    }


// TODO: check all member functions
public:
    Q vvalsmooth(int, double, double, double, int, char);
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, double, double, double, int);
    void K1_setvert(int, int, int, Q);
    void K2_setvert(int, int, int, int, Q);
    void K3_setvert(int, int, int, int, int, Q);
    Q K1_vval(int, int, int);
    Q K2_vval(int, int, int, int);
    Q K3_vval(int, int, int, int, int);
    Q K1_vvalsmooth(int, double, int);
    Q K2_vvalsmooth(int, double, double, int);
    Q K3_vvalsmooth(int, double, double, double, int);

    friend tvert operator*(double alpha, const tvert& vertex);
    friend tvert operator*(const tvert& vertex, double alpha);
    friend tvert operator+(const tvert& vertex1, const tvert& vertex2);
    /*
    Q vvalsmooth(int, int, int, double, double, double, char);
    Q vvalsmooth(int, int, int, double, double, double, char, int, char);//second to last arguement: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    Q vvalsmooth(int, int, int, int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
    Q vvalsmooth(int, int, int, double, double, double);
    void K1_setvert( int, int, int, int, Q);
    void K2_setvert( int, int, int, int, int, Q);
    void K3_setvert( int, int, int, int, int, int, Q);
    Q K1_vval(int, int, int, int);
    Q K2_vval(int, int, int, int, int);
    Q K3_vval(int, int, int, int, int, int);
    Q K1_vvalsmooth(int, int, int, double);
    Q K2_vvalsmooth(int, int, int, double, double);
    Q K3_vvalsmooth(int, int, int, double, double, double);
    */
};


template <class Q>
class irreducible{
public:
    Q U_bare;

// TODO: check all member functions
public:
    Q vval();
    Q vvalsmooth();
    Q vvalsmooth(double,double,double,char,int,char);

    void setvert(Q);

    friend irreducible operator*(double alpha, const irreducible & vertex);
    friend irreducible  operator*(const irreducible & vertex, double alpha);
    friend irreducible  operator+(const irreducible & vertex1, const irreducible & vertex2);
    /*
    Q vval(int, int, int);
    Q vvalsmooth(int, int, int);
    Q vvalsmooth(int, int, int,double,double,double,char,int,char);

    void setvert(int,int,int,double);
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


// TODO: check all member functions
public:
    Q vvalsmooth(int,int,int,double,double,double,char);
    Q vvalsmooth(int,int,int,double,double,double,char,int,char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    Q vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map
    friend fullvert operator*(double alpha, const fullvert& vertex);
    friend fullvert operator+(const fullvert& vertex1,const fullvert& vertex2);
    /*
    double vvalsmooth(int,int,int,double,double,double,char);
    double vvalsmooth(int,int,int,double,double,double,char,int,char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map
    */
};




//define parvert as tuple of spin and density vertex
template <class T>
class parvert{//define a tuple for parametrized vertices that contains one spin vertex and one density vertex
public:
    T spinvertex;
    T densvertex;


};

// TODO: check all those functions
//define operators for parvert
template <typename Q> parvert<avert<Q> > operator+(parvert<avert<Q> > ,parvert<avert<Q> >);
template <typename Q> parvert<pvert<Q> > operator+(parvert<pvert<Q> > ,parvert<pvert<Q> >);
template <typename Q> parvert<tvert<Q> > operator+(parvert<tvert<Q> > ,parvert<tvert<Q> >);
template <typename Q> parvert<irreducible<Q> > operator+(parvert<irreducible<Q> > ,parvert<irreducible<Q> >);
template <typename Q> parvert<avert<Q> > operator+=(parvert<avert<Q> > ,parvert<avert<Q> >);
template <typename Q> parvert<pvert<Q> > operator+=(parvert<pvert<Q> > ,parvert<pvert<Q> >);
template <typename Q> parvert<tvert<Q> > operator+=(parvert<tvert<Q> > ,parvert<tvert<Q> >);
template <typename Q> parvert<irreducible<Q> > operator+=(parvert<irreducible<Q> > ,parvert<irreducible<Q> >);
template <typename Q> parvert<avert<Q> > operator*(double  ,parvert<avert<Q> >&);
template <typename Q> parvert<avert<Q> > operator*(parvert<avert<Q> >& ,double );
template <typename Q> parvert<pvert<Q> > operator*(double  ,parvert<pvert<Q> >&);
template <typename Q> parvert<pvert<Q> > operator*(parvert<pvert<Q> >& ,double );
template <typename Q> parvert<tvert<Q> > operator*(double  ,parvert<tvert<Q> >&);
template <typename Q> parvert<tvert<Q> > operator*(parvert<tvert<Q> >& ,double );
template <typename Q> parvert<irreducible<Q> > operator*(double  ,parvert<irreducible<Q> >&);
template <typename Q> parvert<irreducible<Q> > operator*(parvert<irreducible<Q> >& ,double );



#endif //KELDYSH_MFRG_VERTEX_H
