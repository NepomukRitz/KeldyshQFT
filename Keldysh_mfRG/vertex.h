//
// Created by Sa.Aguirre on 7/18/19.
//

#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <vector>

#include "parameters.h"
#include "data_structures.h"

using namespace std;

/*******************************CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX*********************************/
template <typename Q>
class avert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wa * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wa * nw2_nua * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in);

public:
    /*This function returns the value of the full vertex (i.e. the sum of the diagrammatic classes) for a given
     * combination of Keldysh (first int) and internal structure (second int, set to 0 if no extra structure).*/
    Q vvalsmooth(int, double, double, double, int, char);

    /*Same idea as function above, but is oriented towards the multi-loop implementation */
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)

    /*No clue, suspect is unnecessary fro us since we do not need map or red_side operations*/
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map

    /*This function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45*/
    Q vvalsmooth(int, double, double, double, int);


    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);


    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    Q K1_vval(int, int, int);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    Q K2_vval(int, int, int, int);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    Q K3_vval(int, int, int, int, int);


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    Q K1_vvalsmooth(int, double, int);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    Q K2_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);

//    /*Define the operator of multiplying an a-vertex with a number.    TODO generalize it to multiply by type Q ?*/
//    friend avert operator*(double alpha, const avert& vertex);
//    friend avert operator*(const avert& vertex, double alpha);
//
//    /*Define the addition operation of two a-vertices*/
//    friend avert operator+(const avert& vertex1, const avert& vertex2);

    avert<Q> friend operator*(double alpha, const avert<Q> vertex) {
        avert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    avert<Q> friend operator*(const avert<Q>& vertex, double alpha){
        avert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    avert<Q> friend operator+(const avert<Q>& vertex1, const avert<Q>& vertex2){
        avert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

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
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wp * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);

public:
    /*This function returns the value of the full vertex (i.e. the sum of the diagrammatic classes) for a given
     * combination of Keldysh (first int) and internal structure (second int, set to 0 if no extra structure).*/
    Q vvalsmooth(int, double, double, double, int, char);

    /*Same idea as function above, but is oriented towards the multi-loop implementation */
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)

    /*No clue, suspect is unnecessary fro us since we do not need map or red_side operations*/
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map

    /*This function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45*/
    Q vvalsmooth(int, double, double, double, int);


    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);


    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    Q K1_vval(int, int, int);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    Q K2_vval(int, int, int, int);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    Q K3_vval(int, int, int, int, int);


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    Q K1_vvalsmooth(int, double, int);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    Q K2_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);


//    /*Define the operator of multiplying a p-vertex with a number.    TODO generalize it to multiply by type Q ?*/
//    friend pvert operator*(double alpha, const pvert& vertex);
//    friend pvert operator*(const pvert& vertex, double alpha);
//
//    /*Define the addition operation of two p-vertices*/
//    friend pvert operator+(const pvert& vertex1, const pvert& vertex2);

    pvert<Q> friend operator*(double alpha, const pvert<Q> vertex) {
        pvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    pvert<Q> friend operator*(const pvert<Q>& vertex, double alpha){
        pvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    pvert<Q> friend operator+(const pvert<Q>& vertex1, const pvert<Q>& vertex2){
        pvert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

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
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wt * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);

public:
/*This function returns the value of the full vertex (i.e. the sum of the diagrammatic classes) for a given
     * combination of Keldysh (first int) and internal structure (second int, set to 0 if no extra structure).*/
    Q vvalsmooth(int, double, double, double, int, char);

    /*Same idea as function above, but is oriented towards the multi-loop implementation */
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)

    /*No clue, suspect is unnecessary fro us since we do not need map or red_side operations*/
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map

    /*This function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45*/
    Q vvalsmooth(int, double, double, double, int);


    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);


    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    Q K1_vval(int, int, int);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    Q K2_vval(int, int, int, int);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    Q K3_vval(int, int, int, int, int);


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    Q K1_vvalsmooth(int, double, int);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    Q K2_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);


//    /*Define the operator of multiplying a p-vertex with a number.    TODO generalize it to multiply by type Q ?*/
//    friend tvert operator*(double alpha, const tvert& vertex);
//    friend tvert operator*(const tvert& vertex, double alpha);
//
//    /*Define the addition operation of two t-vertices*/
//    friend tvert operator+(const tvert& vertex1, const tvert& vertex2);

    tvert<Q> friend operator*(double alpha, const tvert<Q> vertex) {
        tvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    tvert<Q> friend operator*(const tvert<Q>& vertex, double alpha){
        tvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    tvert<Q> friend operator+(const tvert<Q>& vertex1, const tvert<Q>& vertex2){
        tvert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }


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

public:

    /*All three functions return the value of the bare vertex. Since this value is, this far, independent of everything,
     * the third function stays the same. However, should count on having to adapt it if an internal structure should
     * come forth where the bare interaction does not remain invariant throughout the system.*/
    Q vval();
    Q vvalsmooth();
    Q vvalsmooth(double,double,double,char,int,char);

    /*Sets the value of the bare interaction to Q TODO are we always going to be working in the PA?*/
    void setvert(Q);

    /*Multiplies the value of the irreducible vertex by alpha*/
    irreducible<Q> friend operator*(double alpha, const irreducible<Q>& vertex) {
        irreducible<Q> result;
        result.U_bare = alpha * vertex.U_bare;
        return result;

    }
    irreducible<Q> friend operator*(const irreducible<Q>& vertex,double alpha) {
        irreducible<Q> result;
        result.U_bare = alpha * vertex.U_bare;
        return result;
    }

    /*Defines the addition operation for two irreducible vertices*/
    irreducible<Q> friend operator+(const irreducible<Q>& vertex1,const irreducible<Q>& vertex2) {
        irreducible<Q> result;
        result.U_bare = vertex1.U_bare + vertex2.U_bare;
        return result;
    }

    /*
    Q vval(int, int, int);
    Q vvalsmooth(int, int, int);
    Q vvalsmooth(int, int, int,double,double,double,char,int,char);

    void setvert(int,int,int,double);
    */
};
/***************************************************************************************************************************/

template <class Q>
class fullvert{
public:
    irreducible<Q> irred;
    avert<Q> avertex;
    pvert<Q> pvertex;
    tvert<Q> tvertex;


    /*Returns the value of the full vertex (i.e. irreducible + diagrammatic classes) for the given channel (char),
     * Keldysh index (1st int), internal structure index (2nd int) and the three frequencies.*/
    Q vvalsmooth(int,double,double,double,int,char);

    /*Same as above but with more toys*/
    Q vvalsmooth(int,double,double,double,int,char,int,char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)

    /*We might not need map, but the red_side idea is very nice for multi-loop implementation */
    Q vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map

    /*Defines multiplication of a full vertex by a number*/
    fullvert<Q> friend operator*(double alpha, const fullvert<Q>& vertex){
        fullvert<Q> result;
        result.irred = alpha * vertex.irred;
        result.pvertex = alpha *vertex.pvertex;
        result.tvertex = alpha * vertex.tvertex;
        result.avertex = alpha * vertex.avertex;
        return result;
    }

    /*Defines addition operation of two full vertices*/
    fullvert<Q> friend operator+( const fullvert<Q>& vertex1, const fullvert<Q>& vertex2) {
        fullvert<Q> result;
        result.irred = vertex1.irred + vertex2.irred;
        result.pvertex = vertex1.pvertex + vertex2.pvertex;
        result.tvertex = vertex1.tvertex + vertex2.tvertex;
        result.avertex = vertex1.avertex + vertex2.avertex;
        return result;
    }

    /*
    double vvalsmooth(int,int,int,double,double,double,char);
    double vvalsmooth(int,int,int,double,double,double,char,int,char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
    double vvalsmooth(int,int,int,int,int,double,double,double,char,int,char);//first two arguments: red_side, map
    */
    };




//define Vertex as tuple of spin and density vertex
//Note that T should be fullvert, so that every (i.e. both) spin component of the vertex inherits a Keldysh substructure
template <class T>
class Vertex{//define a tuple for parametrized vertices that contains one spin vertex and one density vertex
public:
    T spinvertex;
    T densvertex;
};

//define operators for Vertex
template <typename Q> Vertex<avert<Q> > operator+(Vertex<avert<Q> > ,Vertex<avert<Q> >);
template <typename Q> Vertex<pvert<Q> > operator+(Vertex<pvert<Q> > ,Vertex<pvert<Q> >);
template <typename Q> Vertex<tvert<Q> > operator+(Vertex<tvert<Q> > ,Vertex<tvert<Q> >);
template <typename Q> Vertex<irreducible<Q> > operator+(Vertex<irreducible<Q> > ,Vertex<irreducible<Q> >);
template <typename Q> Vertex<avert<Q> > operator+=(Vertex<avert<Q> > ,Vertex<avert<Q> >);
template <typename Q> Vertex<pvert<Q> > operator+=(Vertex<pvert<Q> > ,Vertex<pvert<Q> >);
template <typename Q> Vertex<tvert<Q> > operator+=(Vertex<tvert<Q> > ,Vertex<tvert<Q> >);
template <typename Q> Vertex<irreducible<Q> > operator+=(Vertex<irreducible<Q> > ,Vertex<irreducible<Q> >);
template <typename Q> Vertex<avert<Q> > operator*(double  ,Vertex<avert<Q> >&);
template <typename Q> Vertex<avert<Q> > operator*(Vertex<avert<Q> >& ,double );
template <typename Q> Vertex<pvert<Q> > operator*(double  ,Vertex<pvert<Q> >&);
template <typename Q> Vertex<pvert<Q> > operator*(Vertex<pvert<Q> >& ,double );
template <typename Q> Vertex<tvert<Q> > operator*(double  ,Vertex<tvert<Q> >&);
template <typename Q> Vertex<tvert<Q> > operator*(Vertex<tvert<Q> >& ,double );
template <typename Q> Vertex<irreducible<Q> > operator*(double  ,Vertex<irreducible<Q> >&);
template <typename Q> Vertex<irreducible<Q> > operator*(Vertex<irreducible<Q> >& ,double );



#endif //KELDYSH_MFRG_VERTEX_H
