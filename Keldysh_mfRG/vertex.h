//
// Created by Sa.Aguirre on 7/18/19.
//

#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <vector>
#include <array>
#include "parameters.h"
#include "frequency_grid.h"
#include "Keldysh_symmetries.h"

using namespace std;

/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
template <typename Q>
class avert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wa * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wa * nw2_nua * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in);

    /*Lists of the Keldysh components of K1a relating the respective component to the independent ones through the marked
     * trafo*/
    array<int,4>  list_K1_T0_comp1 = {1, 7, 8, 14};     //In the vertex, comp1 will be iK=0
    array<int,4>  list_K1_T1_comp1 = {2, 4, 11, 13};
    array<int,4>  list_K1_T0_comp3 = {3, 5, 10, 12};    //In the vertex, comp2 will be iK=1

    /*Lists of the Keldysh components of K2a relating the respective component to the independent ones through the marked
    * trafo*/
    array<int,2> list_K2_T0_comp1 = {1, 7};
    array<int,2> list_K2_T1_comp1 = {2, 4};
    array<int,2> list_K2_T2_comp1 = {11, 13};
    array<int,2> list_K2_T3_comp1 = {8, 14};
    array<int,2> list_K2_T0_comp3 = {3, 5};
    array<int,2> list_K2_T3_comp3 = {10, 12};


public:
    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
     * of the necessary indices convertions  and what not...
     * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
     * i.e. only complex numbers
     *
     * This function aims to be the sole function one needs to call to read the full vertex*/
    Q value (int, double, double, double, int, char);

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    Q value (int, double, double, double, int);

    /*This function returns the value of the full vertex (i.e. the sum of the diagrammatic classes) for a given
     * combination of Keldysh (first int) and internal structure (second int, set to 0 if no extra structure).*/
    Q vvalsmooth(int, double, double, double, int, char);

    /*Same idea as function above, but is oriented towards the multi-loop implementation */
    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)

    /*No clue, suspect is unnecessary for us since we do not need map or red_side operations*/
//    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map

    /*This function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45*/
    Q vvalsmooth(int, double, double, double, int);


    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);


    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);


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

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    Q K2b_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    tuple<double, double, double> transfToA(double, double, double, char);

//    /*Overload of previous function to single out the transfer from 3-fermionic frequencies*/
//    tuple<double, double, double> transfToA(double, double, double);
//
//    /*This function transforms the frequency arguments from the a-channel convention to the standard 3-fermionic freqs. input
//     * I.e. is the inverse of the function above*/
//    tuple<double, double, double> transfBackA(double, double, double);


    /*The following three functions return a tuple consisting of the new Keldysh index of the overall vertex (given that
     * legs are switched and the three corresponding frequency inputs for the diagrammatic class*/
    tuple<int, double, double, double, int> indices_T1(int, double, double, double, int);

    /*Symmetry which interchanges the outgoing legs*/
    tuple<int, double, double, double, int> indices_T2(int, double, double, double, int);

    /*Symmetry which interchanges both incoming and outgoing legs*/
    tuple<int, double, double, double, int> indices_T3(int, double, double, double, int);

    /*Symmetry which interchanges both incoming with outgoing legs*/
    tuple<int, double, double, double, int> indices_TC(int, double, double, double, int);

    /*Symmetry transformations for the different diagrammatic classes.
     * K1 */
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K1(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K1(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K1(int, double, double, double, int);

    /*K2*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K2(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K2(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K2(int, double, double, double, int);

    /*K3*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K3(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K3(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K3(int, double, double, double, int);

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the a-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    tuple<int, int> indices_sum(int i0, int i2);


    /*Define the operator of multiplying an a-vertex with a number.*/
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

    /*Define the addition operation of two a-vertices*/
    avert<Q> friend operator+(const avert<Q>& vertex1, const avert<Q>& vertex2){
        avert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

};

template <class Q>
class pvert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wp * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);


    /*Lists of the Keldysh components of K1p relating the respective component to the independent ones through the marked
 * trafo*/
    array<int,4> list_K1_T0_comp1 = {1, 2, 13, 14};
    array<int,4> list_K1_TC_comp1 = {4, 7, 8, 11};
    array<int,4> list_K1_T0_comp5 = {5, 6, 9, 10};

    /*Lists of the Keldysh components of K2p relating the respective component to the independent ones through the marked
    * trafo*/
    array<int,2> list_K2_T0_comp1 = {1, 2};
    array<int,2> list_K2_T1_comp1 = {13, 14};
    array<int,2> list_K2_TC_comp1 = {4, 7};
    array<int,2> list_K2_T2TC_comp1 = {8, 11};
    array<int,2> list_K2_T0_comp5 = {5, 6};
    array<int,2> list_K2_T3_comp5 = {9, 10};


public:
    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
 * of the necessary indices convertions  and what not...
 * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
 * i.e. only complex numbers
 *
 * This function aims to be the sole function one needs to call to read the full vertex*/
    Q value (int, double, double, double, int, char);

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    Q value (int, double, double, double, int);

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


    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);


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

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    Q K2b_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
 * only have the values 'a', 'p', or 't'.*/
    tuple<double, double, double> transfToP(double, double, double, char);

//    /*Overload of previous function to single out the transfer from 3-fermionic frequencies*/
//    tuple<double, double, double> transfToP(double, double, double);
//
//    /*This function transforms the frequency arguments from the a-channel convention to the standard 3-fermionic freqs. input
// * I.e. is the inverse of the function above*/
//    tuple<double, double, double> transfBackP(double, double, double);


    /*The following three functions return a tuple consisting of the new Keldysh index of the overall vertex (given that
     * legs are switched and the three corresponding frequency inputs for the diagrammatic class*/
    tuple<int, double, double, double, int> indices_T1(int, double, double, double, int);

    /*Symmetry which interchanges the outgoing legs*/
    tuple<int, double, double, double, int> indices_T2(int, double, double, double, int);

    /*Symmetry which interchanges both incoming and outgoing legs*/
    tuple<int, double, double, double, int> indices_T3(int, double, double, double, int);

    /*Symmetry which interchanges both incoming with outgoing legs*/
    tuple<int, double, double, double, int> indices_TC(int, double, double, double, int);

    /*Symmetry transformations for the different diagrammatic classes.
     * K1 */
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K1(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K1(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K1(int, double, double, double, int);

    /*K2*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K2(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K2(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K2(int, double, double, double, int);

    /*K3*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K3(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K3(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K3(int, double, double, double, int);

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the p-channel*/
    tuple<int, int> indices_sum(int i0, int i2);


    /*Define the operator of multiplying a p-vertex with a number.*/
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

    /*Define the addition operation of two t-vertices*/
    pvert<Q> friend operator+(const pvert<Q>& vertex1, const pvert<Q>& vertex2){
        pvert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

};

template <class Q>
class tvert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wt * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);


    /*Lists of the Keldysh components of K1t relating the respective component to the independent ones through the marked
 * trafo*/
    array<int,4> list_K1_T0_comp1 = {1, 4, 11, 14};
    array<int,4> list_K1_T1_comp1 = {2, 7, 8, 13};
    array<int,4> list_K1_T0_comp3 = {3, 6, 9, 12};

    /*Lists of the Keldysh components of K2t relating the respective component to the independent ones through the marked
    * trafo*/
    array<int,2> list_K2_T0_comp1 = {1, 11};
    array<int,2> list_K2_T1_comp1 = {2, 8};
    array<int,2> list_K2_T2_comp1 = {7, 13};
    array<int,2> list_K2_T3_comp1 = {4, 14};
    array<int,2> list_K2_T0_comp3 = {3, 9};
    array<int,2> list_K2_T3_comp3 = {6, 12};


public:
    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
 * of the necessary indices convertions  and what not...
 * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
 * i.e. only complex numbers
 *
 * This function aims to be the sole function one needs to call to read the full vertex*/
    Q value (int iK, double, double, double, int, char);

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    Q value (int, double, double, double, int);

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


    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);


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

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    Q K2b_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);


    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    tuple<double, double, double> transfToT(double, double, double, char);

//    /*Overload of previous function to single out the transfer from 3-fermionic frequencies*/
//    tuple<double, double, double> transfToT(double, double, double);
//
//    /*This function transforms the frequency arguments from the a-channel convention to the standard 3-fermionic freqs. input
//     * I.e. is the inverse of the function above*/
//    tuple<double, double, double> transfBackT(double, double, double);


    /*The following three functions return a tuple consisting of the new Keldysh index of the overall vertex (given that
     * legs are switched and the three corresponding frequency inputs for the diagrammatic class*/
    tuple<int, double, double, double, int> indices_T1(int, double, double, double, int);

    /*Symmetry which interchanges the outgoing legs*/
    tuple<int, double, double, double, int> indices_T2(int, double, double, double, int);

    /*Symmetry which interchanges both incoming and outgoing legs*/
    tuple<int, double, double, double, int> indices_T3(int, double, double, double, int);

    /*Symmetry which interchanges both incoming with outgoing legs*/
    tuple<int, double, double, double, int> indices_TC(int, double, double, double, int);

    /*Symmetry transformations for the different diagrammatic classes.
     * K1 */
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K1(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K1(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K1(int, double, double, double, int);

    /*K2*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K2(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K2(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K2(int, double, double, double, int);

    /*K3*/
    /*Symmetry which interchanges the incoming legs*/
    Q T1_K3(int, double, double, double, int);
    /*Symmetry which interchanges the outgoing legs*/
    Q T2_K3(int, double, double, double, int);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    Q T3_K3(int, double, double, double, int);

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the lower(0) and upper(1) vertices of a
    * buuble in the t-channel*/
    tuple<int, int> indices_sum(int i0, int i2);


    /*Define the operator of multiplying a p-vertex with a number.*/
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

    /*Define the addition operation of two t-vertices*/
    tvert<Q> friend operator+(const tvert<Q>& vertex1, const tvert<Q>& vertex2){
        tvert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

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

};

/**********************************************************************************************************************/

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
class Vertex{
public:
    T spinvertex;
    T densvertex;
};


/************************ MEMBER FUNCTIONS OF THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX **************************/

/****************************************** MEMBER FUNCTIONS OF THE A-VERTEX ******************************************/
//Here iK is in 0...15 already. Only need to check to what component to transfer to.
template <typename Q> Q avert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
     * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_a=0., v1_a=0., v2_a=0.;
    tie(w_a, v1_a, v2_a) = transfToA(w,v1,v2,channel);
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_a,i_in);

    }
    else if(isInList(iK,list_K1_T1_comp1, list_K1_T1_comp1.size()))
    {
        tie(iK1, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK1 = 0;
        if(fabs(w_a)>w_lower_b)
            valueK1=0.;
        else
            valueK1 = -K1_vvalsmooth(iK1,w_a,i_in);

    }
    else if(isInList(iK, list_K1_T0_comp3, list_K1_T0_comp3.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_a,i_in);
    }
    else
    {
        valueK1=0.;
    }


    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T2_comp1,list_K2_T2_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp1,list_K2_T3_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T0_comp3,list_K2_T0_comp3.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp3,list_K2_T3_comp3.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}

template <typename Q> Q avert<Q>::value(int iK, double w, double v1, double v2, int i_in){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
     * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_a=w, v1_a=v1, v2_a=v2;
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_a,i_in);

    }
    else if(isInList(iK,list_K1_T1_comp1, list_K1_T1_comp1.size()))
    {
        tie(iK1, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK1 = 0;
        if(fabs(w_a)>w_lower_b)
            valueK1=0.;
        else
            valueK1 = -K1_vvalsmooth(iK1,w_a,i_in);

    }
    else if(isInList(iK, list_K1_T0_comp3, list_K1_T0_comp3.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_a,i_in);
    }
    else
    {
        valueK1=0.;
    }


    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T2_comp1,list_K2_T2_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp1,list_K2_T3_comp1.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 0;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T0_comp3,list_K2_T0_comp3.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp3,list_K2_T3_comp3.size())){
        tie(iK2, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK2 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_a, v1_a, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T1(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 1;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T3(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 5;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_TC(iK, w_a, v1_a, v2_a, i_in);
        tie(iK3, w_a, v1_a, v2_a, i_in) = indices_T2(iK, w_a, v1_a, v2_a, i_in);
        iK3 = 7;
        if(fabs(w_a)>w_upper_b || fabs(v1_a)>w_upper_f || fabs(v2_a>w_upper_f))
            valueK3=0.;
        else
            valueK3 = conj(K3_vvalsmooth(iK3, w_a, v1_a, v2_a, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q avert<Q>::vvalsmooth(int iK, double w, double v1, double v2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double w_a=0., v1_a=0., v2_a=0.;
    tie(w_a, v1_a, v2_a) = transfToA(w,v1,v2,channel);

    Q value;

    value += K1_vvalsmooth(iK, w_a, i_in) + K2_vvalsmooth(iK,w_a,v1_a,i_in) + K3_vvalsmooth(iK, w_a, v1_a, v2_a, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,u,w2_u,i_in)

    return value;
}
template <typename Q> Q avert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in,  char channel, int p, char f) {//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class

    double u = 0., w1_u = 0., w2_u = 0.;
    if (channel == 's') {
        u = -w2 - w1;
        w1_u = (w1 - w2 + q) / 2;
        w2_u = (-w1 + w2 + q) / 2;
    } else if (channel ==
               't') {//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2 - w1;
        w1_u = (w1 + w2 - q) / 2;
        w2_u = (w1 + w2 + q) / 2;
    } else if (channel == 'u') {
        u = q;
        w1_u = w1;
        w2_u = w2;
    } else if (channel == 'v') {//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1 - w2;
        w1_u = q + (w1 - w2) / 2;
        w2_u = (w1 + w2) / 2;
    }
    Q value;

//        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    if (p == 1) {
        if (channel == 'u') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
            else if (f == 'K' || f == 'L') {
                value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in);
            };//if outer legs are conntected to same bare vertex
        } else if (channel == 's' || channel == 't') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                         K2_vvalsmooth(iK, u, w1_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}
        }
    } else if (p == 2) {
        if (channel == 'u') {
            if (f == 'R' || f == 'L') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K2_vvalsmooth(iK, u, w2_u, i_in);
            }//if outer legs are conntected to different bare vertex
            else if (f == 'K' || f == 'M') {
                value += K1_vvalsmooth(iK, u,
                                       i_in);; // + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex
            } else if (channel == 's' || channel == 't') {
                if (f == 'R' || f == 'L') {
                    value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                             K2_vvalsmooth(iK, u, w1_u, i_in);  //+ K2b_vvalsmooth(a,b,c,u,w2_u);
                }
            }
        }
        return value;

    }
}

/*overload of previous function         => I'm pretty sure we won't be needing this function couldn't*/
//template <typename Q> Q avert<Q>::vvalsmooth(int red_side, int map, double q, double w1, double w2, char channel, int p, char f){
//    return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
//}

template <typename Q> Q avert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    double u,w1_u,w2_u;
    u = q;
    w1_u = w1;
    w2_u = w2;
    Q value;
//      if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//      if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//      if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};
    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in) + K3_vvalsmooth(iK,u,w1_u,w2_u, i_in);   // +  K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
    return value;
}

template <typename Q> void avert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wa*n_in + i*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in] = value;
}

template <typename Q> void avert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wa*n_in + i*n_in + i_in] += value;
}
template <typename Q> void avert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in] += value;
}
template <typename Q> void avert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in] += value;
}

template <typename Q> Q avert<Q>::K1_vval(int iK, int i, int i_in){
    return K1[iK*nw1_wa*n_in + i*n_in + i_in];
}
template <typename Q> Q avert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in];
}
template <typename Q> Q avert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in];
}

template <typename Q> Q avert<Q>::K1_vvalsmooth(int iK, double w_a, int i_in){

    int index = fconv_K1_a(w_a);

    double x1 = freqs_a[index];
    double x2 = freqs_a[index+1];
    double xd = (w_a-x1)/(x2-x1);

    Q f1 = K1_vval(iK, index, i_in);
    Q f2 = K1_vval(iK, index+1, i_in);

    return (1.-xd)*f1 + xd*f2;
}
template <typename Q> Q avert<Q>::K2_vvalsmooth(int iK, double w_a, double v1_a, int i_in){

    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_a(w_a, v1_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    double xd = (w_a-x1)/(x2-x1);
    double yd = (v1_a-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q avert<Q>::K2b_vvalsmooth(int iK, double w_a, double v1_a, int i_in){
    //TODO implement this method correctly!
    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_a(w_a, v1_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    double xd = (w_a-x1)/(x2-x1);
    double yd = (v1_a-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q avert<Q>::K3_vvalsmooth(int iK, double w_a, double v1_a, double v2_a, int i_in){

    int index_b, index_f, index_fp;
    tie(index_b,index_f, index_fp) = fconv_K3_a(w_a, v1_a, v2_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];
    double z1 = freqs_a[index_fp];
    double z2 = freqs_a[index_fp+1];

    double xd = (w_a-x1)/(x2-x1);
    double yd = (v1_a-y1)/(y2-y1);
    double zd = (v2_a-z1)/(z2-z1);

    Q f111 = K3_vval(iK, index_b, index_f, index_fp, i_in);
    Q f112 = K3_vval(iK, index_b, index_f, index_fp+1, i_in);
    Q f121 = K3_vval(iK, index_b, index_f+1, index_fp, i_in);
    Q f122 = K3_vval(iK, index_b, index_f+1, index_fp+1, i_in);
    Q f211 = K3_vval(iK, index_b+1, index_f, index_fp, i_in);
    Q f212 = K3_vval(iK, index_b+1, index_f, index_fp+1, i_in);
    Q f221 = K3_vval(iK, index_b+1, index_f+1, index_fp, i_in);
    Q f222 = K3_vval(iK, index_b+1, index_f+1, index_fp+1, i_in);

    Q c00 = f111*(1.-xd) + f211*xd;
    Q c01 = f112*(1.-xd) + f212*xd;
    Q c10 = f121*(1.-xd) + f221*xd;
    Q c11 = f122*(1.-xd) + f222*xd;

    Q c0 = c00*(1.-yd) + c10*yd;
    Q c1 = c01*(1.-yd) + c11*yd;

    return c0*(1.-zd) + c1*zd;
}

template<typename Q> tuple<double, double, double> avert<Q>::transfToA(double w, double v1, double v2, char channel)
{
    double w_a=0., v1_a=0., v2_a=0.;
    if(channel == 'a') {
        w_a = w;
        v1_a = v1;
        v2_a = v2;}
    else if(channel == 'p'){
        w_a = -v2-v2;
        v1_a = 0.5*(w+v1-v2);
        v2_a = 0.5*(w-v1+v2);}
    else if(channel == 't'){
        w_a = v1-v2;
        v1_a = 0.5*( w+v1+v2);
        v2_a = 0.5*(-w+v1+v2);}
    return make_tuple(w_a, v1_a, v2_a);
}

template<typename Q> tuple<int, double, double, double, int> avert<Q>::indices_T1(int iK, double w_a, double v1_a, double v2_a, int i_in)
{
    int iKp = T_1_Keldysh(iK);
    double trans_w_a, trans_v1_a, trans_v2_a;

//    tie(ferm1p, ferm2p, ferm1) = transfBackA(w_a,v1_a, v2_a);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    tie(trans_w_a, trans_v1_a, trans_v2_a) = transfToA(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = v2_a - v1_a;
    trans_v1_a = 0.5*(v1_a+v2_a-w_a);
    trans_v2_a = 0.5*(v1_a+v2_a+w_a);

    return make_tuple(iKp, trans_w_a, trans_v1_a, trans_v2_a, i_in);
}
template<typename Q> tuple<int, double, double, double, int> avert<Q>::indices_T2(int iK, double w_a, double v1_a, double v2_a, int i_in)
{
    int iKp = T_2_Keldysh(iK);
    double trans_w_a, trans_v1_a, trans_v2_a;

//    tie(ferm1p, ferm2p, ferm1) = transfBackA(w_a,v1_a, v2_a);
//    /*This is the flipping stage!*/
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_a, trans_v1_a, trans_v2_a) = transfToA(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = v1_a - v2_a;
    trans_v1_a = 0.5*(v1_a+v2_a+w_a);
    trans_v2_a = 0.5*(v1_a+v2_a-w_a);

    return make_tuple(iKp, trans_w_a, trans_v1_a, trans_v2_a, i_in);
}
template<typename Q> tuple<int, double, double, double, int> avert<Q>::indices_T3(int iK, double w_a, double v1_a, double v2_a, int i_in)
{
    int iKp = T_3_Keldysh(iK);
    double trans_w_a, trans_v1_a, trans_v2_a;

//    tie(ferm1p, ferm2p, ferm1) = transfBackA(w_a,v1_a, v2_a);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_a, trans_v1_a, trans_v2_a) = transfToA(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    trans_v1_a = v2_a;
    trans_v2_a = v1_a;

    return make_tuple(iKp, trans_w_a, trans_v1_a, trans_v2_a, i_in);
}
template<typename Q> tuple<int, double, double, double, int> avert<Q>::indices_TC(int iK, double w_a, double v1_a, double v2_a, int i_in)
{
    int iKp = T_C_Keldysh(iK);
    double trans_w_a, trans_v1_a, trans_v2_a;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackA(w_a, v1_a, v2_a);
//    //This is the flipping stage
//    double ferm2 = ferm1p+ferm2p-ferm1;
//    double temp1 = ferm1, temp2 = ferm2;
//    ferm1 = ferm1p;
//    ferm2 = ferm2p;
//    ferm1p = temp1;
//    ferm2p = temp2;
//    tie(trans_w_a, trans_v1_a, trans_v2_a) = transfToA(ferm1p, ferm2p, ferm1);

    trans_w_a = w_a;
    trans_v1_a = v2_a;
    trans_v2_a = v1_a;

    return make_tuple(iKp, trans_w_a, trans_v1_a, trans_v2_a, i_in);
}

template<typename Q> Q avert<Q>::T1_K1(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T1(iK, w_a, v1_a, v2_a, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T2_K1(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T2(iK, w_a, v1_a, v2_a, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T3_K1(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T3(iK, w_a, v1_a, v2_a, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T1_K2(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T1(iK, w_a, v1_a, v2_a, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T2_K2(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T2(iK, w_a, v1_a, v2_a, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T3_K2(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T3(iK, w_a, v1_a, v2_a, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T1_K3(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T1(iK, w_a, v1_a, v2_a, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T2_K3(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T2(iK, w_a, v1_a, v2_a, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q avert<Q>::T3_K3(int iK, double w_a, double v1_a, double v2_a, int i_in){
    auto indices = indices_T3(iK, w_a, v1_a, v2_a, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}

template<typename Q> tuple<int, int> avert<Q>::indices_sum(int i0, int i2)
{
    int a1pi0, a2pi0, a1i0, a2i0, a1pi2, a2pi2, a1i2, a2i2;

    tie(a1pi0, a2pi0, a1i0, a2i0) = alphas(i0);
    tie(a1pi2, a2pi2, a1i2, a2i2) = alphas(i2);

    return make_tuple(
            8*(a1pi0-1) + 4*(a2i2-1) + 2*(a1pi2-1) + 1*(a2i0-1),
            8*(a1i2-1) + 4*(a2pi0-1) + 2*(a1i0-1) + 1*(a2pi2-1));
}


/****************************************** MEMBER FUNCTIONS OF THE P-VERTEX ******************************************/
//Here iK is in 0...15 already. Only need to check to what component to transfer to.
template <typename Q> Q pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_p=0., v1_p=0., v2_p=0.;
    tie(w_p, v1_p, v2_p) = transfToP(w,v1,v2,channel);
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_p,i_in);

    }
    else if(isInList(iK,list_K1_TC_comp1, list_K1_TC_comp1.size()))
    {
        tie(iK1, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK1 = 0;
        valueK1 = conj(K1_vvalsmooth(iK1,w_p,i_in));

    }
    else if(isInList(iK, list_K1_T0_comp5, list_K1_T0_comp5.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_p,i_in);
    }
    else
    {
        valueK1=0.;
    }



    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_TC_comp1,list_K2_TC_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        valueK2 = conj(K2_vvalsmooth(iK2, w_p, v1_p, i_in));
    }
    else if(isInList(iK,list_K2_T2TC_comp1,list_K2_T2TC_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = -conj(K2_vvalsmooth(iK2, w_p, v1_p, i_in));
    }
    else if(isInList(iK,list_K2_T0_comp5,list_K2_T0_comp5.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp5,list_K2_T3_comp5.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T3(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T3(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        valueK3 = conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}

template <typename Q> Q pvert<Q>::value(int iK, double w, double v1, double v2, int i_in){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_p=w, v1_p=v1, v2_p=v2;
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_p,i_in);

    }
    else if(isInList(iK,list_K1_TC_comp1, list_K1_TC_comp1.size()))
    {
        tie(iK1, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK1 = 0;
        valueK1 = conj(K1_vvalsmooth(iK1,w_p,i_in));

    }
    else if(isInList(iK, list_K1_T0_comp5, list_K1_T0_comp5.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_p,i_in);
    }
    else
    {
        valueK1=0.;
    }



    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_TC_comp1,list_K2_TC_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        valueK2 = conj(K2_vvalsmooth(iK2, w_p, v1_p, i_in));
    }
    else if(isInList(iK,list_K2_T2TC_comp1,list_K2_T2TC_comp1.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 0;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = -conj(K2_vvalsmooth(iK2, w_p, v1_p, i_in));
    }
    else if(isInList(iK,list_K2_T0_comp5,list_K2_T0_comp5.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp5,list_K2_T3_comp5.size())){
        tie(iK2, w_p, v1_p, v2_p, i_in) = indices_T3(iK, w_p, v1_p, v2_p, i_in);
        iK2 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f)
            valueK2 = 0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_p, v1_p, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T1(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T3(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 5;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        valueK3 = conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_TC(iK, w_p, v1_p, v2_p, i_in);
        tie(iK3, w_p, v1_p, v2_p, i_in) = indices_T2(iK, w_p, v1_p, v2_p, i_in);
        iK3 = 7;
        if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
            valueK3 = 0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_p, v1_p, v2_p, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q pvert<Q>::vvalsmooth(int iK, double w, double v1, double v2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double w_p=0., v1_p=0., v2_p=0.;
    tie (w_p, v1_p, v2_p) = transfToP(w,v1,v2,channel);

    Q value;

    value += K1_vvalsmooth(iK, w_p, i_in) + K2_vvalsmooth(iK,w_p,v1_p,i_in) + K3_vvalsmooth(iK, w, v1_p, v2_p, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,u,w2_u,i_in)

    return value;
}
template <typename Q> Q pvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in,  char channel, int p, char f) {//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class

    double u = 0., w1_u = 0., w2_u = 0.;
    if (channel == 's') {
        u = -w2 - w1;
        w1_u = (w1 - w2 + q) / 2;
        w2_u = (-w1 + w2 + q) / 2;
    } else if (channel ==
               't') {//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2 - w1;
        w1_u = (w1 + w2 - q) / 2;
        w2_u = (w1 + w2 + q) / 2;
    } else if (channel == 'u') {
        u = q;
        w1_u = w1;
        w2_u = w2;
    } else if (channel == 'v') {//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1 - w2;
        w1_u = q + (w1 - w2) / 2;
        w2_u = (w1 + w2) / 2;
    }
    Q value;

//        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    if (p == 1) {
        if (channel == 'u') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
            else if (f == 'K' || f == 'L') {
                value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in);
            };//if outer legs are conntected to same bare vertex
        } else if (channel == 's' || channel == 't') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                         K2_vvalsmooth(iK, u, w1_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}
        }
    } else if (p == 2) {
        if (channel == 'u') {
            if (f == 'R' || f == 'L') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K2_vvalsmooth(iK, u, w2_u, i_in);
            }//if outer legs are conntected to different bare vertex
            else if (f == 'K' || f == 'M') {
                value += K1_vvalsmooth(iK, u,
                                       i_in);; // + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex
            } else if (channel == 's' || channel == 't') {
                if (f == 'R' || f == 'L') {
                    value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                             K2_vvalsmooth(iK, u, w1_u, i_in);  //+ K2b_vvalsmooth(a,b,c,u,w2_u);
                }
            }
        }
        return value;

    }
}

/*overload of previous function         => I'm pretty sure we won't be needing this function*/
//template <typename Q> Q avert<Q>::vvalsmooth(int red_side, int map, double q, double w1, double w2, char channel, int p, char f){
//    return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
//}

template <typename Q> Q pvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    double u,w1_u,w2_u;
    u = q;
    w1_u = w1;
    w2_u = w2;
    Q value;
//      if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//      if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//      if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};
    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in) + K3_vvalsmooth(iK,u,w1_u,w2_u, i_in);   // +  K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
    return value;
}

template <typename Q> void pvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] += value;
}
template <typename Q> void pvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] += value;
}
template <typename Q> void pvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] += value;
}

template <typename Q> Q pvert<Q>::K1_vval(int iK, int i, int i_in){
    return K1[iK*nw1_wp*n_in + i*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in];
}

template <typename Q> Q pvert<Q>::K1_vvalsmooth(int iK, double w_p, int i_in){

    int index = fconv_K1_p(w_p);

    double x1 = freqs_p[index];
    double x2 = freqs_p[index+1];

    double xd = (w_p-x1)/(x2-x1);

    Q f1 = K1_vval(iK, index, i_in);
    Q f2 = K1_vval(iK, index+1, i_in);

    return (1.-xd)*f1 + xd*f2;
}
template <typename Q> Q pvert<Q>::K2_vvalsmooth(int iK, double w_p, double v1_p, int i_in){

    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_p(w_p, v1_p);

    double x1 = freqs_p[index_b];
    double x2 = freqs_p[index_b+1];
    double y1 = freqs_p[index_f];
    double y2 = freqs_p[index_f+1];

    double xd = (w_p-x1)/(x2-x1);
    double yd = (v1_p-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q pvert<Q>::K2b_vvalsmooth(int iK, double w_p, double v1_p, int i_in)
{
    //TODO implement this method correctly!
    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_a(w_p, v1_p);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    double xd = (w_p-x1)/(x2-x1);
    double yd = (v1_p-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q pvert<Q>::K3_vvalsmooth(int iK, double w_p, double v1_p, double v2_p, int i_in){

    int index_b, index_f, index_fp;
    tie(index_b,index_f, index_fp) = fconv_K3_p(w_p, v1_p, v2_p);

    double x1 = freqs_p[index_b];
    double x2 = freqs_p[index_b+1];
    double y1 = freqs_p[index_f];
    double y2 = freqs_p[index_f+1];
    double z1 = freqs_p[index_fp];
    double z2 = freqs_p[index_fp+1];

    Q f111 = K3_vval(iK, index_b, index_f, index_fp, i_in);
    Q f112 = K3_vval(iK, index_b, index_f, index_fp+1, i_in);
    Q f121 = K3_vval(iK, index_b, index_f+1, index_fp, i_in);
    Q f122 = K3_vval(iK, index_b, index_f+1, index_fp+1, i_in);
    Q f211 = K3_vval(iK, index_b+1, index_f, index_fp, i_in);
    Q f212 = K3_vval(iK, index_b+1, index_f, index_fp+1, i_in);
    Q f221 = K3_vval(iK, index_b+1, index_f+1, index_fp, i_in);
    Q f222 = K3_vval(iK, index_b+1, index_f+1, index_fp+1, i_in);

    double xd = (w_p-x1)/(x2-x1);
    double yd = (v1_p-y1)/(y2-y1);
    double zd = (v2_p-z1)/(z2-z1);

    Q c00 = f111*(1.-xd) + f211*xd;
    Q c01 = f112*(1.-xd) + f212*xd;
    Q c10 = f121*(1.-xd) + f221*xd;
    Q c11 = f122*(1.-xd) + f222*xd;

    Q c0 = c00*(1.-yd) + c10*yd;
    Q c1 = c01*(1.-yd) + c11*yd;

    return c0*(1.-zd) + c1*zd;
}


template<typename Q> tuple<double, double, double> pvert<Q>::transfToP(double w, double v1, double v2, char channel) {
    double w_p=0., v1_p = 0., v2_p=0.;
    if(channel == 'a'){
        w_p = v2+v1;
        v1_p = 0.5*(-w+v1-v2);
        v2_p = 0.5*(-w-v1+v2);}
    else if(channel == 'p'){
        w_p = w;
        v1_p = v1;
        v2_p = v2;}
    else if(channel == 't'){
        w_p = v1+v2;
        v1_p = 0.5*( w-v1+v2);
        v2_p = 0.5*(-w-v1+v2);}
    return make_tuple(w_p, v1_p, v2_p);
}

template<typename Q> tuple<int, double, double, double, int> pvert<Q>::indices_T1(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    int iKp = T_1_Keldysh(iK);
    double trans_w_p, trans_v1_p, trans_v2_p;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackP(w_p,v1_p, v2_p);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    tie(trans_w_p, trans_v1_p, trans_v2_p) = transfToP(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = v1_p;
    trans_v2_p = -v2_p;

    return make_tuple(iKp, trans_w_p, trans_v1_p, trans_v2_p, i_in);
}
template<typename Q> tuple<int, double, double, double, int> pvert<Q>::indices_T2(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    int iKp = T_2_Keldysh(iK);
    double trans_w_p, trans_v1_p, trans_v2_p;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackP(w_p,v1_p, v2_p);
//    /*This is the flipping stage!*/
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_p, trans_v1_p, trans_v2_p) = transfToP(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    trans_v2_p = v2_p;

    return make_tuple(iKp, trans_w_p, trans_v1_p, trans_v2_p, i_in);
}
template<typename Q> tuple<int, double, double, double, int> pvert<Q>::indices_T3(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    int iKp = T_3_Keldysh(iK);
    double trans_w_p, trans_v1_p, trans_v2_p;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackP(w_p,v1_p, v2_p);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_p, trans_v1_p, trans_v2_p) = transfToP(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    trans_v2_p = -v2_p;

    return make_tuple(iKp, trans_w_p, trans_v1_p, trans_v2_p, i_in);
}
template<typename Q> tuple<int, double, double, double, int> pvert<Q>::indices_TC(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    int iKp = T_C_Keldysh(iK);
    double trans_w_p, trans_v1_p, trans_v2_p;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackP(w_p, v1_p, v2_p);
//    //This is the flipping stage
//    double ferm2 = ferm1p+ferm2p-ferm1;
//    double temp1 = ferm1, temp2 = ferm2;
//    ferm1 = ferm1p;
//    ferm2 = ferm2p;
//    ferm1p = temp1;
//    ferm2p = temp2;
//    tie(trans_w_p, trans_v1_p, trans_v2_p) = transfToP(ferm1p, ferm2p, ferm1);

    trans_w_p = w_p;
    trans_v1_p = v2_p;
    trans_v2_p = v1_p;

    return make_tuple(iKp, trans_w_p, trans_v1_p, trans_v2_p, i_in);
}

template<typename Q> Q pvert<Q>::T1_K1(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T1(iK, w_p, v1_p, v2_p, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T2_K1(int iK, double w_p, double v1_p, double v2_p, int i_in){
    auto indices = indices_T2(iK, w_p, v1_p, v2_p, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T3_K1(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T3(iK, w_p, v1_p, v2_p, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T1_K2(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T1(iK, w_p, v1_p, v2_p, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T2_K2(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T2(iK, w_p, v1_p, v2_p, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T3_K2(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T3(iK, w_p, v1_p, v2_p, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T1_K3(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T1(iK, w_p, v1_p, v2_p, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T2_K3(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T2(iK, w_p, v1_p, v2_p, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q pvert<Q>::T3_K3(int iK, double w_p, double v1_p, double v2_p, int i_in)
{
    auto indices = indices_T3(iK, w_p, v1_p, v2_p, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}

template<typename Q> tuple<int, int> pvert<Q>::indices_sum(int i0, int i2)
{
    int a1pi0, a2pi0, a1i0, a2i0, a1pi2, a2pi2, a1i2, a2i2;

    tie(a1pi0, a2pi0, a1i0, a2i0) = alphas(i0);
    tie(a1pi2, a2pi2, a1i2, a2i2) = alphas(i2);

    return make_tuple(
            8*(a1pi0-1) + 4*(a2pi0-1) + 2*(a1pi2-1) + 1*(a2pi2-1),
            8*(a1i2-1) + 4*(a2i2-1) + 2*(a1i0-1) + 1*(a2i0));
}

/****************************************** MEMBER FUNCTIONS OF THE T-VERTEX ******************************************/

template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_t=0., v1_t=0., v2_t=0.;
    tie(w_t, v1_t, v2_t) = transfToT(w,v1,v2,channel);
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_t,i_in);

    }
    else if(isInList(iK,list_K1_T1_comp1, list_K1_T1_comp1.size()))
    {
        tie(iK1, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK1 = 0;
        if(fabs(w_t)>w_upper_b)
            valueK1 = 0.;
        else
            valueK1 = -K1_vvalsmooth(iK1,w_t,i_in);

    }
    else if(isInList(iK, list_K1_T0_comp3, list_K1_T0_comp3.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_t,i_in);
    }
    else
    {
        valueK1=0.;
    }


    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T2_comp1,list_K2_T2_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp1,list_K2_T3_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T0_comp3,list_K2_T0_comp3.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp3,list_K2_T3_comp3.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        valueK3 = conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}

template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/

    double w_t=w, v1_t=v1, v2_t=v2;
    int iK1, iK2, iK3;
    Q valueK1, valueK2, valueK3;

    if(isInList(iK,list_K1_T0_comp1, list_K1_T0_comp1.size()))
    {
        iK1=0;
        valueK1 = K1_vvalsmooth(iK1,w_t,i_in);

    }
    else if(isInList(iK,list_K1_T1_comp1, list_K1_T1_comp1.size()))
    {
        tie(iK1, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK1 = 0;
        if(fabs(w_t)>w_upper_b)
            valueK1 = 0.;
        else
            valueK1 = -K1_vvalsmooth(iK1,w_t,i_in);

    }
    else if(isInList(iK, list_K1_T0_comp3, list_K1_T0_comp3.size()))
    {
        iK1 = 1;
        valueK1 = K1_vvalsmooth(iK1,w_t,i_in);
    }
    else
    {
        valueK1=0.;
    }


    if(isInList(iK,list_K2_T0_comp1,list_K2_T0_comp1.size()))
    {
        iK2 = 0;
        valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T1_comp1,list_K2_T1_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T2_comp1,list_K2_T2_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = -K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp1,list_K2_T3_comp1.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 0;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T0_comp3,list_K2_T0_comp3.size())){
        iK2 = 1;
        valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else if(isInList(iK,list_K2_T3_comp3,list_K2_T3_comp3.size())){
        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK2 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
            valueK2=0.;
        else
            valueK2 = K2_vvalsmooth(iK2, w_t, v1_t, i_in);
    }
    else
    {
        valueK2 = 0.;
    }


    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7)
    {
        iK3=iK;
        valueK3 = K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK == 2)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK==4)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        valueK3 = conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }
    else if(iK == 6)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }
    else if(iK == 8)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }
    else if(iK == 9)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 10)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 11)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
    }

    else if(iK == 12)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    else if(iK == 13)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        valueK3 = conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    else if(iK == 14)
    {
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
        iK3 = 7;
        if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
            valueK3=0.;
        else
            valueK3 = -conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
    }

    return valueK1 + valueK2 + conj(valueK2) + valueK3;
}


//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double w, double v1, double v2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double w_t, v1_t, v2_t;
    tie(w_t, v1_t, v2_t) = transfToT(w, v1,v2,channel);

    Q value;

    value += K1_vvalsmooth(iK, w_t, i_in) + K2_vvalsmooth(iK,w_t,v1_t,i_in) + K3_vvalsmooth(iK, w_t, v1_t, v2_t, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,w_p,w2_u,i_in)

    return value;
}
template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in,  char channel, int p, char f) {//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class

    double u = 0., w1_u = 0., w2_u = 0.;
    if (channel == 's') {
        u = -w2 - w1;
        w1_u = (w1 - w2 + q) / 2;
        w2_u = (-w1 + w2 + q) / 2;
    } else if (channel ==
               't') {//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2 - w1;
        w1_u = (w1 + w2 - q) / 2;
        w2_u = (w1 + w2 + q) / 2;
    } else if (channel == 'u') {
        u = q;
        w1_u = w1;
        w2_u = w2;
    } else if (channel == 'v') {//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1 - w2;
        w1_u = q + (w1 - w2) / 2;
        w2_u = (w1 + w2) / 2;
    }
    Q value;

//        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    if (p == 1) {
        if (channel == 'u') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
            else if (f == 'K' || f == 'L') {
                value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in);
            };//if outer legs are conntected to same bare vertex
        } else if (channel == 's' || channel == 't') {
            if (f == 'R' || f == 'M') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                         K2_vvalsmooth(iK, u, w1_u, i_in);
            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}
        }
    } else if (p == 2) {
        if (channel == 'u') {
            if (f == 'R' || f == 'L') {
                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K2_vvalsmooth(iK, u, w2_u, i_in);
            }//if outer legs are conntected to different bare vertex
            else if (f == 'K' || f == 'M') {
                value += K1_vvalsmooth(iK, u,
                                       i_in);; // + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex
            } else if (channel == 's' || channel == 't') {
                if (f == 'R' || f == 'L') {
                    value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
                             K2_vvalsmooth(iK, u, w1_u, i_in);  //+ K2b_vvalsmooth(a,b,c,u,w2_u);
                }
            }
        }
        return value;

    }
}

/*overload of previous function         => I'm pretty sure we won't be needing this function*/
//template <typename Q> Q avert<Q>::vvalsmooth(int red_side, int map, double q, double w1, double w2, char channel, int p, char f){
//    return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
//}

template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    double u,w1_u,w2_u;
    u = q;
    w1_u = w1;
    w2_u = w2;
    Q value;
//      if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//      if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//      if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};
    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in) + K3_vvalsmooth(iK,u,w1_u,w2_u, i_in);   // +  K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
    return value;
}

template <typename Q> void tvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] += value;
}
template <typename Q> void tvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] += value;
}
template <typename Q> void tvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] += value;
}

template <typename Q> Q tvert<Q>::K1_vval(int iK, int i, int i_in){
    return K1[iK*nw1_wt*n_in + i*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in];
}

template <typename Q> Q tvert<Q>::K1_vvalsmooth(int iK, double w_t, int i_in){

    int index = fconv_K1_t(w_t);

    double x1 = freqs_t[index];
    double x2 = freqs_t[index+1];

    double xd = (w_t-x1)/(x2-x1);

    Q f1 = K1_vval(iK, index, i_in);
    Q f2 = K1_vval(iK, index+1, i_in);

    return (1.-xd)*f1 + xd*f2;
}
template <typename Q> Q tvert<Q>::K2_vvalsmooth(int iK, double w_t, double v1_t, int i_in){

    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_t(w_t, v1_t);

    double x1 = freqs_t[index_b];
    double x2 = freqs_t[index_b+1];
    double y1 = freqs_t[index_f];
    double y2 = freqs_t[index_f+1];

    double xd = (w_t-x1)/(x2-x1);
    double yd = (v1_t-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q tvert<Q>::K2b_vvalsmooth(int iK, double w_t, double v1_t, int i_in){
    //TODO implement this method correctly!
    int index_b, index_f;
    tie(index_b, index_f) = fconv_K2_a(w_t, v1_t);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    double xd = (w_t-x1)/(x2-x1);
    double yd = (v1_t-y1)/(y2-y1);

    Q f11 = K2_vval(iK, index_b, index_f, i_in);
    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);

    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
}
template <typename Q> Q tvert<Q>::K3_vvalsmooth(int iK, double w_t, double v1_t, double v2_t, int i_in){


    int index_b, index_f, index_fp;
    tie(index_b,index_f, index_fp) = fconv_K3_t(w_t, v1_t, v2_t);

    double x1 = freqs_t[index_b];
    double x2 = freqs_t[index_b+1];
    double y1 = freqs_t[index_f];
    double y2 = freqs_t[index_f+1];
    double z1 = freqs_t[index_fp];
    double z2 = freqs_t[index_fp+1];

    Q f111 = K3_vval(iK, index_b, index_f, index_fp, i_in);
    Q f112 = K3_vval(iK, index_b, index_f, index_fp+1, i_in);
    Q f121 = K3_vval(iK, index_b, index_f+1, index_fp, i_in);
    Q f122 = K3_vval(iK, index_b, index_f+1, index_fp+1, i_in);
    Q f211 = K3_vval(iK, index_b+1, index_f, index_fp, i_in);
    Q f212 = K3_vval(iK, index_b+1, index_f, index_fp+1, i_in);
    Q f221 = K3_vval(iK, index_b+1, index_f+1, index_fp, i_in);
    Q f222 = K3_vval(iK, index_b+1, index_f+1, index_fp+1, i_in);

    double xd = (w_t-x1)/(x2-x1);
    double yd = (v1_t-y1)/(y2-y1);
    double zd = (v2_t-z1)/(z2-z1);

    Q c00 = f111*(1.-xd) + f211*xd;
    Q c01 = f112*(1.-xd) + f212*xd;
    Q c10 = f121*(1.-xd) + f221*xd;
    Q c11 = f122*(1.-xd) + f222*xd;

    Q c0 = c00*(1.-yd) + c10*yd;
    Q c1 = c01*(1.-yd) + c11*yd;

    return c0*(1.-zd) + c1*zd;
}

template<typename Q> tuple<double, double, double> tvert<Q>::transfToT(double w, double v1, double v2, char channel) {
    double w_t=0.,v1_t=0., v2_t=0.;
    if(channel == 'a'){
        w_t = v1-v2;
        v1_t = 0.5*( w+v1+v2);
        v2_t = 0.5*(-w+v1+v2);}
    else if(channel == 'p'){
        w_t = v1-v2;
        v1_t = 0.5*(w-v1-v2);
        v2_t = 0.5*(w+v1+v2);}
    else if(channel == 't'){
        w_t = w;
        v1_t = v1;
        v2_t = v2;}
    return make_tuple(w_t, v1_t, v2_t);
}

template<typename Q> tuple<int, double, double, double, int> tvert<Q>::indices_T1(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    int iKp = T_1_Keldysh(iK);
    double trans_w_t, trans_v1_t, trans_v2_t;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t, v1_t, v2_t);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = v2_t - v1_t;
    trans_v1_t = 0.5*(v1_t+v2_t-w_t);
    trans_v2_t = 0.5*(v1_t+v2_t+w_t);

    return make_tuple(iKp, trans_w_t, trans_v1_t, trans_v2_t, i_in);
}
template<typename Q> tuple<int, double, double, double, int> tvert<Q>::indices_T2(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    int iKp = T_2_Keldysh(iK);
    double trans_w_t, trans_v1_t, trans_v2_t;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t,v1_t, v2_t);
//    /*This is the flipping stage!*/
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = v1_t-v2_t;
    trans_v1_t = 0.5*(v1_t+v2_t+w_t);
    trans_v2_t = 0.5*(v1_t+v2_t+w_t);

    return make_tuple(iKp, trans_w_t, trans_v1_t, trans_v2_t, i_in);
}
template<typename Q> tuple<int, double, double, double, int> tvert<Q>::indices_T3(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    int iKp = T_3_Keldysh(iK);
    double trans_w_t, trans_v1_t, trans_v2_t;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t,v1_t, v2_t);
//    /*This is the flipping stage!*/
//    ferm1 = ferm1p+ferm2p-ferm1;
//    double temp = ferm1p;
//    ferm1p = ferm2p;
//    ferm2p = temp;
//    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;

    return make_tuple(iKp, trans_w_t, trans_v1_t, trans_v2_t, i_in);
}
template<typename Q> tuple<int, double, double, double, int> tvert<Q>::indices_TC(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    int iKp = T_C_Keldysh(iK);
    double trans_w_t, trans_v1_t, trans_v2_t;
    double ferm1p, ferm2p, ferm1;

//    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t, v1_t, v2_t);
//    //This is the flipping stage
//    double ferm2 = ferm1p+ferm2p-ferm1;
//    double temp1 = ferm1, temp2 = ferm2;
//    ferm1 = ferm1p;
//    ferm2 = ferm2p;
//    ferm1p = temp1;
//    ferm2p = temp2;
//    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v1_t;
    trans_v2_t = v2_t;

    return make_tuple(iKp, trans_w_t, trans_v1_t, trans_v2_t, i_in);
}


template<typename Q> Q tvert<Q>::T1_K1(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T1(iK, w_t, v1_t, v2_t, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T2_K1(int iK, double w_t, double v1_t, double v2_t, int i_in) {
    auto indices = indices_T2(iK, w_t, v1_t, v2_t, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T3_K1(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T3(iK, w_t, v1_t, v2_t, i_in);
    return K1_vvalsmooth(get<0>(indices), get<1>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T1_K2(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T1(iK, w_t, v1_t, v2_t, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T2_K2(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T2(iK, w_t, v1_t, v2_t, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T3_K2(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T3(iK, w_t, v1_t, v2_t, i_in);
    return K2_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T1_K3(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T1(iK, w_t, v1_t, v2_t, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T2_K3(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T2(iK, w_t, v1_t, v2_t, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}
template<typename Q> Q tvert<Q>::T3_K3(int iK, double w_t, double v1_t, double v2_t, int i_in)
{
    auto indices = indices_T3(iK, w_t, v1_t, v2_t, i_in);
    return K3_vvalsmooth(get<0>(indices), get<1>(indices), get<2>(indices), get<3>(indices), get<4>(indices));
}

template<typename Q> tuple<int, int> tvert<Q>::indices_sum(int i0, int i2)
{
    int a1pi0, a2pi0, a1i0, a2i0, a1pi2, a2pi2, a1i2, a2i2;

    tie(a1pi0, a2pi0, a1i0, a2i0) = alphas(i0);
    tie(a1pi2, a2pi2, a1i2, a2i2) = alphas(i2);

    return make_tuple(
            8*(a1pi0-1) + 4*(a1i2-1) + 2*(a1i0-1) + 1*(a1pi2),
            8*(a2i2-1) + 4*(a2pi0-1) + 2*(a2pi2-1) + 1*(a2i0-1));
}

/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> Q irreducible<Q>::vval() {
    return U_bare;
}
template <typename Q> Q irreducible<Q>::vvalsmooth() {
    return U_bare;
}
template <typename Q> Q irreducible<Q>::vvalsmooth(double q, double w1, double w2, char channel, int par, char f){
    return U_bare;
}
template <typename Q> void irreducible<Q>::setvert( Q value){
    U_bare = value;}


/************************************* MEMBER FUNCTIONS OF THE VERTEX "fullvertex" ************************************/
//arguments are equivalent to those in the simple vertex functions

template <typename Q> Q fullvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel){
    Q result = irred.vvalsmooth() + pvertex.vvalsmooth(iK,q,w1,w2,i_in,channel) + tvertex.vvalsmooth(iK,q,w1,w2,i_in,channel) + avertex.vvalsmooth(iK,q,w1,w2,i_in,channel);
    if(abs(result)>1e-16){
        return  result;}
    else{return 0.;}
}
template <typename Q> Q fullvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel, int p, char f){
    Q result=0;

    if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
        result += irred.vvalsmooth();

    }
    else if( p==2 && (f=='K' || f== 'M')){
        result += irred.vvalsmooth();

    }

    result +=  pvertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f) + tvertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f) + avertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f);
    //if(p==2 && f=='L'){cout <<  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << " " << tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f)  << " " <<  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << endl;};
    if(abs(result)<1e-16){
        result  =0;}
    return result;
}
template <typename Q> Q fullvert<Q>::vvalsmooth(int red_side, int map, int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){// red_side: if only complementary channel in one vertex, which one is reduced? (0/1/2), p: is this the left/upper (1) or the right/lower (2) vertex of the bubble?, f: diagrammatic class that is computed
    Q result=0;

    if(red_side != p){
        if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
            result = irred.vvalsmooth(a,b,c);
        }
        else if( p==2 && (f=='K' || f== 'M')){
            result = irred.vvalsmooth(a,b,c);
        };
        result +=  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
    }


    else if(red_side == p){

        if(map==0){
            if(channel == 's'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 'u'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}}
        else if(map==1){
            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
            else if(channel == 'u'){  result =  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
        }

    };

    if(abs(result)<1e-16){
        result  =0;};
    return result;
}

//template <typename Q> fullvert<Q> operator*(double alpha, const fullvert<Q>& vertex){
//    fullvert<Q> result;
//    result.irred = alpha * vertex.irred;
//    result.pvertex = alpha *vertex.pvertex;
//    result.tvertex = alpha * vertex.tvertex;
//    result.avertex = alpha * vertex.avertex;
//    return result;
//}
//template <typename Q> fullvert<Q> operator+( const fullvert<Q>& vertex1, const fullvert<Q>& vertex2){
//    fullvert<Q> result;
//    result.irred = vertex1.irred + vertex2.irred ;
//    result.pvertex = vertex1.pvertex + vertex2.pvertex ;
//    result.tvertex = vertex1.tvertex + vertex2.tvertex ;
//    result.avertex = vertex1.avertex + vertex2.avertex ;
//    return result;
//}

/**************************************** operators concerning Vertex objects *****************************************/

template <typename Q> Vertex<avert<Q> > operator+(Vertex<avert<Q> > vertex1, Vertex<avert<Q> > vertex2){
    Vertex<avert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<pvert<Q> > operator+(Vertex<pvert<Q> > vertex1, Vertex<pvert<Q> > vertex2){
    Vertex<pvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<tvert<Q> > operator+(Vertex<tvert<Q> > vertex1, Vertex<tvert<Q> > vertex2){
    Vertex<tvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<irreducible<Q> > operator+(Vertex<irreducible<Q> > vertex1,Vertex<irreducible<Q> > vertex2){
    Vertex<irreducible<Q> > result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<avert<Q> > operator+=(Vertex<avert<Q> > vertex1, Vertex<avert<Q> > vertex2){
    Vertex<avert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<pvert<Q> > operator+=(Vertex<pvert<Q> > vertex1, Vertex<pvert<Q> > vertex2){
    Vertex<pvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<tvert<Q> > operator+=(Vertex<tvert<Q> > vertex1, Vertex<tvert<Q> > vertex2){
    Vertex<tvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<irreducible<Q> > operator+=(Vertex<irreducible<Q> > vertex1,Vertex<irreducible<Q> > vertex2){
    Vertex<irreducible<Q> > result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> Vertex<avert<Q> > operator*(double alpha, Vertex<avert<Q> > &vertex){
    Vertex<avert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<avert<Q> > operator*(Vertex<avert<Q> > &vertex,double alpha){
    Vertex<avert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<pvert<Q> > operator*(double alpha, Vertex<pvert<Q> > &vertex){
    Vertex<pvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<pvert<Q> > operator*(Vertex<pvert<Q> > &vertex,double alpha){
    Vertex<pvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<tvert<Q> > operator*(double alpha, Vertex<tvert<Q> > &vertex){
    Vertex<tvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<tvert<Q> > operator*(Vertex<tvert<Q> > &vertex,double alpha){
    Vertex<tvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<irreducible<Q> > operator*(double alpha, Vertex<irreducible<Q> > &vertex){
    Vertex<irreducible<Q> > result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> Vertex<irreducible<Q> > operator*(Vertex<irreducible<Q> > &vertex, double alpha){
    Vertex<irreducible<Q> > result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}

/**********************************************************************************************************************/


#endif //KELDYSH_MFRG_VERTEX_H
