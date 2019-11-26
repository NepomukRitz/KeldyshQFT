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
#include "a_vertex.h"
#include "p_vertex.h"
#include "t_vertex.h"

using namespace std;

/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
//Irreducible
//The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
template <class Q>
class irreducible{
public:
    vec<Q> bare = vec<Q>(16);

    irreducible() = default;;

    /*All three functions return the value of the bare vertex. Since this value is, this far, independent of everything,
     * the third function stays the same. However, should count on having to adapt it if an internal structure should
     * come forth where the bare interaction does not remain invariant throughout the system.*/
    Q vval(int iK);

    /*Sets the value of the bare interaction to Q*/
    void setvert(int iK, Q);

    /*Multiplies the value of the irreducible vertex by alpha anf the addition operation for two irreducible vertices*/
    irreducible<Q> operator+ (const irreducible<Q>& vertex)
    {
        this->bare + vertex.bare;
        return *this;
    }
    irreducible<Q> operator+= (const irreducible<Q>& vertex)
    {
        this->bare +=vertex.bare;
        return *this;
    }
    irreducible<Q> operator-= (const irreducible<Q>& vertex)
    {
        this->bare -=vertex.bare;
        return *this;
    }
    irreducible<Q> operator* (double alpha)
    {
        this->bare *alpha;
        return *this;
    }

    irreducible<Q> operator*= (double alpha)
    {
        this->bare *=alpha;
        return *this;
    }
};

/**********************************************************************************************************************/
//The class fullvert
//The class defining a vertex with full channel decomposition i.e. irreducible (bare) a, p and t channels
template <class Q>
class fullvert {
public:
    irreducible<Q> irred;
    avert<Q> avertex;
    pvert<Q> pvertex;
    tvert<Q> tvertex;

    fullvert() = default;;


    /*Returns the value of the full vertex (i.e. irreducible + diagrammatic classes) for the given channel (char),
     * Keldysh index (1st int), internal structure index (2nd int) and the three frequencies.*/
    Q value(int, double, double, double, int, char);

    /* Returns the sum of the contributions of the diagrammatic classes r' =/= r */
    Q gammaRb(int, double, double, double, int, char);


    fullvert<Q> operator+ (const fullvert<Q>& vertex1) {
        this->irred   + vertex1.irred;
        this->pvertex + vertex1.pvertex;
        this->tvertex + vertex1.tvertex;
        this->avertex + vertex1.avertex;
        return *this;
    }
    fullvert<Q> operator+= (const fullvert<Q>& vertex1) {
        this->irred   += vertex1.irred;
        this->pvertex += vertex1.pvertex;
        this->tvertex += vertex1.tvertex;
        this->avertex += vertex1.avertex;
        return *this;
    }
    fullvert<Q> operator* (double alpha){
        this->irred   *alpha;
        this->pvertex *alpha;
        this->tvertex *alpha;
        this->avertex *alpha;
        return *this;
    }

    fullvert<Q> operator*= (double alpha){
        this->irred   *=alpha;
        this->pvertex *=alpha;
        this->tvertex *=alpha;
        this->avertex *=alpha;
        return *this;
    }
    fullvert<Q> operator-= (const fullvert<Q>& vertex1) {
        this->irred   -= vertex1.irred;
        this->pvertex -= vertex1.pvertex;
        this->tvertex -= vertex1.tvertex;
        this->avertex -= vertex1.avertex;
        return *this;
    }
};


//define Vertex as tuple of spin and density vertex
//Note that T should be fullvert, so that every (i.e. both) spin component of the vertex inherits a Keldysh substructure
template <class T>
class Vertex{
public:
    T spinvertex;
    T densvertex;

    Vertex() = default;;

    Vertex<T> operator+ (const Vertex<T>& vertex1)
    {
        this->densvertex + vertex1.densvertex;
//        this->spinvertex + vertex1.spinvertex;
        return *this;
    }
    Vertex<T> operator+= (const Vertex<T>& vertex1)
    {
        this->densvertex += vertex1.densvertex;
//        this->spinvertex += vertex1.spinvertex;
        return *this;
    }
    Vertex<T> operator* (double alpha)
    {
        this->densvertex*alpha;
//        this->spinvertex*alpha;
        return *this;
    }
    Vertex<T> operator*= (double alpha)
    {
        this->densvertex*=alpha;
//        this->spinvertex*=alpha;
        return *this;
    }
    Vertex<T> operator-= (const Vertex<T>& vertex1)
    {
        this->densvertex -= vertex1.densvertex;
//        this->spinvertex -= vertex1.spinvertex;
        return *this;
    }

///************************************* MEMBER FUNCTIONS OF THE VERTEX "Vertex" ****************************************/
//
//    template <typename Q> void set_densvertex(Vertex<avert<Q> >& avertex, Vertex<pvert<Q> >& pvertex, Vertex<tvert<Q> >& tvertex) {
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wa * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
//            int iwa = (i % (nw1_wa * n_in)) / n_in;
//            int i_in = i % n_in;
//            densvertex.avertex.K1_setvert(iK, iwa, i_in, avertex.densvertex.K1_vval(iK, iwa, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wp * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
//            int iwp = (i % (nw1_wp * n_in)) / n_in;
//            int i_in = i % n_in;
//            densvertex.pvertex.K1_setvert(iK, iwp, i_in, pvertex.densvertex.K1_vval(iK, iwp, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wt * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
//            int iwt = (i % (nw1_wt * n_in)) / n_in;
//            int i_in = i % n_in;
//            densvertex.tvertex.K1_setvert(iK, iwt, i_in, tvertex.densvertex.K1_vval(iK, iwt, i_in));
//        }
//
//
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wa*nw2_nua*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wa*nw2_nua*n_in))/(nw2_wa*nw2_nua*n_in);
//            int iwa = (i%(nw2_wa*nw2_nua*n_in))/(nw2_nua*n_in);
//            int iva = (i%(nw2_nua*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.avertex.K2_setvert(iK, iwa, iva, i_in, avertex.densvertex.K2_vval(iK, iwa, iva, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wp*nw2_nup*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wp*nw2_nup*n_in))/(nw2_wp*nw2_nup*n_in);
//            int iwp = (i%(nw2_wp*nw2_nup*n_in))/(nw2_nup*n_in);
//            int ivp = (i%(nw2_nup*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.pvertex.K2_setvert(iK, iwp, ivp, i_in, pvertex.densvertex.K2_vval(iK, iwp, ivp, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wt*nw2_nut*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wt*nw2_nut*n_in))/(nw2_wt*nw2_nut*n_in);
//            int iwt = (i%(nw2_wt*nw2_nut*n_in))/(nw2_nut*n_in);
//            int ivt = (i%(nw2_nut*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.tvertex.K2_setvert(iK, iwt, ivt, i_in, tvertex.densvertex.K2_vval(iK, iwt, ivt, i_in));
//        }
//
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wa*nw3_nua*nw3_nuap*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wa*nw3_nua*nw3_nuap*n_in))/(nw3_wa*nw3_nua*nw3_nuap*n_in);
//            int iwa = (i%(nw3_wa*nw3_nua*nw3_nuap*n_in))/(nw3_nua*nw3_nuap*n_in);
//            int iva = (i%(nw3_nua*nw3_nuap*n_in))/(nw3_nuap*n_in);
//            int ivap= (i%(nw3_nuap*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.avertex.K3_setvert(iK, iwa, iva, ivap, i_in, avertex.densvertex.K3_vval(iK, iwa, iva, ivap, i_in));
//        }
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wp*nw3_nup*nw3_nupp*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wp*nw3_nup*nw3_nupp*n_in))/(nw3_wp*nw3_nup*nw3_nupp*n_in);
//            int iwp = (i%(nw3_wp*nw3_nup*nw3_nupp*n_in))/(nw3_nup*nw3_nupp*n_in);
//            int ivp = (i%(nw3_nup*nw3_nupp*n_in))/(nw3_nupp*n_in);
//            int ivpp= (i%(nw3_nupp*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.pvertex.K3_setvert(iK, iwp, ivp, ivpp, i_in, pvertex.densvertex.K3_vval(iK, iwp, ivp, ivpp, i_in));
//        }
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wt*nw3_nut*nw3_nutp*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wt*nw3_nut*nw3_nutp*n_in))/(nw3_wt*nw3_nut*nw3_nutp*n_in);
//            int iwt = (i%(nw3_wt*nw3_nut*nw3_nutp*n_in))/(nw3_nut*nw3_nutp*n_in);
//            int ivt = (i%(nw3_nut*nw3_nutp*n_in))/(nw3_nutp*n_in);
//            int ivtp= (i%(nw3_nutp*n_in))/(n_in);
//            int i_in = i%n_in;
//            densvertex.tvertex.K3_setvert(iK, iwt, ivt, ivtp, i_in, tvertex.densvertex.K3_vval(iK, iwt, ivt, ivtp, i_in));
//        }
//    }
//
//    template <typename Q> void set_spinvertex(Vertex<avert<Q> >& avertex, Vertex<pvert<Q> >& pvertex, Vertex<tvert<Q> >& tvertex) {
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wa * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
//            int iwa = (i % (nw1_wa * n_in)) / n_in;
//            int i_in = i % n_in;
//            spinvertex.avertex.K1_setvert(iK, iwa, i_in, avertex.spinvertex.K1_vval(iK, iwa, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wp * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
//            int iwp = (i % (nw1_wp * n_in)) / n_in;
//            int i_in = i % n_in;
//            spinvertex.pvertex.K1_setvert(iK, iwp, i_in, pvertex.spinvertex.K1_vval(iK, iwp, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K1 * nw1_wt * n_in; ++i) {
//            int iK = (i % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
//            int iwt = (i % (nw1_wt * n_in)) / n_in;
//            int i_in = i % n_in;
//            spinvertex.tvertex.K1_setvert(iK, iwt, i_in, tvertex.spinvertex.K1_vval(iK, iwt, i_in));
//        }
//
//
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wa*nw2_nua*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wa*nw2_nua*n_in))/(nw2_wa*nw2_nua*n_in);
//            int iwa = (i%(nw2_wa*nw2_nua*n_in))/(nw2_nua*n_in);
//            int iva = (i%(nw2_nua*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.avertex.K2_setvert(iK, iwa, iva, i_in, avertex.spinvertex.K2_vval(iK, iwa, iva, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wp*nw2_nup*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wp*nw2_nup*n_in))/(nw2_wp*nw2_nup*n_in);
//            int iwp = (i%(nw2_wp*nw2_nup*n_in))/(nw2_nup*n_in);
//            int ivp = (i%(nw2_nup*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.pvertex.K2_setvert(iK, iwp, ivp, i_in, pvertex.spinvertex.K2_vval(iK, iwp, ivp, i_in));
//        }
//#pragma omp parallel for
//        for (int i = 0; i < nK_K2*nw2_wt*nw2_nut*n_in; ++i){
//            int iK = (i%(nK_K2*nw2_wt*nw2_nut*n_in))/(nw2_wt*nw2_nut*n_in);
//            int iwt = (i%(nw2_wt*nw2_nut*n_in))/(nw2_nut*n_in);
//            int ivt = (i%(nw2_nut*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.tvertex.K2_setvert(iK, iwt, ivt, i_in, tvertex.spinvertex.K2_vval(iK, iwt, ivt, i_in));
//        }
//
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wa*nw3_nua*nw3_nuap*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wa*nw3_nua*nw3_nuap*n_in))/(nw3_wa*nw3_nua*nw3_nuap*n_in);
//            int iwa = (i%(nw3_wa*nw3_nua*nw3_nuap*n_in))/(nw3_nua*nw3_nuap*n_in);
//            int iva = (i%(nw3_nua*nw3_nuap*n_in))/(nw3_nuap*n_in);
//            int ivap= (i%(nw3_nuap*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.avertex.K3_setvert(iK, iwa, iva, ivap, i_in, avertex.spinvertex.K3_vval(iK, iwa, iva, ivap, i_in));
//        }
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wp*nw3_nup*nw3_nupp*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wp*nw3_nup*nw3_nupp*n_in))/(nw3_wp*nw3_nup*nw3_nupp*n_in);
//            int iwp = (i%(nw3_wp*nw3_nup*nw3_nupp*n_in))/(nw3_nup*nw3_nupp*n_in);
//            int ivp = (i%(nw3_nup*nw3_nupp*n_in))/(nw3_nupp*n_in);
//            int ivpp= (i%(nw3_nupp*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.pvertex.K3_setvert(iK, iwp, ivp, ivpp, i_in, pvertex.spinvertex.K3_vval(iK, iwp, ivp, ivpp, i_in));
//        }
//#pragma omp parallel for
//        for(int i=0; i<nK_K3*nw3_wt*nw3_nut*nw3_nutp*n_in; ++i){
//            int iK = (i%(nK_K3*nw3_wt*nw3_nut*nw3_nutp*n_in))/(nw3_wt*nw3_nut*nw3_nutp*n_in);
//            int iwt = (i%(nw3_wt*nw3_nut*nw3_nutp*n_in))/(nw3_nut*nw3_nutp*n_in);
//            int ivt = (i%(nw3_nut*nw3_nutp*n_in))/(nw3_nutp*n_in);
//            int ivtp= (i%(nw3_nutp*n_in))/(n_in);
//            int i_in = i%n_in;
//            spinvertex.tvertex.K3_setvert(iK, iwt, ivt, ivtp, i_in, tvertex.spinvertex.K3_vval(iK, iwt, ivt, ivtp, i_in));
//        }
//    }
};


/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> Q irreducible<Q>::vval(int iK) {
    return bare[iK];
}
template <typename Q> void irreducible<Q>::setvert(int iK, Q value){
    bare[iK] = value;
}


/************************************* MEMBER FUNCTIONS OF THE VERTEX "fullvertex" ************************************/
//arguments are equivalent to those in the simple vertex functions

//template <typename Q> Q fullvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel){
//    Q result = irred.vvalsmooth() + pvertex.vvalsmooth(iK,q,w1,w2,i_in,channel) + tvertex.vvalsmooth(iK,q,w1,w2,i_in,channel) + avertex.vvalsmooth(iK,q,w1,w2,i_in,channel);
//    if(abs(result)>1e-16){
//        return  result;}
//    else{return 0.;}
//}
//template <typename Q> Q fullvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel, int p, char f){
//    Q result=0;
//
//    if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
//        result += irred.vvalsmooth();
//
//    }
//    else if( p==2 && (f=='K' || f== 'M')){
//        result += irred.vvalsmooth();
//
//    }
//
//    result +=  pvertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f) + tvertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f) + avertex.vvalsmooth(iK,q,w1,w2,i_in,channel,p,f);
//    //if(p==2 && f=='L'){cout <<  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << " " << tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f)  << " " <<  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << endl;};
//    if(abs(result)<1e-16){
//        result  =0;}
//    return result;
//}
//template <typename Q> Q fullvert<Q>::vvalsmooth(int red_side, int map, int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){// red_side: if only complementary channel in one vertex, which one is reduced? (0/1/2), p: is this the left/upper (1) or the right/lower (2) vertex of the bubble?, f: diagrammatic class that is computed
//    Q result=0;
//
//    if(red_side != p){
//        if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
//            result = irred.vvalsmooth(a,b,c);
//        }
//        else if( p==2 && (f=='K' || f== 'M')){
//            result = irred.vvalsmooth(a,b,c);
//        };
//        result +=  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
//    }
//
//
//    else if(red_side == p){
//
//        if(map==0){
//            if(channel == 's'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
//            else if(channel == 't'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
//            else if(channel == 'u'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}}
//        else if(map==1){
//            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
//            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
//            else if(channel == 'u'){  result =  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
//        }
//
//    };
//
//    if(abs(result)<1e-16){
//        result  =0;};
//    return result;
//}

template <typename Q> Q fullvert<Q>::value (int iK, double w, double v1, double v2, int i_in, char channel)
{
    Q result = irred.vval(iK) + avertex.value(iK, w, v1, v2, i_in, channel, tvertex) + pvertex.value(iK, w, v1, v2, i_in, channel) + tvertex.value(iK, w, v1, v2, i_in, channel, avertex);
    return result;
}

template <typename Q> Q fullvert<Q>::gammaRb (int iK, double w, double v1, double v2, int i_in, char r)
{
   Q resp;
    switch(r){
        case 'a':
            resp = pvertex.value(iK, w, v1, v2, i_in, r)          + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
            break;
        case 'p':
            resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
            break;
        case 't':
            resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + pvertex.value(iK, w, v1, v2, i_in, r);
            break;
        default :
            resp = 0.;
            cout << "Something's going wrong with gammaRb"<< endl;
    }

//    if(r == 'a')
//        resp = pvertex.value(iK, w, v1, v2, i_in, r)          + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
//    else if(r == 'p')
//        resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
//    else if(r == 't')
//        resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + pvertex.value(iK, w, v1, v2, i_in, r);
//    else
//        resp = 0.;

    return resp;
}

#endif //KELDYSH_MFRG_VERTEX_H
