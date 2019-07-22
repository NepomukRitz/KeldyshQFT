//
// Created by Sa.Aguirre on 7/19/19.
//

#include "vertex.h"


/*****************FUNCTIONS FOR THE A-VERTEX********************************************/

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q avert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double u=0.,w1_u=0., w2_u=0.;
    if(channel == 's'){
        u = -w2-w1;
        w1_u = (w1-w2+q)/2;
        w2_u = (-w1+w2+q)/2;}
    else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2-w1;
        w1_u = (w1+w2-q)/2;
        w2_u = (w1+w2+q)/2;}
    else if(channel == 'u'){
        u = q;
        w1_u = w1;
        w2_u = w2;}
    else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1-w2;
        w1_u = q + (w1-w2)/2;
        w2_u = (w1+w2)/2;};

    Q value;

//    if(abs(u) < freqs_a[nw1_wa/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//    if(abs(w1_u) < freqs_a[nw1_wa/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//    if(abs(w2_u) < freqs_a[nw1_wa/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK,u,w1_u,i_in) + K3_vvalsmooth(iK, u, w1_u, w2_u, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,u,w2_u,i_in)

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
    };
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

//overload of previous function         => I'm pretty sure we don't andwon't be needing this function
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
    vec_K1[iK*nw1_wa*n_in + i*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    vec_K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in] = value;
}

template <typename Q> Q avert<Q>::K1_vval(int iK, int i, int i_in){
    return vec_K1[iK*nw1_wa*n_in + i*n_in + i_in];
}
template <typename Q> Q avert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return vec_K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in];
}
template <typename Q> Q avert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in];
}

template <typename Q> Q avert<Q>::K1_vvalsmooth(int iK, double u, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    /*I assume a linear grid for this implementation. TODO: define a set of functions (one per type of grid to use) to convert double in range to index on grid*/
    /*TODO: since we might want to interpolate according to channel and diagrammatic class, rethink previous TODO */

    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw1_wa-1));     //nw1_wa because we're interpolating for K1 in channel a
    auto index = (int)((u-w_lower_b)/dw_b_a);

    double x1 = freqs_a[index];
    double x2 = freqs_a[index+1];

    Q f1 = vec_K1[iK*nw1_wa*n_in + (index)*n_in + i_in];
    Q f2 = vec_K1[iK*nw1_wa*n_in + (index+1)*n_in + i_in];

    return f1 + (u-x1)*(f2-f1)/(x2-x1);
}
template <typename Q> Q avert<Q>::K2_vvalsmooth(int iK, double u, double w1, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw2_wa-1));     //nw2_wa because we're interpolating for bosonic freq in K2 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw2_nua-1));    //nw2_nua because we're interpolating for fermionic freq in K2 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    Q f11 = vec_K2[iK*nw2_wa*nw2_nua*n_in + (index_b)*nw2_nua*n_in + (index_f)*n_in + i_in];
    Q f12 = vec_K2[iK*nw2_wa*nw2_nua*n_in + (index_b)*nw2_nua*n_in + (index_f+1)*n_in + i_in];
    Q f21 = vec_K2[iK*nw2_wa*nw2_nua*n_in + (index_b+1)*nw2_nua*n_in + (index_f)*n_in + i_in];
    Q f22 = vec_K2[iK*nw2_wa*nw2_nua*n_in + (index_b+1)*nw2_nua*n_in + (index_f+1)*n_in + i_in];

    return (y2-w1)/(y2-y1)*((x2-u)/(x2-x1)*f11 + (u-x1)/(x2-x1)*f21) + (w1-y1)/(y2-y1)*((x2-u)/(x2-x1)*f12 + (u-x1)/(x2-x1)*f22);
}
template <typename Q> Q avert<Q>::K3_vvalsmooth(int iK, double u, double w1, double w2, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw3_wa-1));     //nw3_wa because we're interpolating for bosonic freq in K3 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nua-1));    //nw3_nua because we're interpolating for fermionic freq in K3 in channel a
    double dwp_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nuap-1));  //nw3_nuap because we're interpolating for fermionic freq in K3 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);
    auto index_fp = (int)((w2-w_lower_f)/dwp_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];
    double z1 = freqs_a[index_fp];
    double z2 = freqs_a[index_fp+1];

    Q f111 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b)*nw3_nua*nw3_nuap*n_in + (index_f)*nw3_nuap*n_in + (index_fp)*n_in + i_in];
    Q f112 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b)*nw3_nua*nw3_nuap*n_in + (index_f)*nw3_nuap*n_in + (index_fp+1)*n_in + i_in];
    Q f121 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b)*nw3_nua*nw3_nuap*n_in + (index_f+1)*nw3_nuap*n_in + (index_fp)*n_in + i_in];
    Q f122 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b)*nw3_nua*nw3_nuap*n_in + (index_f+1)*nw3_nuap*n_in + (index_fp+1)*n_in + i_in];
    Q f211 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b+1)*nw3_nua*nw3_nuap*n_in + (index_f)*nw3_nuap*n_in + (index_fp)*n_in + i_in];
    Q f212 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b+1)*nw3_nua*nw3_nuap*n_in + (index_f)*nw3_nuap*n_in + (index_fp+1)*n_in + i_in];
    Q f221 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b+1)*nw3_nua*nw3_nuap*n_in + (index_f+1)*nw3_nuap*n_in + (index_fp)*n_in + i_in];
    Q f222 = vec_K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + (index_b+1)*nw3_nua*nw3_nuap*n_in + (index_f+1)*nw3_nuap*n_in + (index_fp+1)*n_in + i_in];

    double xd = (u-x1)/(x2-x1);
    double yd = (w1-y1)/(y2-y1);
    double zd = (w2-z1)/(z2-z1);

    Q c00 = f111*(1-xd) + f211*xd;
    Q c01 = f112*(1-xd) + f212*xd;
    Q c10 = f121*(1-xd) + f221*xd;
    Q c11 = f122*(1-xd) + f222*xd;

    Q c0 = c00*(1-yd) + c10*yd;
    Q c1 = c01*(1-yd) + c11*yd;

    return c0*(1-zd) + c1*zd;
}
//non-member functions
template <typename Q> avert<Q> operator*(double alpha, const avert<Q>& vertex){
    avert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wa; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wa*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wa*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wa; ++i){
                for(int j=0; j<nw2_nua; ++j){
                    vertex2.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wa; ++i){
                for(int j=0; j<nw3_nua; j++){
                    for(int k=0; k<nw3_nuap; ++k){
                        vertex2.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> avert<Q> operator*(const avert<Q>& vertex,double alpha){
    avert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wa; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wa*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wa*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wa; ++i){
                for(int j=0; j<nw2_nua; ++j){
                    vertex2.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wa; ++i){
                for(int j=0; j<nw3_nua; j++){
                    for(int k=0; k<nw3_nuap; ++k){
                        vertex2.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> avert<Q> operator+(const avert<Q>& vertex1,const avert<Q>& vertex2){
    avert<Q> vertex3;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wa; ++i){
                //K1 contributions
                vertex3.vec_K1[nk*nw1_wa*n_in + i*n_in + nin] = vertex1.vec_K1[nk*nw1_wa*n_in + i*n_in + nin] + vertex2.vec_K1[nk*nw1_wa*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wa; ++i){
                for(int j=0; j<nw2_nua; ++j){
                    vertex3.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin] = vertex1.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin] + vertex2.vec_K2[nk*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wa; ++i){
                for(int j=0; j<nw3_nua; j++){
                    for(int k=0; k<nw3_nuap; ++k){
                        vertex3.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin] =
                                vertex1.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin] + vertex2.vec_K3[nk*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}


/*

template <typename Q> Q avert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double u,w1_u,w2_u;
        if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2;
            w2_u = (-w1+w2+q)/2;}
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2;
            w2_u = (w1+w2+q)/2;}
        else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2;
            w2_u = (w1+w2)/2;};

        Q value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};



//        if(u > bfreqs[nw-1] && abs(u - bfreqs[nw-1]) < 1e-12){u = bfreqs[nw-1];}
//        else if(u<bfreqs[0] && abs(u - bfreqs[0]) < 1e-12){u = bfreqs[0];};


//        if(w1_u > ffreqs[nw-1] && abs(w1_u - ffreqs[nw-1]) < 1e-12){w1_u = ffreqs[nw-1];}
//        else if(w1_u<ffreqs[0] && abs(w1_u - ffreqs[0]) < 1e-12){w1_u = ffreqs[0];};


//        if(w2_u > ffreqs[nw-1] && abs(w2_u - ffreqs[nw-1]) < 1e-12){w2_u = ffreqs[nw-1];}
//        else if(w2_u<ffreqs[0] && abs(w2_u - ffreqs[0]) < 1e-12){w2_u = ffreqs[0];};


        value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u)  + K2b_vvalsmooth(a,b,c,u,w2_u)  ;//K2b is extracted from K2 by the symmetry relations

        return value;  }
    else{return 0;}



}
template <typename Q> Q avert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
    if(distance(a,b,c) <= d_c){//cutoff distance



        double u,w1_u,w2_u;
        if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2;
            w2_u = (-w1+w2+q)/2;}
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2;
            w2_u = (w1+w2+q)/2;}
        else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2;
            w2_u = (w1+w2)/2;};

        Q value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};



//        if(u > bfreqs[nw-1] && abs(u - bfreqs[nw-1]) < 1e-12){u = bfreqs[nw-1];}
//        else if(u<bfreqs[0] && abs(u - bfreqs[0]) < 1e-12){u = bfreqs[0];};


//        if(w1_u > ffreqs[nw-1] && abs(w1_u - ffreqs[nw-1]) < 1e-12){w1_u = ffreqs[nw-1];}
//        else if(w1_u<ffreqs[0] && abs(w1_u - ffreqs[0]) < 1e-12){w1_u = ffreqs[0];};


//        if(w2_u > ffreqs[nw-1] && abs(w2_u - ffreqs[nw-1]) < 1e-12){w2_u = ffreqs[nw-1];}
//        else if(w2_u<ffreqs[0] && abs(w2_u - ffreqs[0]) < 1e-12){w2_u = ffreqs[0];};

        if(p==1){
            if(channel=='u'){
                if(f == 'R' || f == 'M'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='t'){
                if(f == 'R' || f== 'M'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);};
            };
        }

        else if(p==2){
            if(channel=='u'){
                if(f == 'R' || f == 'L'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,u) + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='t'){
                if(f == 'R' || f== 'L'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);
                };
            };
        };




        return value;  }
    else{return 0;}



}
//overload of previous function
template <typename Q> Q avert<Q>::vvalsmooth(int red_side, int map,int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
    return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
}
template <typename Q> Q avert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance


        double u,w1_u,w2_u;

        u = q;
        w1_u = w1;
        w2_u = w2;

        Q value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

        value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
        return value;  }
    else{return 0;}



}
template <typename Q> void avert<Q>::K1_setvert(int a, int b, int c, int i, Q value){

    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>= 0 ){
            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> void avert<Q>::K2_setvert(int a, int b, int c,int i, int j,Q value){


    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0 ){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> void avert<Q>::K3_setvert(int a, int b, int c, int i, int j, int k, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){


            K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> Q avert<Q>::K1_vval(int a_raw, int b_raw, int c_raw,  int i){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw1 && i >= 0 ){
        if(sym==0){

            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 || sym==2){
            if (i>=nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);
                i=nw1-1-i;
                value += K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)];};
        };};
    if(abs(value)<1e-100){value =0;};

    return value;
}
template <typename Q> Q avert<Q>::K2_vval(int a_raw, int b_raw, int c_raw,  int i,int j){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0 ){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j >= nw2/2){


                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
            else if (j<nw2/2){

                site x = site_switch(a,b,c);

                j = nw2-1-j;
                value += conj(K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)]);};
        }
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j;x = site_switch(x.a,x.b,x.c);};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };};
    if(abs(value)<1e-100){value =0;};



    return value;
}
template <typename Q> Q avert<Q>::K3_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0 ){

            value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}


        else if(sym == 1 ){
            if(i>=nw3/2){
                if(j >=nw3/2){

                    value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

                else if(j < nw3/2){
                    site x = site_switch(a,b,c);

                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)]);};}

            else if (i<nw3/2 ){
                if(k >= nw3/2){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)];}//note that k and j are interchanged

                else if ( k < nw3/2){


                    i = nw3-1-i;
                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);};//note that k and j are interchanged

            };}
        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;x = site_switch(x.a,x.b,x.c);};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};
            value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};};
    if(abs(value)<1e-100){value=0;};



    return value;
}
template <typename Q> Q avert<Q>::K1_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(u)+ 1e-6<= ffreqs[(nw1+nw)/2-1]){
        if(sym==0){

            int U = fconv_n(u,nw1);

            int U_vert = fconv_n(u,nw1)-(nw1-nw1_q);

            value += ((K1[a+(nuc_eff-1)/2][b][c-1][U_vert ]*(bfreqs[(nw-nw1)/2+U+1]-u)
                       +K1[a+(nuc_eff-1)/2][b][c-1][U_vert +1]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));}
        else if(sym==1 ||sym==2){
            if (u>0){


                int U = fconv_n(u,nw1);

                int U_vert = fconv_n(u,nw1)-(nw1-nw1_q);

                value += ((K1[a+(nuc_eff-1)/2][b][c-1][U_vert]*(bfreqs[(nw-nw1)/2+U+1]-u)
                           +K1[a+(nuc_eff-1)/2][b][c-1][U_vert+1]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));}
            else if (u<0){

                site x = site_switch(a,b,c);
                int U = fconv_n(-u,nw1);

                int U_vert = fconv_n(-u,nw1)-(nw1-nw1_q);
                u = -u;

                value += (K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert]*(bfreqs[(nw-nw1)/2+U+1]-u)
                          +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1]*(-bfreqs[(nw-nw1)/2+U]+u))/(bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U]);};};}
    else{
        int i;
        if(abs(u)<= bfreqs[nw1-1] ){
            i= fconv_n(u,nw1);

            value += K1_vval(a,b,c,i);};    };


    return value;
}
template <typename Q> Q avert<Q>::K2_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    Q value=0;
    if(abs(u)+ 1e-6<= bfreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1,nw2);

            int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){

            if(w1 > 0){

                int U = fconv_n(u,nw2);
                int W1 = fconv_n(w1,nw2);

                int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
            }
            else if (w1<0){

                site x = site_switch(a,b,c);

                int U = fconv_n(u,nw2);
                int W1 = fconv_n(-w1,nw2);

                int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(-w1,nw2)-(nw2-nw2_w1);
                double w1_eff = -w1 ;
                value += conj((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert ][W1_vert ]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                               +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert +1][W1_vert ]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                               +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert ][W1_vert +1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                               +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert +1][W1_vert +1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

            };
        }
        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1;
            if(u<0){u = -u;x = site_switch(x.a,x.b,x.c); };
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);
            value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        };}
    else{ int i,j;
        if(abs(u)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(u,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};   };



    return value;
}
template <typename Q> Q avert<Q>::K3_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(u)+ 1e-6< ffreqs[(nw3+nw)/2-1] && abs(w1)+ 1e-6<ffreqs[(nw3+nw)/2-1]&& abs(w2)+ 1e-6<=ffreqs[(nw3+nw)/2-1]){
        if (sym==0){

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);


            int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

            value += ((K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                      ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){
            if(u > 0 ){
                if( (w1)>0){

                    int U = fconv_n(u,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);


                    int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);
                    value += ((K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                              ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));}
                else if((w1)<0){
                    site x = site_switch(a,b,c);

                    int U = fconv_n(u,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);

                    int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;
                    value += conj((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                  ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}

            else if (u < 0 ){
                if((w2)>0){
                    site x = site_switch(a,b,c);

                    int U = fconv_n(-u,nw3);
                    int W1 = fconv_n(w2 ,nw3);
                    int W2 = fconv_n(w1 ,nw3);


                    int U_vert = fconv_n(-u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2 ,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1 ,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;
                    u = -u;

                    value += ((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                              ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }


                else if ((w2)<0){


                    int U = fconv_n(-u,nw3);
                    int W1 = fconv_n(-w2 ,nw3);
                    int W2 = fconv_n(-w1 ,nw3);

                    int U_vert = fconv_n(-u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2 ,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1 ,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w2;
                    double w2_eff = -w1;
                    u = -u;

                    value += conj((K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                  ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));};

            };}

        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, u_eff =u;
            if(u<0){u_eff = -u; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if(w1_eff < 0){w1_eff = -w1_eff; w2_eff = -w2_eff;x = site_switch(x.a,x.b,x.c);};

            u = u_eff;

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);
            value = ((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }


        ;}
    else{
        int i,j,k;
        if(abs(u)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(u,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value = K3_vval(a,b,c,i,j,k);};
    };




    return value;
}

*/

/*****************FUNCTIONS FOR THE P-VERTEX********************************************/

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q pvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double u=0.,w1_u=0., w2_u=0.;
    if(channel == 's'){
        u = -w2-w1;
        w1_u = (w1-w2+q)/2;
        w2_u = (-w1+w2+q)/2;}
    else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2-w1;
        w1_u = (w1+w2-q)/2;
        w2_u = (w1+w2+q)/2;}
    else if(channel == 'u'){
        u = q;
        w1_u = w1;
        w2_u = w2;}
    else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1-w2;
        w1_u = q + (w1-w2)/2;
        w2_u = (w1+w2)/2;};

    Q value;

//    if(abs(u) < freqs_a[nw1_wa/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//    if(abs(w1_u) < freqs_a[nw1_wa/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//    if(abs(w2_u) < freqs_a[nw1_wa/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK,u,w1_u,i_in) + K3_vvalsmooth(iK, u, w1_u, w2_u, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,u,w2_u,i_in)

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
    };
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

//overload of previous function         => I'm pretty sure we don't andwon't be needing this function
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
    vec_K1[iK*nw1_wp*n_in + i*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    vec_K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] = value;
}

template <typename Q> Q pvert<Q>::K1_vval(int iK, int i, int i_in){
    return vec_K1[iK*nw1_wp*n_in + i*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return vec_K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in];
}

template <typename Q> Q pvert<Q>::K1_vvalsmooth(int iK, double u, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    /*I assume a linear grid for this implementation. TODO: define a set of functions (one per type of grid to use) to convert double in range to index on grid*/
    /*TODO: since we might want to interpolate according to channel and diagrammatic class, rethink previous TODO */

    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw1_wa-1));     //nw1_wa because we're interpolating for K1 in channel a
    auto index = (int)((u-w_lower_b)/dw_b_a);

    double x1 = freqs_a[index];
    double x2 = freqs_a[index+1];

    Q f1 = vec_K1[iK*nw1_wp*n_in + (index)*n_in + i_in];
    Q f2 = vec_K1[iK*nw1_wp*n_in + (index+1)*n_in + i_in];

    return f1 + (u-x1)*(f2-f1)/(x2-x1);
}
template <typename Q> Q pvert<Q>::K2_vvalsmooth(int iK, double u, double w1, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw2_wa-1));     //nw2_wa because we're interpolating for bosonic freq in K2 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw2_nua-1));    //nw2_nua because we're interpolating for fermionic freq in K2 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    Q f11 = vec_K2[iK*nw2_wp*nw2_nup*n_in + (index_b)*nw2_nup*n_in + (index_f)*n_in + i_in];
    Q f12 = vec_K2[iK*nw2_wp*nw2_nup*n_in + (index_b)*nw2_nup*n_in + (index_f+1)*n_in + i_in];
    Q f21 = vec_K2[iK*nw2_wp*nw2_nup*n_in + (index_b+1)*nw2_nup*n_in + (index_f)*n_in + i_in];
    Q f22 = vec_K2[iK*nw2_wp*nw2_nup*n_in + (index_b+1)*nw2_nup*n_in + (index_f+1)*n_in + i_in];

    return (y2-w1)/(y2-y1)*((x2-u)/(x2-x1)*f11 + (u-x1)/(x2-x1)*f21) + (w1-y1)/(y2-y1)*((x2-u)/(x2-x1)*f12 + (u-x1)/(x2-x1)*f22);
}
template <typename Q> Q pvert<Q>::K3_vvalsmooth(int iK, double u, double w1, double w2, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw3_wa-1));     //nw3_wa because we're interpolating for bosonic freq in K3 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nua-1));    //nw3_nua because we're interpolating for fermionic freq in K3 in channel a
    double dwp_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nuap-1));  //nw3_nuap because we're interpolating for fermionic freq in K3 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);
    auto index_fp = (int)((w2-w_lower_f)/dwp_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];
    double z1 = freqs_a[index_fp];
    double z2 = freqs_a[index_fp+1];

    Q f111 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b)*nw3_nup*nw3_nupp*n_in + (index_f)*nw3_nupp*n_in + (index_fp)*n_in + i_in];
    Q f112 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b)*nw3_nup*nw3_nupp*n_in + (index_f)*nw3_nupp*n_in + (index_fp+1)*n_in + i_in];
    Q f121 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b)*nw3_nup*nw3_nupp*n_in + (index_f+1)*nw3_nupp*n_in + (index_fp)*n_in + i_in];
    Q f122 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b)*nw3_nup*nw3_nupp*n_in + (index_f+1)*nw3_nupp*n_in + (index_fp+1)*n_in + i_in];
    Q f211 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b+1)*nw3_nup*nw3_nupp*n_in + (index_f)*nw3_nupp*n_in + (index_fp)*n_in + i_in];
    Q f212 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b+1)*nw3_nup*nw3_nupp*n_in + (index_f)*nw3_nupp*n_in + (index_fp+1)*n_in + i_in];
    Q f221 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b+1)*nw3_nup*nw3_nupp*n_in + (index_f+1)*nw3_nupp*n_in + (index_fp)*n_in + i_in];
    Q f222 = vec_K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + (index_b+1)*nw3_nup*nw3_nupp*n_in + (index_f+1)*nw3_nupp*n_in + (index_fp+1)*n_in + i_in];

    double xd = (u-x1)/(x2-x1);
    double yd = (w1-y1)/(y2-y1);
    double zd = (w2-z1)/(z2-z1);

    Q c00 = f111*(1-xd) + f211*xd;
    Q c01 = f112*(1-xd) + f212*xd;
    Q c10 = f121*(1-xd) + f221*xd;
    Q c11 = f122*(1-xd) + f222*xd;

    Q c0 = c00*(1-yd) + c10*yd;
    Q c1 = c01*(1-yd) + c11*yd;

    return c0*(1-zd) + c1*zd;
}
//non-member functions
template <typename Q> pvert<Q> operator*(double alpha, const pvert<Q>& vertex){
    pvert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wp; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wp*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wp*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wp; ++i){
                for(int j=0; j<nw2_nup; ++j){
                    vertex2.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wp; ++i){
                for(int j=0; j<nw3_nup; j++){
                    for(int k=0; k<nw3_nupp; ++k){
                        vertex2.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> pvert<Q> operator*(const pvert<Q>& vertex,double alpha){
    pvert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wp; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wp*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wp*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wp; ++i){
                for(int j=0; j<nw2_nup; ++j){
                    vertex2.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wp; ++i){
                for(int j=0; j<nw3_nup; j++){
                    for(int k=0; k<nw3_nupp; ++k){
                        vertex2.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> pvert<Q> operator+(const pvert<Q>& vertex1,const pvert<Q>& vertex2){
    pvert<Q> vertex3;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wp; ++i){
                //K1 contributions
                vertex3.vec_K1[nk*nw1_wp*n_in + i*n_in + nin] = vertex1.vec_K1[nk*nw1_wp*n_in + i*n_in + nin] + vertex2.vec_K1[nk*nw1_wp*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wp; ++i){
                for(int j=0; j<nw2_nup; ++j){
                    vertex3.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin] = vertex1.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin] + vertex2.vec_K2[nk*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wp; ++i){
                for(int j=0; j<nw3_nup; j++){
                    for(int k=0; k<nw3_nupp; ++k){
                        vertex3.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin] =
                                vertex1.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin] + vertex2.vec_K3[nk*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}

/*

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q pvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (w1-w2+q)/2;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (-w1+w2-q)/2;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2;
            w2_s = w2-(q+w1)/2;};
        Q value=0;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};


//        if(s > bfreqs[nw-1] && abs(s - bfreqs[nw-1]) < 1e-12){s = bfreqs[nw-1];}
//        else if(s<bfreqs[0] && abs(s - bfreqs[0]) < 1e-12){s = bfreqs[0];};


//        if(w1_s > ffreqs[nw-1] && abs(w1_s - ffreqs[nw-1]) < 1e-12){w1_s = ffreqs[nw-1];}
//        else if(w1_s<ffreqs[0] && abs(w1_s - ffreqs[0]) < 1e-12){w1_s = ffreqs[0];};


//        if(w2_s > ffreqs[nw-1] && abs(w2_s - ffreqs[nw-1]) < 1e-12){w2_s = ffreqs[nw-1];}
//        else if(w2_s<ffreqs[0] && abs(w2_s - ffreqs[0]) < 1e-12){w2_s = ffreqs[0];};

        value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations


        return value;}
    else{return 0;}

}
//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
template <typename Q> Q pvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (w1-w2+q)/2;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (-w1+w2-q)/2;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2;
            w2_s = w2-(q+w1)/2;};
        Q value=0;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};


//        if(s > bfreqs[nw-1] && abs(s - bfreqs[nw-1]) < 1e-12){s = bfreqs[nw-1];}
//        else if(s<bfreqs[0] && abs(s - bfreqs[0]) < 1e-12){s = bfreqs[0];};


//        if(w1_s > ffreqs[nw-1] && abs(w1_s - ffreqs[nw-1]) < 1e-12){w1_s = ffreqs[nw-1];}
//        else if(w1_s<ffreqs[0] && abs(w1_s - ffreqs[0]) < 1e-12){w1_s = ffreqs[0];};


//        if(w2_s > ffreqs[nw-1] && abs(w2_s - ffreqs[nw-1]) < 1e-12){w2_s = ffreqs[nw-1];}
//        else if(w2_s<ffreqs[0] && abs(w2_s - ffreqs[0]) < 1e-12){w2_s = ffreqs[0];};

        if(p==1){
            if(channel=='s'){
                if(f == 'R' || f == 'M'){value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2b_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'M'){value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        }

        else if(p==2){
            if(channel=='s'){
                if(f == 'R' || f == 'L'){value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to differentbare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,s) + K2b_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'L'){value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        };



        return value;}
    else{return 0;}

}
//overload with first extra firs two arguments: red_side = vertex number with only complementary channels (0,1,2), map = determines if channel mapping is turned on for this vertex (0=off/1=on)
template <typename Q> Q pvert<Q>::vvalsmooth(int red_side,int map, int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){

    //THIS FUNCTION IS NEEDED ONLY WHEN BUBBLE FUNCTIONS ARE USED WITH VERTEX OF TYPE "PARVERT" INSTEAD OF "FULLVERT". LEADS TO PREVIOUS FUNCTION.
    return vvalsmooth(a, b, c, q, w1,w2, channel, p, f);
}
template <typename Q> Q pvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;
        s = q;
        w1_s = w1;
        w2_s = w2;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};



        Q value=0;
        value += K3_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}

}
template <typename Q> void pvert<Q>::K1_setvert(int a, int b, int c,int i, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>=0 ){

            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};

    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
template <typename Q> void pvert<Q>::K2_setvert(int a, int b, int c, int i, int j, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value ;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> void pvert<Q>::K3_setvert(int a, int b, int c, int i, int j, int k, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){

            K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};
    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> Q pvert<Q>::K1_vval(int a_raw, int b_raw, int c_raw,  int i){



    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i<nw1 && i>= 0 ){
        if (sym==0){

            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 || sym==2){
            if (i>= nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);
                i = nw1-1-i;
                value += conj(K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)]);};
        }

        ;};
    if(abs(value)<1e-100){value=0;};



    return value;
}
template <typename Q> Q pvert<Q>::K2_vval(int a_raw, int b_raw, int c_raw,  int i, int j ){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j >= nw2/2){

                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
            else if (j< nw2/2){


                site x = site_switch(a,b,c);
                j = nw2 -1 - j;

                value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];

            };

        }
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i; x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j; x = site_switch(x.a,x.b,x.c);};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };
    };

    if(abs(value)<1e-100){value=0;};


    return value;
}
template <typename Q> Q pvert<Q>::K3_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    Q value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0){

            value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

        else if(sym == 1 ){

            if(i>=nw3/2){
                if(j>=nw3/2){
                    value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}
                else if(j<nw3/2){
                    site x = site_switch(a,b,c);
                    j=nw3-1-j;
                    k = nw3-1-k;
                    value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];
                };

            }
            else if(i<nw3/2){
                if(k>=nw3/2){
                    site x = site_switch(a,b,c);
                    i=nw3-1-i;

                    value += conj(K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);//note that j and k are interchanged
                }
                else if(k<nw3/2){
                    i=nw3-1-i;
                    j=nw3-1-j;
                    k = nw3-1-k;
                    value +=conj(K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);//note that j and k are interchanged
                }
            };
        }
        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){x = site_switch(x.a,x.b,x.c);i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};
            value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};


    };
    if(abs(value)<1e-100){value=0;};


    return value;
}
template <typename Q> Q pvert<Q>::K1_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(s) + 1e-6 < ffreqs[(nw1+nw)/2-1]){
        if(sym==0){
            int S = fconv_n(s,nw1);
            int S_vert = fconv_n(s,nw1)-(nw1-nw1_q);
            value = ((K1[a+(nuc_eff-1)/2][b][c-1][S_vert]*(bfreqs[(nw-nw1)/2+S+1]-s)
                      +K1[a+(nuc_eff-1)/2][b][c-1][S_vert+1]*(-bfreqs[(nw-nw1)/2+S]+s))/(bfreqs[(nw-nw1)/2+S+1]-bfreqs[(nw-nw1)/2+S]));}
        else if(sym==1 || sym==2){


            if (s > 0){

                int S = fconv_n(s,nw1);
                int S_vert = fconv_n(s,nw1)-(nw1-nw1_q);
                value += ((K1[a+(nuc_eff-1)/2][b][c-1][S_vert]*(bfreqs[S+1]-s)
                           +K1[a+(nuc_eff-1)/2][b][c-1][S_vert+1]*(-bfreqs[S]+s))/(bfreqs[S+1]-bfreqs[S]));}
            else if (s<0){
                site x = site_switch(a,b,c);

                double s_eff = -s;
                int S = fconv_n(s_eff,nw1);
                int S_vert = fconv_n(s_eff,nw1)-(nw1-nw1_q);

                if(sym==1){

                    value = conj((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert]*(bfreqs[S+1]-s_eff)
                                  +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));}
                else{
                    value = ((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert]*(bfreqs[S+1]-s_eff)
                              +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));};
            };};

    }
    else{
        int i;
        if(abs(s)<= bfreqs[nw1-1]){
            i= fconv_n(s,nw1);

            value += K1_vval(a,b,c,i);}

        ;};




    return value;

}
template <typename Q> Q pvert<Q>::K2_vvalsmooth(int a_raw, int b_raw, int c_raw,  double s, double w1){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    Q value=0;

    if(abs(s) + 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){
            int S = fconv_n(s,nw2);
            int W1 = fconv_n(w1,nw2);

            int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){
            if(w1 > 0){


                int S = fconv_n(s,nw2);
                int W1 = fconv_n(w1,nw2);


                int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
            else if(w1 < 0){

                site x = site_switch(a,b,c);

                int S = fconv_n(s,nw2);
                int W1 = fconv_n(-w1,nw2);

                int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(-w1,nw2)-(nw2-nw2_w1);
                double w1_eff = -w1 ;
                value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                           +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                           +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                           +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));};
        }
        else if(sym==2){

            site x(a,b,c);
            double w1_eff = w1;
            double s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};


            int S = fconv_n(s_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int S_vert = fconv_n(s_eff,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);

            value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                       +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

        };}
    else{
        if(abs(s)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            int i,j;

            i= fconv_n(s,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};
    };




    return value;
}
template <typename Q> Q pvert<Q>::K3_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(s)+1e-6 < ffreqs[(nw3+nw)/2-1] && abs(w1)+1e-6<ffreqs[(nw3+nw)/2-1] && abs(w2)+1e-6 <ffreqs[(nw3+nw)/2-1]){//if frequency arguments are out of range, vertex vanishes
        if (sym==0){

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);
            int S_vert = fconv_n(s,nw3);
            int W1_vert = fconv_n(w1,nw3);
            int W2_vert = fconv_n(w2,nw3);


            value += ((K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                      ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){

            if(s > 0){
                if((w1)>0 ){
                    int S = fconv_n(s,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);

                    int S_vert = fconv_n(s,nw3)-(nw3-nw3_w1);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

                    value += ((K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                               +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                              ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }
                else if((w1)<0 ){
                    site x = site_switch(a,b,c);

                    int S = fconv_n(s,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);
                    int S_vert = fconv_n(s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;

                    value += ((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                               +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                              ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}

            else if (s < 0){
                if((w2)>0){

                    site x = site_switch(a,b,c);

                    int S = fconv_n(-s,nw3);
                    int W1 = fconv_n(w2,nw3);
                    int W2 = fconv_n(w1,nw3);
                    int S_vert = fconv_n(-s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;
                    s = -s;


                    value += conj((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                  ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }

                else if ((w2)<0){
                    int S = fconv_n(-s,nw3);
                    int W1 = fconv_n(-w2,nw3);
                    int W2 = fconv_n(-w1,nw3);

                    int S_vert = fconv_n(-s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w2;
                    double w2_eff = -w1;
                    s = -s;


                    value += conj((K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                   +K3[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                  ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

                };};

        }
        else if(sym==2){
            site x(a,b,c);

            double w1_eff = w1, w2_eff = w2, s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if((w1_eff) < 0){w1_eff = -w1_eff; w2_eff = -w2_eff; x = site_switch(x.a,x.b,x.c);};
            s = s_eff;

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int S_vert = fconv_n(s,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);

            value = ((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))
                     / ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));



        };}
    else{
        int i,j,k;
        if(abs(s)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(s,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value += K3_vval(a,b,c,i,j,k);};};





    return value;
}

*/

/*****************FUNCTIONS FOR THE T-VERTEX********************************************/

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45

    double u=0.,w1_u=0., w2_u=0.;
    if(channel == 's'){
        u = -w2-w1;
        w1_u = (w1-w2+q)/2;
        w2_u = (-w1+w2+q)/2;}
    else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
        u = w2-w1;
        w1_u = (w1+w2-q)/2;
        w2_u = (w1+w2+q)/2;}
    else if(channel == 'u'){
        u = q;
        w1_u = w1;
        w2_u = w2;}
    else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
        u = w1-w2;
        w1_u = q + (w1-w2)/2;
        w2_u = (w1+w2)/2;};

    Q value;

//    if(abs(u) < freqs_a[nw1_wa/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
//    if(abs(w1_u) < freqs_a[nw1_wa/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
//    if(abs(w2_u) < freqs_a[nw1_wa/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK,u,w1_u,i_in) + K3_vvalsmooth(iK, u, w1_u, w2_u, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,u,w2_u,i_in)

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
    };
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

//overload of previous function         => I'm pretty sure we don't andwon't be needing this function
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
    vec_K1[iK*nw1_wt*n_in + i*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    vec_K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] = value;
}

template <typename Q> Q tvert<Q>::K1_vval(int iK, int i, int i_in){
    return vec_K1[iK*nw1_wt*n_in + i*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return vec_K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in];
}

template <typename Q> Q tvert<Q>::K1_vvalsmooth(int iK, double u, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    /*I assume a linear grid for this implementation. TODO: define a set of functions (one per type of grid to use) to convert double in range to index on grid*/
    /*TODO: since we might want to interpolate according to channel and diagrammatic class, rethink previous TODO */

    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw1_wa-1));     //nw1_wa because we're interpolating for K1 in channel a
    auto index = (int)((u-w_lower_b)/dw_b_a);

    double x1 = freqs_a[index];
    double x2 = freqs_a[index+1];

    Q f1 = vec_K1[iK*nw1_wt*n_in + (index)*n_in + i_in];
    Q f2 = vec_K1[iK*nw1_wt*n_in + (index+1)*n_in + i_in];

    return f1 + (u-x1)*(f2-f1)/(x2-x1);
}
template <typename Q> Q tvert<Q>::K2_vvalsmooth(int iK, double u, double w1, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw2_wa-1));     //nw2_wa because we're interpolating for bosonic freq in K2 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw2_nua-1));    //nw2_nua because we're interpolating for fermionic freq in K2 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];

    Q f11 = vec_K2[iK*nw2_wt*nw2_nut*n_in + (index_b)*nw2_nut*n_in + (index_f)*n_in + i_in];
    Q f12 = vec_K2[iK*nw2_wt*nw2_nut*n_in + (index_b)*nw2_nut*n_in + (index_f+1)*n_in + i_in];
    Q f21 = vec_K2[iK*nw2_wt*nw2_nut*n_in + (index_b+1)*nw2_nut*n_in + (index_f)*n_in + i_in];
    Q f22 = vec_K2[iK*nw2_wt*nw2_nut*n_in + (index_b+1)*nw2_nut*n_in + (index_f+1)*n_in + i_in];

    return (y2-w1)/(y2-y1)*((x2-u)/(x2-x1)*f11 + (u-x1)/(x2-x1)*f21) + (w1-y1)/(y2-y1)*((x2-u)/(x2-x1)*f12 + (u-x1)/(x2-x1)*f22);
}
template <typename Q> Q tvert<Q>::K3_vvalsmooth(int iK, double u, double w1, double w2, int i_in){

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw3_wa-1));     //nw3_wa because we're interpolating for bosonic freq in K3 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nua-1));    //nw3_nua because we're interpolating for fermionic freq in K3 in channel a
    double dwp_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nuap-1));  //nw3_nuap because we're interpolating for fermionic freq in K3 in channel a

    auto index_b = (int)((u-w_lower_b)/dw_b_a);
    auto index_f = (int)((w1-w_lower_f)/dw_f_a);
    auto index_fp = (int)((w2-w_lower_f)/dwp_f_a);

    double x1 = freqs_a[index_b];
    double x2 = freqs_a[index_b+1];
    double y1 = freqs_a[index_f];
    double y2 = freqs_a[index_f+1];
    double z1 = freqs_a[index_fp];
    double z2 = freqs_a[index_fp+1];

    Q f111 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b)*nw3_nut*nw3_nutp*n_in + (index_f)*nw3_nutp*n_in + (index_fp)*n_in + i_in];
    Q f112 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b)*nw3_nut*nw3_nutp*n_in + (index_f)*nw3_nutp*n_in + (index_fp+1)*n_in + i_in];
    Q f121 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b)*nw3_nut*nw3_nutp*n_in + (index_f+1)*nw3_nutp*n_in + (index_fp)*n_in + i_in];
    Q f122 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b)*nw3_nut*nw3_nutp*n_in + (index_f+1)*nw3_nutp*n_in + (index_fp+1)*n_in + i_in];
    Q f211 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b+1)*nw3_nut*nw3_nutp*n_in + (index_f)*nw3_nutp*n_in + (index_fp)*n_in + i_in];
    Q f212 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b+1)*nw3_nut*nw3_nutp*n_in + (index_f)*nw3_nutp*n_in + (index_fp+1)*n_in + i_in];
    Q f221 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b+1)*nw3_nut*nw3_nutp*n_in + (index_f+1)*nw3_nutp*n_in + (index_fp)*n_in + i_in];
    Q f222 = vec_K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + (index_b+1)*nw3_nut*nw3_nutp*n_in + (index_f+1)*nw3_nutp*n_in + (index_fp+1)*n_in + i_in];

    double xd = (u-x1)/(x2-x1);
    double yd = (w1-y1)/(y2-y1);
    double zd = (w2-z1)/(z2-z1);

    Q c00 = f111*(1-xd) + f211*xd;
    Q c01 = f112*(1-xd) + f212*xd;
    Q c10 = f121*(1-xd) + f221*xd;
    Q c11 = f122*(1-xd) + f222*xd;

    Q c0 = c00*(1-yd) + c10*yd;
    Q c1 = c01*(1-yd) + c11*yd;

    return c0*(1-zd) + c1*zd;
}
//non-member functions
template <typename Q> tvert<Q> operator*(double alpha, const tvert<Q>& vertex){
    tvert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wt; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wt*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wt*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wt; ++i){
                for(int j=0; j<nw2_nut; ++j){
                    vertex2.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wt; ++i){
                for(int j=0; j<nw3_nut; j++){
                    for(int k=0; k<nw3_nutp; ++k){
                        vertex2.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> tvert<Q> operator*(const tvert<Q>& vertex,double alpha){
    tvert<Q> vertex2;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wt; ++i){
                //K1 contributions
                vertex2.vec_K1[nk*nw1_wt*n_in + i*n_in + nin] = alpha * vertex.vec_K1[nk*nw1_wt*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wt; ++i){
                for(int j=0; j<nw2_nut; ++j){
                    vertex2.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin] = alpha*vertex.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wt; ++i){
                for(int j=0; j<nw3_nut; j++){
                    for(int k=0; k<nw3_nutp; ++k){
                        vertex2.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin] = alpha*vertex.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}
template <typename Q> tvert<Q> operator+(const tvert<Q>& vertex1,const tvert<Q>& vertex2){
    tvert<Q> vertex3;
#pragma omp for collapse(3)
    for(int nk=0; nk<nK_K1; ++nk){
        for(int nin=0; nin<n_in; ++nin){
            for(int i=0; i<nw1_wt; ++i){
                //K1 contributions
                vertex3.vec_K1[nk*nw1_wt*n_in + i*n_in + nin] = vertex1.vec_K1[nk*nw1_wt*n_in + i*n_in + nin] + vertex2.vec_K1[nk*nw1_wt*n_in + i*n_in + nin];
            }

            for(int i=0; i<nw2_wt; ++i){
                for(int j=0; j<nw2_nut; ++j){
                    vertex3.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin] = vertex1.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin] + vertex2.vec_K2[nk*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + nin];
                }
            }

            for(int i=0; i<nw3_wt; ++i){
                for(int j=0; j<nw3_nut; j++){
                    for(int k=0; k<nw3_nutp; ++k){
                        vertex3.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin] =
                                vertex1.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin] + vertex2.vec_K3[nk*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + nin];
                    }
                }
            }



        }

    }
    return vertex2;
}

/*

template <typename Q> Q tvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2;
            w2_t = (-w1-w2+q)/2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
        else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2;
            w2_t =(w1+w2+q)/2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2;
            w2_t = w1+(q-w2)/2;

        };


        Q value=0;

        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};



//        if(t > bfreqs[nw-1] && abs(t - bfreqs[nw-1]) < 1e-12){t = bfreqs[nw-1];}
//        else if(t<bfreqs[0] && abs(t - bfreqs[0]) < 1e-12){t = bfreqs[0];};


//        if(w1_t > ffreqs[nw-1] && abs(w1_t - ffreqs[nw-1]) < 1e-12){w1_t = ffreqs[nw-1];}
//        else if(w1_t<ffreqs[0] && abs(w1_t - ffreqs[0]) < 1e-12){w1_t = ffreqs[0];};


//        if(w2_t > ffreqs[nw-1] && abs(w2_t - ffreqs[nw-1]) < 1e-12){w2_t = ffreqs[nw-1];}
//        else if(w2_t<ffreqs[0] && abs(w2_t - ffreqs[0]) < 1e-12){w2_t = ffreqs[0];};


        value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t)+K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}


}
template <typename Q> Q tvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2;
            w2_t = (-w1-w2+q)/2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
        else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2;
            w2_t =(w1+w2+q)/2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2;
            w2_t = w1+(q-w2)/2;

        };


        Q value=0;

        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};



//        if(t > bfreqs[nw-1] && abs(t - bfreqs[nw-1]) < 1e-12){t = bfreqs[nw-1];}
//        else if(t<bfreqs[0] && abs(t - bfreqs[0]) < 1e-12){t = bfreqs[0];};


//        if(w1_t > ffreqs[nw-1] && abs(w1_t - ffreqs[nw-1]) < 1e-12){w1_t = ffreqs[nw-1];}
//        else if(w1_t<ffreqs[0] && abs(w1_t - ffreqs[0]) < 1e-12){w1_t = ffreqs[0];};


//        if(w2_t > ffreqs[nw-1] && abs(w2_t - ffreqs[nw-1]) < 1e-12){w2_t = ffreqs[nw-1];}
//        else if(w2_t<ffreqs[0] && abs(w2_t - ffreqs[0]) < 1e-12){w2_t = ffreqs[0];};

        if(p==1){
            if(channel=='t'){
                if(f == 'R' || f == 'M'){value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2b_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='u'){
                if(f == 'R' || f== 'M'){value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);};
            };
        }

        else if(p==2){
            if(channel=='t'){
                if(f == 'R' || f == 'L'){value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,t) + K2b_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='u'){
                if(f == 'R' || f== 'L'){value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);

                };
            };
        };



        return value;}
    else{return 0;}


}
//overload of previous function
template <typename Q> Q tvert<Q>::vvalsmooth(int red_side, int map,int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
    return vvalsmooth( a, b, c, q, w1,  w2,  channel,p,  f);

}
template <typename Q> Q tvert<Q>::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see K3euther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        t = q;
        w1_t = w1;
        w2_t = w2;
        Q value=0;


        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};

        value += K3_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}


}
template <typename Q> void tvert<Q>::K1_setvert(int a, int b, int c,int i, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>=0 ){

            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};

    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
template <typename Q> void tvert<Q>::K2_setvert(int a, int b, int c, int i, int j, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value ;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> void tvert<Q>::K3_setvert(int a, int b, int c, int i, int j, int k, Q value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){

            K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};
    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
template <typename Q> Q tvert<Q>::K1_vval(int a_raw, int b_raw,int c_raw,  int i){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw1 && i >= 0){
        if(sym==0){

            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 ||sym==2){
            if (i>=nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);

                i = nw1-1-i;
                value += K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)];};
        };};
    if(abs(value)<1e-100){value=0;};



    return value;
}
template <typename Q> Q tvert<Q>::K2_vval(int a_raw, int b_raw,int c_raw,  int i, int j){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp

    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    Q value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >=0){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j>=nw2/2){

                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}

            else if(j < nw2/2){


                j = nw2-1-j;
                value += conj(K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)]);
            };}
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;};
            if(j < nw2/2){j = nw2-1-j;};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };};
    if(abs(value)<1e-100){value=0;};


    return value;
}
template <typename Q> Q tvert<Q>::K3_vval(int a_raw, int b_raw,int c_raw,  int i, int j, int k){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0 ){

            value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}
        else if(sym == 1 ){

            if(i >= nw3/2){
                if(j >= nw3/2 ){

                    value += K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

                else if(i >= nw3/2  && j < nw3/2 ){

                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(K3[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)]);};}


            else if(i < nw3/2){
                if(k >= nw3/2 ){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)];}//note that k and j are interchanged

                else if( k < nw3/2 ){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);};//note that k and j are interchanged
            };}

        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;x = site_switch(x.a,x.b,x.c);};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;};
            value += K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};
    };
    if(abs(value)<1e-100){value=0;};



    return value;
}
template <typename Q> Q tvert<Q>::K1_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t){



    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(t)+ 1e-6< bfreqs[(nw1+nw)/2-1]){
        if(sym==0){

            int T = fconv_n(t,nw1);

            int T_vert = fconv_n(t,nw1)-(nw1-nw1_q);
            value = ((K1[a+(nuc_eff-1)/2][b][c-1][T_vert ]*(ffreqs[(nw-nw1)/2+T+1]-t)
                      +K1[a+(nuc_eff-1)/2][b][c-1][T_vert +1]*(-ffreqs[(nw-nw1)/2+T]+t))/((ffreqs[(nw-nw1)/2+T+1]-ffreqs[(nw-nw1)/2+T])));}
        else if(sym==1 || sym==2){
            if (t>0){
                int T = fconv_n(t,nw1);
                int T_vert = fconv_n(t,nw1)-(nw1-nw1_q);
                value = ((K1[a+(nuc_eff-1)/2][b][c-1][T_vert ]*(bfreqs[T+1]-t)
                          +K1[a+(nuc_eff-1)/2][b][c-1][T_vert +1]*(-bfreqs[T]+t))/((bfreqs[T+1]-bfreqs[T])));}
            else if (t<0){

                site x = site_switch(a,b,c);
                double t_eff = -t;

                int T = fconv_n(t_eff,nw1);
                int T_vert = fconv_n(t_eff,nw1)-(nw1-nw1_q);

                value = ((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert ]*(bfreqs[T+1]-t_eff)
                          +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert +1]*(-bfreqs[T]+t_eff))/(bfreqs[T+1]-bfreqs[T]));
            };};}
    else{ int i;
        if(abs(t)<= bfreqs[nw1-1] ){
            i= fconv_n(t,nw1);

            value = K1_vval(a,b,c,i);};
    };



    return value;
}
template <typename Q> Q tvert<Q>::K2_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w1){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(t)+ 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int T = fconv_n(t,nw2);
            int W1 = fconv_n(w1,nw2);

            int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                       +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                      ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){
            if(w1 > 0){

                int T = fconv_n(t,nw2);
                int W1 = fconv_n(w1,nw2);

                int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                           +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                          ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
            else if (w1<0){

                double w1_eff = - w1;
                int T = fconv_n(t,nw2);
                int W1 =  fconv_n(w1_eff,nw2);

                int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
                int W1_vert =  fconv_n(w1_eff,nw2)-(nw2-nw2_w1);



                value += conj((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                               +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                               +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                               +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                              ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));};}
        else if(sym==2){
            double w1_eff = w1;
            double t_eff = t;
            if(t<0){t_eff = -t; };
            if(w1_eff<0){w1_eff = -w1_eff;};

            int T = fconv_n(t_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int T_vert = fconv_n(t_eff,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);
            value = ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                      +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                     ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        };}
    else{ int i,j;
        if(abs(t)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(t,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};
    };
    if(abs(value)<1e-100){value=0;};



    return value;
}
template <typename Q> Q tvert<Q>::K3_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    Q value=0;
    if(abs(t)+ 1e-6< bfreqs[(nw3+nw)/2-1] && abs(w1)+ 1e-6<ffreqs[(nw3+nw)/2-1]&& abs(w2)+ 1e-6<ffreqs[(nw3+nw)/2-1]){
        if (sym==0){


            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);

            int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

            value = ((K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){
            if(t > 0){
                if((w1)>0 ){

                    int T = fconv_n(t,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);

                    int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

                    value = ((K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                              +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                             ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }

                else if( (w1)<0){

                    int T = fconv_n(t,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);

                    int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;
                    value = conj((K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                 ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}


            else if (t < 0  ){
                if((w2)>0 ){
                    site x = site_switch(a,b,c);


                    int T = fconv_n(-t,nw3);
                    int W1 = fconv_n(w2,nw3);
                    int W2 = fconv_n(w1,nw3);

                    int T_vert = fconv_n(-t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;

                    t = -t;

                    value = (K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                             +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2]));
                }



                else if ((w2)<0 ){
                    site x = site_switch(a,b,c);


                    int T = fconv_n(-t,nw3);
                    int W1 = fconv_n(-w2,nw3);
                    int W2 = fconv_n(-w1,nw3);

                    int T_vert = fconv_n(-t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1,nw3)-(nw3-nw3_w2);

                    double w1_eff = -w2;
                    double w2_eff = -w1;

                    t = -t;

                    value = conj((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                                  +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                                 ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };};}
        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, t_eff =t;
            if(t<0){t_eff = -t; };
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1; x = site_switch(x.a,x.b,x.c);};
            if((w1_eff)<0){w1_eff = -w1_eff; w2_eff = -w2_eff;};
            t = t_eff;


            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);

            value = ((K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

        };

    }
    else{
        int i,j,k;
        if(abs(t)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(t,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value += K3_vval(a,b,c,i,j,k);};};
    if(abs(value)<1e-100){value=0;};



    return value;
}

*/

/*****************FUNCTIONS FOR THE IRREDUCIBLE VERTEX********************************************/
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
//operators for irreducible vertex
template <typename Q> irreducible<Q> operator*(double alpha, const irreducible<Q>& vertex) {
    irreducible<Q> result;
    result.U_bare = alpha * vertex.U_bare;
    return result;

}
template <typename Q> irreducible<Q> operator*(const irreducible<Q>& vertex,double alpha) {
    irreducible<Q> result;
    result.U_bare = alpha * vertex.U_bare;
    return result;
}
template <typename Q> irreducible<Q> operator+(const irreducible<Q>& vertex1,const irreducible<Q>& vertex2) {
    irreducible<Q> result;
    result.U_bare = vertex1.U_bare + vertex2.U_bare;
    return result;
}

/*

template <typename Q> Q irreducible<Q>::vval(int a_raw, int b_raw, int c_raw){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
template <typename Q> Q irreducible<Q>::vvalsmooth(int a_raw, int b_raw, int c_raw){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        //cout << a << " " << b << " " << c << endl;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
template <typename Q> Q irreducible<Q>::vvalsmooth(int a_raw, int b_raw, int c_raw,double q, double w1, double w2, char channel, int par, char f){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
template <typename Q> void irreducible<Q>::setvert(int a_raw, int b_raw, int c_raw,  Q value){
    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    U_bare[a+(nuc_eff-1)/2][b][c-1] = value;}

*/


/*****************************************operators concerning parvert objects********************************************************/

template <typename Q> parvert<avert<Q> > operator+(parvert<avert<Q> >  vertex1,parvert<avert<Q> >  vertex2){
    parvert<avert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<pvert<Q> >  operator+(parvert<pvert<Q> >  vertex1,parvert<pvert<Q> >  vertex2){
    parvert<pvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<tvert<Q> >  operator+(parvert<tvert<Q> >  vertex1,parvert<tvert<Q> >  vertex2){
    parvert<tvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<irreducible<Q> > operator+(parvert<irreducible<Q> > vertex1,parvert<irreducible<Q> > vertex2){
    parvert<irreducible<Q> > result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<avert<Q> >  operator+=(parvert<avert<Q> >  vertex1,parvert<avert<Q> >  vertex2){
    parvert<avert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<pvert<Q> >  operator+=(parvert<pvert<Q> >  vertex1,parvert<pvert<Q> >  vertex2){
    parvert<pvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<tvert<Q> >  operator+=(parvert<tvert<Q> >  vertex1,parvert<tvert<Q> >  vertex2){
    parvert<tvert<Q> >  result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<irreducible<Q> > operator+=(parvert<irreducible<Q> > vertex1,parvert<irreducible<Q> > vertex2){
    parvert<irreducible<Q> > result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
template <typename Q> parvert<avert<Q> >  operator*(double alpha ,parvert<avert<Q> > & vertex){
    parvert<avert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<avert<Q> >  operator*(parvert<avert<Q> > & vertex,double alpha){
    parvert<avert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<pvert<Q> >  operator*(double alpha ,parvert<pvert<Q> > & vertex){
    parvert<pvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<pvert<Q> >  operator*(parvert<pvert<Q> > & vertex,double alpha){
    parvert<pvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<tvert<Q> >  operator*(double alpha ,parvert<tvert<Q> > & vertex){
    parvert<tvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<tvert<Q> >  operator*(parvert<tvert<Q> > & vertex,double alpha){
    parvert<tvert<Q> >  result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<irreducible<Q> > operator*(double alpha ,parvert<irreducible<Q> >& vertex){
    parvert<irreducible<Q> > result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
template <typename Q> parvert<irreducible<Q> > operator*(parvert<irreducible<Q> >& vertex,double alpha){
    parvert<irreducible<Q> > result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
/*************************************************************************************************************/






/*****************************************FUNCTIONS FOR FULL VERTEX "FULLVERT"********************************************************/
//arguments are equivalent to those in the simple vertex functions

template <typename Q> Q fullvert<Q>::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel){
    Q result = irred.vvalsmooth(a,b,c) + pvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel);
    if(abs(result)>1e-16){
        return  result;}
    else{return 0;};
}
template <typename Q> Q fullvert<Q>::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){
    Q result=0;

    if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
        result += irred.vvalsmooth(a,b,c);

    }
    else if( p==2 && (f=='K' || f== 'M')){
        result += irred.vvalsmooth(a,b,c);

    };



    result +=  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
    //if(p==2 && f=='L'){cout <<  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << " " << tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f)  << " " <<  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << endl;};
    if(abs(result)<1e-16){
        result  =0;};
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
            else if(channel == 'u'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);};}
        else if(map==1){
            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
            else if(channel == 'u'){  result =  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
        };

    };

    if(abs(result)<1e-16){
        result  =0;};
    return result;
}
template <typename Q> fullvert<Q> operator*(double alpha, const fullvert<Q>& vertex){
    fullvert<Q> result;
    result.irred = alpha * vertex.irred;
    result.pvertex = alpha *vertex.pvertex;
    result.tvertex = alpha * vertex.tvertex;
    result.avertex = alpha * vertex.avertex;
    return result;
}
template <typename Q> fullvert<Q> operator+( const fullvert<Q>& vertex1, const fullvert<Q>& vertex2){
    fullvert<Q> result;
    result.irred = vertex1.irred + vertex2.irred ;
    result.pvertex = vertex1.pvertex + vertex2.pvertex ;
    result.tvertex = vertex1.tvertex + vertex2.tvertex ;
    result.avertex = vertex1.avertex + vertex2.avertex ;
    return result;
}

/*

template <typename Q> Q fullvert<Q>::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel){
    Q result = irred.vvalsmooth(a,b,c) + pvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel);
    if(abs(result)>1e-16){
        return  result;}
    else{return 0;};
}
template <typename Q> Q fullvert<Q>::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){
    Q result=0;

    if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
        result += irred.vvalsmooth(a,b,c);

    }
    else if( p==2 && (f=='K' || f== 'M')){
        result += irred.vvalsmooth(a,b,c);

    };



    result +=  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
    //if(p==2 && f=='L'){cout <<  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << " " << tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f)  << " " <<  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << endl;};
    if(abs(result)<1e-16){
        result  =0;};
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
            else if(channel == 'u'){  result =  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);};}
        else if(map==1){
            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
            else if(channel == 'u'){  result =  avertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  pvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
        };

    };

    if(abs(result)<1e-16){
        result  =0;};
    return result;
}

*/