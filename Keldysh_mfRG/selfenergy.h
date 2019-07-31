//
// Created by E.Walter on 7/31/19.
//

#ifndef KELDYSH_MFRG_SELFENERGY_H
#define KELDYSH_MFRG_SELFENERGY_H

/******************CLASS FOR SELF ENERGY *************/
template <typename Q>
class self{
    vec<Q> selfenergy =  vec<Q> (nSE);
public:
    void setself(int, Q);
    Q sval(int);
    Q svalsmooth(double);
    friend self operator+(const self& self1, const self& self2);
    friend self operator+=(const self& self1,const self& self2);
    friend self operator*(Q alpha, const self& self1);
    friend self operator*(const self& self1, Q alpha);
};



/*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
template <typename Q>
Q self<Q>::sval(int i){
    return selfenergy[i];
}

template <typename Q>
Q self<Q>::svalsmooth(double w){//smoothly interpolates for values between discrete frequency values of mesh
    Q value;
    int W = fconv(w); // TODO: define fconv
    value += ((selfenergy[W]*(ffreqs[W+1]-w)+selfenergy[W+1]*(-ffreqs[W]+w))/(ffreqs[W+1]-ffreqs[W])); // TODO: make it similar to the version for vertex
    return value;
}

template <typename Q>
void self<Q>::setself(int i, Q val){
    selfenergy[i] = val;
}

//operators for self energy
self operator*(double alpha, const self& self1){//product operator overloading
    self self2;
    self2.selfenergy = self1.selfenergy * alpha;
    return self2;
}
self operator*(const self& self1, double alpha){//product operator overloading
    self self2;
    self2.selfenergy = self1.selfenergy * alpha;
    return self2;
}
self operator+(const self& self1, const self& self2){//sum operator overloading
    self self3;
    self3.selfenergy = self1.selfenergy + self2.selfenergy;
    return self3;
}
self operator+=(const self& self1, const self& self2){//sum operator overloading
    self self3;
    self3.selfenergy = self1.selfenergy + self2.selfenergy;
    return self3;
}


//TODO: check this below (and define state first)
//operators containing state objects
state operator+(state state1, state state2){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred + state2.vertex.spinvertex.irred;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex + state2.vertex.spinvertex.svertex;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex + state2.vertex.spinvertex.tvertex;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex + state2.vertex.spinvertex.uvertex;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred + state2.vertex.densvertex.irred;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex + state2.vertex.densvertex.svertex;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex + state2.vertex.densvertex.tvertex;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex + state2.vertex.densvertex.uvertex;
    result.selfenergy = state1.selfenergy + state2.selfenergy;
    return result;
}
state operator*(double alpha, state state1){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex * alpha;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex * alpha;

    result.selfenergy = alpha * state1.selfenergy;
    return result;
}
state operator*(state state1, double alpha){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex * alpha;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex * alpha;

    result.selfenergy = alpha * state1.selfenergy;
    return result;
}



#endif //KELDYSH_MFRG_SELFENERGY_H
