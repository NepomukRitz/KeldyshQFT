//
// Created by Sa.Aguirre on 6/27/19.
//

#ifndef KELDYSH_SELFENERGY_SOPT_COMPARING_PROPAGATOR_H
#define KELDYSH_SELFENERGY_SOPT_COMPARING_PROPAGATOR_H

#include <complex>
#include<vector>

using namespace std;
typedef complex<double> comp;
typedef vector<comp> cvec;
typedef vector<cvec> cmat;


class Propagator : public cvec
{
public:
    cvec vec;
    explicit Propagator (int size)
            :vec(cvec(size)) {}

    Propagator operator* (double dt){
        for (int i = 0; i<vec.size(); i++)
            vec[i] *= dt;
        return *this;
    }

    Propagator operator* (comp a){
        for (int i = 0; i<vec.size(); i++)
            vec[i] *= a;
        return *this;
    }

    Propagator operator* (Propagator b){
        for(int i=0; i<vec.size(); i++){
            vec[i] *= b[i];
        }
        return *this;
    }

    Propagator operator+ (Propagator& b)
    {
        for (int i=0; i<vec.size(); i++){
            vec[i] += b[i];
        }
        return *this;
    }

    Propagator operator- ()
    {
        for (int i=0; i<vec.size(); i++){
            vec[i] *=-1.;
        }
        return *this;
    }

    Propagator operator- (Propagator& b)
    {
        for (int i=0; i<vec.size(); i++){
            vec[i] -= b[i];
        }
        return *this;
    }

    comp& operator[] (int i)
    {
        return vec[i];
    }

//    Propagator operator= (Propagator rhs){
//        for(int i=0;i<vec.size();i++)
//        {
//            vec[i] = rhs[i];
//        }
//        return *this;
//    }

    int size(){
        return vec.size();
    }

    void setToZero(){
        for (int i=0; i<vec.size(); i++)
            vec[i] = 0.;
    }

    Propagator conjugate()
    {
        for(int i=0; i<vec.size(); i++)
            vec[i] = vec[i].real() - (comp)1.i*vec[i].imag();
        return *this;
    }
};

#endif //KELDYSH_SELFENERGY_SOPT_COMPARING_PROPAGATOR_H
