//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <iostream>

#include "selfenergy.h"
#include "data_structures.h"
#include "parameters.h"

using namespace std;

//TODO: potentially this could also be a template type (don't know if it's necessary though)
//TODO: add Keldysh component!!

/*******PROPAGATOR FUNCTION***********/
comp propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type);


/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/
comp propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type) {

    if(w!=0){
        complex<double> value;
        complex<double> iw(0.0,w);
        if(type == 'g'){//regular undifferentiated greensfunction
            complex<double> g0(0.0,0.0);
            if(reg==1){

                if(fabs(w) > Lambda){
                    g0 =1./iw;
                    value += 1./(1./g0-selfenergy.svalsmooth(w));

                }
                else if (fabs(w) == Lambda){
                    g0 = 1./iw;
                    value += 0.5/(1./g0-selfenergy.svalsmooth(w));
                    //this is the implementation of Morris Lemma, see SB.II, p. 9
                }
            }
            else if(reg==2){

                g0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw;

                value = 1./(1./g0-selfenergy.svalsmooth(w));
            };

        }
        else if(type == 's'){//single scale propagator

            if(reg==1){
                if(fabs(w) == Lambda){

                    value += (-1./(iw-selfenergy.svalsmooth(w)));}
                else{return 0.;};
            }
            else if(reg==2){

                complex<double> G0 = (1.-exp(-pow(fabs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += pow((1. + G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * 1./iw;//see page 20 in SB II

            };
        }
        else if(type == 'k'){//katanin substitution
            if(reg==1){
                complex<double> S,G;
                if(fabs(w) == Lambda){ S = - 1./(iw-selfenergy.svalsmooth(w));};//single scale propagator
                if(fabs(w) > Lambda){ G = 1./(iw-selfenergy.svalsmooth(w));}
                else if(fabs(w) == Lambda) {G = 0.5/(iw-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9
                value = (S + G * G * diffselfenergy.svalsmooth(w));
            }
            else if(reg==2){
                complex<double> G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += pow((1. + G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * 1./iw
                         + G * G * diffselfenergy.svalsmooth(w);//see page 20 in SB II

            };

        }
        else if(type == 'e'){//only self energy extension of katanin substitution ( needed for sharp regulator)
            if(reg==1){
                complex<double> G;
                if(fabs(w) > Lambda){ G = 1./(iw-selfenergy.svalsmooth(w));}
                else if(fabs(w) == Lambda) {G = 0.5/(iw-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9
                value += (G * G * diffselfenergy.svalsmooth(w));
            }
            else if (reg==2){
                complex<double> G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += G * G * diffselfenergy.svalsmooth(w);
            };

        }
        else{cout << "could not resolve propagator type" << endl;};

        return value;
    }
    else{return 0;};}



#endif //KELDYSH_MFRG_PROPAGATOR_H
