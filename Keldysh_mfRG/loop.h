//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include "selfenergy.h"
#include "vertex.h"
#include "propagator.h"


SelfEnergy<comp> loop(Vertex<fullvert<comp> >& fullvertex, Propagator& prop)
{
    cvec integrandR(nPROP);
    cvec integrandR2(nPROP);
    cvec integrandK(nPROP);
    cvec integrandK2(nPROP);
    SelfEnergy<comp> resp = SelfEnergy<comp> ();
    for (int i = 0; i<nPROP; ++i){
        for(int j = 0; j<nPROP; ++j){
            double w = ffreqs[i];
            double wp= ffreqs[j];
            comp pR = prop.pvalsmooth(0, wp);
            comp pA = conj(pR);
            comp pK = prop.pvalsmooth(1, wp);

            integrandR[j] = -(fullvertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pR +    //In this formulas, the frequencies
                             fullvertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pA +    //have already been converted to the
                             fullvertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pK +    //respective channel convention!

                             fullvertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
                             fullvertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
                             fullvertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +

                             fullvertex.densvertex.tvertex.value(3, 0., wp, w, 0)*pR +
                             fullvertex.densvertex.tvertex.value(6, 0., wp, w, 0)*pA +
                             fullvertex.densvertex.tvertex.value(7, 0., wp, w, 0)*pK +

                             fullvertex.densvertex.irred.vval(3)*pR +
                             fullvertex.densvertex.irred.vval(6)*pA +
                             fullvertex.densvertex.irred.vval(7)*pK);

            integrandR2[j] = -0.5*U*SK(1.0, wp, 0.5, 0.0, 0.5);

//TODO finish debugging and get rid of debugging stuff
            comp var1 = fullvertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pR;
            comp var2 = fullvertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pA;
            comp var3 = fullvertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pK;

            comp var4 = fullvertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR;
            comp var5 = fullvertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA;
            comp var6 = fullvertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK;

            comp var7 = fullvertex.densvertex.tvertex.value(1, 0., wp, w, 0)*pR;
            comp var8 = fullvertex.densvertex.tvertex.value(4, 0., wp, w, 0)*pA;
            comp var9 = fullvertex.densvertex.tvertex.value(5, 0., wp, w, 0)*pK;

            comp var10 = fullvertex.densvertex.irred.vval(1)*pR;
            comp var11 = fullvertex.densvertex.irred.vval(4)*pA;
            comp var12 = fullvertex.densvertex.irred.vval(5)*pK;


            integrandK[j] = -( var1 + var2+ var3+ var4+ var5+ var6+ var7+ var8+ var9+ var10+ var11+ var12);

            integrandK2[j] = -0.5*U*(SR(1.0, wp, 0.5) + conj(SR(1.0, wp, 0.5)));

//            integrandK[j] =  (  fullvertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pR +
//                                fullvertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pA +
//                                fullvertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pK +
//
//                                fullvertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
//                                fullvertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
//                                fullvertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +
//
//                                fullvertex.densvertex.tvertex.value(1, 0., wp, w, 0)*pR +
//                                fullvertex.densvertex.tvertex.value(4, 0., wp, w, 0)*pA +
//                                fullvertex.densvertex.tvertex.value(5, 0., wp, w, 0)*pK +
//
//                                fullvertex.densvertex.irred.vval(1)*pR +
//                                fullvertex.densvertex.irred.vval(4)*pA +
//                                fullvertex.densvertex.irred.vval(5)*pK);
//
            //TODO include spinvertex contributions!!

        }

        comp integratedR = integrator(integrandR);
        comp integratedK = integrator(integrandK);
        comp integratedR2 = integrator(integrandR2);
        comp integratedK2 = integrator(integrandK2);


        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);

    }

    return resp;
}


#endif //KELDYSH_MFRG_LOOP_H
