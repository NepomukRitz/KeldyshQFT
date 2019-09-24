//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include "selfenergy.h"
#include "vertex.h"
#include "propagator.h"

//TODO One would actually want to have integrands to feed an integrand, not vectors
SelfEnergy<comp> loop(Vertex<fullvert<comp> >& fullvertex, Propagator& prop)
{
    cvec integrandR(nPROP);
    cvec integrandK(nPROP);
    SelfEnergy<comp> resp = SelfEnergy<comp> ();
    for (int i = 0; i<nPROP; i++){
        for(int j = 0; j<nPROP; j++){
            double w = ffreqs[i];
            double wp= ffreqs[j];
            comp pR = prop.pval(0, j);
            comp pA = conj(pR);
            comp pK = prop.pval(1, j);

            integrandR[j] = -(fullvertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pR +    //In this formulas, the frequencies
                              fullvertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pA +    //have already been converted to the
                              fullvertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pK +    //respective channel convention!

                              fullvertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
                              fullvertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
                              fullvertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +

                              fullvertex.densvertex.tvertex.value(3, 0., wp, w, 0, fullvertex.densvertex.avertex)*pR +
                              fullvertex.densvertex.tvertex.value(6, 0., wp, w, 0, fullvertex.densvertex.avertex)*pA +
                              fullvertex.densvertex.tvertex.value(7, 0., wp, w, 0, fullvertex.densvertex.avertex)*pK +

                              fullvertex.densvertex.irred.vval(3)*pR +
                              fullvertex.densvertex.irred.vval(6)*pA +
                              fullvertex.densvertex.irred.vval(7)*pK);


            integrandK[j] = -(fullvertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pR +
                              fullvertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pA +
                              fullvertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, fullvertex.densvertex.tvertex)*pK +

                              fullvertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
                              fullvertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
                              fullvertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +

                              fullvertex.densvertex.tvertex.value(1, 0., wp, w, 0, fullvertex.densvertex.avertex)*pR +
                              fullvertex.densvertex.tvertex.value(4, 0., wp, w, 0, fullvertex.densvertex.avertex)*pA +
                              fullvertex.densvertex.tvertex.value(5, 0., wp, w, 0, fullvertex.densvertex.avertex)*pK +

                              fullvertex.densvertex.irred.vval(1)*pR +
                              fullvertex.densvertex.irred.vval(4)*pA +
                              fullvertex.densvertex.irred.vval(5)*pK);

            //TODO include spinvertex contributions!!

        }

        comp integratedR = integrator(integrandR);
        comp integratedK = integrator(integrandK);


        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);

    }

    return resp;
}


#endif //KELDYSH_MFRG_LOOP_H
