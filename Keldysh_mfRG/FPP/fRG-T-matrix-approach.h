//
// Created by Marcel on 11.08.2021.
//

#ifndef MAIN_CPP_FRG_T_MATRIX_APPROACH_H
#define MAIN_CPP_FRG_T_MATRIX_APPROACH_H

#include "Ladder-approximation.h"
#include "Momentum-integral-Bubble.h"
#include "../util.h"
#include "../solvers.h"

comp rhs_test3(const comp& y, double Lambda) {
    //comp y;
    return SimpleBubble(0.03, -2.0, 0.2, 0.4,Lambda, 'c', 'c');
}



/*
template <typename Q>
class rhs_1lfRG {
private:
    double w, q;
    const double Lambda_i, Lambda_f;
    char i, j, r;

public:
    /**
     * Constructor:
     */ /*
    rhs_1lfRG(double w_in, double q_in, const double Lambda_i_in, const double Lambda_f_in, char i_in, char j_in,
              char r_in)
            : w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), i(i_in), j(j_in), r(r_in) {
    };

    auto operator() (double Lambda) const -> Q {
        return  sharp_frequency_exact_bare_bubble (w, Lambda, q, i, j, r);
    };
};

template <typename Q> auto set_rhs_1lfRG(const rhs_1lfRG<comp>& rhs1LfRG_p, double Lambda) -> rhs_1lfRG<comp>{
     flow_rhs(Lambda);
}; */

comp rhs_test2(const comp& y, double Lambda) {
    //comp y;
    return 2.*pow(-gint(1e4,1e-10,1)+y,2)*sharp_frequency_exact_bare_bubble(0.0,Lambda,0.0,'c','d','p');
}

/*
comp solve_1lfRG_nsc {
    //rhs_1lfRG<comp> rhs_p(w,q,Lambda_i,Lambda_f,'c','d','p');
    comp y_fin;
    const comp y_ini = 0.0;
    ODE_solver_RK4(y_fin, 1e4, 0.0, 1e-10, rhs_test2, sq_substitution, sq_resubstitution, nODE);
    return y_fin;
};

*/

#endif //MAIN_CPP_FRG_T_MATRIX_APPROACH_H
