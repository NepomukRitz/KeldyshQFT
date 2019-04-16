#ifndef X_BUBBLE_RPA_ZERO_MAG_15032017
#define X_BUBBLE_RPA_ZERO_MAG_15032017

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Syma_Matrix.h"

template <int mode> class Integrand_X_bubble_rpa_zero_mag{
 public:
 double external_freq;
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 double measure_flow;
 int l, k; // l muss stets groesser als k sein, wenn INTEGRAND_WITHOUT_SYMA 0!
 int jmin, jmax, imin, imax; 
 syma<complex<double> > Gu_diff_syma, Gu_sum_syma, Su_syma; 

#if INTEGRAND_WITHOUT_SYMA
 Syma_Matrix<complex<double> > trafo;
 matrix<complex<double> > Gu_diff_matrix, Gu_sum_matrix, Su_matrix;
#endif

 Integrand_X_bubble_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l, int k);
 matrix<complex<double> > operator()(double internal);
 matrix<double> select(matrix<complex<double> > &M);
};









#endif
