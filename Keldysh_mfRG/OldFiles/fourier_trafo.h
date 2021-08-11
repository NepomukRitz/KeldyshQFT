#ifndef FOURIER_TRAFO_H
#define FOURIER_TRAFO_H

#include <cmath>                    // exp (to adjust Fourier convention), M_PI = 3.1415...
#include <fftw3.h>                  // Fast Fourier transform library
#include "../data_structures.h"        // real/complex vector classes, imag. unit
#include "../utilities/write_data2file.h"        // writing data into text or hdf5 files
#include <iostream>                 // text input/output
#include "../propagator.h"             // propagator to perform second-order perturbation theory (SOPT)
#include "../selfenergy.h"             // self-energy filled in SOPT
#include "../r_vertex.h"               // reducible vertex in channel r (for K1a, K1p in SOPT)
#include "../utilities/util.h"                   // sign function

using namespace std;

double Theta(double x) { // Heaviside step function for analytical Fourier transform
    if (x > 0) return 1.;
    else if (x < 0) return 0.;
    else return 1./2.;
}



// auxiliary function to flip and adjust an array
void flip_adjust(cvec& Gout, const cvec& Gin) {
    const int N = Gin.size(); // vector size
    Gout[0] = Gin[0]; // adjust boundary values
    for(int i=1; i<N; ++i) // flip array with shift by one
        Gout[i] = Gin[N-i];
}

// linear interpolation from equidistant grid as used in FFT
// function given as G(x_i) where x_i = x_0 + i * dx, output G(y_i)
void interp1_FFT(cvec& Gy, const rvec& y, const cvec& Gx, const rvec& x) {
    // vector sizes
    const int Nx = Gx.size();
    const int Ny = Gy.size();
    if(Nx != x.size() || Ny != y.size())
        cout << "Error in interp1_FFT: Vector sizes must be equal. "
                "Got " << Nx << ", " << x.size() << "; " << Ny << ", " << y.size() << "." << endl;
    // properties of x grid
    const double x0 = x[0];
    const double x_max = x[Nx-1];
    const double dx = x[1] - x[0];
    double temp_diff; // temporary difference
    for (int i = 1; i < Nx; ++i) {
        temp_diff = x[i] - x[i-1];
        if(abs(temp_diff - dx) > 1e-12) {
            cout << "Error in interp1_FFT: x grid must be linear. Deviation in difference is: " << (temp_diff-dx) << "." << endl;
            break;
        }
    }
    // iterate over all y points
    for(int i=0; i<Ny; ++i) {
        if(y[i]<x0 || y[i]>x_max) // outside grid: 0
            Gy[i] = 0.;
        else if(y[i]==x_max) // end point of grid
            Gy[i] = Gx[Nx-1];
        else { // in between
            int fl_y_idx = (int)((y[i] - x0)/dx); // get integer corresponding to y rounded down
            double x1 = x[fl_y_idx]; // lower x value
            double x2 = x[fl_y_idx+1]; // upper x value
            double xd = (y[i] - x1) / (x2 - x1); // distance between grid points
            Gy[i] = (1. - xd) * Gx[fl_y_idx] + xd * Gx[fl_y_idx+1]; // linearly interpolated value
        }
    }
}
void test_interp1_FFT() {
    rvec x{0., 1., 2.};
    cvec fx{0., 1., 4.};
    rvec y{-0.5, 0.5, 1.5, 4};
    cvec fy(4);
    interp1_FFT(fy, y, fx, x);
    cout << "Testing interp1_FFT. Differences from exact results are: " << fy[0] - 0. << " " << fy[1] - 0.5 << " " << fy[2] - 2.5 << " "
         << fy[3] - 0. << "." << endl;
}

// Fourier trafo from real frequency (v, Gv) to real time (t, Gt) using fftw
// G(t) = G(t0 + nt Dt) = 1/2pi \int dv e^{-ivt} G(v) = Dv/2pi \sum_{nv} e^{-i(v0 + nv Dv)(t0 + nt Dt) G(v0 + nv Dv)
//      = e^{-i v0 t} 1/T \sum_{nv} e^{-2pi i nv nt} e^{-i (v-v0) t0} G(v) = e^{-i v0 t} 1/T DFT(e^{-i (v-v0) t0} G(v))
// Dt = T/N, Dv = V/N = 2pi/T
// hyb encodes analytic high-frequency decay tail_coef/(v+i tail_hyb), setting coef=0 ignores this
void ft_v2t(rvec& t, cvec& Gt, const rvec& v, const cvec& Gv, const comp tail_coef, const double tail_hyb) {
    const int N = v.size(); // vector size, all vector sizes must be equal
    if(N != t.size() || N != Gt.size() || N != Gv.size())
        cout << "Error in ft_v2t: Vector sizes must be equal." << endl;
    const double Vtot = (v[1]-v[0])*(double)N; // frequency span
    const double Ttot = (2.*M_PI*(double)N)/Vtot; // time span
    for(int i=0; i<N; ++i) // fill time grid
        t[i] = -Ttot/2. + (Ttot/(double)N)*(double)i;
    fftw_complex in[N], out[N]; // input and output array in fftw format
    comp temp; // temporary variable to adjust Fourier convention (see above)
    if(abs(tail_coef)>1e-10) { // include high-frequency decay tail_coef/(v+i*tail_hyb) explicitly
        for (int i = 0; i < N; ++i) { // fill input array
            temp = ( Gv[i] - tail_coef / (v[i] + glb_i * tail_hyb) )
                   * exp(-glb_i * (v[i] - v[0]) * t[0]); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    else { // standard procedure
        for(int i=0; i<N; ++i) { // fill input array
            temp = Gv[i] * exp(-glb_i*(v[i]-v[0])*t[0]); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    fftw_plan p; // fftw forward Fourier transform, result saved in 'out'
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    if(abs(tail_coef)>1e-10) { // FT tail_coef/(v+i*tail_hyb) = -i tail_coef * sign(tail_hyb) Theta(sign(tail_hyb)*t) e^{-tail_hyb*t}
        for (int i = 0; i < N; ++i) // fill output
            Gt[i] = ( -glb_i * tail_coef * sign(tail_hyb) * Theta(tail_hyb*t[i] + 1e-16) * exp(-abs(tail_hyb*t[i])) ) +
                    comp(out[i][0], out[i][1]) * exp(-glb_i * v[0] * t[i]) / Ttot; // see Fourier convention above
    }
    else { // standard procedure
        for (int i = 0; i < N; ++i) // fill output
            Gt[i] = comp(out[i][0], out[i][1]) * exp(-glb_i * v[0] * t[i]) / Ttot; // see Fourier convention above
    }
    fftw_cleanup(); // clear memory from fftw variables
}
void ft_v2t(rvec& t, cvec& Gt, const rvec& v, const cvec& Gv) {
    ft_v2t(t, Gt, v, Gv, 0., 0.);
}

// Fourier trafo from real time (t, Gt) to real frequency (v, Gv) using fftw
// G(v) = G(v0 + nv Dv) = \int dt e^{ivt} G(t) = Dt \sum_{nt} e^{i(v0 + nv Dv)(t0 + nt Dt) G(t0 + nt Dt)
//      = e^{i v t0} T/N \sum_{nt} e^{2pi i nv nt} e^{i v0 (t-t0)} G(t) = e^{i v0 t} T/N IDFT(e^{i v0 (t-t0)} G(t))
// Dt = T/N, Dv = V/N = 2pi/T
// hyb encodes analytic high-frequency decay tail_coef/(v+i tail_hyb), setting coef=0 ignores this
void ft_t2v(rvec& v, cvec& Gv, const rvec& t, const cvec& Gt, const comp tail_coef, const double tail_hyb) {
    const int N = t.size(); // vector size, all vector sizes must be equal
    if(N != v.size() || N != Gv.size() || N != Gt.size())
        cout << "Error in ft_t2v: Vector sizes must be equal." << endl;
    const double Ttot = (t[1]-t[0])*(double)N; // time span
    const double Vtot = (2.*M_PI*(double)N)/Ttot; // frequency span
    for(int i=0; i<N; ++i) // fill frequency grid
        v[i] = -Vtot/2. + (Vtot/(double)N)*(double)i;
    fftw_complex in[N], out[N]; // input and output array in fftw format
    comp temp; // temporary variable to adjust Fourier convention (see above)
    if(abs(tail_coef)>1e-10) { // include long-time decay -i sgn(tail_hyb) tail_coef \Theta(sgn(tail_hyb)*t) e^{-tail_hyb t} explicitly, here enforce Theta(0)=1 via +1e-16
        for (int i = 0; i < N; ++i) { // fill input array
            temp = ( Gt[i] + glb_i * tail_coef * sign(tail_hyb) * Theta(tail_hyb*t[i] + 1e-16) * exp(-abs(tail_hyb*t[i])) )
                   * exp(glb_i * v[0] * (t[i] - t[0])); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    else { // standard procedure
        for (int i = 0; i < N; ++i) { // fill input array
            temp = Gt[i] * exp(glb_i * v[0] * (t[i] - t[0])); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    fftw_plan p; // fftw forward Fourier transform, result saved in 'out'
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    if(abs(tail_coef)>1e-10) { // FT tail_coef/(v+i*tail_hyb) = -i tail_coef * sign(tail_hyb) Theta(sign(tail_hyb)*t) e^{-tail_hyb*t}
        for (int i = 0; i < N; ++i) // fill output
            Gv[i] = tail_coef / (v[i] + glb_i * tail_hyb) +
                    comp(out[i][0], out[i][1]) * exp(glb_i * v[i] * t[0]) * (Ttot / (double) N); // see Fourier convention above
    }
    else { // standard procedure
        for (int i = 0; i < N; ++i) // fill output
            Gv[i] = comp(out[i][0], out[i][1]) * exp(glb_i * v[i] * t[0]) * (Ttot / (double) N); // see Fourier convention above
    }
    fftw_cleanup(); // clear memory from fftw variables
}
void ft_t2v(rvec& v, cvec& Gv, const rvec& t, const cvec& Gt) {
    ft_t2v(v, Gv, t, Gt, 0., 0.);
}

// Fourier trafo from imaginary frequency (vn, Gvn) to imaginary time (tau, Gtau) using fftw
// possibly including the tail analytically (incl_tail=true): Fourier trafo of 1/(iv) is (-1/2 + Theta(-tau))
// G(tau) = G(tau0 + ntau Dtau) = 1/beta \sum_{nv} e^{-i vn t} G(vn)
//        = 1/beta \sum_{nv} e^{-i(v0 + nv Dv)(tau0 + ntau Dtau) G(v0 + nv Dv)
//        = e^{-i v0 tau} 1/beta \sum_{nv} e^{-2pi i nv ntau} e^{-i (vn-v0) tau0} G(vn)
//        = e^{-i v0 tau} 1/beta DFT(e^{-i (vn-v0) t0} G(vn))
// Dt = T/N, Dv = V/N = 2pi/T
void ft_vn2tau(rvec& tau, cvec& Gtau, const rvec& vn, const cvec& Gvn, const bool incl_tail) {
    const int N = vn.size(); // vector size, all vector sizes must be equal
    if(N != tau.size() || N != Gtau.size() || N != Gvn.size())
        cout << "Error in ft_vn2tau: Vector sizes must be equal." << endl;
    const double Vtot = (vn[1]-vn[0])*(double)N; // frequency span = 2pi/beta*N
    const double beta = (2.*M_PI*(double)N)/Vtot; // time span = beta = inverse temperature
    for(int i=0; i<N; ++i) // fill time grid
        tau[i] = (beta/(double)N)*(double)i;
    fftw_complex in[N], out[N]; // input and output array in fftw format
    comp temp; // temporary variable to adjust Fourier convention (see above)
    if(incl_tail) { // incorporate FT(1/iv) = -1/2 * Theta(-tau) exactly
        for(int i=0; i<N; ++i) { // fill input array
            temp = (Gvn[i] + glb_i/vn[i])
                   * exp(-glb_i*(vn[i]-vn[0])*tau[0]); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    else {
        for(int i=0; i<N; ++i) { // fill input array
            temp = Gvn[i] * exp(-glb_i*(vn[i]-vn[0])*tau[0]); // see Fourier convention above
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    fftw_plan p; // fftw forward Fourier transform, result saved in 'out'
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    if(incl_tail) { // incorporate FT(iv) = -1/2 + Theta(-tau) exactly
        for (int i = 0; i < N; ++i) // fill output
            Gtau[i] = -0.5 + Theta(-tau[i]) +
                      comp(out[i][0], out[i][1]) * exp(-glb_i*vn[0]*tau[i]) / beta; // see Fourier convention above
    }
    else {
        for(int i=0; i<N; ++i) // fill output
            Gtau[i] = comp(out[i][0], out[i][1]) * exp(-glb_i*vn[0]*tau[i]) / beta; // see Fourier convention above
    }
    fftw_cleanup(); // clear up memory from fftw variables
}

// Fourier trafo from imaginary time (tau, Gtau) to imaginary frequency (vn, Gvn) using fftw,
// for fermions (type='f') or bosons (type='b')
// possibly including the tail analytically (incl_tail=true): Fourier trafo of 1/(iv) is (-1/2 + Theta(-tau))
// G(vn) = G(v0 + nv Dv) = \int dtau_0^\beta e^{i vn tau} G(tau)
//       = Dtau \sum_{ntau} e^{i(v0 + nv Dv)(tau0 + ntau Dtau) G(tau0 + ntau Dtau)
//       = e^{i vn tau0} \beta/N \sum_{ntau} e^{2pi i nv ntau} e^{i v0 (tau-tau0)} G(tau)
//       = e^{i vn tau0} \beta/N IDFT(e^{i v0 (tau-tau0)} G(tau))
// Dtau = beta/N, Dv = V/N = 2pi/beta
void ft_tau2vn(rvec& vn, cvec& Gvn, const rvec& tau, const cvec& Gtau, const char type, const bool incl_tail) {
    const int N = tau.size(); // vector size, all vector sizes must be equal
    if(N != vn.size() || N != Gvn.size() || N != Gtau.size())
        cout << "Error in ft_tau2vn: Vector sizes must be equal." << endl;
    const double beta = (tau[1]-tau[0])*(double)N; // time span = inverse temperature
    const int nmin = (int)(-floor(0.5*(double)N)); // lowest Matsubara-frequency index
    if(type=='f') { // fermionic Matsubara frequencies
        for (int i = 0; i < N; ++i) // fill frequency grid
            vn[i] = (2. * (double)(nmin + i) + 1.) / (M_PI * beta);
    }
    else if(type=='b') { // bosonic Matsubara frequencies
        for(int i=0; i<N; ++i) // fill frequency grid
            vn[i] = (2. * (double)(nmin + i)) / (M_PI*beta);
    }
    else
        cout << "Error in ft_tau2vn: Type must be f or b." << endl;
    fftw_complex in[N], out[N]; // input and output array in fftw format
    comp temp; // temporary variable to adjust Fourier convention (see above)
    if(incl_tail) { // incorporate FT(1/iv) = -1/2 + Theta(-tau) exactly
        for (int i = 0; i < N; ++i) { // fill input array
            temp = ( Gtau[i] - (-0.5 + Theta(-tau[0])) ) // subtract FT(1/iv)
                   * exp(glb_i * vn[0] * (tau[i] - tau[0])); // see Fourier convention above, tau0=0
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    else {
        for (int i = 0; i < N; ++i) { // fill input array
            temp = Gtau[i] * exp(glb_i * vn[0] * (tau[i] - tau[0])); // see Fourier convention above, tau0=0
            in[i][0] = temp.real(); // real part
            in[i][1] = temp.imag(); // imaginary part
        }
    }
    fftw_plan p; // fftw forward Fourier transform, result saved in 'out'
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    if(incl_tail) { // incorporate FT(iv) = -1/2 * Theta(-tau) exactly
        for (int i = 0; i < N; ++i) // fill output
            Gvn[i] = -glb_i / vn[i] + // explicit 1/iv
                     comp(out[i][0], out[i][1]) * exp(glb_i * vn[i] * tau[0]) * (beta / (double)N); // see Fourier convention above
    }
    else {
        for (int i = 0; i < N; ++i) // fill output
            Gvn[i] = comp(out[i][0], out[i][1]) * exp(glb_i * vn[i] * tau[0]) * (beta / (double)N); // see Fourier convention above
    }
    fftw_cleanup(); // clear up memory from fftw variables
}

/// Perform second-order perturbation theory (SOPT) using fast Fourier transform (FFT)
/// Output arguments: Fill results for self-energy (SE), K1 in a channel (K1a), K1 in p channel (K1p) in cvecs
/// Input arguments: Bare or full propagator Gin, value of interaction Uin (Gamma0 = -Uin/2), FFT parameters
/// FFT arguments: number of grid points (nFFT), frequency range (-V_FFT/2 ... V_FFT/2), suitable choice: 10000, 80.
/// Flags: whether to compute SE, K1a, K1p; whether to perform ladder summation
void SOPT_FFT(cvec& SEout_R, cvec& SEout_K, cvec& K1aout_R, cvec& K1aout_K, cvec& K1pout_R, cvec& K1pout_K, const Propagator<comp>& Gin,
        const double Uin, const int nFFT, const double V_FFT, const bool SEflag_FFT, const bool K1aflag_FFT, const bool K1pflag_FFT, const bool ladflag_FFT) {
    bool writeflag_FFT = true; // whether to write data into file
    bool tailflag_FFT = true; // whether to include analytic treatment of tail
    double hyb_FFT = 0; // hybridization for high-freuqency decay, e.g., G^R(v)=1/(v+i*hyb)
    if(tailflag_FFT) hyb_FFT = (glb_Gamma+Gin.Lambda)/2.;
    /// ------------------- READ INPUT --------------------- ///
    rvec vFFT(nFFT); // allocate equidistant frequency grid for Fourier trafo
    const double dvFFT = V_FFT / (double)nFFT; // frequency spacing
    if(dvFFT > glb_T) cout << "Warning in SOPT_FFT: frequency spacing dv=" << dvFFT << " is greater than temperature T=" << glb_T << "." << endl;
    for (int i = 0; i < nFFT; ++i) // fill frequency grid
        vFFT[i] = -V_FFT / 2. + dvFFT * (double)i;
    cvec GvR (nFFT), GvK (nFFT); // separate vectors for retarded and Keldysh component of input GF
    for (int i = 0; i < nFFT; ++i) {
        GvR[i] = Gin.valsmooth(0, vFFT[i], 0); // retarded component
        GvK[i] = Gin.valsmooth(1, vFFT[i], 0); // Keldysh component
    }
    /// ------------------- TRANSFORM INPUT --------------------- ///
    rvec tFFT (nFFT); // vector for real time
    cvec GtR (nFFT), GtK (nFFT), GtAc (nFFT), GtK_(nFFT); // vectors for G^R(t), G^K(t), G^A(t)*=G^R(-t), G^K(-t)
    // Fourier trafo, retarded component
    if(tailflag_FFT) ft_v2t(tFFT, GtR, vFFT, GvR, 1., hyb_FFT);
    else ft_v2t(tFFT, GtR, vFFT, GvR);
    ft_v2t(tFFT, GtK, vFFT, GvK); // Fourier trafo, Keldysh component
    flip_adjust(GtAc, GtR); // store G^A(t)*=G^R(-t) for convenience
    flip_adjust(GtK_, GtK); // store G^K(-t) for convenience
    if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_G.h5",
            {"t","Gt_R_R","Gt_R_I","Gt_K_R","Gt_K_I","v","Gv_R_R","Gv_R_I","Gv_K_R","Gv_K_I"},
            {tFFT,GtR.real(),GtR.imag(),GtK.real(),GtK.imag(),vFFT,GvR.real(),GvR.imag(),GvK.real(),GvK.imag()});
    /// ------------------- ALLOCATE K1 OUTPUT ON FFT GRID --------------------- ///
    cvec K1av_R(nFFT), K1av_K(nFFT), K1pv_R(nFFT), K1pv_K(nFFT); // these are needed to compute ladder contributions to self-energy
    if(K1aflag_FFT) {
        /// ------------------- K1a RETARDED COMPONENT --------------------- ///
        // sum of indices to the left of bare vertex is even, right sum is odd: retarded component of bosonic self-energy
        // vectors for K1a(t) terms arising in SOPT for retarded component
        cvec Piat_RK(nFFT), Piat_KA(nFFT);
        // individual terms
        Piat_RK = GtR * GtK_; // G^R(t) G^K(-t)
        Piat_KA = GtK * GtR.conj(); // G^K(t) G^A(-t)
        // resulting Pia retarded component in time and frequency
        cvec Piat_R(nFFT), Piav_R(nFFT);
        Piat_R = (Piat_RK + Piat_KA) * (-glb_i); // prefactor: 1/i (1/(2pi) is in Fourier trafo)
        ft_t2v(vFFT, Piav_R, tFFT, Piat_R); // Fourier trafo, retarded component
        K1av_R = Piav_R * (Uin*Uin/4.); // prefactor: \Gamma_0^2 = (-U/2)^2
        if(ladflag_FFT) // ladder summation: \Gamma_0^2 ( Pia + Pia \Gamma_0 Pia + Pia \Gamma_0 Pia \Gamma_0 Pia + ...)
            K1av_R *= (Piav_R*(Uin/2.) + 1.).inv();
        /// ------------------- K1a KELDYSH COMPONENT --------------------- ///
        // sum of indices to the left of bare vertex is odd, right sum is odd: Keldysh component of bosonic self-energy
        // vectors for K1a(t) terms arising in SOPT for Keldysh component
        cvec Piat_KK(nFFT), Piat_AR(nFFT), Piat_RA(nFFT);
        // individual terms
        Piat_KK = GtK * GtK_; // G^K(t) G^K(-t)
        Piat_AR = GtAc.conj() * GtAc; // G^A(t) G^R(-t)
        Piat_RA = GtR * GtR.conj(); // G^R(t) G^A(-t)
        if(tailflag_FFT) {
            cvec Piav_KK(nFFT), Piav_AR(nFFT), Piav_RA(nFFT);
            // Fourier trafos
            ft_t2v(vFFT, Piav_KK, tFFT, Piat_KK);
            ft_t2v(vFFT, Piav_AR, tFFT, Piat_AR, -glb_i*1., -2.*hyb_FFT); // G^A(t) G^R(-t) = Theta(-t)e^{2ht} = (-i) i Theta(-t)e^{-(-2h)t}
            ft_t2v(vFFT, Piav_RA, tFFT, Piat_RA, glb_i*1., 2.*hyb_FFT); // G^R(t) G^A(-t) = Theta(t)e^{-2ht} = i (-i) Theta(t)e^{-(2h)t}
            K1av_K = (Piav_KK + Piav_AR + Piav_RA) * (-glb_i*Uin*Uin/4.); // prefactor: (-U/2)^2 1/i
        }
        else {
            // resulting Pia Keldysh component in time and frequency
            cvec Piat_K(nFFT), Piav_K(nFFT);
            Piat_K = Piat_KK + Piat_AR + Piat_RA; // combined Pia Keldysh component
            ft_t2v(vFFT, Piav_K, tFFT, Piat_K); // Fourier trafo, Keldysh component
            K1av_K = Piav_K * (-glb_i*Uin*Uin/4.); // prefactor: (-U/2)^2 1/i
        }
        if(ladflag_FFT) // ladder summation: \Gamma_0^2 ( Pia + Pia \Gamma_0 Pia + Pia \Gamma_0 Pia \Gamma_0 Pia + ...)
            K1av_K *= ( (Piav_R*(Uin/2.) + 1.) * (Piav_R.conj()*(Uin/2.) + 1.) ).inv();
            //K1av_K_lad = K1av_K * ( (Piav_R*(Uin/2.) + 1.) * (Piav_R.conj()*(Uin/2.) + 1.) ).inv();
        /// ------------------- K1a OUTPUT --------------------- ///
        // resulting K1a components on external grid (nw1_a elements for frequencies bfreqs)
        interp1_FFT(K1aout_R, bfreqs, K1av_R, vFFT); // linear interpolation for retarded component
        interp1_FFT(K1aout_K, bfreqs, K1av_K, vFFT); // linear interpolation for Keldysh component
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_K1a.h5", {"v","K1a_R_R","K1a_R_I","K1a_K_R","K1a_K_I"}, {vFFT,K1av_R.real(),K1av_R.imag(),K1av_K.real(),K1av_K.imag()});
    }
    if(K1pflag_FFT) {
        /// ------------------- K1p RETARDED COMPONENT --------------------- ///
        // sum of indices to the left of bare vertex is even, right sum is odd: retarded component of bosonic self-energy
        // vectors for K1p(t) terms arising in SOPT for retarded component
        cvec Pipt_RK(nFFT);
        // only one term but appears twice
        Pipt_RK = GtR * GtK; // G^R(t) G^K(t)
        // resulting Pip retarded component in time and frequency and K1p retarded component in frequency
        cvec Pipt_R(nFFT), Pipv_R(nFFT);
        Pipt_R = Pipt_RK * (-glb_i*2.); // prefactor: 1/i (1/(2pi) is in Fourier trafo) twice
        ft_t2v(vFFT, Pipv_R, tFFT, Pipt_R); // Fourier trafo, retarded component
        K1pv_R = Pipv_R * (Uin*Uin/4.); // prefactor: \Gamma_0^2 = (-U/2)^2
        if(ladflag_FFT) // ladder summation: \Gamma_0^2 ( Pip + Pip \Gamma_0 Pip + Pip \Gamma_0 Pip \Gamma_0 Pip + ...)
            K1pv_R *= (Pipv_R*(Uin/2.) + 1.).inv();
        /// ------------------- K1p KELDYSH COMPONENT --------------------- ///
        // sum of indices to the left of bare vertex is odd, right sum is odd: Keldysh component of bosonic self-energy
        // vectors for K1p(t) terms arising in SOPT for Keldysh component
        cvec Pipt_KK(nFFT), Pipt_RR(nFFT), Pipt_AA(nFFT);
        // individual terms
        Pipt_KK = GtK * GtK; // G^K(t) G^K(t)
        Pipt_RR = GtR * GtR; // G^R(t) G^R(t)
        Pipt_AA = GtAc.conj() * GtAc.conj(); // G^A(t) G^A(t)
        if(tailflag_FFT) {
            cvec Pipv_KK(nFFT), Pipv_RR(nFFT), Pipv_AA(nFFT);
            // Fourier trafos
            ft_t2v(vFFT, Pipv_KK, tFFT, Pipt_KK);
            ft_t2v(vFFT, Pipv_RR, tFFT, Pipt_RR, -glb_i*1., 2.*hyb_FFT); // G^R(t)^2 = -Theta(t)e^{-2ht} = (-i) (-i) Theta(t)e^{-(2h)t}
            ft_t2v(vFFT, Pipv_AA, tFFT, Pipt_AA, glb_i*1., -2.*hyb_FFT); // G^A(t)^2 = -Theta(-t)e^{2ht} = (i) i Theta(-t)e^{-(-2h)t}
            K1pv_K = (Pipv_KK + Pipv_RR + Pipv_AA) * (-glb_i*Uin*Uin/4.); // prefactor: (-U/2)^2 1/i
        }
        else {
            // resulting Pip Keldysh component in time and frequency
            cvec Pipt_K(nFFT), Pipv_K(nFFT);
            Pipt_K = Pipt_KK + Pipt_RR + Pipt_AA; // combined Pip Keldysh component
            ft_t2v(vFFT, Pipv_K, tFFT, Pipt_K); // Fourier trafo, Keldysh component
            K1pv_K = Pipv_K * (-glb_i*Uin*Uin/4.); // prefactor: (-U/2)^2 1/i
        }
        if(ladflag_FFT) // ladder summation: \Gamma_0^2 ( Pip + Pip \Gamma_0 Pip + Pip \Gamma_0 Pip \Gamma_0 Pip + ...)
            K1pv_K *= ( (Pipv_R*(Uin/2.) + 1.) * (Pipv_R.conj()*(Uin/2.) + 1.) ).inv();
        /// ------------------- K1p OUTPUT --------------------- ///
        // resulting K1p components on external grid (nw1_p elements for frequencies bfreqs)
        interp1_FFT(K1pout_R, bfreqs, K1pv_R, vFFT); // linear interpolation for component 1
        interp1_FFT(K1pout_K, bfreqs, K1pv_K, vFFT); // linear interpolation for component 3
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_K1p_lad.h5", {"v","K1p_R_R","K1p_R_I","K1p_K_R","K1p_K_I"}, {vFFT,K1pv_R.real(),K1pv_R.imag(),K1pv_K.real(),K1pv_K.imag()});
    }
    if(SEflag_FFT) {
        /// ------------------- SELF-ENERGY RETARDED COMPONENT --------------------- ///
        // vectors for \Sigma(t) terms arising in SOPT for retarded component
        cvec SEt_RKK(nFFT), SEt_KKA(nFFT), SEt_KRK(nFFT), SEt_RRA(nFFT);
        // individual terms
        SEt_RKK = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
        SEt_KKA = GtK * GtK * GtR.conj(); // G^K(t) G^K(t) G^A(-t)
        SEt_KRK = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
        SEt_RRA = GtR * GtR * GtR.conj(); // G^R(t) G^R(t) G^A(-t)
        // resulting retarded self-energy in time and frequency
        cvec SEt_R(nFFT), SEv_R(nFFT);
        if(tailflag_FFT) {
            cvec SEv_RRA(nFFT);
            SEt_R = SEt_RKK + SEt_KKA + SEt_KRK; // without RAR component
            // Fourier trafos
            ft_t2v(vFFT, SEv_R, tFFT, SEt_R); // Fourier trafo, without RAR component
            ft_t2v(vFFT, SEv_RRA, tFFT, SEt_RRA, 1., 3.*hyb_FFT); // G^R(t)^2 G^A(-t) = -i Theta(t)e^{-3ht}
            SEv_R = (SEv_R + SEv_RRA) * (Uin*Uin/4.); // prefactor: -(-U/2/i)^2
        }
        else {
            SEt_R = (SEt_RKK + SEt_KKA + SEt_KRK + SEt_RRA) * (Uin*Uin/4.); // prefactor: -(-U/2/i)^2
            ft_t2v(vFFT, SEv_R, tFFT, SEt_R); // Fourier trafo, retarded component
        }
        /// ------------------- SELF-ENERGY KELDYSH COMPONENT --------------------- ///
        // vectors for \Sigma(t) terms arising in SOPT for Keldysh component
        cvec SEt_RKA(nFFT), SEt_RRK(nFFT), SEt_AKR(nFFT), SEt_AAK(nFFT), SEt_KRA(nFFT), SEt_KAR(nFFT), SEt_KKK(nFFT);
        // individual terms
        SEt_RKA = GtR * GtK * GtR.conj(); // G^R(t) G^K(t) G^A(-t)
        SEt_RRK = GtR * GtR * GtK_; // G^R(t) G^R(t) G^K(-t)
        SEt_AKR = GtAc.conj() * GtK * GtAc; // G^A(t) G^K(t) G^R(-t)
        SEt_AAK = GtAc.conj() * GtAc.conj() * GtK_; // G^A(t) G^A(t) G^K(-t)
        SEt_KRA = GtK * GtR * GtR.conj(); // G^K(t) G^R(t) G^A(-t)
        SEt_KAR = GtK * GtAc.conj() * GtAc; // G^K(t) G^A(t) G^R(-t)
        SEt_KKK = GtK * GtK * GtK_; // G^K(t) G^K(t) G^K(-t)
        // resulting Keldysh self-energy in time and frequency
        cvec SEt_K(nFFT), SEv_K(nFFT);
        SEt_K = (SEt_RKA + SEt_RRK + SEt_AKR + SEt_AAK + SEt_KRA + SEt_KAR + SEt_KKK) * (Uin*Uin/4.); // prefactor: -(-U/2/i)^2
        ft_t2v(vFFT, SEv_K, tFFT, SEt_K); // Fourier trafo, Keldysh component
        if(ladflag_FFT) { // use K1a and K1p ladders for self-energy correction
            /// PREPARE K1 LADDERS ///
            cvec K1at_lad_R(nFFT), K1at_lad_Ac(nFFT), K1at_lad_K(nFFT); // K1a ladders in time
            cvec K1pt_lad_R(nFFT), K1pt_lad_Ac(nFFT), K1pt_lad_K(nFFT); // K1p ladders in time
            ft_v2t(tFFT, K1at_lad_R, vFFT, K1av_R); // Fourier trafo, retarded component, K1av_R was updated to ladder
            flip_adjust(K1at_lad_Ac, K1at_lad_R); // store K1a^A(t)* = K1a^R(-t) for convenience
            ft_v2t(tFFT, K1at_lad_K, vFFT, K1av_K); // Fourier trafo, Keldysh component, K1av_K was updated to ladder
            ft_v2t(tFFT, K1pt_lad_R, vFFT, K1pv_R); // Fourier trafo, retarded component, K1pv_R was updated to ladder
            flip_adjust(K1pt_lad_Ac, K1pt_lad_R); // store K1p^A(t)* = K1p^R(-t) for convenience
            ft_v2t(tFFT, K1pt_lad_K, vFFT, K1pv_K); // Fourier trafo, Keldysh component, K1pv_K was updated to ladder
            /// INSERT LADDERS INTO SELF-ENERGY ///
            cvec SEt_alad_R(nFFT), SEt_plad_R(nFFT), SEt_alad_K(nFFT), SEt_plad_K(nFFT); // self-energy ladders in time
            // SE = - contracted vertex -> prefactor -(1/i) (1/(2pi) is in Fourier trafo)
            SEt_alad_R = (K1at_lad_R * GtK + K1at_lad_K * GtR)*(glb_i); // K1a^R(t) * G^K(t) + K1a^K(t) * G^R(t)
            SEt_plad_R = (K1pt_lad_R * GtK_ + K1pt_lad_K * GtR.conj())*(glb_i); // K1p^R(t) * G^K(-t) + K1p^K(t) * G^A(-t)
            SEt_alad_K = (K1at_lad_R * GtR + K1at_lad_Ac.conj() * GtAc.conj() + K1at_lad_K * GtK)*(glb_i); // K1a^R(t) * G^R(t) + K1a^A(t) * G^A(t) + K1a^K(t) * G^K(t)
            SEt_plad_K = (K1pt_lad_R * GtR.conj() + K1pt_lad_Ac.conj() * GtAc + K1pt_lad_K * GtK_)*(glb_i); // K1p^R(t) * G^A(-t) + K1p^A(t) * G^R(-t) + K1a^K(t) * G^K(-t)
            cvec SEv_alad_R(nFFT), SEv_plad_R(nFFT), SEv_alad_K(nFFT), SEv_plad_K(nFFT); // self-energy ladders in frequency
            ft_t2v(vFFT, SEv_alad_R, tFFT, SEt_alad_R); // Fourier trafo, retarded component
            ft_t2v(vFFT, SEv_plad_R, tFFT, SEt_plad_R); // Fourier trafo, retarded component
            ft_t2v(vFFT, SEv_alad_K, tFFT, SEt_alad_K); // Fourier trafo, Keldysh component
            ft_t2v(vFFT, SEv_plad_K, tFFT, SEt_plad_K); // Fourier trafo, Keldysh component
            SEv_R = SEv_alad_R + SEv_plad_R - SEv_R; // combine and subtract SOPT which is contained in both ladders
            SEv_K = SEv_alad_K + SEv_plad_K - SEv_K; // combine and subtract SOPT which is contained in both ladders
        }
        /// ------------------- SELF-ENERGY OUTPUT --------------------- ///
        // resulting self-energy components on external grid (nSE elements for frequencies ffreqs)
        interp1_FFT(SEout_R, ffreqs, SEv_R, vFFT); // linear interpolation for retarded component
        interp1_FFT(SEout_K, ffreqs, SEv_K, vFFT); // linear interpolation for Keldysh component
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_SE.h5", {"v","SE_R_R","SE_R_I","SE_K_R","SE_K_I"}, {vFFT,SEv_R.real(),SEv_R.imag(),SEv_K.real(),SEv_K.imag()});
    }
}

void SOPT_FFT_SE_R(cvec& SEout_R, const Propagator<comp>& Gin, const double Uin, const int nFFT, const double V_FFT) {
    cvec K1a_R, K1a_K, K1p_R, K1p_K; // dummy cvecs for output that is not computed at all
    cvec SE_K(nSE); // dummy cvec for output that is not returned
    SOPT_FFT(SEout_R, SE_K, K1a_R, K1a_K, K1p_R, K1p_K, Gin, Uin, nFFT, V_FFT, true, false, false, false);
    for (int iv = 0; iv < nSE; ++iv) // add first-order perturbation theory: Hartree term
        SEout_R[iv] += glb_U / 2.;
}

void SOPT_FFT_K1a_R(cvec& K1aout_R, const Propagator<comp>& Gin, const double Uin, const int nFFT, const double V_FFT) {
    cvec SE_R, SE_K, K1p_R, K1p_K; // dummy cvecs for output that is not computed at all
    cvec K1a_K(nw1_a); // dummy cvec for output that is not returned
    SOPT_FFT(SE_R, SE_K, K1aout_R, K1a_K, K1p_R, K1p_K, Gin, Uin, nFFT, V_FFT, false, true, false, false);
}

void SOPT_FFT_K1p_R(cvec& K1pout_R, const Propagator<comp>& Gin, const double Uin, const int nFFT, const double V_FFT) {
    cvec SE_R, SE_K, K1a_R, K1a_K; // dummy cvecs for output that is not computed at all
    cvec K1p_K(nw1_p); // dummy cvec for output that is not returned
    SOPT_FFT(SE_R, SE_K, K1a_R, K1a_K, K1pout_R, K1p_K, Gin, Uin, nFFT, V_FFT, false, false, true, false);
}

void SOPT_FFT(SelfEnergy<comp>& SEout, rvert<comp>& gammaAout, rvert<comp>& gammaPout, const Propagator<comp>& Gin,
        const double Uin, const int nFFT, const double V_FFT, const bool SEflag_FFT, const bool K1aflag_FFT, const bool K1pflag_FFT, const bool ladflag_FFT) {
    /// prepare output ///
    cvec SE_R(nSE), SE_K(nSE); // cvecs for SOPT self-energy
    cvec K1a_R(nw1_a), K1a_K(nw1_a); // cvecs for SOPT K1a
    cvec K1p_R(nw1_p), K1p_K(nw1_p); // cvecs for SOPT K1p
    /// perform SOPT ///
    SOPT_FFT(SE_R, SE_K, K1a_R, K1a_K, K1p_R, K1p_K, Gin, Uin, nFFT, V_FFT, SEflag_FFT, K1aflag_FFT, K1pflag_FFT, ladflag_FFT);
    /// fill output ///
    if(SEflag_FFT) {
        SEout.initialize(glb_U/2., 0.); // self-energy in first-order perturbation theory: Hartree term
        for (int iv = 0; iv < nSE; ++iv) { // self-energy
            SEout.addself(0, iv, 0, SE_R[iv]); // retarded component
            SEout.addself(1, iv, 0, SE_K[iv]); // Keldyish component
        }
    }
    if(K1aflag_FFT) {
        for (int iw = 0; iw < nw1_a; ++iw) { // K1a
            gammaAout.K1_setvert(0, iw, 0, K1a_R[iw]); // retarded component
            gammaAout.K1_setvert(1, iw, 0, K1a_K[iw]); // Keldysh component
        }
    }
    if(K1pflag_FFT) {
        for (int iw = 0; iw < nw1_p; ++iw) { // K1p
            gammaPout.K1_setvert(0, iw, 0, K1p_R[iw]); // retarded component
            gammaPout.K1_setvert(1, iw, 0, K1p_K[iw]); // Keldysh component
        }
    }
}

void SOPT_FFT(SelfEnergy<comp>& SEout, const Propagator<comp>& Gin, const double Uin, const int nFFT, const double V_FFT, const bool ladflag_FFT) {
    rvert<comp> gammaA ('a'); rvert<comp> gammaP ('p'); // dummy objects for output that is not computed at all
    SOPT_FFT(SEout, gammaA, gammaP, Gin, Uin, nFFT, V_FFT, true, false, false, ladflag_FFT);
}

/// Compute scale-differentiated second-order perturbation theory (SOPT) using fast Fourier transform (FFT), arguments as above
void diffSOPT_FFT(cvec& dSEout_R, cvec& dK1aout_R, const Propagator<comp>& Gin, const Propagator<comp>& Sin,
              const double Uin, const int nFFT, const double V_FFT, const bool SEflag_FFT, const bool K1aflag_FFT) {
    bool writeflag_FFT = false; // whether to write data into file
    /// ------------------- READ INPUT --------------------- ///
    rvec vFFT(nFFT); // allocate equidistant frequency grid for Fourier trafo
    const double dvFFT = V_FFT / (double)nFFT; // frequency spacing
    if(dvFFT > glb_T) cout << "Warning in SOPT_FFT: frequency spacing dv=" << dvFFT << " is greater than temperature T=" << glb_T << "." << endl;
    for (int i = 0; i < nFFT; ++i) // fill frequency grid
        vFFT[i] = -V_FFT / 2. + dvFFT * (double)i;
    cvec GvR (nFFT), GvK (nFFT), SvR (nFFT), SvK (nFFT); // separate vectors for retarded and Keldysh component of input GF
    for (int i = 0; i < nFFT; ++i) {
        GvR[i] = Gin.valsmooth(0, vFFT[i], 0); // retarded component
        GvK[i] = Gin.valsmooth(1, vFFT[i], 0); // Keldysh component
        SvR[i] = Sin.valsmooth(0, vFFT[i], 0); // retarded component
        SvK[i] = Sin.valsmooth(1, vFFT[i], 0); // Keldysh component
    }
    /// ------------------- TRANSFORM INPUT --------------------- ///
    rvec tFFT (nFFT); // vector for real time
    cvec GtR (nFFT), GtK (nFFT), GtAc (nFFT), GtK_(nFFT); // vectors for G^R(t), G^K(t), G^A(t)*=G^R(-t), G^K(-t)
    cvec StR (nFFT), StK (nFFT), StAc (nFFT), StK_(nFFT); // vectors for G^R(t), G^K(t), G^A(t)*=G^R(-t), G^K(-t)
    // Fourier trafo, retarded component
    ft_v2t(tFFT, GtR, vFFT, GvR);
    ft_v2t(tFFT, StR, vFFT, SvR);
    // Fourier trafo, Keldysh component
    ft_v2t(tFFT, GtK, vFFT, GvK);
    ft_v2t(tFFT, StK, vFFT, SvK);
    flip_adjust(GtAc, GtR); // store G^A(t)*=G^R(-t) for convenience
    flip_adjust(GtK_, GtK); // store G^K(-t) for convenience
    flip_adjust(StAc, StR); // store S^A(t)*=S^R(-t) for convenience
    flip_adjust(StK_, StK); // store S^K(-t) for convenience
    if(SEflag_FFT) {
        /// ------------------- SELF-ENERGY RETARDED COMPONENT --------------------- ///
        // vectors for \Sigma(t) terms arising in SOPT for retarded component
        cvec SEt_RKK_SGG(nFFT), SEt_RKK_GSG(nFFT), SEt_RKK_GGS(nFFT), SEt_KKA_SGG(nFFT), SEt_KKA_GGS(nFFT), SEt_RRA_SGG(nFFT), SEt_RRA_GGS(nFFT);
        // individual terms
        SEt_RKK_SGG = StR * GtK * GtK_; // S^R(t) G^K(t) G^K(-t)
        SEt_RKK_GSG = GtR * StK * GtK_; // G^R(t) S^K(t) G^K(-t)
        SEt_RKK_GGS = GtR * GtK * StK_; // G^R(t) G^K(t) S^K(-t)

        SEt_KKA_SGG = StK * GtK * GtR.conj(); // S^K(t) G^K(t) G^A(-t)
        SEt_KKA_GGS = GtK * GtK * StR.conj(); // G^K(t) G^K(t) S^A(-t)

        SEt_RRA_SGG = StR * GtR * GtR.conj(); // S^R(t) G^R(t) G^A(-t)
        SEt_RRA_GGS = GtR * GtR * StR.conj(); // G^R(t) G^R(t) S^A(-t)
        // resulting retarded self-energy in time and frequency
        cvec SEtR(nFFT), SEvR(nFFT);
        SEtR = ( (SEt_RKK_SGG+SEt_RKK_GSG+SEt_RKK_GGS)*2. + (SEt_KKA_SGG*2.+SEt_KKA_GGS) + (SEt_RRA_SGG*2.+SEt_RRA_GGS) ) * (Uin*Uin/4.); // prefactor: -(-U/2/i)^2
        ft_t2v(vFFT, SEvR, tFFT, SEtR); // Fourier trafo, retarded component
        /// ------------------- SELF-ENERGY OUTPUT --------------------- ///
        interp1_FFT(dSEout_R, ffreqs, SEvR, vFFT); // linear interpolation for retarded component
        if(writeflag_FFT) write_h5_rvecs("dSOPT_FFT_SE.h5", {"v","SE_R_R","SE_R_I"}, {vFFT,SEvR.real(),SEvR.imag()});
    }
    if(K1aflag_FFT) {
        /// ------------------- K1a RETARDED COMPONENT --------------------- ///
        // sum of indices to the left of bare vertex is even, right sum is odd: retarded component of bosonic self-energy
        // vectors for K1a(t) terms arising in SOPT for retarded component
        cvec K1at_RK_GS(nFFT), K1at_KA_GS(nFFT), K1at_RK_SG(nFFT), K1at_KA_SG(nFFT);
        // individual terms
        K1at_RK_GS = GtR * StK_; // G^R(t) S^K(-t)
        K1at_KA_GS = GtK * StR.conj(); // G^K(t) S^A(-t)
        K1at_RK_SG = StR * GtK_; // S^R(t) G^K(-t)
        K1at_KA_SG = StK * GtR.conj(); // S^K(t) G^A(-t)
        // resulting K1a retarded component in time and frequency
        cvec K1at_R(nFFT), K1av_R(nFFT);
        K1at_R = (K1at_RK_GS + K1at_KA_GS + K1at_RK_SG + K1at_KA_SG) * (-glb_i*Uin*Uin/4.); // prefactor: (-U/2)^2 1/i
        ft_t2v(vFFT, K1av_R, tFFT, K1at_R); // Fourier trafo, retarded component
        /// ------------------- K1a OUTPUT --------------------- ///
        interp1_FFT(dK1aout_R, bfreqs, K1av_R, vFFT); // linear interpolation for component 1
        if(writeflag_FFT) write_h5_rvecs("dSOPT_FFT_K1a.h5", {"v","K1a_R_R","K1a_R_I"}, {vFFT,K1av_R.real(),K1av_R.imag()});
    }
}

void diffSOPT_FFT_SE_R(cvec& dSEout_R, const Propagator<comp>& Gin, const Propagator<comp>& Sin, const double Uin, const int nFFT, const double V_FFT) {
    cvec dK1a_R; // dummy cvecs for output that is not computed at all
    diffSOPT_FFT(dSEout_R, dK1a_R, Gin, Sin, Uin, nFFT, V_FFT, true, false);
}

void diffSOPT_FFT_K1a_R(cvec& dK1aout_R, const Propagator<comp>& Gin, const Propagator<comp>& Sin, const double Uin, const int nFFT, const double V_FFT) {
    cvec dSE_R; // dummy cvecs for output that is not computed at all
    diffSOPT_FFT(dSE_R, dK1aout_R, Gin, Sin, Uin, nFFT, V_FFT, false, true);
}

#endif // FOURIER_TRAFO_H