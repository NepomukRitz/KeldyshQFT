#ifndef FOURIER_TRAFO_H
#define FOURIER_TRAFO_H

#include <cmath> // needed for exponential function to adjust Fourier convention
#include <fftw3.h> // fftw library
#include "write_data2file.h"

using namespace std;

double Theta(double x) {
    // Heaviside step function
    if (x > 0) return 1.;
    else if (x < 0) return 0.;
    else return 1./2.;
}

double sign(double x) {
    // sign function
    if (x > 0) return 1.;
    else if (x < 0) return -1.;
    else return 0.;
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
        cout << "Error in interp1_FFT: Vector sizes must be equal." << endl;
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
    cout << "Differences from exact results are: " << fy[0] - 0. << " " << fy[1] - 0.5 << " " << fy[2] - 2.5 << " "
         << fy[3] - 0. << "." << endl;
}

// Fourier trafo from real frequency (v, Gv) to real time (t, Gt) using fftw
// G(t) = G(t0 + nt Dt) = 1/2pi \int dv e^{-ivt} G(v) = Dv/2pi \sum_{nv} e^{-i(v0 + nv Dv)(t0 + nt Dt) G(v0 + nv Dv)
//      = e^{-i v0 t} 1/T \sum_{nv} e^{-2pi i nv nt} e^{-i (v-v0) t0} G(v) = e^{-i v0 t} 1/T DFT(e^{-i (v-v0) t0} G(v))
// Dt = T/N, Dv = V/N = 2pi/T
// hyb encodes analytic high-frequency decay tail_coef/(v+i tail_hyb), setting coef=0 ignores this
void ft_v2t(rvec& t, cvec& Gt, const rvec& v, const cvec& Gv, const comp tail_coef, const double tail_hyb)
{
    const int N = v.size(); // vector size, all vector sizes must be equal
    if(N != t.size() || N != Gt.size() || N != Gv.size())
        cout << "Error in ft_v2t: Vector sizes must be equal." << endl;
    const double Vtot = (v[1]-v[0])*(double)N; // frequency span
    const double Ttot = (2.*pi*(double)N)/Vtot; // time span
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
            Gt[i] = ( -glb_i * tail_coef * sign(tail_hyb) * Theta(sign(tail_hyb)*t[i] + 1e-16) * exp(-tail_hyb*t[i]) ) +
                    comp(out[i][0], out[i][1]) * exp(-glb_i * v[0] * t[i]) / Ttot; // see Fourier convention above
    }
    else { // standard procedure
        for (int i = 0; i < N; ++i) // fill output
            Gt[i] = comp(out[i][0], out[i][1]) * exp(-glb_i * v[0] * t[i]) / Ttot; // see Fourier convention above
    }
    fftw_cleanup(); // clear up memory from fftw variables
}
void ft_v2t(rvec& t, cvec& Gt, const rvec& v, const cvec& Gv) {
    ft_v2t(t, Gt, v, Gv, 0., 0.);
}

// Fourier trafo from real time (t, Gt) to real frequency (v, Gv) using fftw
// G(v) = G(v0 + nv Dv) = \int dt e^{ivt} G(t) = Dt \sum_{nt} e^{i(v0 + nv Dv)(t0 + nt Dt) G(t0 + nt Dt)
//      = e^{i v t0} T/N \sum_{nt} e^{2pi i nv nt} e^{i v0 (t-t0)} G(t) = e^{i v0 t} T/N IDFT(e^{i v0 (t-t0)} G(t))
// Dt = T/N, Dv = V/N = 2pi/T
// hyb encodes analytic high-frequency decay tail_coef/(v+i tail_hyb), setting coef=0 ignores this
void ft_t2v(rvec& v, cvec& Gv, const rvec& t, const cvec& Gt, const comp tail_coef, const double tail_hyb)
{
    const int N = t.size(); // vector size, all vector sizes must be equal
    if(N != v.size() || N != Gv.size() || N != Gt.size())
        cout << "Error in ft_t2v: Vector sizes must be equal." << endl;
    const double Ttot = (t[1]-t[0])*(double)N; // time span
    const double Vtot = (2.*pi*(double)N)/Ttot; // frequency span
    for(int i=0; i<N; ++i) // fill frequency grid
        v[i] = -Vtot/2. + (Vtot/(double)N)*(double)i;
    fftw_complex in[N], out[N]; // input and output array in fftw format
    comp temp; // temporary variable to adjust Fourier convention (see above)
    if(abs(tail_coef)>1e-10) { // include long-time decay -i sgn(tail_hyb) tail_coef \Theta(sgn(tail_hyb)*t) e^{-tail_hyb t} explicitly, here enforce Theta(0)=1 via +1e-16
        for (int i = 0; i < N; ++i) { // fill input array
            temp = ( Gt[i] + glb_i * tail_coef * sign(tail_hyb) * Theta(sign(tail_hyb)*t[i] + 1e-16) * exp(-tail_hyb*t[i]) )
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
    fftw_cleanup(); // clear up memory from fftw variables
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
void ft_vn2tau(rvec& tau, cvec& Gtau, const rvec& vn, const cvec& Gvn, const bool incl_tail)
{
    const int N = vn.size(); // vector size, all vector sizes must be equal
    if(N != tau.size() || N != Gtau.size() || N != Gvn.size())
        cout << "Error in ft_vn2tau: Vector sizes must be equal." << endl;
    const double Vtot = (vn[1]-vn[0])*(double)N; // frequency span = 2pi/beta*N
    const double beta = (2.*pi*(double)N)/Vtot; // time span = beta = inverse temperature
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
void ft_tau2vn(rvec& vn, cvec& Gvn, const rvec& tau, const cvec& Gtau, const char type, const bool incl_tail)
{
    const int N = tau.size(); // vector size, all vector sizes must be equal
    if(N != vn.size() || N != Gvn.size() || N != Gtau.size())
        cout << "Error in ft_tau2vn: Vector sizes must be equal." << endl;
    const double beta = (tau[1]-tau[0])*(double)N; // time span = inverse temperature
    const int nmin = (int)(-floor(0.5*(double)N)); // lowest Matsubara-frequency index
    if(type=='f') { // fermionic Matsubara frequencies
        for (int i = 0; i < N; ++i) // fill frequency grid
            vn[i] = (2. * (double)(nmin + i) + 1.) / (pi * beta);
    }
    else if(type=='b') { // bosonic Matsubara frequencies
            for(int i=0; i<N; ++i) // fill frequency grid
                vn[i] = (2. * (double)(nmin + i)) / (pi*beta);
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


// TODO: overload = operator for basic_vec
void SOPTbare_FFT_SelfEnergy(SelfEnergy<comp>& SEout, const double Lambda, const double Uin, const int nFFT, const double V_FFT) {
    // whether to compute self-energy (SE), antiparallel/particle-hole bubble (PiA), parallel/particle-particle bubble (PiP)
    bool SEflag_FFT = true;
    bool PiAflag_FFT = true;
    bool PiPflag_FFT = true;
    bool writeflag_FFT = false; // whether to write data into file
    bool tailflag_FFT = true; // whether to include analytic treatment of tail
    double hyb_FFT = 0; // hybridization for high-freuqency decay, e.g., G^R(v)=1/(v+i*hyb)
    if(tailflag_FFT) hyb_FFT = (glb_Gamma_REG+Lambda)/2.;
    /// ------------------- READ INPUT --------------------- ///
    rvec vFFT(nFFT); // allocate equidistant frequency grid for Fourier trafo
    const double dvFFT = V_FFT / (double)nFFT; // frequency spacing
    if(dvFFT > glb_T) cout << "Warning in SOPT_FFT: frequency spacing dv=" << dvFFT << " is greater than temperature T=" << glb_T << "." << endl;
    for (int i = 0; i < nFFT; ++i) // fill frequency grid
        vFFT[i] = -V_FFT / 2. + dvFFT * (double)i;
    cvec GvR (nFFT), GvK (nFFT); // separate vectors for retarded and Keldysh component of input GF
    for (int i = 0; i < nFFT; ++i) {
        GvR[i] = gR(Lambda, vFFT[i]); // retarded component
        GvK[i] = gK(Lambda, vFFT[i]); // Keldysh component
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
    if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_G.h5", {"t","G_R_R","G_R_I","G_K_R","G_K_I"}, {tFFT,GtR.real(),GtR.imag(),GtK.real(),GtK.imag()});
    if(SEflag_FFT) {
        /// ------------------- SELF-ENERGY RETARDED COMPONENT --------------------- ///
        // vectors for \Sigma(t) terms arising in SOPT for retarded component
        cvec SEt_RKK(nFFT), SEt_KAK(nFFT), SEt_KKR(nFFT), SEt_RAR(nFFT);
        // individual terms
        SEt_RKK = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
        SEt_KAK = GtK * GtK * GtR.conj(); // G^K(t) G^K(t) G^A(-t)
        SEt_KKR = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
        SEt_RAR = GtR * GtR * GtR.conj(); // G^R(t) G^R(t) G^A(-t) //!
        // resulting retarded self-energy in time and frequency
        cvec SEtR(nFFT), SEvR(nFFT);
        if(tailflag_FFT) {
            cvec SEv_RAR(nFFT);
            SEtR = SEt_RKK + SEt_KAK + SEt_KKR; // without RAR component
            // Fourier trafos
            ft_t2v(vFFT, SEvR, tFFT, SEtR); // Fourier trafo, without RAR component
            ft_t2v(vFFT, SEv_RAR, tFFT, SEt_RAR, 1., 3.*hyb_FFT); // G^R(t)^2 G^A(-t) = -i Theta(t)e^{-3ht}
            SEvR = (SEvR + SEv_RAR) * (Uin*Uin/4.); // prefactor: -(\pm U/2/i)^2
        }
        else {
            SEtR = (SEt_RKK + SEt_KAK + SEt_KKR + SEt_RAR) * (Uin*Uin/4.); // prefactor: -(\pm U/2/i)^2
            ft_t2v(vFFT, SEvR, tFFT, SEtR); // Fourier trafo, retarded component
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
        cvec SEtK(nFFT), SEvK(nFFT);
        SEtK = (SEt_RKA + SEt_RRK + SEt_AKR + SEt_AAK + SEt_KRA + SEt_KAR + SEt_KKK) * (Uin*Uin/4.); // prefactor: -(\pm U/2/i)^2
        ft_t2v(vFFT, SEvK, tFFT, SEtK); // Fourier trafo, Keldysh component
        /// ------------------- SELF-ENERGY OUTPUT --------------------- ///
        cvec SEoutR(nSE), SEoutK(nSE); // resulting self-energy components on external grid
        interp1_FFT(SEoutR, ffreqs, SEvR, vFFT); // linear interpolation for retarded component
        interp1_FFT(SEoutK, ffreqs, SEvK, vFFT); // linear interpolation for Keldysh component
        for (int i = 0; i < nSE; ++i) { // fill results in output self-energy
            SEout.setself(0, i, SEoutR[i]); // retarded component
            SEout.setself(1, i, SEoutK[i]); // Keldysh component
        }
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_SE.h5", {"v","SE_R_R","SE_R_I","SE_K_R","SE_K_I"}, {vFFT,SEvR.real(),SEvR.imag(),SEvK.real(),SEvK.imag()});
    }
    if(PiAflag_FFT) {
        /// ------------------- K1a COMPONENT 1 --------------------- ///
        // left sum of indices is even, right sum of indices is odd, retarded component of bosonic self-energy
        // vectors for \K1a(t) terms arising in SOPT for component 1
        cvec K1at_RK(nFFT), K1at_KA(nFFT);
        // individual terms
        K1at_RK = GtR * GtK_; // G^R(t) G^K(-t)
        K1at_KA = GtK * GtR.conj(); // G^K(t) G^A(-t)
        // resulting K1a component 1 in time and frequency
        cvec K1at_1(nFFT), K1av_1(nFFT);
        K1at_1 = (K1at_RK + K1at_KA) * (-glb_i*Uin*Uin/4.); // prefactor: (\pm U/2)^2 1/i
        ft_t2v(vFFT, K1av_1, tFFT, K1at_1); // Fourier trafo, component 3
        /// ------------------- K1a COMPONENT 3 --------------------- ///
        // left sum of indices is odd, right sum of indices is odd, Keldysh component of bosonic self-energy
        // vectors for \K1a(t) terms arising in SOPT for component 3
        cvec K1at_KK(nFFT), K1at_AR(nFFT), K1at_RA(nFFT);
        // individual terms
        K1at_KK = GtK * GtK_; // G^K(t) G^K(-t)
        K1at_AR = GtAc.conj() * GtAc; // G^A(t) G^R(-t)
        K1at_RA = GtR * GtR.conj(); // G^R(t) G^A(-t)
        // resulting K1a component 3 in time and frequency
        cvec K1at_3(nFFT), K1av_3(nFFT);
        if(tailflag_FFT) {
            cvec K1av_KK(nFFT), K1av_AR(nFFT), K1av_RA(nFFT);
            // Fourier trafos
            ft_t2v(vFFT, K1av_KK, tFFT, K1at_KK);
            ft_t2v(vFFT, K1av_AR, tFFT, K1at_AR, -glb_i*1., -2.*hyb_FFT); // G^A(t) G^R(-t) = Theta(-t)e^{2ht} = (-i) i Theta(-t)e^{-(-2h)t}
            ft_t2v(vFFT, K1av_RA, tFFT, K1at_RA, glb_i*1., 2.*hyb_FFT); // G^R(t) G^A(-t) = Theta(t)e^{-2ht} = i (-i) Theta(t)e^{-(2h)t}
            K1av_3 = (K1av_KK + K1av_AR + K1av_RA) * (-glb_i*Uin*Uin/4.); // prefactor: (\pm U/2)^2 1/i
        }
        else {
            K1at_3 = (K1at_KK + K1at_AR + K1at_RA) * (-glb_i*Uin*Uin/4.); // prefactor: (\pm U/2)^2 1/i
            ft_t2v(vFFT, K1av_3, tFFT, K1at_3); // Fourier trafo, component 3
        }
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_K1a.h5", {"v","K1a_1_R","K1a_1_I","K1a_3_R","K1a_3_I"}, {vFFT,K1av_1.real(),K1av_1.imag(),K1av_3.real(),K1av_3.imag()});
    }
    if(PiPflag_FFT) {
        /// ------------------- K1p COMPONENT 1 --------------------- ///
        // left sum of indices is even, right sum of indices is odd, retarded component of bosonic self-energy
        // vectors for \K1p(t) terms arising in SOPT for component 1
        cvec K1pt_RK(nFFT);
        // only one term but appears twice
        K1pt_RK = GtR * GtK; // G^R(t) G^K(t)
        // resulting K1p component 1 in time and frequency
        cvec K1pt_1(nFFT), K1pv_1(nFFT);
        K1pt_1 = K1pt_RK * (-glb_i*Uin*Uin/2.); // prefactor: (\pm U/2)^2 1/i twice
        ft_t2v(vFFT, K1pv_1, tFFT, K1pt_1); // Fourier trafo, component 3
        /// ------------------- K1p COMPONENT 3 --------------------- ///
        // left sum of indices is odd, right sum of indices is odd, Keldysh component of bosonic self-energy
        // vectors for \K1a(t) terms arising in SOPT for component 3
        cvec K1pt_KK(nFFT), K1pt_RR(nFFT), K1pt_AA(nFFT);
        // individual terms
        K1pt_KK = GtK * GtK; // G^K(t) G^K(t)
        K1pt_RR = GtR * GtR; // G^R(t) G^R(t)
        K1pt_AA = GtAc.conj() * GtAc.conj(); // G^A(t) G^A(t)
        // resulting K1p component 3 in time and frequency
        cvec K1pt_3(nFFT), K1pv_3(nFFT);
        if(tailflag_FFT) {
            cvec K1pv_KK(nFFT), K1pv_RR(nFFT), K1pv_AA(nFFT);
            // Fourier trafos
            ft_t2v(vFFT, K1pv_KK, tFFT, K1pt_KK);
            ft_t2v(vFFT, K1pv_RR, tFFT, K1pt_RR, -glb_i*1., 2.*hyb_FFT); // G^R(t)^2 = -Theta(t)e^{-2ht} = (-i) (-i) Theta(t)e^{-(2h)t}
            ft_t2v(vFFT, K1pv_AA, tFFT, K1pt_AA, glb_i*1., -2.*hyb_FFT); // G^A(t)^2 = -Theta(-t)e^{2ht} = (i) i Theta(-t)e^{-(-2h)t}
            K1pv_3 = (K1pv_KK + K1pv_RR + K1pv_AA) * (-glb_i*Uin*Uin/4.); // prefactor: (\pm U/2)^2 1/i
        }
        else {
            K1pt_3 = (K1pt_KK + K1pt_RR + K1pt_AA) * (-glb_i*Uin*Uin/4.); // prefactor: (\pm U/2)^2 1/i
            ft_t2v(vFFT, K1pv_3, tFFT, K1pt_3); // Fourier trafo, component 3
        }
        if(writeflag_FFT) write_h5_rvecs("SOPT_FFT_K1p.h5", {"v","K1p_1_R","K1p_1_I","K1p_3_R","K1p_3_I"}, {vFFT,K1pv_1.real(),K1pv_1.imag(),K1pv_3.real(),K1pv_3.imag()});
    }
}


void SOPT_FFT_SelfEnergy(SelfEnergy<comp>& SEout, rvec& v, SelfEnergy<comp>& Gin, const double Uin) {
    /// ------------------- READ INPUT --------------------- ///
    const int N = v.size(); // vector size, all vector sizes must be equal
    // linear grid is required for Fourier trafo
    const double diff = v[1] - v[0];
    double temp_diff; // temporary difference
    for (int i = 1; i < N; ++i) {
        temp_diff = v[i] - v[i-1];
        if(abs(temp_diff - diff) > 1e-14) {
            cout << "Error in SOPT_SelfEnergy: Frequency grid must be linear." << endl;
            break;
        }
    }
    cvec GvR (N), GvK (N); // separate vectors for retarded and Keldysh component of input GF
    for (int i = 0; i < N; ++i) {
        GvR[i] = Gin.sval(0, i); // retarded component
        GvK[i] = Gin.sval(1, i); // Keldysh component
    }
    /// ------------------- TRANSFORM INPUT --------------------- ///
    rvec t (N); // vector for real time
    cvec GtR (N), GtK (N), GtAc (N), GtK_(N); // vectors for G^R(t), G^K(t), G^A(t)*=G^R(-t), G^K(-t)
    ft_v2t(t, GtR, v, GvR); // Fourier trafo, retarded component
    ft_v2t(t, GtK, v, GvK); // Fourier trafo, Keldysh component
    flip_adjust(GtAc, GtR); // store G^A(t)*=G^R(-t) for convenience
    flip_adjust(GtK_, GtK); // store G^K(-t) for convenience
    /// ------------------- SELF-ENERGY RETARDED COMPONENT --------------------- ///
    // vectors for \Sigma(t) terms arising in SOPT for retarded component
    cvec SEt_RKK (N), SEt_KAK (N), SEt_KKR (N), SEt_RAR (N);
    // individual terms
    SEt_RKK = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
    SEt_KAK = GtK * GtK * GtR.conj(); // G^K(t) G^K(t) G^A(-t)
    SEt_KKR = GtR * GtK * GtK_; // G^R(t) G^K(t) G^K(-t)
    SEt_RAR = GtR * GtR * GtR.conj(); // G^R(t) G^R(t) G^A(-t)
    // resulting retarded self-energy in time and frequency
    cvec SEtR (N), SEvR (N);
    SEtR = (SEt_RKK + SEt_KAK + SEt_KKR + SEt_RAR) * (Uin*Uin/4.); // prefactor: -(\pm U/2/i)^2
    ft_t2v(v, SEvR, t, SEtR); // Fourier trafo, retarded component
    /// ------------------- SELF-ENERGY KELDYSH COMPONENT --------------------- ///
    // vectors for \Sigma(t) terms arising in SOPT for Keldysh component
    cvec SEt_RKA (N), SEt_RRK (N), SEt_AKR (N), SEt_AAK (N), SEt_KRA (N), SEt_KAR (N), SEt_KKK (N);
    // individual terms
    SEt_RKA = GtR * GtK * GtR.conj(); // G^R(t) G^K(t) G^A(-t)
    SEt_RRK = GtR * GtR * GtK_; // G^R(t) G^R(t) G^K(-t)
    SEt_AKR = GtAc.conj() * GtK * GtAc; // G^A(t) G^K(t) G^R(-t)
    SEt_AAK = GtAc.conj() * GtAc.conj() * GtK_; // G^A(t) G^A(t) G^K(-t)
    SEt_KRA = GtK * GtR * GtR.conj(); // G^K(t) G^R(t) G^A(-t)
    SEt_KAR = GtK * GtAc.conj() * GtAc; // G^K(t) G^A(t) G^R(-t)
    SEt_KKK = GtK * GtK * GtK_; // G^K(t) G^K(t) G^K(-t)
    // resulting Keldysh self-energy in time and frequency
    cvec SEtK (N), SEvK (N);
    SEtK = (SEt_RKA + SEt_RRK + SEt_AKR + SEt_AAK + SEt_KRA + SEt_KAR + SEt_KKK) * (Uin*Uin/4.); // prefactor: -(\pm U/2/i)^2
    ft_t2v(v, SEvK, t, SEtK); // Fourier trafo, Keldysh component
    /// ------------------- SELF-ENERGY OUTPUT --------------------- ///
    for (int i = 0; i < N; ++i) { // fill results in output self-energy
        SEout.setself(0, i, SEvR[i]); // retarded component
        SEout.setself(1, i, SEvK[i]); // Keldysh component
    }
}

// Comment for other function
// Dyson eq.: G^K = G^R G^A ( \Sigma^K + \Delta^K)
// equilibrium: \Sigma^K = (1+2n_F)(\Sigma^R-\Sigma^A), acc. for \Delta^K
// \Rightarrow G^K = (1-2n_F) G^R G^A [ (\Sigma+\Delta)^R - (\Sigma+\Delta)^A ]
//                 = (1-2n_F) G^R G^A (G_0^A - G_0^R) = (1+2n_F) (G^R-G^A)



#endif // FOURIER_TRAFO_H