#ifndef __PHYSICAL_BASIC__
#define __PHYSICAL_BASIC__

#include <complex>
#include <matrix.h>
#include <approxxpp.h>

double fermi (double x, double mu, double T);

//Compute the retarded and keldysh component of the green's function and the single scale propagator. Returns them as a matrix
//returns (GR, GK, SR, SK)
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//            - distribution: distribution function in the artificial leads
matrix<matrix<std::complex<double> > > green_and_single_vary_distribution_with_Lambda(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, linear_ipol_bin<matrix<complex<double> > >  & distribution, linear_ipol_bin<matrix<complex<double> > > & ddistribution, double taul, double Lambda);

//Compute the retarded and keldysh component of the green's function and the single scale propagator. Returns them as a matrix
//returns (GR, GK, SR, SK)
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//            - distribution: distribution function in the artificial leads
matrix<matrix<std::complex<double> > > green_and_single_various_distribution_functions(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, linear_ipol_bin<matrix<complex<double> > > distribution, double taul, double Lambda);

//Compute the retarded component of the green's function and the single scale propagator. Returns them as a matrix
//returns (GR, SR)
//parameters:
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
matrix<syma<std::complex<double> > > green_and_single_eq_non_ps(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda);


//Compute the retarded component of the green's function. Returns them as a matrix
//returns (GR, SR)
//parameters:
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
syma<std::complex<double> > green_eq_non_ps(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda);

//Compute the retarded and keldysh component of the green's function and the single scale propagator. Returns them as a matrix
//returns (GR, GK, SR, SK)
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//            - shiftL: shift of the left lead
//shiftL, shiftR: due to different densities, there is a shift in the left and right lead, with some interpolation in between, in the electrostatic energy. The shift of the leads is determined by shiftL, shiftR (where a POSITIVE value shifts the lead UP). 
//CAVEAT: The interpolation between these shifts is NOT accounted for automatically and has to be taken care of in the self-energy (which also contains the hamiltonian).
matrix<matrix<std::complex<double> > > green_and_single_electrostatic(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, matrix<double> mus, matrix<double> Ts, double taul, double shiftL, double shiftR, double Lambda);

//Compute the retarded component of the green's function and the single scale propagator. Returns them as a matrix s.t. ret(0) = green, ret(1) = single_scale;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium
matrix<syma<std::complex<double> > > green_and_single_Req_ps(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda);
matrix<syma<std::complex<double> > > green_and_single_Req(double w, syma<std::complex<double> > &ERet, syma<std::complex<double> > &dH0, double h, double taul, double Lambda);
matrix<matrix<std::complex<double> > > green_and_single_R(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, matrix<double> mus, matrix<double> Ts, double taul, double Lambda);
//Compute the retarded component of the green's function;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - self-energy (retarded) ERet (has absorbed the hamiltonian)
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium
syma<std::complex<double> > greenReq(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda);
//The same as previously, but parity symmetry is used and the result is stored in the matrices Ge and Go.
syma<std::complex<double> > greenReq_ps_decomp(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda, syma<std::complex<double> > &Ge, syma<std::complex<double> > &Go);

//Compute the keldysh component of the green's function;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - Retarded green's function GR
//            - self-energy (keldysh) EKel
//            - magnetic field h
//            - coupling to the leads taul
//            - chemical potentials mu
//            - temperatures T
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium
syma<std::complex<double> > greenKeq(double w, syma<std::complex<double> > &GRet, double h, double taul, double mu, double T, double Lambda);

//Use the following flow-parameter: at every site a lead is attached, these leads are slowly turned off -> Sigma gets modified. Use mean temperature and chemical potential as offsets for the additional leads.
//Compute the retarded component of the single-scale propagator;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - Retarded Green's function GRet
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium
syma<std::complex<double> > single_scaleReq(double w, syma<std::complex<double> > &GRet, double taul, double Lambda);
//Use the following flow-parameter: at every site a lead is attached, these leads are slowly turned off -> Sigma gets modified. Use mean temperature and chemical potential as offsets for the additional leads.
//Compute the retarded component of the single-scale propagator;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - Retarded Green's function GRet
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium, explicitely uses parity
syma<std::complex<double> > single_scaleReq_ps(double w, syma<std::complex<double> > &H, double h, double taul, double Lambda);
//The same as previously, but the Green's functions are already known and stored as even and odd components in Ge and Go.
syma<std::complex<double> > single_scaleReq_ps_from_decomp(double w, int N, syma<std::complex<double> > &Ge, syma<std::complex<double> > &Go, double h, double taul, double Lambda);

//Use the following flow-parameter: at every site a lead is attached, these leads are slowly turned off -> Sigma gets modified. Use mean temperature and chemical potential as offsets for the additional leads.
//Compute the keldysh component of the single-scale propagator;
//parameters:
//            - Components i1, i2; 1 corresponds to q and 0 corresponds to cl
//            - frequency w
//            - Hamiltionian H,
//            - retarded Green's function GRet
//            - magnetic field h
//            - coupling to the leads taul
//            - chemical potentials mu
//            - temperatures T
//            - value of the flow parameter Lambda
//            - Band-widths of the artificial leads at each site ArtLeads
//in equilibrium
syma<std::complex<double> > single_scaleKeq(double w, syma<std::complex<double> > &SRet, double h, double taul, double mu, double T);

//Map the above computations for convenient usage later on
syma<std::complex<double> > greenEq       (int i2, int i4, double w, syma<std::complex<double> > &GRet, double h, double taul, double muL, double TL, double Lambda); 
syma<std::complex<double> > single_scaleEq(int i1, int i3, double w, syma<std::complex<double> > &SRet, double h, double taul, double muL, double TL, double Lambda); 

//line to circle
double subst_integral(double x);

//circle to line
double resu_integral(double x);

//measure on circle
double weight_integral(double x);

//lead contribution to self-energy
std::complex<double> SigmaRLead(std::complex<double> w, double taul);
//flow-param-derivative of lead contribution to self-energy
std::complex<double> dSigmaRLead(std::complex<double> w, double taul);

//line to finite line segments
double subst_concatenated (double omega);

//finite line segments to line
double resu_concatenated (double x);

//measure on finite line segments
double weight_concatenated (double x);

//line to finite line segments
double subst_concatenated (double omega, double taul);

//finite line segments to line
double resu_concatenated (double x, double taul);

//measure on finite line segments
double weight_concatenated (double x, double taul);

//line to finite line segments
double subst_concatenated (double omega, double taul, double Lambda);

//finite line segments to line
double resu_concatenated (double y, double taul, double Lambda);

//measure on finite line segments
double weight_concatenated (double y, double taul, double Lambda);

//line to finite line segments
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double subst_concatenated (double omega, matrix<double> &breaks, const double Lambda);

//finite line segments to line
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double resu_concatenated (double x, matrix<double> &breaks, const double Lambda);

//measure on finite line segments
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double weight_concatenated (double x, matrix<double> breaks, double Lambda);

//line to finite line segments
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double subst_concatenated_fast_decay (double omega, matrix<double> &breaks, const double Lambda);

//finite line segments to line
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double resu_concatenated_fast_decay (double x, matrix<double> &breaks, const double Lambda);

//measure on finite line segments
double weight_concatenated_fast_decay (double x, matrix<double> breaks, double Lambda);

#endif
