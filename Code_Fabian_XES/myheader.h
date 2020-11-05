//
//  header.h
//  Multiloop fRG for X-ray-edge singularity
//
//  Created by Fabian Kugler on 10.01.17.
//  Copyright © 2016 Fabian Kugler. All rights reserved.
//

typedef complex<double> Cp;
const static Cp I = Cp(0,1);

typedef vector<Cp> CVec;
typedef vector<CVec> CMat;
typedef vector<CMat> CTen;
typedef vector<double> DVec;
typedef vector<DVec> DMat;

#define round2int(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))
#define round2intUp(f) ((int)(f + 0.9999999999))
#define round2intDown(f) ((int)(f))

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//basic mathematical relations and global variables PiOverBeta, mu_integer

const double PiOverBeta = M_PI/BETA;

double sign(const int x) {
    // sign function
	if (x>=0) return 1.;
    else if (x<0) return -1.;
    return 0.;
}
double Theta(double x) {
    // Heaviside step function
	if (x > 0) return 1.;
    else if (x < 0) return 0.;
    else return 1./2.;
}

int sumindex_dataindex_fer(const int i) {
	// get odd sumindex (fermion) from dataindex
	return 2*i+1;
}
int sumindex_dataindex_bos(const int i) {
	// get even sumindex (boson) from dataindex
	return 2*i;
}
int sumindex_dataindex_ferMat(const int i, const int N_freq) {
	// get odd sumindex (fermion) from dataindex of CMat
	// fermionic frequencies: i=0, 1, ..., 2*N_freq-1 relates to sumindex=-2*N_freq+1, -2*N_freq+3, ..., 2*N_freq-1
	return 2*i+1-2*N_freq;
}

int triangularIndex(int i, int j, int N) {
	// consider symmetric matrix of size N with indices i and j, only store upper triangular part
	// wlog i <= j: l = j + i*N - i*(i+1)/2
	// i=0, j=0,...,N-1: l=0,...,N-1
	// i=1, j=1,...,N-1: l=N,...,2N-2
	// i=2, j=2,...,N-1: l=2N-1,...,3N-3
	// ...
	// i=N-1, j=N-1: l=(N-1)(N+2)/2=(N^2+N)/2-1
	// in total: (N^2+N)/2 values stored
	int i1 = min(i, j); int j1 = max(i, j);
	int l = j1 + N*i1 - (i1*(i1+1))/2;
	return l;
}

// get Matsubara frequency index corresponding to positive value Lambda, rounded up or down or in closest way
int getOddIndexUp(const double Lambda) {
	int oddIndex = round2intUp(Lambda/PiOverBeta);
	if (oddIndex%2 == 0) oddIndex += 1;
	return oddIndex;
}
int getOddIndexDown(const double Lambda) {
	int oddIndex = round2intDown(Lambda/PiOverBeta);
	if (oddIndex%2 == 0) oddIndex -= 1;
	return oddIndex;
}
int getOddIndex(const double Lambda) {
	int index = round2intDown(Lambda/PiOverBeta);
	if (index%2 == 0) index += 1;
	return index;
}

// Matsubara frequency index corresponding to the chemical potential
const int mu_integer = getOddIndexDown(MU);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functors for fermionic propagators

// functor for single-scale propagator depending on Lambda
struct get_SDL {
#if REGULATOR==0
	// for delta regulator, -(beta/2pi)/(iw-xiD) can be precomputed
	CVec SD;
	get_SDL() : SD(2*NAUX) {
		for (unsigned int i = 0; i < SD.size(); ++i) {
			int w = 2*i+1 - 2*NAUX;
			SD[i] = -1./(I*PiOverBeta*(double)w - XID)/(2.*PiOverBeta);
		}
	}
#endif
#if REGULATOR==2 or REGULATOR==3
	// for other regulators, Gd = 1/(iw-xiD) can be precomputed
	CVec GD;
	get_SDL() : GD(2*NAUX) {
		for (unsigned int i = 0; i < GD.size(); ++i) {
			int w = 2*i+1 - 2*NAUX;
			GD[i] = 1./(I*PiOverBeta*(double)w - XID);
		}
	}
#endif
	
	// return single-scale propagator at given frequency and Lambda, for form of functions confer formulae.pdf
	// the effect of delta or Theta functions is (possibly) commented out as it is incorporated by the boundary conditions of Matsubara sums implicitly
	// computation without precomputed values is commented out
	Cp operator()(const int freq_integer, const double Lambda) const {
#if REGULATOR==0
		double factor1 = Theta(Lambda - PiOverBeta*(abs(freq_integer)-1));
		double factor2 = Theta(PiOverBeta*(abs(freq_integer)+1) - Lambda);
		//if(factor1*factor2 != 1.) cout << "OOHHH" << endl;
		return SD[(freq_integer-1)/2 + NAUX] * factor1*factor2;
		//return -1./(I * sign(freq_integer) * Lambda - XID) * factor1*factor2 / (2.*PiOverBeta);
#elif REGULATOR==1
		double factor = Theta(Lambda - PiOverBeta*abs(freq_integer));
		//if(factor != 1.) cout << "OOHHH" << endl;
		const Cp temp = I * sign(freq_integer) * Lambda - XID;
		return I * ( -sign(freq_integer)/( temp * temp ) ) * factor; // Theta function also incorporated by loop boundaries
#elif REGULATOR==2
		const double freq_value = PiOverBeta * (double)freq_integer;
		const double x = pow(abs(freq_value/Lambda), REGEXP); 
		//return - par.regExp * x/Lambda * exp(-x) / (I * freq_value - XID);
		return (- REGEXP * x/Lambda * exp(-x) ) * GD[(freq_integer-1)/2 + NAUX];
#elif REGULATOR==3
		const double freq_value = PiOverBeta * (double)freq_integer;
		const Cp x = 1.-I*OSCFACTOR*sign(freq_integer);
		const double y = Lambda / freq_value;
		//return -2.*y/freq_value*x * exp( -y*y*x ) / (I * freq_value - XID);
		return ( -2.*y/freq_value*x * exp( -y*y*x ) ) * GD[(freq_integer-1)/2 + NAUX];
#else
		return 0.;
#endif
	}
};

// functor for propagator depending on Lambda
struct get_GDL {
#if REGULATOR!=1
	// for regulator except Litim, 1/(iw-xiD) can be precomputed
	CVec GD;
	get_GDL() : GD(2*NAUX) {
		for (unsigned int i = 0; i < GD.size(); ++i) {
			int w = 2*i+1 - 2*NAUX;
			GD[i] = 1./(I*PiOverBeta*(double)w - XID);
		}
	}
#endif
	
	
	Cp operator()(const int freq_integer, const double Lambda) const {
//#if REGULATOR!=0
		const double freq_value = PiOverBeta * (double)freq_integer;
//#endif
// the effect of delta or Theta functions is (possibly) commented out as it is incorporated by the boundary conditions of Matsubara sums implicitly
// computation without precomputed values is commented out
#if REGULATOR==0
		double factor = Theta(abs(freq_value) - Lambda); // incorporated by if conditions
		//if(factor != 1.) cout << "OOHHH" << endl;
		//return 1. / (I*freq_value - par.xiD);// * factor
		return GD[(freq_integer-1)/2 + NAUX] * factor;
#elif REGULATOR==1
		return 1./(I*sign(freq_integer)*max(abs(freq_value),Lambda) - XID);
#elif REGULATOR==2
		const double x = pow(abs(freq_value/Lambda), REGEXP);
		//return ( 1. - exp(-x) ) / (I * freq_value - XID);
		return ( 1. - exp(-x) ) * GD[(freq_integer-1)/2 + NAUX];
#elif REGULATOR==3
		const Cp x = 1.-I*OSCFACTOR*sign(freq_integer);
		double y = Lambda / freq_value;
		//return exp( -y*y*x ) / (I * freq_value - XID);
		return exp( -y*y*x ) * GD[(freq_integer-1)/2 + NAUX];
#else 
		return 0.;
#endif
	}
};

// functor for scale-independent propagator, fully precomputed
struct get_GD {	
	CVec GD;
	get_GD() : GD(2*NAUX) {
		for (unsigned int i = 0; i < GD.size(); ++i) {
			int w = 2*i+1 - 2*NAUX;
			GD[i] = 1./(I*PiOverBeta*(double)w - XID);
		}
	}
	
	Cp operator()(const int freq_integer) const {
		return GD[(freq_integer-1)/2 + NAUX];
	}
};

// create (global) instances of functors
get_GDL getGDL;
get_SDL getSDL;
get_GD getGD;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// initialization and operators for new types CVec, etc.

CVec initCVec(const int N_freq) {
	// zero array of sive N_freq
    CVec Vec(N_freq); return Vec;
}
CMat initCMat(const int N1, const int N2) {
	// zero matrix of size N1*N2
    CMat Mat(N1, CVec(N2)); return Mat;
}
CMat initCMat(const int N1, const int N2, Cp a) {
	// matrix of size N1*N2 with constant value a
    CMat Mat(N1, CVec(N2, a)); return Mat;
}
CTen initCTen(const int N1, const int N2, const int N3) {
	// zero tensor of size N1*N2*N3
	CTen Ten(N1, CMat(N2, CVec(N3))); return Ten;
}

CVec operator+(const CVec& lhs, const CVec& rhs) {
	unsigned int N = lhs.size(); CVec out(N);
	if ( rhs.size() != N ) return out;
	for (unsigned int i = 0; i < N; ++i) out[i] = lhs[i] + rhs[i];
	return out;
}
CMat operator+(const CMat& lhs, const CMat& rhs) {
	unsigned int N1 = lhs.size(); if (N1==0) {CMat out; return out;}
	unsigned int N2 = lhs[0].size(); CMat out(N1, CVec(N2));
	if ( rhs.size() != N1 or rhs[0].size() != N2 ) return out;
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) out[i][j] = lhs[i][j] + rhs[i][j];
	}
	return out;
}
CTen operator+(const CTen& lhs, const CTen& rhs) {
	unsigned int N1 = lhs.size(); if (N1==0) {CTen out; return out;}
	unsigned int N2 = lhs[0].size(); if (N2==0) {CTen out; return out;}
	unsigned int N3 = lhs[0][0].size(); CTen out(N1, CMat(N2, CVec(N3)));
	if ( rhs.size() != N1 or rhs[0].size() != N2 or rhs[0][0].size() != N3 ) return out;
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) {
			for (unsigned int k = 0; k < N3; ++k) out[i][j][k] = lhs[i][j][k] + rhs[i][j][k];
		}
	}
	return out;
}

CVec& operator+=(CVec& lhs, const CVec& rhs) {
	unsigned int N = lhs.size();
	for (unsigned int i = 0; i < N; ++i) lhs[i] += rhs[i];
	return lhs;
}
CMat& operator+=(CMat& lhs, const CMat& rhs) {
	unsigned int N1 = lhs.size(); if (N1==0) return lhs;
	unsigned int N2 = lhs[0].size(); 
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) lhs[i][j] += rhs[i][j];
	}
	return lhs;
}
CTen& operator+=(CTen& lhs, const CTen& rhs) {
	unsigned int N1 = lhs.size(); if (N1==0) return lhs;
	unsigned int N2 = lhs[0].size(); if (N2==0) return lhs;
	unsigned int N3 = lhs[0][0].size(); 
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) {
			for (unsigned int k = 0; k < N3; ++k) lhs[i][j][k] += rhs[i][j][k];
		}
	}
	return lhs;
}

CVec operator*(const CVec& lhs, const double factor) {
	unsigned int N = lhs.size(); CVec out(N);
	for (unsigned int i = 0; i < N; ++i) out[i] = lhs[i] * factor;
	return out;
}
CMat operator*(const CMat& lhs, const double factor) {
	unsigned int N1 = lhs.size(); if (N1==0) {CMat out; return out;}
	unsigned int N2 = lhs[0].size(); CMat out(N1, CVec(N2));
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) out[i][j] = lhs[i][j] * factor;
	}
	return out;
}
CTen operator*(const CTen& lhs, const double factor) {
	unsigned int N1 = lhs.size(); if (N1==0) {CTen out; return out;}
	unsigned int N2 = lhs[0].size(); if (N2==0) {CTen out; return out;}
	unsigned int N3 = lhs[0][0].size(); CTen out(N1, CMat(N2, CVec(N3)));
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) {
			for (unsigned int k = 0; k < N3; ++k) out[i][j][k] = lhs[i][j][k] * factor;
		}
	}
	return out;
}

CVec& operator*=(CVec& lhs, const double factor) {
	unsigned int N = lhs.size();
	for (unsigned int i = 0; i < N; ++i) lhs[i] = lhs[i] * factor;
	return lhs;
}
CMat& operator*=(CMat& lhs, const double factor) {
	unsigned int N1 = lhs.size(); if (N1==0) return lhs;
	unsigned int N2 = lhs[0].size(); 
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) lhs[i][j] = lhs[i][j] * factor;
	}
	return lhs;
}
CTen& operator*=(CTen& lhs, const double factor) {
	unsigned int N1 = lhs.size(); if (N1==0) return lhs;
	unsigned int N2 = lhs[0].size(); if (N2==0) return lhs;
	unsigned int N3 = lhs[0][0].size(); 
	for (unsigned int i = 0; i < N1; ++i) {
		for (unsigned int j = 0; j < N2; ++j) {
			for (unsigned int k = 0; k < N3; ++k) lhs[i][j][k] = lhs[i][j][k] * factor;
		}
	}
	return lhs;
}

// different ways to compare arrays using the first few elements
double compareCVec(const CVec& PiX1, const CVec& PiX2, const int N) {
	double diff = 0.; 
	for (int i = 0; i < N; ++i) diff += abs(PiX1[i] - PiX2[i]);
	return diff;
}
double compareCVec2(const CVec& PiX1, const CVec& PiX2, const int N) {
	double diff = 0.;
	for (int i = 0; i < N; ++i) diff = max(diff, abs(PiX1[i] - PiX2[i]) );
	return diff;
}
double compareCMat(const CMat& Gamma1, const CMat& Gamma2, const int N) {
	double diff = 0.; double mean = N*N;//0.;
	for (int i = 0; i < N; ++i) {
		for (unsigned int j = Gamma2[0].size()-N/2; j < Gamma2[0].size()+N/2; ++j) {
			diff += abs(Gamma1[i][j] - Gamma2[i][j]);
			//mean += abs(Gamma1[i][j] + Gamma2[i][j]) / 2.;
		}
	}
	return diff/mean;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CVertex class for parametrization of a 4-point vertex as suggested by Toschi et al.
// I(\omega, \nu, \bar{\omega}) = const. + Theta(Omega1 - |\bar{\omega}|) K1(\bar{\omega})
// 									+ Theta(Omega2 - |\bar{\omega}|) Theta(Omega2 - |\omega|) K2(\bar{\omega}, \omega)
// 									+ Theta(Omega2 - |\bar{\omega}|) Theta(Omega2 - |\nu|) K2ns(\bar{\omega}, \nu)
// 									+ Theta(Omega2 - |\bar{\omega}|) Theta(Omega2 - |\omega|) Theta(Omega2 - |\nu|) K3ns(\bar{\omega}, \omega, \nu)
// for a symmetric vertex: K2 = K2ns; K3ns symmetric in \omega, \nu -> K3

// move semantics possible in c++11 and more efficient. can be turned off by defining NOMOVE

class CVertex {
public:	
	CVec K1;
	CMat K2, K3;
	CMat K2ns; CTen K3ns; //only for non-symmetric vertex
	
	CVertex(); ~CVertex();
	CVertex(const CVertex&);
	CVertex& operator=(const CVertex&);
    CVertex operator+(const CVertex&) const&;
	CVertex operator*(const double) const&;
#ifndef NOMOVE
	CVertex& operator=(CVertex&&);
	CVertex& operator+(const CVertex&) &&;
	CVertex& operator+(CVertex&&) const&;
	CVertex& operator+(CVertex&&) &&;
	CVertex& operator*(const double) &&;
#endif
	Cp getVal(const int, const int, const int) const;
	Cp getValns(const int, const int, const int) const;
};

CVertex::CVertex() {
	K1 = initCVec(NUMFREQ1);
	K2 = initCMat(NUMFREQ2, 2*NUMFREQ2);
	K3 = initCMat(NUMFREQ3, (2*NUMFREQ3*(2*NUMFREQ3+1))/2);
}
CVertex::CVertex(const CVertex& rhs) {
	K1 = rhs.K1; K2 = rhs.K2; K3 = rhs.K3;
}
CVertex::~CVertex() {}

CVertex& CVertex::operator=(const CVertex& rhs) {
	K1 = rhs.K1; K2 = rhs.K2; K3 = rhs.K3;
	return *this;
}
CVertex CVertex::operator+(const CVertex& rhs) const& {
	CVertex out;
	out.K1 = K1 + rhs.K1; out.K2 = K2 + rhs.K2; out.K3 = K3 + rhs.K3;
    return out;
}
CVertex CVertex::operator*(const double factor) const& {
    CVertex out;
	out.K1 = K1 * factor; out.K2 = K2 * factor; out.K3 = K3 * factor;
    return out;
}
#ifndef NOMOVE
CVertex& CVertex::operator=(CVertex&& rhs) {
	K1 = move(rhs.K1); K2 = move(rhs.K2); K3 = move(rhs.K3);
	return *this;
}
CVertex& CVertex::operator+(CVertex&& rhs) const& {
	rhs.K1 += K1; rhs.K2 += K2; rhs.K3 += K3;
	return rhs;
}
CVertex& CVertex::operator+(CVertex&& rhs) && {
	K1 += rhs.K1; K2 += rhs.K2; K3 += rhs.K3; 
    return *this;
}
CVertex& CVertex::operator+(const CVertex& rhs) && {
	K1 += rhs.K1; K2 += rhs.K2; K3 += rhs.K3; 
    return *this;
}
CVertex& CVertex::operator*(const double factor) && {
	K1 *= factor; K2 *= factor; K3 *= factor; 
    return *this;
}
#endif

// get value, straightforward sum of Theta functions, using symmetry
Cp CVertex::getVal(const int w, const int v, const int w_bar) const {
	Cp out(-UVAL, 0.);
	int k = abs(w_bar)/2;
	
	if ( k >= NUMFREQ1 ) return out;
	
	if (w_bar >= 0) out += K1[k];
	else out += conj(K1[k]);
	if ( k >= NUMFREQ2 ) return out;
	
	int j = (v-1)/2;
	if ( j < NUMFREQ2 and j >= -NUMFREQ2 ) {
		if (w_bar==0) out += ( K2[0][j + NUMFREQ2] + conj(K2[0][NUMFREQ2-j-1]) )/2.;
		else if (w_bar > 0) out += K2[k][j + NUMFREQ2];
		else out += conj(K2[k][NUMFREQ2-j-1]);
	}
	
	int i = (w-1)/2; 
	if ( i < NUMFREQ2 and i >= -NUMFREQ2 ) {
		if (w_bar==0) out += ( K2[0][i + NUMFREQ2] + conj(K2[0][NUMFREQ2-i-1]) )/2.;
		else if (w_bar > 0) out += K2[k][i + NUMFREQ2];
		else out += conj(K2[k][NUMFREQ2-i-1]);
	}

	if ( k >= NUMFREQ3 or j >= NUMFREQ3 or j < -NUMFREQ3 or i >= NUMFREQ3 or i < -NUMFREQ3) return out;
	if (w_bar==0) out += ( K3[0][triangularIndex(i+NUMFREQ3, j+NUMFREQ3, 2*NUMFREQ3)] + conj(K3[0][triangularIndex(NUMFREQ3-i-1, NUMFREQ3-j-1, 2*NUMFREQ3)]) )/2.;
	else if (w_bar > 0) out += K3[k][triangularIndex(i+NUMFREQ3, j+NUMFREQ3, 2*NUMFREQ3)];
	else out += conj(K3[k][triangularIndex(NUMFREQ3-i-1, NUMFREQ3-j-1, 2*NUMFREQ3)]); 
	
	return out;
}

// get value, straightforward sum of Theta functions, nonsymmetric vertex
Cp CVertex::getValns(const int w, const int v, const int w_bar) const {
	Cp out;
	int k = abs(w_bar)/2;
	
	if ( k >= NUMFREQ1 ) return out;
	
	if (w_bar >= 0) out += K1[k];
	else out += conj(K1[k]);
	if ( k >= NUMFREQ2 ) return out;
	
	int j = (v-1)/2;
	if ( j < NUMFREQ2 and j >= -NUMFREQ2 ) {
		if (w_bar==0) out += ( K2ns[0][j + NUMFREQ2] + conj(K2ns[0][NUMFREQ2-j-1]) )/2.;
		else if (w_bar > 0) out += K2ns[k][j + NUMFREQ2];
		else out += conj(K2ns[k][NUMFREQ2-j-1]);
	}
	
	int i = (w-1)/2; 
	if ( i < NUMFREQ2 and i >= -NUMFREQ2 ) {
		if (w_bar==0) out += ( K2[0][i + NUMFREQ2] + conj(K2[0][NUMFREQ2-i-1]) )/2.;
		else if (w_bar > 0) out += K2[k][i + NUMFREQ2];
		else out += conj(K2[k][NUMFREQ2-i-1]);
	}

	if ( k >= NUMFREQ3 or j >= NUMFREQ3 or j < -NUMFREQ3 or i >= NUMFREQ3 or i < -NUMFREQ3) return out;
	if (w_bar==0) out += ( K3ns[0][i+NUMFREQ3][j+NUMFREQ3] + conj(K3ns[0][NUMFREQ3-i-1][NUMFREQ3-j-1]) )/2.;
	else if (w_bar > 0) out += K3ns[k][i+NUMFREQ3][j+NUMFREQ3];
	else out += conj(K3ns[k][NUMFREQ3-i-1][NUMFREQ3-j-1]); 
	
	return out;
}

// get (symmetric) IT from (nonsymmetric) IL by IT(w, v, wbar) = IL(w, v, wbar) + IL(v, w, wbar)
void symmetrizeCVertex(const CVertex& IL, CVertex& IT) {
	IT.K1 = IL.K1*2.;
	IT.K2 = IL.K2 + IL.K2ns;
	for (int k = 0; k < NUMFREQ3; ++k) {
		for (int i = 0; i < 2*NUMFREQ3; ++i) {
			for (int j = i; j < 2*NUMFREQ3; ++j) IT.K3[k][triangularIndex(i, j, 2*NUMFREQ3)] = IL.K3ns[k][i][j] + IL.K3ns[k][j][i];
		}
	}
}

// get (symmetric) IT from (nonsymmetric) IL and (symmetric) IC by IT(w, v, wbar) = IL(w, v, wbar) + IL(v, w, wbar) + IC(w, v, wbar)
void symmetrizeCVertex_withC(const CVertex& IL, const CVertex& IC, CVertex& IT) {
	IT.K1 = IL.K1*2. + IC.K1;
	IT.K2 = IL.K2 + IL.K2ns + IC.K2;
	for (int k = 0; k < NUMFREQ3; ++k) {
		for (int i = 0; i < 2*NUMFREQ3; ++i) {
			for (int j = i; j < 2*NUMFREQ3; ++j) IT.K3[k][triangularIndex(i, j, 2*NUMFREQ3)] = IL.K3ns[k][i][j] + IL.K3ns[k][j][i] + IC.K3[k][triangularIndex(i, j, 2*NUMFREQ3)];
		}
	}
}

// get IT = IL + IR, only needed for check
void combineCVertex(const CVertex& IL, const CVertex& IR, CVertex& IT) {
	IT.K1 = IL.K1 + IR.K1;
	IT.K2 = IL.K2 + IR.K2;
	for (int k = 0; k < NUMFREQ3; ++k) {
		for (int i = 0; i < 2*NUMFREQ3; ++i) {
			for (int j = i; j < 2*NUMFREQ3; ++j) IT.K3[k][triangularIndex(i, j, 2*NUMFREQ3)] = IL.K3ns[k][i][j] + IR.K3ns[k][i][j];
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////

// old write functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void writeVector(string path, const CVec& Vec, const int N) {
	// write real part of matrix to file
    ofstream Data;
    Data.open(path);
    for (int i = 0; i < N; i++) {
        Data << Vec[i].real() << " " << Vec[i].imag() << " ";
        Data << endl;
    }
    Data.close();
}

void writeMatReal(string path, const CMat& Mat, const int N1, const int N2) {
	// write real part of matrix to file
    ofstream Data;
    Data.open(path);
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) Data << Mat[i][j].real() << " ";
        Data << endl;
    }
    Data.close();
}

void writeMatImag(string path, const CMat& Mat, const int N1, const int N2) {
	// write imaginary part of matrix to file
    ofstream Data;
    Data.open(path);
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) Data << Mat[i][j].imag() << " ";
        Data << endl;
    }
    Data.close();
}

void readMatReal(string path, CMat& Mat, const int N1, const int N2) {
	// read real part of matrix from file
	ifstream Data;
    Data.open(path);
	for (int i = 0; i < N1; i++) {
		for (int j = 0; j < N2; j++) {
			double x;
			Data >> x;
			Mat[i][j] += x;
		}
	}
    Data.close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// state class, which contains all vertices
// Ip: vertex irreducible in parallel lines, corresponds to gamma_a, vertex reducible in antiparallel lines
// Iap: vertex irreducible in antiparallel lines, corresponds to gamma_p, vertex reducible in parallel lines
// move semantics as above, for different MODEs confer main

class State {
public:
	CVertex Ip, Iap;   

	CVec PiX, PiY; // irrel. for purely fermionic system 
	CMat GammaCDX, GammaCDY; // irrel. for purely fermionic system 

	State(); ~State();
	State(const State&);
	
	State& operator=(const State&);
	State operator+(const State&) const&;
	State operator*(const double) const&;
#ifndef NOMOVE 
	State& operator=(State&&);
	State& operator+(const State&) &&;
	State& operator+(State&&) const&;
	State& operator+(State&&) &&;
	State& operator*(const double) &&;
#endif

	Cp getGX(const int) const; // irrel. for purely fermionic system 
	Cp getGY(const int) const; // irrel. for purely fermionic system 
	Cp getGammaCDX(const int, const int) const; // irrel. for purely fermionic system 
	Cp getGammaCDY(const int, const int) const; // irrel. for purely fermionic system 
	
	void getCorr(CVec&) const; // get particle-hole susceptibility ("correlator") from vertices 
	void getCorr(CVec&, const double) const;
	void getCorrSD(CVec&) const; // irrel. for purely fermionic system 
	void getCorrGamma3(CVec&) const; // irrel. for purely fermionic system 
	void writeVec(string, const CVec&) const;
    void write(string) const;
};

State::State() : Ip(), Iap() {
#if MODE>1
	PiX = initCVec(NFREQ1); PiY = initCVec(NFREQ1); // irrel. for purely fermionic system 
	GammaCDX = initCMat(NFREQ2, 2*NFREQ2, 1.); GammaCDY = initCMat(NFREQ2, 2*NFREQ2, 1.); // irrel. for purely fermionic system 
#endif
}

State::State(const State& rhs) {
	Ip = CVertex(rhs.Ip); Iap = CVertex(rhs.Iap);
#if MODE>1 // irrel. for purely fermionic system 
	PiX = rhs.PiX; PiY = rhs.PiY;
	GammaCDX = rhs.GammaCDX; GammaCDY = rhs.GammaCDY;
#endif
}
State::~State() {}
State& State::operator=(const State& rhs) {
	Ip = rhs.Ip; Iap = rhs.Iap;
#if MODE>1 // irrel. for purely fermionic system 
	PiX = rhs.PiX; PiY = rhs.PiY;
	GammaCDX = rhs.GammaCDX; GammaCDY = rhs.GammaCDY;
#endif
	return *this;
}
State State::operator+(const State& rhs) const& {
	State out(*this);
	out.Ip = Ip + rhs.Ip; out.Iap = Iap + rhs.Iap;    
#if MODE>1
	out.PiX = PiX + rhs.PiX; out.PiY = PiY + rhs.PiY;
	out.GammaCDX = GammaCDX + rhs.GammaCDX; out.GammaCDY = GammaCDY + rhs.GammaCDY;
#endif
    return out;
}
State State::operator*(const double factor) const& {
    State out(*this);
	out.Ip = Ip * factor; out.Iap = Iap * factor;    
#if MODE>1
	out.PiX = PiX * factor; out.PiY = PiY * factor;
	out.GammaCDX = GammaCDX * factor; out.GammaCDY = GammaCDY * factor;
#endif
    return out;
}

#ifndef NOMOVE
State& State::operator=(State&& rhs) {
	Ip = move(rhs.Ip); Iap = move(rhs.Iap);
#if MODE>1
	PiX = move(rhs.PiX); PiY = move(rhs.PiY);
	GammaCDX = move(rhs.GammaCDX); GammaCDY = move(rhs.GammaCDY);
#endif
	return *this;
}

State& State::operator+(State&& rhs) const& {
	rhs.Ip = Ip + rhs.Ip; rhs.Iap = Iap + rhs.Iap;
#if MODE>1
	rhs.PiX = PiX + rhs.PiX; rhs.PiY = PiY + rhs.PiY;
	rhs.GammaCDX = GammaCDX + rhs.GammaCDX; rhs.GammaCDY = GammaCDY + rhs.GammaCDY;
#endif
	return rhs;
}
	
State& State::operator+(State&& rhs) && {
	Ip = Ip + rhs.Ip; Iap = Iap + rhs.Iap;    
#if MODE>1
	PiX = PiX + rhs.PiX; PiY = PiY + rhs.PiY;
	GammaCDX = GammaCDX + rhs.GammaCDX; GammaCDY = GammaCDY + rhs.GammaCDY;
#endif
    return *this;
}
	
State& State::operator+(const State& rhs) && {
	Ip = Ip + rhs.Ip; Iap = Iap + rhs.Iap;    
#if MODE>1
	PiX = PiX + rhs.PiX; PiY = PiY + rhs.PiY;
	GammaCDX = GammaCDX + rhs.GammaCDX; GammaCDY = GammaCDY + rhs.GammaCDY;
#endif
    return *this;
}

State& State::operator*(const double factor) && {
	Ip = Ip * factor; Iap = Iap * factor;    
#if MODE>1
	PiX = PiX * factor; PiY = PiY * factor;
	GammaCDX = GammaCDX * factor; GammaCDY = GammaCDY * factor;
#endif
    return *this;
}
#endif

// self-energy stored as CVec, approximate as constant for frequencies exceeding stored ones, use \Pi(-w) = \Pi(w)^*, dressed propagator via Dyson equation
Cp State::getGX(const int w_bar) const {
	if (not WITH_BOSONIC_SELFENERGIES) return -UXVAL;
	else {
		Cp Pi;
		int i = min( abs(w_bar)/2, NFREQ1-1 );    
		if (w_bar >= 0) Pi = PiX[i];
		else Pi = conj(PiX[i]);
		return -UXVAL / ( 1. + Pi * UXVAL );
	}
}
Cp State::getGY(const int w_bar) const {
	if (not WITH_BOSONIC_SELFENERGIES) return -UYVAL;
	else {
		Cp Pi;
		int i = min( abs(w_bar)/2, NFREQ1-1 );    
		if (w_bar >= 0) Pi = PiY[i];
		else Pi = conj(PiY[i]);
		return -UYVAL / ( 1. + Pi * UYVAL );
	}
}
// 3-point vertex stored as CMat, approx. as constant for..., use \Gamma(w,-w_bar)=\Gamma(-w,w_bar)^*, note: for this to hold also for \Gamma^{\bar{c} \bar{d} \psi}, one needs to include the factor 1/i
// freq_value_C=1,3,...		freq_index_C=0,1,...,N_freq_C-1		freq_index_C+N_freq_C=N_freq_C,...,2N_freq_C-1		-freq_index_C-1+N_freq_C=N_freq_C-1,...,0
// freq_value_C=-1,-3,...	freq_index_C=-1,...,-N_freq_C		freq_index_C+N_freq_C=N_freq_C-1,...,0				-freq_index_C-1+N_freq_C=N_freq_C,...,2N_freq_C-1 
Cp State::getGammaCDX(const int w, const int w_bar) const {  
		int i = min( abs(w_bar)/2, NFREQ2-1 );   
	int j = (w-1)/2;
	if ( j >= NFREQ2 ) j = NFREQ2-1;    
	else if ( j < -NFREQ2 ) j = -NFREQ2; 
	if (w_bar >= 0) return GammaCDX[i][j + NFREQ2];
	else return conj(GammaCDX[i][NFREQ2 - j - 1]);   
}
Cp State::getGammaCDY(const int w, const int w_bar) const {  
	int i = min( abs(w_bar)/2, NFREQ2-1 );   
	int j = (w-1)/2;
	if ( j >= NFREQ2 ) j = NFREQ2-1;    
	else if ( j < -NFREQ2 ) j = -NFREQ2; 
	if (w_bar >= 0) return GammaCDY[i][j + NFREQ2];
	else return conj(GammaCDY[i][NFREQ2 - j - 1]);   
}

void State::getCorr(CVec& PiX) const {
// basic relation between four-point correlator and four-point vertex in purely fermionic theory
// \Pi_{\bar{\omega}} = \DimInt{\omega} G^d_{\omega-\bar{\omega}} G^c_{\omega}
// + \DDimInt{\omega \nu} G^d_{\omega-\bar{\omega}} G^d_{\nu-\bar{\omega}} G^c_{\omega} G^c_{\nu} \tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\omega-\bar{\omega}, \omega, \nu, \nu-\bar{\omega}}
	#pragma omp parallel for
	for (unsigned int i = 0; i < PiX.size(); ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		Cp sum1;
		Cp sum2;
		for (int w = -mu_integer; w <= mu_integer; w+=2) {
			sum1 += getGD(w-w_bar) * sign(w);
			for (int v = -mu_integer; v <= mu_integer; v+=2) {
				Cp gdccd = Ip.getVal(w-w_bar, v-w_bar, w_bar) + Iap.getVal(w-w_bar, v-w_bar, w+v-w_bar) + UVAL; //d_freq corrected
				sum2 += getGD(w-w_bar) * getGD(v-w_bar) * sign(w) * sign(v) * gdccd;
			}
		}
		PiX[i] = sum1 * I*(-PiOverBeta) - sum2 * (PiOverBeta*PiOverBeta) ;
	}
}

void State::getCorr(CVec& PiX, const double Lambda) const {
// basic relation between four-point correlator and four-point vertex in purely fermionic theory
// \Pi_{\bar{\omega}} = \DimInt{\omega} G^d_{\omega-\bar{\omega}} G^c_{\omega}
// + \DDimInt{\omega \nu} G^d_{\omega-\bar{\omega}} G^d_{\nu-\bar{\omega}} G^c_{\omega} G^c_{\nu} \tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\omega-\bar{\omega}, \omega, \nu, \nu-\bar{\omega}}
	#pragma omp parallel for
	for (unsigned int i = 0; i < PiX.size(); ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		Cp sum1;
		Cp sum2;
		for (int w = -mu_integer; w <= mu_integer; w+=2) {
			sum1 += getGDL(w-w_bar, Lambda) * sign(w);
			for (int v = -mu_integer; v <= mu_integer; v+=2) {
				Cp gdccd = Ip.getVal(w-w_bar, v-w_bar, w_bar) + Iap.getVal(w-w_bar, v-w_bar, w+v-w_bar) + UVAL; //d_freq corrected
				sum2 += getGDL(w-w_bar, Lambda) * getGDL(v-w_bar, Lambda) * sign(w) * sign(v) * gdccd;
			}
		}
		PiX[i] = sum1 * I*(-PiOverBeta) - sum2 * (PiOverBeta*PiOverBeta) ;
	}
}

void State::getCorrSD(CVec& PiX) const {
// correlator from Schwinger-Dyson equation
// 	\Pi_{\bar{\omega}} = \DimInt{\omega} G^d_{\omega-\bar{\omega}} G^c_{\omega} \tilde{\Gamma}^{\bar{c} d \gamma}_{\omega, \omega-\bar{\omega}, \bar{\omega}}
	#pragma omp parallel for
	for (unsigned int i = 0; i < PiX.size(); ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		Cp sum;
		for (int w = -mu_integer; w <= mu_integer; w+=2) sum += getGD(w-w_bar) * sign(w) * getGammaCDX(w-w_bar, w_bar); //d_freq corrected
		PiX[i] = sum * I*(-PiOverBeta);
	}
}

void State::getCorrGamma3(CVec& PiX) const {
// correlator-vertex relation after HS trafo, where the fermionic four-point vertex is neglected
// \Pi_{\bar{\omega}} = \DimInt{\omega} G^d_{\omega-\bar{\omega}} G^c_{\omega}
// + \DDimInt{\omega \nu} G^d_{\omega-\bar{\omega}} G^d_{\nu-\bar{\omega}} G^c_{\omega} G^c_{\nu} * (
// \tilde{\Gamma}^{\bar{c} d \chi}_{\omega, \omega-\bar{\omega}, \bar{\omega}} \tilde{\Gamma}^{\bar{c} d \chi}_{\nu, \nu-\bar{\omega}, \bar{\omega}} G^{\chi}_{\bar{\omega}}
// \tilde{\Gamma}^{\bar{c} \bar{d} \psi}_{\omega, \nu-\bar{\omega}, \omega+\nu-\bar{\omega}} \tilde{\Gamma}^{\bar{c} \bar{d} \psi}_{\nu, \omega\bar{\omega}, \omega+\nu-\bar{\omega}} G^{\psi}_{\omega+\nu-\bar{\omega}} )

// second summand = G^{\chi}_{\bar{\omega}} ( \DimInt{\omega} G^d_{\omega-\bar{\omega}} G^c_{\omega}  \tilde{\Gamma}^{\bar{c} d \chi}_{\omega, \omega-\bar{\omega}, \bar{\omega}} )^2
	bool simplify_integral = true;
	if (simplify_integral) {
		#pragma omp parallel for
		for (unsigned int i = 0; i < PiX.size(); i++) {
			const int w_bar = sumindex_dataindex_bos(i);
			Cp sum0; // = 0.;
			Cp sum1; // = 0.;
			Cp sum2; // = 0.;

			Cp gx = getGX(w_bar);
			for (int w = -mu_integer; w <= mu_integer; w+=2) {
				sum0 += getGD(w - w_bar) * sign(w);
				sum1 += getGD(w - w_bar) * sign(w) * getGammaCDX(w-w_bar, w_bar); //d_freq corrected
				for (int v = -mu_integer; v <= mu_integer; v+=2) {
					Cp gy = getGY(w + v - w_bar);
					Cp gcy1 = getGammaCDY(w-w_bar, w + v - w_bar); //d_freq corrected
					Cp gcy2 = getGammaCDY(v-w_bar, w + v - w_bar); //d_freq corrected
					sum2 += getGD(w - w_bar) * getGD(v - w_bar) * sign(w) * sign(v) * gcy1 * gy * gcy2;
				}
			}
			PiX[i] = sum0 * I*(-PiOverBeta) - gx*sum1*sum1 * (PiOverBeta*PiOverBeta) - sum2 * (PiOverBeta*PiOverBeta);
		}
	}
	else {
		#pragma omp parallel for
		for (unsigned int i = 0; i < PiX.size(); ++i) {
			const int w_bar = sumindex_dataindex_bos(i);
			Cp sum1; // = 0.;
			Cp sum2; // = 0.;
			
			for (int w = -mu_integer; w <= mu_integer; w+=2) {
				sum1 += getGD(w - w_bar) * sign(w);
				for (int v = -mu_integer; v <= mu_integer; v+=2) {
					Cp gdccdHS = getGX(w_bar)*getGammaCDX(w-w_bar, w_bar)*getGammaCDX(v-w_bar, w_bar) + getGY(w + v - w_bar)*getGammaCDY(w-w_bar, w + v - w_bar)*getGammaCDY(v-w_bar, w + v - w_bar); //d_freq corrected
					sum2 += getGD(w - w_bar) * sign(w) * gdccdHS * getGD(v - w_bar) * sign(v);
				}
			}
			PiX[i] = sum1 * I*(-PiOverBeta) - sum2 * (PiOverBeta*PiOverBeta);
		}
	}
}

// write complex data of vector into file, top lines contain most important parameters
void State::writeVec(string path, const CVec& Vec) const {
    ofstream Data; Data.open(path);
	Data << NUMFREQ1 << " " << BETA << endl << MU << " " << XID << endl << UVAL << " " << UVAL << endl << UVAL << " " << REGULATOR << endl << REGULATOR << " " << REGULATOR << endl;
    for (unsigned int i = 0; i < Vec.size(); ++i) Data << Vec[i].real() << " " << Vec[i].imag() << endl;
    Data.close();
}

// write vertex functions of state into files
void State::write(string suffix) const {
#if MODE<3
	CVec Pi = initCVec(400); getCorr(Pi);
	writeVec("C_" + suffix, Pi);
	writeVector("IpK1_" + suffix, Ip.K1, NUMFREQ1);
#endif
#if MODE>0	
	//writeVector("IpK1_" + suffix, Ip.K1, NUMFREQ1);
	writeVector("IapK1_" + suffix, Iap.K1, NUMFREQ1);
	writeMatReal("IpK2_" + suffix, Ip.K2, NUMFREQ2/2, NUMFREQ2);
	writeMatReal("IapK2_" + suffix, Iap.K2, NUMFREQ2/2, NUMFREQ2);
	writeMatReal("IpK3_" + suffix, Ip.K3, NUMFREQ3/2, NUMFREQ3*NUMFREQ3);
	writeMatReal("IapK3_" + suffix, Iap.K3, NUMFREQ3/2, NUMFREQ3*NUMFREQ3);
#endif
#if MODE>1
	writeVec("Pi_" + suffix, PiX);
#endif
}

double getError(const State& s1, const State& s2) {
	//deviation between different states, needed for adaptive ODE solver
	const int N = 20; 
#if MODE>2
	return compareCVec2(s1.PiX, s2.PiX, N);
#else
	CVec PiX1(N); CVec PiX2(N); 
	s1.getCorr(PiX1); s2.getCorr(PiX2);
	return compareCVec2(PiX1, PiX2, N);
	/*
	double err1 = max( compareCVec(s1.Ip.K1, s2.Ip.K1, N), compareCVec(s1.Iap.K1, s2.Iap.K1, N) );
	double err2 = max( compareCMat(s1.Ip.K2, s2.Ip.K2, N), compareCMat(s1.Iap.K2, s2.Iap.K2, N) );
	double err3 = max( compareCMat(s1.Ip.K3, s2.Ip.K3, N), compareCMat(s1.Iap.K3, s2.Iap.K3, N) );
	return max(err1, max(err2, err3) );*/
#endif
}

void writeFlow(string suffix, const DVec& l0, const DVec& L0, const CVec& C0, const CVec& C1, const CVec& C2, const CVec& C3, const CVec& C4) {
    ofstream Data; Data.open("Flow_" + suffix);
    for (unsigned int i = 0; i < l0.size(); ++i) Data << l0[i] << " " << L0[i] << " " << C0[i].real() << " " << C1[i].real() << " " << C2[i].real() << " " << C3[i].real() << C4[i].real() <<  endl;
    Data.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// coefficients for ODE solvers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fehlberg(DVec& c, DMat& a, DVec& b, DVec& d) {
	c[0] = 0.; c[1] = 1./4.; c[2] = 3./8.; c[3] = 12./13.; c[4] = 1.; c[5] = 1./2.;
	b[0] = 16./135.; b[1] = 0.; b[2] = 6656./12825.; b[3] = 28561./56430.; b[4] = -9./50.; b[5] = 2./55.;
	d[0] = 25./216.; d[1] = 0.; d[2] = 1408./2565.; d[3] = 2197./4104.; d[4] = -1./5.; d[5] = 0.;
	
	a[1][0] = 1./4.;
	a[2][0] = 3./32.; a[2][1] = 9./32.;
	a[3][0] = 1932./2197.; a[3][1] = -7200./2197.; a[3][2] = 7296./2197.;
	a[4][0] = 439./216.; a[4][1] = -8.; a[4][2] = 3680./513.; a[4][3] = -845./4104.;
	a[5][0] = -8./27.; a[5][1] = 2.; a[5][2] = -3544./2565.; a[5][3] = 1859./4104.; a[5][4] = -11./40.;
}

void cashKarp(DVec& c, DMat& a, DVec& b, DVec& d) {
	c[0] = 0.; c[1] = 1./5.; c[2] = 3./10.; c[3] = 3./5.; c[4] = 1.; c[5] = 7./8.;
	b[0] = 37./378.; b[1] = 0.; b[2] = 250./621.; b[3] = 125./594.; b[4] = 0.; b[5] = 512./1771.;
	d[0] = 2825./27648.; d[1] = 0.; d[2] = 18575./48384.; d[3] = 13525./55296.; d[4] = 277./14336.; d[5] = 1./4.;
	
	a[1][0] = 1./5.;
	a[2][0] = 3./40.; a[2][1] = 9./40.;
	a[3][0] = 3./10.; a[3][1] = -9./10.; a[3][2] = 6./5.;
	a[4][0] = -11./54.; a[4][1] = 5./2.; a[4][2] = -70./27.; a[4][3] = 35./27.;
	a[5][0] = 1631./55296.; a[5][1] = 175./512.; a[5][2] = 575./13824.; a[5][3] = 44275./110592.; a[5][4] = 253./4096.;
}

void adamsB(DVec& c, const int order) {
	if (order == 1) {
		c[0] = 1.;
	}
	else if (order == 2) {
		c[0] = -1./2.; c[1] = 3./2.; 
	}
	else if (order == 3) {
		c[0] = 5./12.; c[1] = -4./3.; c[2] = 23./12.; 
	}
	else if (order == 4) {
		c[0] = -3./8.; c[1] = 37./24.; c[2] = -59./24.; c[3] = 55./24.; 
	}
	else if (order == 5) {
		c[0] = 251./720.; c[1] = -637./360.; c[2] = 109./30.; c[3] = -1387./360.; c[4] = 1901./720.; 
	}
	else if (order == 6) {
		 c[0] = -475./1440.; c[1] = 2877./1440.; c[2] = -7298/1440.; c[3] = 9982./1440.; c[4] = -7923./1440.; c[5] = 4277./1440.;
	}
	else if (order == 7) {
		c[0] = 19087./60480.; c[1] = -134472./60480.; c[2] = 407139./60480.; c[3] = -688256./60480.; c[4] = 705549./60480.; c[5] = -447288./60480.; c[6] = 198721./60480.;
	}
	else if (order == 8) {
		c[0] = -36799./120960.; c[1] = 295767./120960.; c[2] = -1041723./120960.; c[3] = 2102243./120960.; c[4] = -2664477./120960.; c[5] = 2183877./120960.; c[6] = -1152169./120960.; c[7] = 434241./120960.;
	}
	else if (order == 9) {
		c[0] = 1070017./3628800.; c[1] = -9664106./3628800.; c[2] = 38833486./3628800.; c[3] = -91172642./3628800.; c[4] = 137968480./3628800.; c[5] = -139855262./3628800.; c[6] = 95476786./3628800.; c[7] = -43125206./3628800.; c[8] = 14097247./3628800.;
	}
	else if (order == 10 or order == 11) {
		c[0] = -2082753./7257600.; c[1] = 20884811./7257600.; c[2] = -94307320./7257600.; c[3] = 252618224./7257600.; c[4] = -444772162./7257600.; c[5] = 538363838./7257600.; c[6] = -454661776./7257600.; c[7] = 265932680./7257600.; c[8] = -104995189./7257600.; c[9] = 30277247./7257600.;
	}
	else if (order == 12) {
		c[0] = -262747265./958003200.; c[1] = 3158642445./958003200.; c[2] = -17410248271./958003200.; c[3] = 58189107627./958003200.; c[4] = -131365867290./958003200.; c[5] = 211103573298./958003200.; c[6] = -247741639374./958003200.; c[7] = 214139355366./958003200.; c[8] = -135579356757./958003200.; c[9] = 61633227185./958003200.; c[10] = -19433810163./958003200.; c[11] = 4527766399./958003200.;
	}
}

void adamsM(DVec& d, const int order) {
	if (order == 1) {
		d[0] = 1.;
	}
	else if (order == 2) {
		d[0] = 1./2.; d[1] = 1./2.; 
	}
	else if (order == 3) {
		d[0] = -1./12.; d[1] = 2./3.; d[2] = 5./12.; 
	}
	else if (order == 4) {
		d[0] = 1./24.; d[1] = -5./24.; d[2] = 19./24.; d[3] = 3./8.; 
	}
	else if (order == 5) {
		d[0] = -19./720.; d[1] = 106./720.; d[2] = -264./720.; d[3] = 646./720.; d[4] = 251./720.; 
	}
}

void rungeKutta(DVec& c, DVec& b, DMat& a, const int order, const int variation) {
	if (order == 1) {//Euler
		c[0] = 0.;
		b[0] = 1.;
	}
	else if (order == 2 and variation == 0) {//Heun
		c[0] = 0.; c[1] = 1.;
		b[0] = 0.5; b[1] = 0.5;
		a[1][0] = 1.;
	}
	else if (order == 2 and variation == 1) {//Ralston
		c[0] = 0.; c[1] = 2./3.;
		b[0] = 1./4.; b[1] = 3./4.;
		a[1][0] = 2./3.;
	}
	else if (order == 3) {//Kutta 3
		c[0] = 0.; c[1] = 0.5; c[2] = 1.;
		b[0] = 1./6.; b[1] = 2./3.;  b[2] = 1./6.;
		a[1][0] = 0.5;
		a[2][0] = -1.; a[2][1] = 2.;
	}
	else if (order == 4 and variation == 0) {//RK4
		c[0] = 0.; c[1] = 0.5; c[2] = 0.5; c[3] = 1.;
		b[0] = 1./6.; b[1] = 1./3.;  b[2] = 1./3.;  b[3] = 1./6.;
		a[1][0] = 0.5;
		a[2][0] = 0.; a[2][1] = 0.5;
		a[3][0] = 0.; a[3][1] = 0.; a[3][2] = 1.;
	}
	else if (order == 4 and variation == 1) {//Kutta 3/8
		c[0] = 0.; c[1] = 1./3.; c[2] = 2./3.; c[3] = 1.;
		b[0] = 1./8.; b[1] = 3./8.;  b[2] = 3./8.;  b[3] = 1./8.;
		a[1][0] = 1./3.;
		a[2][0] = -1./3.; a[2][1] = 1.;
		a[3][0] = 1.; a[3][1] = -1.; a[3][2] = 1.;
	}
	else if (order == 6 and variation == 0) {//Fehlberg, 5th order solution
		c[0] = 0.; c[1] = 1./4.; c[2] = 3./8.; c[3] = 12./13.; c[4] = 1.; c[5] = 1./2.;
		b[0] = 16./135.; b[1] = 0.; b[2] = 6656./12825.; b[3] = 28561./56430.; b[4] = -9./50.; b[5] = 2./55.;
		a[1][0] = 1./4.;
		a[2][0] = 3./32.; a[2][1] = 9./32.;
		a[3][0] = 1932./2197.; a[3][1] = -7200./2197.; a[3][2] = 7296./2197.;
		a[4][0] = 439./216.; a[4][1] = -8.; a[4][2] = 3680./513.; a[4][3] = -845./4104.;
		a[5][0] = -8./27.; a[5][1] = 2.; a[5][2] = -3544./2565.; a[5][3] = 1859./4104.; a[5][4] = -11./40.;
	}
	else if (order == 6 and variation == 1) {//Cash-Karp, 5th order solution
		c[0] = 0.; c[1] = 1./5.; c[2] = 3./10.; c[3] = 3./5.; c[4] = 1.; c[5] = 7./8.;
		b[0] = 37./378.; b[1] = 0.; b[2] = 250./621.; b[3] = 125./594.; b[4] = 0.; b[5] = 512./1771.;		
		a[1][0] = 1./5.;
		a[2][0] = 3./40.; a[2][1] = 9./40.;
		a[3][0] = 3./10.; a[3][1] = -9./10.; a[3][2] = 6./5.;
		a[4][0] = -11./54.; a[4][1] = 5./2.; a[4][2] = -70./27.; a[4][3] = 35./27.;
		a[5][0] = 1631./55296.; a[5][1] = 175./512.; a[5][2] = 575./13824.; a[5][3] = 44275./110592.; a[5][4] = 253./4096.;
	}

}
