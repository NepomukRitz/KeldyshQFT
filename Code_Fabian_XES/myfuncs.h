//
//  funcs.h
//  Multiloop fRG for X-ray-edge singularity
//
//  Created by Fabian Kugler on 10.01.17.
//  Copyright © 2016 Fabian Kugler. All rights reserved.
//

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// different types of bubbles

void getSumBoundaries(int& lb1, int& ub1, int& lb2, int& ub2, int& lb3, int& ub3, const int w_bar, const int mu_integer, const int Lambda_integer, const int looptype) {
// adjust summation boundaries due to a Theta or delta function in the integrand
// delta regulator: G^d(\omega) \propto \Theta(|\omega|-\Lambda), S^d(\omega) \propto \delta(|\omega|-\Lambda)
// Litim regulator: S^d(\omega) \propto \Theta(\Lambda - |\omega|)
#if REGULATOR==0
	if (looptype == 1) { // delta( |w1-w_bar| - Lambda )
		const int mLw = -Lambda_integer + w_bar; const int pLw = Lambda_integer + w_bar;
		lb1 = max(-mu_integer, mLw); ub1 = min(mLw, -1); lb2 = max(mLw, 1); ub2 = min(mu_integer, mLw); lb3 = Lambda_integer + w_bar; ub3 = min(mu_integer, pLw);
	}
	else if (looptype > 1) { // Theta( |w1-w_bar| - Lambda )
		lb1 = -mu_integer; ub1 = min(-Lambda_integer + w_bar, -1); lb2 = 1; ub2 = min(mu_integer, -Lambda_integer + w_bar); lb3 = Lambda_integer + w_bar; ub3 = mu_integer;
	}
#elif REGULATOR==1
	if (looptype == 1) { // Theta( Lambda - |w1-w_bar| ) // note: for Lambda < PiOverBeta, Lambda_integer = -1, and there will be no summation
		lb1 = max(-mu_integer, -Lambda_integer + w_bar); lb2 = max(1, -Lambda_integer + w_bar); ub2 = min(mu_integer, Lambda_integer + w_bar);
	}
#elif REGULATOR==2
	if (looptype == 1) { // due to decay use Theta( Lambda' - |w1-w_bar| ) // note: for Lambda' < PiOverBeta, Lambda_integer = -1, and there will be no summation
		lb1 = max(-mu_integer, -Lambda_integer + w_bar); lb2 = max(1, -Lambda_integer + w_bar); ub2 = min(mu_integer, Lambda_integer + w_bar);
	}
	else if (looptype > 1 and Lambda_integer > 0) { // due to decay use Theta( |w1-w_bar| - Lambda' )
		lb1 = -mu_integer; ub1 = min(-Lambda_integer + w_bar, -1); lb2 = 1; ub2 = min(mu_integer, -Lambda_integer + w_bar); lb3 = Lambda_integer + w_bar; ub3 = mu_integer;
	}
#elif REGULATOR==3
	// due to decay use Theta( |w1-w_bar| - Lambda' )
	lb1 = -mu_integer; ub1 = min(-Lambda_integer + w_bar, -1); lb2 = 1; ub2 = min(mu_integer, -Lambda_integer + w_bar); lb3 = Lambda_integer + w_bar; ub3 = mu_integer;
#endif
}

// mixed fermionic-bosonic flow
Cp bubble_Pi_summand(const int w_bar, const int w1, const double Lambda, function<Cp(int,int,double)>& retSD, function<Cp(int,int)>& retV) {
	return retSD(w_bar, w1, Lambda) * retV(w_bar, w1);
}
void bubble_Pi(const State& state, CVec& dPi, const double Lambda, const int ph) {
    // calculate rhs of flow eq. for \chi self-energy, \psi self-energy, or \gamma self-energy
	// \partial_{\Lambda} \Pi^{\chi}_{\Lambda, \bar{\omega}} = \DimInt_{\omega} S^d_{\Lambda, \omega-\bar{\omega}} G^c_{\omega} \big(\tilde{\Gamma}^{\bar{c} d \chi}_{\Lambda, \omega, \omega-\bar{\omega}, \bar{\omega}} \big)^2
	// \partial_{\Lambda} \Pi^{\psi}_{\Lambda, \bar{\omega}} = \DimInt_{\omega} S^d_{\Lambda, \bar{\omega}-\omega} G^c_{\omega} \big(\tilde{\Gamma}^{\bar{c} d \psi}_{\Lambda, \omega, \bar{\omega}-\omega, \bar{\omega}} \big)^2
	
	int Lambda_integer = 0; // looptype always 1, i.e., bubble with single-scale propagator
#if REGULATOR==0
	Lambda_integer = getOddIndex(Lambda); 
#elif REGULATOR==1
	Lambda_integer = getOddIndexDown(Lambda);
#elif REGULATOR==2
	Lambda_integer = getOddIndexDown(Lambda*4.); // artificial Theta function, since for frequencies larger than 4*Lambda, single-scale propagator strongly dampened
#elif REGULATOR==3
	Lambda_integer = getOddIndexUp(Lambda/4.); // artificial Theta function, since for frequencies smaller than Lambda/4, single-scale propagator strongly dampened
#endif
	
	function<Cp(int, int, double)> retSD;
	function<Cp(int, int)> retV;
	if (ph == 1 or ph == 3) {// particle-hole bubble, \chi or \gamma self-energy
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w1-w_bar, Lambda); };
		retV = [&](int w_bar, int w1) { return state.getGammaCDX(w1-w_bar, w_bar) * state.getGammaCDX(w1-w_bar, w_bar); };
	}
	else {// particle-particle bubble, \psi self-energy
		if( not WITH_BOSONIC_SELFENERGIES ) return;
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w_bar-w1, Lambda); };
		retV = [&](int w_bar, int w1) { return state.getGammaCDY(w_bar-w1, w_bar) * state.getGammaCDY(w_bar-w1, w_bar); };
	}
	
	#pragma omp parallel for
	for (int i = 0; i < NFREQ1; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1; // definition here because of parallelization
		getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, 1);
		Cp sum; // = 0.;
		for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_Pi_summand(w_bar, w1, Lambda, retSD, retV); // w1 < 0, Gc = i pi rho 
		for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_Pi_summand(w_bar, w1, Lambda, retSD, retV); // w1 > 0, Gc = - i pi rho 
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
		for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_Pi_summand(w_bar, w1, Lambda, retSD, retV); // w1 > 0, Gc = - i pi rho // this extra part of summation is only necessary for certain regulators 
#endif
		dPi[i] = sum * I*(-PiOverBeta);
	}
}

Cp bubble_Gamma3_summand(const int w, const int w_bar, const int w1, const double Lambda, function<Cp(int,int,double)>& retSD, function<Cp(int,int,int)>& retV) {
	return retSD(w_bar, w1, Lambda) * retV(w, w_bar, w1);
}
void bubble_Gamma3(const State& state, CMat& dGamma3, const double Lambda, const int ph) {
// 	 calculate rhs of flow eq. for \chi, \psi, or \gamma 3-point vertex
// 	
// 	\partial_{\Lambda} \tilde{\Gamma}^{\bar{c} d \chi}_{\Lambda, \omega, \omega-\bar{\omega}, \bar{\omega}} = 
// 	\DimInt_{\omega'} S^d_{\Lambda, \omega-\bar{\omega}} G^c_{\omega} 
// 	\Gamma^{\bar{c} d \chi}_{\Lambda, \omega, \omega-\bar{\omega}, \bar{\omega}} G^{\psi}_{\Lambda, \omega+\omega'-\bar{\omega}} 
// 	\tilde{\Gamma}^{\bar{c} \bar{d} \psi }_{\Lambda, \omega, \omega'-\bar{\omega}, \omega+\omega'-\bar{\omega}}
// 	\tilde{\Gamma}^{\bar{c} \bar{d} \psi}_{\Lambda, \omega', \omega-\bar{\omega}, \omega+\omega'-\bar{\omega}}
// 
// 	\partial_{\Lambda} \tilde{Gamma}^{\bar{c} d \gamma}_{\Lambda, \omega, \omega-\bar{\omega}, \bar{\omega}} = 
// 	\DimInt_{\omega'} S^d_{\Lambda, \omega-\bar{\omega}} G^c_{\omega} 
// 	\Gamma^{\bar{c} d \gamma}_{\Lambda, \omega, \omega-\bar{\omega}, \bar{\omega}}
// 	\tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\Lambda, \omega'-\bar{\omega}, \omega', \omega, \omega-\bar{\omega}}
//
// 	\partial_{\Lambda} \tilde{\Gamma}^{\bar{c} d \psi}_{\Lambda, \omega, \bar{\omega}-\omega, \bar{\omega}} = 
// 	\DimInt_{\omega'} S^d_{\Lambda, \bar{\omega}-\omega} G^c_{\omega} 
// 	\Gamma^{\bar{c} d \psi}_{\Lambda, \omega, \bar{\omega}-\omega, \bar{\omega}} G^{\chi}_{\Lambda, \omega+\omega'-\bar{\omega}} 
// 	\tilde{\Gamma}^{\bar{c} \bar{d} \chi}_{\Lambda, \omega, \bar{\omega}-\omega', \omega+\omega'-\bar{\omega}}
// 	\tilde{\Gamma}^{\bar{c} \bar{d} \chi}_{\Lambda, \omega', \bar{\omega}-\omega, \omega+\omega'-\bar{\omega}}

	int Lambda_integer = 0; // looptype always 1
#if REGULATOR==0
	Lambda_integer = getOddIndex(Lambda); 
#elif REGULATOR==1
	Lambda_integer = getOddIndexDown(Lambda);
#elif REGULATOR==2
	Lambda_integer = getOddIndexDown(Lambda*4.); // artificial Theta function, see above
#elif REGULATOR==3
	Lambda_integer = getOddIndexUp(Lambda/4.); // artificial Theta function, see above
#endif
	
	function<Cp(int, int, double)> retSD;
	function<Cp(int, int, int)> retV;
	if (ph == 1) { // \Gamma^{\bar{c} d \chi}
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w1-w_bar, Lambda); };
		retV = [&](int w, int w_bar, int w1) { return state.getGammaCDX(w1-w_bar, w_bar) * state.getGY(w+w1) * state.getGammaCDY(w, w+w1) * state.getGammaCDY(w1-w_bar, w+w1); };
	}
	else if (ph == 3) { // \Gamma^{\bar{c} c \gamma}
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w1-w_bar, Lambda); };
		retV = [&](int w, int w_bar, int w1) { return state.getGammaCDX(w1-w_bar, w_bar) * ( state.Ip.getVal(w1-w_bar, w, w_bar) + state.Iap.getVal(w1-w_bar, w, w+w1) + UVAL ); };
	}
	else if (ph == 5) { // \Gamma^{\bar{c} c \gamma} with four-point vertex constant
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w1-w_bar, Lambda); };
		retV = [&](int w, int w_bar, int w1) { return state.getGammaCDX(w1-w_bar, w_bar) * (-UVAL); };
	}
	else { // \Gamma^{\bar{c} \bar{d} \psi}
		if( UXVAL == 0. ) return; // UXVAL=0 -> GX=0, quit function
		retSD = [](int w_bar, int w1, double Lambda) { return getSDL(w_bar-w1, Lambda); };
		retV = [&](int w, int w_bar, int w1) { return state.getGammaCDY(w_bar-w1, w_bar) * state.getGX(w1-w) * state.getGammaCDX(w, w1-w) * state.getGammaCDX(w_bar-w1, w1-w); };
	}
	
	#pragma omp parallel for
	for (int i = 0; i < NFREQ2; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		for (int j = 0; j < 2*NFREQ2; ++j) {
			const double w = sumindex_dataindex_ferMat(j,NFREQ2);
			int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
			getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, 1);
			Cp sum; // = 0.;
			for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_Gamma3_summand(w, w_bar, w1, Lambda, retSD, retV); // w1 < 0
			for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_Gamma3_summand(w, w_bar, w1, Lambda, retSD, retV); // w1 < 0
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
			for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_Gamma3_summand(w, w_bar, w1, Lambda, retSD, retV); // w1 > 0
#endif
			dGamma3[i][j] = sum * I*(-PiOverBeta);
		}
	}
}

// vertex flow

// integrand of bubbles, w1 is integration frequency
Cp bubble_summand(const int w, const int v, const int w_bar, const int w1, const double Lambda, function<Cp(int,int,double)>& retGD, function<Cp(int,int,int,int)>& retLeftV, function<Cp(int,int,int,int)>& retRightV) {
	return retLeftV(w, v, w_bar, w1) * retRightV(w, v, w_bar, w1) * retGD(w_bar, w1, Lambda); 
}
void bubble(const State& state, CVertex& oldI, CVertex& newI, const double Lambda, const int ph, const int looptype) {
// 	calculate rhs of flow eq. for fermionic 4-point vertex in certain channel and other vertex-connected bubbles
//  note: in the following two equations, frequencies are labeled in the sequence: lower left, upper left, upper right, lower right
//	
// 	\partial_{\Lambda} I_{p, \Lambda, \omega-\bar{\omega}, \omega, \nu, \nu-\bar{\omega}} = 
// 	\DimInt{\omega'} S^d_{\Lambda, \omega'-\bar{\omega}} G^c_{\omega'}
// 	\tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\Lambda, \omega-\bar{\omega}, \omega, \omega', \omega'-\bar{\omega}}
// 	\tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\Lambda, \omega'-\bar{\omega}, \omega', \nu, \nu-\bar{\omega}}
// 	
// 	\partial_{\Lambda} I_{ap, \Lambda, \bar{\omega}-\nu, \omega, \nu, \bar{\omega}-\omega} = 
// 	\DimInt{\omega'} S^d_{\Lambda, \bar{\omega}-\omega'} G^c_{\omega'}
// 	\tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\Lambda, \bar{\omega}-\nu, \omega', \nu, \bar{\omega}-\omega'}
// 	\tilde{\Gamma}^{\bar{d} c \bar{c} d}_{\Lambda, \bar{\omega}-\omega', \omega, \omega', \bar{\omega}-\omega}

	int Lambda_integer = 0;
#if REGULATOR==0
	if (looptype == 1)
		Lambda_integer = getOddIndex(Lambda); 
	else if (looptype > 1)
		Lambda_integer = getOddIndexUp(Lambda); 
#elif REGULATOR==1
	if (looptype == 1)
		Lambda_integer = getOddIndexDown(Lambda);
#elif REGULATOR==2
	if (looptype == 1)
		Lambda_integer = getOddIndexDown(Lambda*4.); // artificial Theta function, see above
	else if (looptype > 1)
		Lambda_integer = getOddIndexDown(Lambda/29.); // artificial Theta function, see above
#elif REGULATOR==3
	Lambda_integer = getOddIndexUp(Lambda/4.); // artificial Theta function, see above
#endif
	
	// allocate functions used to compute bubbles, note: the (constant) c propagator GC(w1) is attributed to the integral measure and not included in the integrand
	function<Cp(int, int, double)> retGD;
	function<Cp(int, int, int, int)> retLeftV;
	function<Cp(int, int, int, int)> retRightV;
	
	// parametrization via three frequencies: lower left fermionic frequency w, lower right fermionic frequency v, channel-dependent bosonic frequency w_bar; translation: w_bar_p = w_bar_a + w + v
	// e.g., BSE: \gamma_a(w,v,w_bar_a) = \sum_w1 I_a(w, w1-w_bar_a, w_bar_p=w_bar_a+(w+w1-w_bar_a)=w+w1) GC(w1) GD(w1-w_bar_a) [ I_p(w1-w_bar_a, v, w_bar_a) +  I_a(w1-w_bar_a, v, w_bar_p=w_bar_a+(w1-w_bar_a+v)=w1+v) - \Gamma_0 ]
	// e.g., BSE: \gamma_p(w,v,w_bar_p) = \sum_w1 I_p(w, w_bar_p-w1, w_bar_a=w_bar_p-(w+w_bar_a-w1)=w1-w) GC(w1) GD(w_bar_p-w1) [ I_a(w_bar_p-w1, v, w_bar_p) +  I_p(w_bar_p-w1, v, w_bar_a=w_bar_p-(w_bar_p-w1+v)=w1-v) - \Gamma_0 ]
	if (ph == 1) { // a channel
		if (looptype == 0) { // Bethe-Salpeter: scale-independent bubble with Iap (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGD(w1-w_bar); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Iap.getVal(w, w1-w_bar, w+w1); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 8) { // Bethe-Salpeter: scale-dependent bubble with Iap (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Iap.getVal(w, w1-w_bar, w+w1); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 1) { // single-scale bubbble with two full vertices
			retGD = [](int w_bar, int w1, double Lambda) { return getSDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w1-w_bar, w_bar) + state.Iap.getVal(w, w1-w_bar, w+w1) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 2) { // standard bubble with I_ap_dot (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return oldI.getVal(w, w1-w_bar, w+w1) + UVAL; }; //gamma_p_dot=I_ap_dot has no constant U background
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 3) { // standard bubble with I_ap_dot (right) and full vertex (left)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w1-w_bar, w_bar) + state.Iap.getVal(w, w1-w_bar, w+w1) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return oldI.getVal(w1-w_bar, v, v+w1) + UVAL; }; //I_ap_dot
		}
		else if (looptype == 4) { // standard bubble with nonsymmetric dIpL (right) and full vertex (left)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w1-w_bar, w_bar) + state.Iap.getVal(w, w1-w_bar, w+w1) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return oldI.getValns(w1-w_bar, v, w_bar); }; //I_p //getValns has no UVAL initializer
		}
		else if (looptype == 5) { // standard bubble with nonsymmetric dIpR (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return oldI.getValns(w, w1-w_bar, w_bar); }; //I_p //getValns has no UVAL initializer
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 6) { // for tests: standard bubble with two full vertices
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w1-w_bar, w_bar) + state.Iap.getVal(w, w1-w_bar, w+w1) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL; };
		}
		else if (looptype == 7) { // average of 4 and 5 encoded in one function, note dIpR(w, w1-w_bar, w_bar) = dIpL(w1-w_bar, w, w_bar)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w1-w_bar, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return (state.Ip.getVal(w, w1-w_bar, w_bar) + state.Iap.getVal(w, w1-w_bar, w+w1) + UVAL) * oldI.getValns(w1-w_bar, v, w_bar) + oldI.getValns(w1-w_bar, w, w_bar) * (state.Ip.getVal(w1-w_bar, v, w_bar) + state.Iap.getVal(w1-w_bar, v, v+w1) + UVAL); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return 0.5; }; 
		}
	}
	else { // p channel
		if (looptype == 0) { // Bethe-Salpeter: scale-independent bubble with Ip (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGD(w_bar-w1); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 8) { // Bethe-Salpeter: scale-dependent bubble with Ip (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 1) { // single-scale bubbble with two full vertices
			retGD = [](int w_bar, int w1, double Lambda) { return getSDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w) + state.Iap.getVal(w, w_bar-w1, w_bar) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 2) { // standard bubble with I_p_dot (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return oldI.getVal(w, w_bar-w1, w1-w) + UVAL; }; //gamma_ap_dot=I_p_dot has no constant U background
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 3) { // standard bubble with I_p_dot (right) and full vertex (left)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w) + state.Iap.getVal(w, w_bar-w1, w_bar) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return oldI.getVal(w_bar-w1, v, w1-v) + UVAL; }; //I_p
		}
		else if (looptype == 4) { // standard bubble with nonsymmetric dIapL (right) and full vertex (left)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w) + state.Iap.getVal(w, w_bar-w1, w_bar) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return oldI.getValns(w_bar-w1, v, w_bar); }; //I_ap //getValns has no UVAL initializer
		}
		else if (looptype == 5) { // standard bubble with nonsymmetric dIapR (left) and full vertex (right)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return oldI.getValns(w, w_bar-w1, w_bar); }; //I_ap //getValns has no UVAL initializer
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 6) { // for tests: standard bubble with two full vertices
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w, w_bar-w1, w1-w) + state.Iap.getVal(w, w_bar-w1, w_bar) + UVAL; };
			retRightV = [&](int w, int v, int w_bar, int w1) { return state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL; };
		}
		else if (looptype == 7) { // average of 4 and 5 encoded in one function, note dIapR(w, w1-w_bar, w_bar) = dIapL(w1-w_bar, w, w_bar)
			retGD = [](int w_bar, int w1, double Lambda) { return getGDL(w_bar-w1, Lambda); };
			retLeftV = [&](int w, int v, int w_bar, int w1) { return (state.Ip.getVal(w, w_bar-w1, w1-w) + state.Iap.getVal(w, w_bar-w1, w_bar) + UVAL) * oldI.getValns(w_bar-w1, v, w_bar) + oldI.getValns(w_bar-w1, w, w_bar) * (state.Ip.getVal(w_bar-w1, v, w1-v) + state.Iap.getVal(w_bar-w1, v, w_bar) + UVAL); };
			retRightV = [&](int w, int v, int w_bar, int w1) { return 0.5; }; //I_ap //getValns has no UVAL initializer
		}
	}
	
	int Omega = 10001;//2*NUMFREQ1-1; // large fermionic frequency to simulate w, v \to \infty in bubble function
	
	// loop over external arguments (k, i, j) and perform summation (w1) for K1, K2, K2ns (for nonsymmetric vertices), K3
	#pragma omp parallel for
	for (int k = 0; k < NUMFREQ1; ++k) { // newI.K1: w, v \to \infty
		int w_bar = sumindex_dataindex_bos(k);
		int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
		getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, looptype);
		Cp sum; 
		for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_summand(Omega+4, Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV); // shift +4 to avoid diagonal structure
		for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_summand(-Omega-4, -Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
		for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_summand(-Omega-4, -Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#endif
		newI.K1[k] = sum * I*(-PiOverBeta);
	}
	
	#pragma omp parallel for
	for (int k = 0; k < NUMFREQ2; ++k) { // newI.K2: v \to \infty
		int w_bar = sumindex_dataindex_bos(k);
		int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
		getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, looptype);
		for (int i = 0; i < 2*NUMFREQ2; ++i) {
			int w = sumindex_dataindex_ferMat(i,NUMFREQ2);
			Cp sum;
			for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_summand(w, Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
			for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_summand(w, -Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
			for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_summand(w, -Omega, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#endif
			newI.K2[k][i] = sum * I*(-PiOverBeta) - newI.K1[k];
		}
	}
	
	if (looptype == 2 or looptype == 3) { // in case of nonsymmetric vertex (e.g., dIpL, dIapL): newI.K2ns: w \to \infty
		#pragma omp parallel for
		for (int k = 0; k < NUMFREQ2; ++k) {
			int w_bar = sumindex_dataindex_bos(k);
			int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
			getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, looptype);
			for (int j = 0; j < 2*NUMFREQ2; ++j) {
				int v = sumindex_dataindex_ferMat(j,NUMFREQ2);
				Cp sum;
				for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_summand(Omega, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
				for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_summand(-Omega, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
				for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_summand(-Omega, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#endif
				newI.K2ns[k][j] = sum * I*(-PiOverBeta) - newI.K1[k];
			}
		}
	}
	
	if (looptype == 2 or looptype == 3) { // in case of nonsymmetric vertex (e.g., dIpL, dIapL): newI.K3ns: w, v arbitrary
		#pragma omp parallel for
		for (int k = 0; k < NUMFREQ3; ++k) {
			int w_bar = sumindex_dataindex_bos(k);
			int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
			getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, looptype);
			for (int i = 0; i < 2*NUMFREQ3; ++i) {
				int w = sumindex_dataindex_ferMat(i,NUMFREQ3);
				for (int j = 0; j < 2*NUMFREQ3; ++j) {
					int v = sumindex_dataindex_ferMat(j,NUMFREQ3);
					Cp sum;
					for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
					for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
					for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#endif
					newI.K3ns[k][i][j] = sum * I*(-PiOverBeta) - newI.K1[k] - newI.K2[k][i+NUMFREQ2-NUMFREQ3] - newI.K2ns[k][j+NUMFREQ2-NUMFREQ3];
				}
			}
		}
	}
	else { // newI.K3: v >= w
		#pragma omp parallel for
		for (int k = 0; k < NUMFREQ3; ++k) {
			int w_bar = sumindex_dataindex_bos(k);
			int lb1 = -mu_integer; int ub1 = -1; int lb2 = 1; int ub2 = mu_integer; int lb3 = 1; int ub3 = -1;
			getSumBoundaries(lb1, ub1, lb2, ub2, lb3, ub3, w_bar, mu_integer, Lambda_integer, looptype);
			for (int i = 0; i < 2*NUMFREQ3; ++i) {
				int w = sumindex_dataindex_ferMat(i,NUMFREQ3);
				for (int j = i; j < 2*NUMFREQ3; ++j) {
					int v = sumindex_dataindex_ferMat(j,NUMFREQ3);
					Cp sum;
					for (int w1 = lb1; w1 <= ub1; w1 += 2) sum -= bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
					for (int w1 = lb2; w1 <= ub2; w1 += 2) sum += bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#if REGULATOR==0 or REGULATOR==2 or REGULATOR==3
					for (int w1 = lb3; w1 <= ub3; w1 += 2) sum += bubble_summand(w, v, w_bar, w1, Lambda, retGD, retLeftV, retRightV);
#endif
					newI.K3[k][triangularIndex(i, j, 2*NUMFREQ3)] = sum * I*(-PiOverBeta) - newI.K1[k] - newI.K2[k][i+NUMFREQ2-NUMFREQ3] - newI.K2[k][j+NUMFREQ2-NUMFREQ3];
				}
			}
		}
	}
}

void bubble(const State& state, CVertex& newI, const double Lambda, const int ph, const int looptype) {
// calculate bubble with full vertex on both sides, such that oldI is not needed and can be set to anything (e.g. newI)
	bubble(state, newI, newI, Lambda, ph, looptype);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parquet solver

void parquetUpdate(State& state, double Lambda, double damping) {
// perform the update of the self-consistent equation with damping and/or decreasing Lambda
	State newState;
	if(Lambda > 0.) {
		bubble(state, newState.Ip, Lambda, 1, 8); // Lambda, a channel, Bethe-Salpeter bubble
    	bubble(state, newState.Iap, Lambda, 2, 8); // Lambda, p channel, Bethe-Salpeter bubble
	}	
	else {
		bubble(state, newState.Ip, 0., 1, 0); // Lambda=0 irrelevant, a channel, Bethe-Salpeter bubble
		bubble(state, newState.Iap, 0., 2, 0); // Lambda=0 irrelevant, p channel, Bethe-Salpeter bubble
	}
	state = state*damping + newState*(1.-damping);
}

void parquetSolver(State& state, int steps) {
// solve parquet equations iteratively until consecutive results for the ph-suscep. are almost equal or start to differ increasingly
    int N = 25; // number of frequencies for output susceptibility
    double diff1 = 1e20; double diff2 = 0.;
    CVec PiX1(N); // all zeros
    CVec PiX2(N);
    cout << "parquet self-consistency loop" << endl;
    for (int m = 1; m < steps; m++) {
        state.getCorr(PiX2); diff2 = compareCVec(PiX1, PiX2, N);
        cout << "step " << m << "\tdiff: " << diff2 << "\t\tval at w=0: " << PiX2[0] << ", w=10: " << PiX2[10] << ", w=20: " << PiX2[20] << endl;
        if ( diff2 < 1e-4 ) {
            cout << "convergence up to 1e-4 achieved" << endl;
            break;
        }
        if ( diff2 > diff1 ) {
            cout << "increasing difference" << endl;
        }
        if ( diff2 > 1e1 and diff2 > diff1 ) {
            cout << "increasing difference, iteration does not converge" << endl;
            break;
        }
        PiX1 = PiX2; diff1 = diff2;
        parquetUpdate(state, 0., 0.8);
    }
}

void parquetSolver2(State& state, int steps) {
// solve parquet equations iteratively while decreasing Lambda
    int N = 25; // number of frequencies for output susceptibility
    double diff1 = 1e20; double diff2 = 0.;
    CVec PiX1(N); // all zeros
    CVec PiX2(N);
    double Lambda = 100.; // initial value of Lambda
    double factor = 0.96; // Lambda = factor * Lambda, exponential decay
    cout << "parquet self-consistency loop" << endl;
    for (int m = 1; m < steps; m++) {
        state.getCorr(PiX2); diff2 = compareCVec(PiX1, PiX2, N);
        //if (diff2 < 1e-3) {
            if (diff2 < 1e-4 and Lambda*(1.-factor) < 2.*PiOverBeta)
                Lambda = Lambda - 2.*PiOverBeta;
            else
                Lambda = Lambda * factor;
        //}
        cout << "step " << m << "\tdiff: " << diff2 << "\t\tval at w=0: " << PiX2[0] << ", w=10: " << PiX2[10] << ", w=20: " << PiX2[20] << endl;
        if ( Lambda < PiOverBeta and diff2 < 1e-4 ) {
            cout << "convergence up to 1e-4 achieved" << endl;
            break;
        }
        if ( diff2 > diff1 ) {
            cout << "increasing difference" << endl;
        }
        if ( Lambda < PiOverBeta and diff2 > 1e1 and diff2 > diff1 ) {
            cout << "increasing difference, iteration does not converge" << endl;
            break;
        }
        if (Lambda > PiOverBeta)
            cout << "\tLambda index: " << getOddIndex(Lambda) << endl;
        PiX1 = PiX2; diff1 = diff2;
        parquetUpdate(state, Lambda, 0.);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
double mean(const CVec& vec) {
	double mean = 0.; int N = vec.size();
	for (int i = 0; i < N; ++i) mean += vec[i].real()/((double)N);
	return mean;
}
double zero(const CVec& vec) {
	return vec[0].real();
}
void monitor_func(const double Lambda, const double delta, const int step, const State& state, const CTen& dIpL, const CTen& dIpR, const CTen& dIpC,  const CTen& dIapL, const CTen& dIapR, const CTen& dIapC, string suffix, string loop) {
	double Lambda_int = state.par.beta*Lambda/M_PI;
	double delta_int = state.par.beta*delta/M_PI;
	if ( step % 50 == 0 ) {
	CVec c = state.getCorr(10);
	State cState(state); cState.Ip = cState.Ip + dIpL * delta;
	CVec cIpL = cState.getCorr(10); 
	cState = state; cState.Ip = cState.Ip + dIpR * delta;
	CVec cIpR = cState.getCorr(10);
	cState = state; cState.Ip = cState.Ip + dIpC * delta;
	CVec cIpC = cState.getCorr(10);
	cState = state; cState.Iap = cState.Iap + dIapL * delta;
	CVec cIapL = cState.getCorr(10);
	cState = state; cState.Iap = cState.Iap + dIapR * delta;
	CVec cIapR = cState.getCorr(10);
	cState = state; cState.Iap = cState.Iap + dIapC * delta;
	CVec cIapC = cState.getCorr(10);
	
    ofstream Data; Data.open("./monitor/monitor_"+suffix+"_"+to_string(step)+"_"+loop); Data << setprecision(8) << fixed;
	Data <<"LambdaInt\t"<< Lambda_int <<"\tdeltaInt\t"<< delta_int <<"\tCorr(0)\t\t";
	Data <<zero(c)<<"\t\t"<<zero(cIpL-c)<<"\t\t"<<zero(cIpR-c)<<"\t\t"<<zero(cIpC-c)<<"\t\t"<<zero(cIapL-c)<<"\t\t"<<zero(cIapR-c)<<"\t\t"<<zero(cIapC-c)<<endl;
	Data <<"LambdaInt\t"<< Lambda_int <<"\tdeltaInt\t"<< delta_int <<"\tCorr.mn\t\t";
	Data <<mean(c)<<"\t\t"<<mean(cIpL-c)<<"\t\t"<<mean(cIpR-c)<<"\t\t"<<mean(cIpC-c)<<"\t\t"<<mean(cIapL-c)<<"\t\t"<<mean(cIapR-c)<<"\t\t"<<mean(cIapC-c)<<endl;
    Data.close();
	}
}
*/
