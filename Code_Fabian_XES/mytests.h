// compare analytical parquet result to fRG flow of susceptibility and three-point vertex with trivial four-point vertex (see paper Fermi-edge singularity and functional renormalization group)
bool test_analyticParquet(bool print) {
	//MODE = 3 //NFREQ arbitrary (1 is sufficient)
#if MODE!=3
	cout << "MODE must be set to 3 for this test.\n";
	return true;
#endif
	double dev; double max_dev = 0.;
	State state;
// 	int steps;
// #if REGULATOR==0
// 	steps = 4000;
// #else
// 	steps = 800;
// #endif
	int steps = 12000; int ODE_solver = 2; double Lambda_0 = 10000.; double Lambda_f_fraction = 0.001;
	
 	double Lambda_f = Lambda_f_fraction * PiOverBeta;
 	const double alpha = 1. / (steps-1.) * log( Lambda_0 / Lambda_f ); const double alphafactor = 1. - exp(-alpha);
 	double delta, Lambda;
 	for (int m = 0; m < steps; m++) {
 		Lambda = Lambda_0 * exp(- alpha * m);
 		delta = - Lambda * alphafactor;
 		State dState;
 		bubble_Pi(state, dState.PiX, Lambda, 1);
 		bubble_Gamma3(state, dState.GammaCDX, Lambda, 5);
 		state = state + dState * delta;
 		if (print and m%500==0) cout << "steps to go: " << steps-1-m << " at Lambda " << Lambda << endl;
 	}
	
//	propagate(state, Lambda_0, Lambda_f_fraction, steps, ODE_solver, false, "");
	
	if (print) cout << "----------------\nTEST: exact parquet result and flow of bosonic self-energy and three-point vertex using constant four-point vertex\nfrequency \t num. flow \t analytical \t deviation\n";
	for (int i = 0; i < min(100, NFREQ1); ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp analyticParquet = 1./(2.*UVAL) * ( 1. - exp( -2.*UVAL* log( (I*freq + XID)/ (-MU) ) ) );
		dev = abs(state.PiX[i] - analyticParquet)/abs(analyticParquet); max_dev = max(max_dev, dev);
		if (print) cout << freq << "\t\t" << state.PiX[i] << "\t" << analyticParquet << "\t\t" << dev << endl;
		//if (print) cout << freq << "\t\t" << state.PiX[i]+9.107564942603041 << endl;
		//if ( dev > 0.1 ) return false;
	}
	if (print) cout << "maximal relative deviation " << max_dev << "\n---------------\n";
	return true;
}

// compare analytical bubbles against fRG flow of single-scale bubble
bool test_singleScaleBubble(bool print) {
	// only NFREQ1 relevant, set MODE=4
#if MODE!=4
	cout << "MODE must be set to 4 for this test.\n";
	return true;
#endif
#if NFREQ1<10
	cout << "Increase NFREQ1 for this test.\n";
	return true;
#endif
	double dev1; double max_dev1 = 0.;
	double dev2; double max_dev2 = 0.;
	int steps = 400; double Lambda_0 = 10000.; double Lambda_f = 0.001 * PiOverBeta;
	const double alpha = 1. / (steps-1.) * log( Lambda_0 / Lambda_f ); const double alphafactor = 1. - exp(-alpha);
	State state;
	
	double delta, Lambda;
	for (int m = 0; m < steps; m++) {
		Lambda = Lambda_0 * exp(- alpha * m);
		delta = - Lambda * alphafactor;
		State dState; //dState.GammaCDX = 1
		bubble_Pi(state, dState.PiX, Lambda, 1);
		bubble_Pi(state, dState.PiY, Lambda, 2);
		state.PiX = state.PiX + dState.PiX * delta;
		state.PiY = state.PiY + dState.PiY * delta;
		if (print and m%50==0) cout << "steps to go: " << steps-1-m << " at Lambda " << Lambda << endl;
	}
	
	if (print) cout << "-----------------\nTEST: bare bubble in Hubbard-Stratonovich self-energies\nfrequency \t ap: num. flow \t analytical \t deviation \t p: num. flow \t analytical \t deviation\n";
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp apBubble = log( (I*freq + XID)/ (-MU) );
		Cp pBubble = - log( (I*freq - XID)/ (MU) );
		dev1 = abs(state.PiX[i] - apBubble)/abs(apBubble); max_dev1 = max(max_dev1, dev1);
		dev2 = abs(state.PiY[i] - pBubble)/abs(pBubble); max_dev2 = max(max_dev2, dev2);
		if (print) cout << freq << "\t\t" << state.PiX[i] << "\t" << apBubble << "\t" << dev1 << "\t\t" << state.PiY[i] << "\t" << pBubble << "\t" << dev2 << endl;
		if ( max(dev1, dev2) > 0.05 ) return false;
	}
	if (print) cout << "maximal relative deviation: \t ap: " << max_dev1 << " \t p: " << max_dev2 << "\n------------------\n";
	return true;
}

// compare analytical bubbles against fRG flow of vertices that result in bubbles
bool test_singleScaleBubble_Ipap(bool print) {
	// NUMFREQ2 and NUMFREQ3 arbitrary
	double dev1; double max_dev1 = 0.;
	double dev2; double max_dev2 = 0.;
	int steps = 800; double Lambda_0 = 10000.; double Lambda_f = 0.001 * PiOverBeta;
	const double alpha = 1. / (steps-1.) * log( Lambda_0 / Lambda_f ); const double alphafactor = 1. - exp(-alpha);
	State state, outstate;
		
	double delta, Lambda;
	for (int m = 0; m < steps; m++) {
		Lambda = Lambda_0 * exp(- alpha * m);
		delta = - Lambda * alphafactor;
		State dState;
		bubble(state, dState.Ip, Lambda, 1, 1);
		bubble(state, dState.Iap, Lambda, 2, 1);
		outstate = outstate + dState * (delta/(UVAL*UVAL));
		if (print and m%50==0) cout << "steps to go: " << steps-1-m << " at Lambda " << Lambda << endl;
	}
	
	if (print) cout << "---------------\nTEST: bare bubble in vertices \nfrequency \t ap: num. flow \t analytical \t deviation \t p: num. flow \t analytical \t deviation \n";
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp apBubble = log( (I*freq + XID)/ (-MU) );
		Cp pBubble = - log( (I*freq - XID)/ (MU) );
		dev1 = abs(outstate.Ip.getVal(0, 0, w_bar)+UVAL - apBubble)/abs(apBubble); max_dev1 = max(max_dev1, dev1);
		dev2 = abs(outstate.Iap.getVal(0, 0, w_bar)+UVAL - pBubble)/abs(pBubble); max_dev2 = max(max_dev2, dev2);
		if (print) cout << freq << "\t\t" << outstate.Ip.getVal(0, 0, w_bar)+UVAL << "\t" << apBubble << "\t" << dev1 << "\t\t" << outstate.Iap.getVal(0, 0, w_bar)+UVAL << "\t" << pBubble << "\t" << dev2 << endl;
		if ( max(dev1, dev2) > 0.05 ) return false;
	}
	if (print) cout << "maximal relative deviation: \t ap: " << max_dev1 << " \t p: " << max_dev2 << "\n----------------\n";
	return true;
}

// compare analytical bubble against numerical bubble without any flow
bool test_bubble(bool print) {
	// NUMFREQ2 and NUMFREQ3 arbitrary
	double dev1; double max_dev1 = 0.;
	double dev2; double max_dev2 = 0.;
	State state, outstate;
	
	bubble(state, outstate.Ip, 0., 1, 0); outstate.Ip = outstate.Ip * ( 1./(UVAL*UVAL) );
	bubble(state, outstate.Iap, 0., 2, 0); outstate.Iap = outstate.Iap * ( 1./(UVAL*UVAL) );
	
	if (print) cout << "----------------\nTEST: bare bubble \nfrequency \t ap: numerical \t analytical \t deviation \t p: numerical \t analytical \t deviation\n";
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp apBubble = log( (I*freq + XID)/ (-MU) );
		Cp pBubble = - log( (I*freq - XID)/ (MU) );
		dev1 = abs(outstate.Ip.getVal(0, 0, w_bar)+UVAL - apBubble)/abs(apBubble); max_dev1 = max(max_dev1, dev1);
		dev2 = abs(outstate.Iap.getVal(0, 0, w_bar)+UVAL - pBubble)/abs(pBubble); max_dev2 = max(max_dev2, dev2);
		if (print) cout << freq << "\t\t" << outstate.Ip.getVal(0, 0, w_bar)+UVAL << "\t" << apBubble << "\t" << dev1 << "\t\t" << outstate.Iap.getVal(0, 0, w_bar)+UVAL << "\t" << pBubble << "\t" << dev2 << endl;
		if ( max(dev1, dev2) > 0.01 ) return false;
	}
	if (print) cout << "maximal relative deviation: \t ap: " << max_dev1 << " \t p: " << max_dev2 << "\n-----------------\n";

	return true;
}

// compare analytical ladder against fRG flow that results in ladder
bool test_ladder(bool print) {
	// NUMFREQ2 and NUMFREQ3 arbitrary
	double dev1; double max_dev1 = 0.;
	double dev2; double max_dev2 = 0.;
	int steps = 2000; double Lambda_0 = 10000.; double Lambda_f = 0.001 * PiOverBeta;
	const double alpha = 1. / (steps-1.) * log( Lambda_0 / Lambda_f ); const double alphafactor = 1. - exp(-alpha);
	State state1, state2;
		
	double delta, Lambda;
	for (int m = 0; m < steps; m++) {
		Lambda = Lambda_0 * exp(- alpha * m);
		delta = - Lambda * alphafactor;
		State dState1, dState2;
		bubble(state1, dState1.Ip, Lambda, 1, 1);
		bubble(state2, dState2.Iap, Lambda, 2, 1);
		state1 = state1 + dState1 * delta;
		state2 = state2 + dState2 * delta;
		if (print and m%100==0) cout << "steps to go: " << steps-1-m << " at Lambda " << Lambda << endl;
	}
	
	if (print) cout << "------------\nTEST: ladder summation \nfrequency \t ap: numerical \t analytical \t deviation \t p: numerical \t analytical \t deviation\n";
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp apBubble = log( (I*freq + XID)/ (-MU) ); Cp apLadder = 1./(-1./UVAL - apBubble);
		Cp pBubble = - log( (I*freq - XID)/ (MU) ); Cp pLadder = 1./(-1./UVAL - pBubble);
		dev1 = abs(state1.Ip.getVal(0, 0, w_bar) - apLadder)/abs(apLadder); max_dev1 = max(max_dev1, dev1);
		dev2 = abs(state2.Iap.getVal(0, 0, w_bar) - pLadder)/abs(pLadder); max_dev2 = max(max_dev2, dev2);
		if (print) cout << freq << "\t\t" << state1.Ip.getVal(0, 0, w_bar) << "\t" << apLadder << "\t" << dev1 << "\t\t" << state2.Iap.getVal(0, 0, w_bar) << "\t" << pLadder << "\t" << dev2 << endl;
		if ( max(dev1, dev2) > 0.15 ) return false;
	}
	if (print) cout << "maximal relative deviation: \t ap: " << max_dev1 << " \t p: " << max_dev2 << "\n-------------------\n";
	return true;
}

// compare single-scale bubble against numerically differentiated bubble
bool test_dBubble(bool print) {
	// NUMFREQ2, NUMFREQ3 arbitrary
	double dev1; double max_dev1 = 0.;
	double dev2; double max_dev2 = 0.;
	State state, outstate1, outstate2;
	
	double Lambda = 0.8; double dL;
#if REGULATOR==0
	dL = PiOverBeta*1.01;
#else
	dL = 0.001;
#endif
	
	bubble(state, outstate1.Ip, Lambda+dL, 1, 6); bubble(state, outstate2.Ip, Lambda-dL, 1, 6); outstate1.Ip = outstate1.Ip * ( 0.5/(dL) ) + outstate2.Ip * ( -0.5/(dL) );
	bubble(state, outstate1.Iap, Lambda+dL, 2, 6); bubble(state, outstate2.Iap, Lambda-dL, 2, 6); outstate1.Iap = outstate1.Iap * ( 0.5/(dL) ) + outstate2.Iap * ( -0.5/(dL) );
	
	bubble(state, outstate2.Ip, Lambda, 1, 1);
	bubble(state, outstate2.Iap, Lambda, 2, 1);
	
	if (print) cout << "--------------\nTEST: single-scale bubble at Lambda=" << Lambda << "\nfrequency \t ap: num. derivative \t direct \t deviation \t p: num. derivative \t direct \t deviation \n";	
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		dev1 = abs(outstate1.Ip.getVal(0, 0, w_bar) - outstate2.Ip.getVal(0, 0, w_bar))/abs(outstate2.Ip.getVal(0, 0, w_bar)); max_dev1 = max(max_dev1, dev1);
		dev2 = abs(outstate1.Iap.getVal(0, 0, w_bar) - outstate2.Iap.getVal(0, 0, w_bar))/abs(outstate2.Iap.getVal(0, 0, w_bar)); max_dev2 = max(max_dev2, dev2);
		if (print) cout << freq << "\t\t" << outstate1.Ip.getVal(0, 0, w_bar) << "\t" << outstate2.Ip.getVal(0, 0, w_bar) << "\t" << dev1 << "\t\t" << outstate1.Iap.getVal(0, 0, w_bar) << "\t" << outstate2.Iap.getVal(0, 0, w_bar) << "\t" << dev2 << endl;
		if ( max(dev1, dev2) > 0.01 ) return false;
	}
	if (print) cout << "maximal relative deviation: \t ap: " << max_dev1 << " \t p: " << max_dev2 << "\n-----------------\n";

	return true;
}

// test analytical zeroth and first order bubble against numerical correlator obtained from constant vertex
bool test_getCorr(bool print) {
	// all parameters arbitrary
	State state;
	double dev; double max_dev = 0.;
	
	CVec PiX = initCVec(10);
	state.getCorr(PiX);
	
	if (print) cout << "-------------\nTEST: particle-hole susceptibility at zeroth and first order\nfrequency \t  analytical \t numerical \t deviation\n";	
	for (int i = 0; i < 10; ++i) {
		const int w_bar = sumindex_dataindex_bos(i);
		const double freq = PiOverBeta*(double)w_bar;
		Cp temp = log( (I*freq + XID)/ (-MU) );
		Cp analyticResult = temp - UVAL * temp*temp;
		dev = abs(PiX[i] - analyticResult)/abs(analyticResult); max_dev = max(max_dev, dev);
		if (print) cout << freq << "\t\t" << PiX[i] << "\t" << analyticResult << "\t" << dev <<  endl;
		if ( dev > 0.01 ) return false;
	}
	if (print) cout << "maximal relative deviation " << max_dev << "\n-----------------\n";
	return true;
}

void test_run(bool print) {
 	bool test1 = test_analyticParquet(print); cout << "Test1 passed? " << test1 << endl;
// 	bool test2 = test_singleScaleBubble(print); cout << "Test2 passed? " << test2 << endl;
// 	bool test3 = test_singleScaleBubble_Ipap(print); cout << "Test3 passed? " << test3 << endl;
// 	bool test4 = test_bubble(print); cout << "Test4 passed? " << test4 << endl;
//	bool test5 = test_ladder(print); cout << "Test5 passed? " << test5 << endl;
// 	bool test6 = test_dBubble(print); cout << "Test6 passed? " << test6 << endl;
//    bool test7 = test_getCorr(print); cout << "Test7 passed? " << test7 << endl;
}
