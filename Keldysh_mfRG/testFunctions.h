//
// Created by Sa.Aguirre on 2/19/20.
//

#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

void testBubbles(Propagator& g1, Propagator& g2){

    vector<comp> Pia9  (nBOS); vector<comp> Pip9  (nBOS); vector<comp> Pit9  (nBOS);
    vector<comp> Pia12 (nBOS); vector<comp> Pip12 (nBOS); vector<comp> Pit12 (nBOS);
    vector<comp> Pia13 (nBOS); vector<comp> Pip13 (nBOS); vector<comp> Pit13 (nBOS);
    vector<comp> Pia14 (nBOS); vector<comp> Pip14 (nBOS); vector<comp> Pit14 (nBOS);
    vector<comp> Pia15 (nBOS); vector<comp> Pip15 (nBOS); vector<comp> Pit15 (nBOS);


    int i = 0;
    for(auto w:bfreqs){
        IntegrandBubble integrandPia9  (g1, g2, false, w, 9,  'a'); IntegrandBubble integrandPip9  (g1, g2, false, w, 9,  'p'); IntegrandBubble integrandPit9  (g1, g2, false, w, 9,  't');
        IntegrandBubble integrandPia12 (g1, g2, false, w, 12, 'a'); IntegrandBubble integrandPip12 (g1, g2, false, w, 12, 'p'); IntegrandBubble integrandPit12 (g1, g2, false, w, 12, 't');
        IntegrandBubble integrandPia13 (g1, g2, false, w, 13, 'a'); IntegrandBubble integrandPip13 (g1, g2, false, w, 13, 'p'); IntegrandBubble integrandPit13 (g1, g2, false, w, 13, 't');
        IntegrandBubble integrandPia14 (g1, g2, false, w, 14, 'a'); IntegrandBubble integrandPip14 (g1, g2, false, w, 14, 'p'); IntegrandBubble integrandPit14 (g1, g2, false, w, 14, 't');
        IntegrandBubble integrandPia15 (g1, g2, false, w, 15, 'a'); IntegrandBubble integrandPip15 (g1, g2, false, w, 15, 'p'); IntegrandBubble integrandPit15 (g1, g2, false, w, 15, 't');


        Pia9 [i] = integrator(integrandPia9 , w_lower_b, w_upper_b); Pip9 [i] = integrator(integrandPip9 , w_lower_b, w_upper_b); Pit9 [i] = integrator(integrandPit9 , w_lower_b, w_upper_b);
        Pia12[i] = integrator(integrandPia12, w_lower_b, w_upper_b); Pip12[i] = integrator(integrandPip12, w_lower_b, w_upper_b); Pit12[i] = integrator(integrandPit12, w_lower_b, w_upper_b);
        Pia13[i] = integrator(integrandPia13, w_lower_b, w_upper_b); Pip13[i] = integrator(integrandPip13, w_lower_b, w_upper_b); Pit13[i] = integrator(integrandPit13, w_lower_b, w_upper_b);
        Pia14[i] = integrator(integrandPia14, w_lower_b, w_upper_b); Pip14[i] = integrator(integrandPip14, w_lower_b, w_upper_b); Pit14[i] = integrator(integrandPit14, w_lower_b, w_upper_b);
        Pia15[i] = integrator(integrandPia15, w_lower_b, w_upper_b); Pip15[i] = integrator(integrandPip15, w_lower_b, w_upper_b); Pit15[i] = integrator(integrandPit15, w_lower_b, w_upper_b);
        i++;
    }


    ostringstream Pia9file;     ostringstream Pip9file;     ostringstream Pit9file;
    ostringstream Pia12file;    ostringstream Pip12file;    ostringstream Pit12file;
    ostringstream Pia13file;    ostringstream Pip13file;    ostringstream Pit13file;
    ostringstream Pia14file;    ostringstream Pip14file;    ostringstream Pit14file;
    ostringstream Pia15file;    ostringstream Pip15file;    ostringstream Pit15file;

    Pia9file << "Output/Pia9.dat";      Pip9file << "Output/Pip9.dat";   Pit9file << "Output/Pit9.dat";
    Pia12file << "Output/Pia12.dat";    Pip12file << "Output/Pip12.dat"; Pit12file << "Output/Pit12.dat";
    Pia13file << "Output/Pia13.dat";    Pip13file << "Output/Pip13.dat"; Pit13file << "Output/Pit13.dat";
    Pia14file << "Output/Pia14.dat";    Pip14file << "Output/Pip14.dat"; Pit14file << "Output/Pit14.dat";
    Pia15file << "Output/Pia15.dat";    Pip15file << "Output/Pip15.dat"; Pit15file << "Output/Pit15.dat";


    ofstream my_file_Pia9 ; ofstream my_file_Pip9 ; ofstream my_file_Pit9 ;
    ofstream my_file_Pia12; ofstream my_file_Pip12; ofstream my_file_Pit12;
    ofstream my_file_Pia13; ofstream my_file_Pip13; ofstream my_file_Pit13;
    ofstream my_file_Pia14; ofstream my_file_Pip14; ofstream my_file_Pit14;
    ofstream my_file_Pia15; ofstream my_file_Pip15; ofstream my_file_Pit15;

    my_file_Pia9.open( Pia9file.str()) ; my_file_Pip9.open( Pip9file.str()) ; my_file_Pit9.open( Pit9file.str()) ;
    my_file_Pia12.open(Pia12file.str()); my_file_Pip12.open(Pip12file.str()); my_file_Pit12.open(Pit12file.str());
    my_file_Pia13.open(Pia13file.str()); my_file_Pip13.open(Pip13file.str()); my_file_Pit13.open(Pit13file.str());
    my_file_Pia14.open(Pia14file.str()); my_file_Pip14.open(Pip14file.str()); my_file_Pit14.open(Pit14file.str());
    my_file_Pia15.open(Pia15file.str()); my_file_Pip15.open(Pip15file.str()); my_file_Pit15.open(Pit15file.str());

    for(int i = 0; i<nBOS; i++){
        my_file_Pia9  << bfreqs[i] << " " << Pia9 [i].real() << " " << Pia9 [i].imag() << "\n"; my_file_Pip9  << bfreqs[i] << " " << Pip9 [i].real() << " " << Pip9 [i].imag() << "\n"; my_file_Pit9  << bfreqs[i] << " " << Pit9 [i].real() << " " << Pit9 [i].imag() << "\n";
        my_file_Pia12 << bfreqs[i] << " " << Pia12[i].real() << " " << Pia12[i].imag() << "\n"; my_file_Pip12 << bfreqs[i] << " " << Pip12[i].real() << " " << Pip12[i].imag() << "\n"; my_file_Pit12 << bfreqs[i] << " " << Pit12[i].real() << " " << Pit12[i].imag() << "\n";
        my_file_Pia13 << bfreqs[i] << " " << Pia13[i].real() << " " << Pia13[i].imag() << "\n"; my_file_Pip13 << bfreqs[i] << " " << Pip13[i].real() << " " << Pip13[i].imag() << "\n"; my_file_Pit13 << bfreqs[i] << " " << Pit13[i].real() << " " << Pit13[i].imag() << "\n";
        my_file_Pia14 << bfreqs[i] << " " << Pia14[i].real() << " " << Pia14[i].imag() << "\n"; my_file_Pip14 << bfreqs[i] << " " << Pip14[i].real() << " " << Pip14[i].imag() << "\n"; my_file_Pit14 << bfreqs[i] << " " << Pit14[i].real() << " " << Pit14[i].imag() << "\n";
        my_file_Pia15 << bfreqs[i] << " " << Pia15[i].real() << " " << Pia15[i].imag() << "\n"; my_file_Pip15 << bfreqs[i] << " " << Pip15[i].real() << " " << Pip15[i].imag() << "\n"; my_file_Pit15 << bfreqs[i] << " " << Pit15[i].real() << " " << Pit15[i].imag() << "\n";
    }

    my_file_Pia9.close() ; my_file_Pip9.close() ; my_file_Pit9.close() ;
    my_file_Pia12.close(); my_file_Pip12.close(); my_file_Pit12.close();
    my_file_Pia13.close(); my_file_Pip13.close(); my_file_Pit13.close();
    my_file_Pia14.close(); my_file_Pip14.close(); my_file_Pit14.close();
    my_file_Pia15.close(); my_file_Pip15.close(); my_file_Pit15.close();
}

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
