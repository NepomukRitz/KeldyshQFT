//#include <omp.h>
//#include <algorithm>
//#include <chrono>
//#include <cmath>
//#include <vector>
//#include <complex>
//#include <functional>
//#include <iostream>
//#include <set>
//#include <thread>
//
//#include "integrand.hpp"
//#include "include/paid.hpp"
//
//#define M_PI 3.14159265358979323846
//
//using namespace std;
//
//typedef complex<double> comp;
//typedef vector<double > vect;
//typedef vector<comp> cvec;
//
//const int n = 2;
//const double omega_a = 0.;
//const double omega_b = 1.;
//
//const int grid = 100;
//const double lim_a = -10.;
//const double lim_b = 10.;
//const double dx = (lim_b-lim_a)/((double)(grid-1));
//
//const double epsilon = 0.;
//const double Gamma = 1.;
//
//vect freq(grid);
//vect omega(n);
//vect weights(grid);
////IntegrandPAID integrand1;
//
//
////comp bubble(double omega, double freq);
////void setUpOmega();
////void setUpFreq();
////void setUpWeights();
////comp integ(Integrand &y);
//
//
//int main() {
//    Domain1D<comp> D(lim_a, lim_b);
//    cvec constructor(grid);
//    vector<cvec> Pi(omega.size(), constructor);
//
//    setUpOmega();
//    setUpFreq();
//    setUpWeights();
//
//    for (int i = 0; i < omega.size(); i++){
//        for (int j = 0; j < freq.size(); j++){
//            Pi[i][j] = bubble(omega[i], freq[j]);
//        }
//    }
//
//
//    for (int i =0; i<omega.size(); i++){
//
//        integrand1.setVector(Pi[i]);
//        integrand1.setOmega(omega[i]);
//
//        comp result1 = integ(integrand1);
//        cout <<omega[i] << " result1 " << result1.real() << " " << result1.imag()  << "\n";
//
//        funct f;
//
////        for (int j=0; j<grid; j++) {
////            cout << f(0) ;
////        }
////        cout << "\n";
//        vector<PAIDInput> inputs;
//        inputs.emplace_back(D, f, 1);
//
//        PAID p(inputs);
//
//        auto result = p.solve();
//
//
//        cout << omega[i] << " result: " << real(result[1]) << " " << result[1].imag() << "\n";
//
//    }
//
//}
//
//void setUpOmega()
//{
//    auto domega = (omega_b-omega_a)/((double)(n-1));
//    for (int i=0; i<omega.size();i++)
//        omega[i] = omega_a + i*domega;
//}
//
//void setUpFreq()
//{
//    for (int i=0; i<freq.size();i++)
//        freq[i] = lim_a + i*dx;
//}
//
//void setUpWeights()
//{
//    for (int i=0;i<weights.size();i++)
//    {
//        weights[i] = 2+2*(i%2);
//    }
//    weights[0]=1.;
//    weights[weights.size()-1]=1.;
//}
//
//comp integ(Integrand &y)
//{
//    comp resp;
//    auto dfreq = (lim_b-lim_a)/((double)(grid-1));
//    for (int i=0; i<weights.size(); i++)
//    {
//        resp += y[i]*weights[i];
//    }
//    return resp*dfreq/3.;
//}