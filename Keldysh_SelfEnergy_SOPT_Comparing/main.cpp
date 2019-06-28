#include <utility>
#include <iostream>
#include <string>
#include <math.h>
#include <complex.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include "Propagator.h"
#include "util.h"
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <set>
#include <thread>

#include "integrate_new.h"
#include "approxxpp.h"
#include "util.h"
#include "tests/integrand.hpp"
#include "include/paid.hpp"

using namespace std;
typedef complex<double> comp;
typedef vector<double> vect;
typedef vector<comp> cvec;
typedef matrix<double> matr;


//Number of integration grid points
const long int n_integration = 21;
//Number of evolution steps
const int n_evolution = 5;

//Frequency range on which it's integrated
const double x_a = -40.0;
const double x_b = 40.0;
const double dx = (x_b-x_a)/((double)(n_integration));

//Values for the start and finish of the fRG flow
const double L_ini = 3.;
const double L_fin = 0.;

//Temperature and chemical potential
const double T = 0.01;
const double mu = 0.0;

vector<double> lambda(n_evolution);     //Declaration of the evolution vector
vector<double> omega(n_integration);  //Declaration of the integration vectors
vector<double > weights(n_integration); //Declaration of the weights vector

const comp epsilon = 0.; // NOLINT(cert-err58-cpp)
const comp Gamma = 1.; // NOLINT(cert-err58-cpp)

const double nu = 0.5;


class Vertex : public vector<cmat>
{
public:
    vector<vector<cvec> > vert;
    int dim1, dim2, dim3;
    explicit Vertex (vector<vector<cvec> > &vertex)
            :vert(vertex), dim1(vertex.size()), dim2(vertex.size()), dim3(vertex.size()) {}

    explicit Vertex(int dim, vector<cvec> &preVert)
            :vert(vector<vector<cvec> > (dim, preVert)), dim1(dim), dim2(dim), dim3(dim) {}

    Vertex operator+ (Vertex b) {
        Vertex resp(vert);
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    resp[i][j][k] += b[i][j][k];
                }
            }
        }
        return resp;
    }

    Vertex operator* (double dt) {
        Vertex resp(vert);
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    resp[i][j][k] *= dt;
                }
            }
        }
        return resp;
    }

    vector<cvec>& operator[] (int i)
    {
        return vert[i];
    }


};

class Integrand : public cvec
{
public:
    cvec vec;
    Integrand (cvec vector)
            : vec(vector) {

    }

    comp operator() (double x){
        //operator () for the Integrand class. Interpolates using underlying freq grid and linear method
        if (x<=x_a)
            return vec[0];
        if (x>=x_b)
            return vec[vec.size()-1];
        else {
            auto lindex = (int)((x-x_a)/dx);
            int uindex = lindex + 1;
            auto x1 = x_a+lindex*dx;
            auto x2 = x_a+uindex*dx;
            auto y1 = vec[lindex];
            auto y2 = vec[uindex];
            return y1+(x-x1)*((y2-y1)/dx);
        }
    }

    template <class yT>
    matr select(yT intvals)
    {
        matr resp (2);
        resp(0) = intvals.real();
        resp(1) = intvals.imag();
        return resp;
    }

};

class IntegrandPAID{
    cvec vec;
    double omega;
public:
    IntegrandPAID()
            : vec(cvec(n_integration)), omega(0.0) {};

    void setOmega(double ome){
        omega = ome;
    }

    void setVector(cvec& input){
        vec = input;
    }

    comp operator[] (int i){
        return vec[i];
    }
};

//--------------------------------------------------------------------------------------------------------------------//
//Function predeclarations
void setUpLambda();
void setUpOmega();
void setUpWeights();
double Fermi_distribution(double omega);
Propagator derivative_SigmaR(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam);
Propagator derivative_SigmaA(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam);
Propagator derivative_SigmaK(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam);
Propagator derivative_Sigma22(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam);

Propagator GR(double lambda, Propagator &SigmaR);
Propagator SR(double lambda, Propagator &SigmaR);
Propagator GA(double lambda, Propagator &SigmaA);
Propagator SA(double lambda, Propagator &SigmaA);
Propagator GK(Propagator &PropR, Propagator &PropA, Propagator &SigmaK);
Propagator SK(Propagator &PropR, Propagator &PropA, Propagator &PropK);

comp gR(double lam, double ome){return 1. / (ome - epsilon + ((comp) 0.5i) * (Gamma + lam));}
comp sR(double lam, double ome){return ((comp)-0.5i)*gR(lam, ome)*gR(lam, ome);}
comp gA(double lam, double ome){return 1. / (ome - epsilon - ((comp) 0.5i) * (Gamma + lam));}
comp sA(double lam, double ome){return ((comp)+0.5i)*gA(lam, ome)*gA(lam, ome);}
comp gK(double lam, double ome){return (1.-2.*Fermi_distribution(ome))*(gR(lam, ome)-gA(lam, ome));}
comp sK(double lam, double ome){return (1.-2.*Fermi_distribution(ome))*(sR(lam, ome)-sA(lam, ome));}

comp integrate(Propagator& integrand);
comp simpson_integrate(Propagator& integrand);
comp dot_product(cvec  &x, vector<double > &y);
comp lukas_integrate(Propagator& integrand);
comp PAID_integrate(Propagator& integrand);

void evolve(vector<Propagator> &SigmaR, vector<Propagator> &SigmaA, vector<Propagator> &SigmaK, vector<Propagator> &Sigma22,
            vector<Propagator> &PropR, vector<Propagator> &PropA, vector<Propagator> &PropK,
            vector<Propagator> &ssPropR, vector<Propagator> &ssPropA, vector<Propagator> &ssPropK);
void writeOutFile(vector<Propagator> &SigmaR, vector<Propagator> &SigmaA, vector<Propagator> &SigmaK, vector<Propagator> &Sigma22, vector<Propagator> &PropR, vector<Propagator> &PropA, vector<Propagator> &PropK);

Vertex Gamma1111(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1112(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1121(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1122(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1211(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1212(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1221(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma1222(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2111(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2112(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2121(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2122(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2211(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2212(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2221(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);
Vertex Gamma2222(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam);

IntegrandPAID integrand1;

class funct {
public:
    funct(){};

    comp operator() (double x){
        auto lindex = (int)((x-x_a)/dx);
        int uindex = lindex+1;
        double x1 = omega[lindex];
        double x2 = omega[uindex];
        comp y1 = integrand1[lindex];
        comp y2 = integrand1[uindex];

        return y1 + ((y2-y1)/dx)*(x-x1);
    }
};

//--------------------------------------------------------------------------------------------------------------------//
//Main
int main() {

    setUpLambda();                          //Initialization of the evolution vector
    setUpOmega();                           //Initialization of the frequency (grid) vector(s)
    setUpWeights();                        // Initialization of the weights vector

    /*Vector to initialize self_energy, propagator and single-scale propagator structures
     * Strictly auxiliary*/
    Propagator s_omega(omega.size());

    //Initialization of the self energy, the propagator and the single scale propagator
    vector<Propagator> SigmaR(lambda.size(), s_omega);
    vector<Propagator> SigmaA(lambda.size(), s_omega);
    vector<Propagator> SigmaK(lambda.size(), s_omega);
    vector<Propagator> Sigma22(lambda.size(), s_omega);

    vector<Propagator> PropR(lambda.size(), s_omega);
    vector<Propagator> PropA(lambda.size(), s_omega);
    vector<Propagator> PropK(lambda.size(), s_omega);

    vector<Propagator> ssPropR(lambda.size(), s_omega);
    vector<Propagator> ssPropA(lambda.size(), s_omega);
    vector<Propagator> ssPropK(lambda.size(), s_omega);


    //Initial conditions for lambda -> infty i.e. L_ini)
    for (int i=0;i<omega.size();i++)
    {
        SigmaR[0][i] = 0.5;                  //Initial condition of the self-energy
        SigmaA[0][i] = 0.5;
        SigmaK[0][i] = 0.0;
        Sigma22[0][i] = 0.0;
    }

    PropR[0] = GR(L_ini,SigmaR[0]);                    //Initial condition of the full propagator
    PropA[0] = GA(L_ini,SigmaA[0]);                    //Initial condition of the full propagator
    PropK[0] = GK(PropR[0],PropA[0],SigmaK[0]);                    //Initial condition of the full propagator

    ssPropR[0] = SR(L_ini, SigmaR[0]);                  //Initial condition of the single-scale propagator
    ssPropA[0] = SA(L_ini, SigmaA[0]);                  //Initial condition of the single-scale propagator
    ssPropK[0] = SK(PropR[0], PropA[0], PropK[0]);                  //Initial condition of the single-scale propagator


    //Computes the fRG flow
    evolve(SigmaR, SigmaA, SigmaK, Sigma22, PropR, PropA, PropK, ssPropR, ssPropA, ssPropK);

    //Writes out the results in .dat files
    writeOutFile(SigmaR, SigmaA, SigmaK, Sigma22, PropR, PropA, PropK);
    return 0;
}

//--------------------------------------------------------------------------------------------------------------------//
//Implementation of the functions

//Sets up the evolution vector
void setUpLambda()
{
    double dt = (L_fin-L_ini)/((double)(n_evolution-1));
    for (int i =0; i<n_evolution; i++)
        lambda[i]= L_ini + i*dt;

}

//Sets up the frequency (grid) vector
void setUpOmega()
{
    double dx = (x_b-x_a)/((double)(n_integration-1));
    for (int i=0; i<omega.size();i++)
        omega[i] = x_a + dx*i;

}

//Sets up Simpson's rule weights
void setUpWeights()
{
    for(int l=0;l<n_integration;l++)
        weights[l] = 2 + 2 * (l % 2);
    weights[0]=1.;
    weights[n_integration-1]=1.;
}

//Self-explanatory
double Fermi_distribution(double ome)
{
    return 1./(exp((ome-mu)/T)+1);
}

//Returns the retarded single-scale propagator at a given scale Lambda
Propagator SR(double l, Propagator &SigmaR)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++)
        resp[i] = -((comp)0.5i)*GR(l,SigmaR)[i]*GR(l,SigmaR)[i];
    return resp;
}

//Returns the advanced single-scale propagator at a given scale Lambda
Propagator SA(double l, Propagator &SigmaA)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++)
        resp[i] = ((comp)0.5i)*GA(l,SigmaA)[i]*GA(l,SigmaA)[i];
    return resp;
}

//Returns the Keldysh component of the single scale propagator for given GR, GA and GK's.
Propagator SK(Propagator& PropR, Propagator& PropA, Propagator& PropK)
{
    Propagator resp(omega.size());
    for(int i=0; i<omega.size(); i++){
        resp[i] = ((comp)-0.5i)*(PropK[i]*PropA[i]-PropR[i]*PropK[i])-((comp)1.i)*(1.-2.*Fermi_distribution(omega[i]))*PropR[i]*PropA[i];
    }
    return resp;
}

//Returns the full retarded propagator at a given scale Lambda
Propagator GR(double l, Propagator &SigmaR)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++) {
        resp[i] = 1. / (omega[i] - epsilon + ((comp) 0.5i) * (Gamma + l) - SigmaR[i]);
    }
    return resp;
}

//Returns the full advanced propagator at a given scale Lambda
Propagator GA(double l, Propagator &SigmaA)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++) {
        resp[i] = 1. / (omega[i] - epsilon - ((comp) 0.5i) * (Gamma + l) - SigmaA[i]);
    }
    return resp;
}

//Returns the full Keldysh propagator at a given scale Lambda
Propagator GK(Propagator &PropR, Propagator &PropA, Propagator &SigmaK)
{
    Propagator resp(omega.size());
    for(int i=0;i<omega.size();i++){
        resp[i]= PropR[i]*SigmaK[i]*PropA[i];
    }
    return resp;
}

//Returns the value of the derivative of SigmaRet w.r.t to Lambda
Propagator derivative_SigmaR(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam)
{
    double t0 = get_time();
    Propagator resp(omega.size()), integrand(omega.size()), summand1(omega.size()), summand2(omega.size()), summand3(omega.size());
    Vertex vertex1 = Gamma1122(PropR,PropA,PropK,lam);
    Vertex vertex2 = Gamma1221(PropR,PropA,PropK,lam);
    Vertex vertex3 = Gamma1222(PropR,PropA,PropK,lam);
    for(int i=0; i<omega.size(); i++){
        for (int j=0; j<omega.size(); j++){
            summand1[j] = vertex1[i][j][i]*sR(lam,omega[j]);
            summand2[j] = vertex2[i][j][i]*sA(lam,omega[j]);
            summand3[j] = vertex3[i][j][i]*sK(lam,omega[j]);
        }
        integrand = summand1+summand2+summand3;
        resp[i] = -integrate(integrand);
    }
    get_time(t0);
    return resp;
}

//Returns the value of the derivative of SigmaAdv w.r.t to Lambda
Propagator derivative_SigmaA(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam)
{
    double t0 = get_time();
    Propagator resp(omega.size()), integrand(omega.size()), summand1(omega.size()), summand2(omega.size()), summand3(omega.size());
    Vertex vertex1 = Gamma2112(PropR,PropA,PropK,lam);
    Vertex vertex2 = Gamma2211(PropR,PropA,PropK,lam);
    Vertex vertex3 = Gamma2212(PropR,PropA,PropK,lam);
    for(int i=0; i<omega.size(); i++){
        for (int j=0; j<omega.size(); j++){
            summand1[j] = vertex1[i][j][i]*sR(lam,omega[j]);
            summand2[j] = vertex2[i][j][i]*sA(lam,omega[j]);
            summand3[j] = vertex3[i][j][i]*sK(lam,omega[j]);
        }
        integrand = summand1+summand2+summand3;
        resp[i] = -integrate(integrand);
    }
    get_time(t0);
    return resp;
}

//Returns the value of the derivative of SigmaKel w.r.t to Lambda
Propagator derivative_SigmaK(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam)
{
    double t0 = get_time();
    Propagator resp(omega.size()), integrand(omega.size()), summand1(omega.size()), summand2(omega.size()), summand3(omega.size());
    Vertex vertex1 = Gamma1112(PropR,PropA,PropK,lam);
    Vertex vertex2 = Gamma1211(PropR,PropA,PropK,lam);
    Vertex vertex3 = Gamma1212(PropR,PropA,PropK,lam);
    for(int i=0; i<omega.size(); i++){
        for (int j=0; j<omega.size(); j++){
            summand1[j] = vertex1[i][j][i]*sR(lam,omega[j]);
            summand2[j] = vertex2[i][j][i]*sA(lam,omega[j]);
            summand3[j] = vertex3[i][j][i]*sK(lam,omega[j]);
        }
        integrand = summand1+summand2+summand3;
        resp[i] = -integrate(integrand);
    }
    get_time(t0);
    return resp;
}

//Just fot he LOLZ. Sigma22 is zero due to causality
Propagator derivative_Sigma22(Propagator &PropR, Propagator &PropA, Propagator &PropK, Propagator &ssPropR, Propagator &ssPropA, Propagator &ssPropK, double lam)
{
    double t0 = get_time();
    Propagator resp(omega.size()), integrand(omega.size()), summand1(omega.size()), summand2(omega.size()), summand3(omega.size());
    Vertex vertex1 = Gamma2122(PropR,PropA,PropK,lam);
    Vertex vertex2 = Gamma2221(PropR,PropA,PropK,lam);
    Vertex vertex3 = Gamma2222(PropR,PropA,PropK,lam);
    for(int i=0; i<omega.size(); i++){
        for (int j=0; j<omega.size(); j++){
            summand1[j] = vertex1[i][j][i]*sR(lam,omega[j]);
            summand2[j] = vertex2[i][j][i]*sA(lam,omega[j]);
            summand3[j] = vertex3[i][j][i]*sK(lam,omega[j]);
        }
        integrand = summand1+summand2+summand3;
        resp[i] = -integrate(integrand);
    }
    get_time(t0);
    return resp;
}

//Chooses the integration method to be used
comp integrate(Propagator& integrand)
{
    return simpson_integrate(integrand);
//    return lukas_integrate(integrand);
//    return PAID_integrate(integrand);
}

//Integrates the Propagator integrand over the omega grid
comp simpson_integrate(Propagator& integrand)
{
    return dx/3.*dot_product(integrand.vec,weights);
}
//Self-explanatory
comp dot_product(cvec  &x, vector<double > &y)
{
    comp result = 0.;

    for (int i=0; i<y.size(); i++)
        result += x[i]*y[i];

    return result;
}

comp lukas_integrate(Propagator& integrand)
{
    comp answ;
    Integrand temp = integrand.vec;
    return intgk(answ, x_a, x_b, 1e-7, 1e-3, 1e-9, temp);
}


comp PAID_integrate(Propagator& integrand)
{
    Domain1D<comp> D(x_a, x_b);
    integrand1.setVector(integrand.vec);
    funct f;
    vector<PAIDInput> inputs;
    inputs.emplace_back(D, f, 1);

    PAID p(inputs);
    auto result_paid = p.solve();
    return result_paid[1];
}

//Each function here returns the corresponding Keldysh component of the vertex, taken in SOPT
Vertex Gamma1111(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gR(lam,a)+gA(lam,p)+gR(lam,t)) + gA(lam,q)*(gA(lam,a)+gR(lam,p)+gA(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1112(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,a) + gA(lam,q)*(gK(lam,p) + gK(lam,t)) + gK(lam,q)*(gA(lam,a)+gA(lam,p)+gR(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1121(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,t) + gA(lam,q)*(gK(lam,a) + gK(lam, p)) + gK(lam,q)*(gR(lam,a)+gA(lam,p)+gA(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1122(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gA(lam,a)+gA(lam,p)+gA(lam,t)) + gA(lam,q)*(gR(lam,a)+gR(lam,p)+gR(lam,t)) + gK(lam,q)*(gK(lam,a)+gK(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1211(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,p) + gA(lam,q)*(gK(lam,a) + gK(lam, t)) + gK(lam,q)*(gR(lam,a)+gR(lam,p)+gR(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1212(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gA(lam,a)+gR(lam,p)+gR(lam,t)) + gA(lam,q)*(gR(lam,a)+gA(lam,p)+gA(lam,t)) + gK(lam,q)*(gK(lam,a)+gK(lam,p));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1221(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gR(lam,a)+gR(lam,p)+gA(lam,t)) + gA(lam,q)*(gA(lam,a)+gA(lam,p)+gR(lam,t)) + gK(lam,q)*(gK(lam, p)+gK(lam,t));
                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma1222(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gK(lam,a)+gK(lam,p)+gK(lam,t)) + gK(lam,q)*(gA(lam,a)+gR(lam,p)+gA(lam,t));
                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2111(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gK(lam,a)+gK(lam,p)+gK(lam,t)) + gK(lam,q)*(gA(lam,a)+gR(lam,p)+gA(lam,t));
                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2112(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gR(lam,a)+gR(lam,p)+gA(lam,t)) + gA(lam,q)*(gA(lam,a)+gA(lam,p)+gR(lam,t)) + gK(lam,q)*(gK(lam,p)+gK(lam,t));
                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2121(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gA(lam,a)+gR(lam,p)+gR(lam,t)) + gA(lam,q)*(gR(lam,a)+gA(lam,p)+gA(lam,t)) + gK(lam,q)*(gK(lam,a)+gK(lam,p));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2122(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,p) + gA(lam,q)*(gK(lam,a) + gK(lam, t)) + gK(lam,q)*(gR(lam,a)+gR(lam,p)+gR(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2211(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gA(lam,a)+gA(lam,p)+gA(lam,t))+ gA(lam,q)*(gR(lam,a)+gR(lam,p)+gR(lam,t)) + gK(lam,q)*(gK(lam,a)+gK(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2212(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,t) + gA(lam,q)*(gK(lam,a) + gK(lam, p)) + gK(lam,q)*(gR(lam,a)+gA(lam,p)+gA(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2221(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*gK(lam,a) + gA(lam,q)*(gK(lam,p) + gK(lam, t)) + gK(lam,q)*(gA(lam,a)+gA(lam,p)+gR(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}
Vertex Gamma2222(Propagator &PropR, Propagator &PropA, Propagator &PropK, double lam)
{
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex resp (omega.size(), gamma_2);
    for(int i=0; i<omega.size(); i++){
        for(int j=0; j<omega.size(); j++){
            for(int k=0; k<omega.size(); k++){
                Propagator integrand(omega.size());
                for (int l =0; l<omega.size(); l++){
                    double q = omega[l];
                    double a = omega[l]+omega[j]-omega[k];
                    double p = omega[i]+omega[j]-omega[l];
                    double t = omega[l]+omega[k]-omega[i];

                    integrand[l] = gR(lam,q)*(gR(lam,a)+gA(lam,p)+gR(lam,t)) + gA(lam,q)*(gA(lam,a)+gR(lam,p)+gA(lam,t));

                }
                resp[i][j][k] = nu*nu*integrate(integrand);
            }
        }
    }
    return resp;
}




//Computes the full fRG flow (here, one-loop fRG)
void evolve(vector<Propagator> &SigmaR, vector<Propagator> &SigmaA, vector<Propagator> &SigmaK, vector<Propagator> &Sigma22,
            vector<Propagator> &PropR, vector<Propagator> &PropA, vector<Propagator> &PropK,
            vector<Propagator> &ssPropR, vector<Propagator> &ssPropA, vector<Propagator> &ssPropK)
{
    double dt = (L_fin-L_ini)/((double)(n_evolution-1));
    for(int i=1; i<lambda.size(); i++) {
        SigmaR[i] = derivative_SigmaR(PropR[i-1], PropA[i-1], PropK[i-1], ssPropR[i - 1],ssPropA[i-1],ssPropK[i-1], lambda[i])*dt + SigmaR[i-1];
        SigmaA[i] = derivative_SigmaA(PropR[i-1], PropA[i-1], PropK[i-1], ssPropR[i - 1],ssPropA[i-1],ssPropK[i-1], lambda[i])*dt + SigmaA[i-1];
        SigmaK[i] = derivative_SigmaK(PropR[i-1], PropA[i-1], PropK[i-1], ssPropR[i - 1],ssPropA[i-1],ssPropK[i-1], lambda[i])*dt + SigmaK[i-1];
        Sigma22[i] = derivative_Sigma22(PropR[i-1], PropA[i-1], PropK[i-1], ssPropR[i - 1],ssPropA[i-1],ssPropK[i-1], lambda[i])*dt + Sigma22[i-1];

        PropR[i] = GR(lambda[i], SigmaR[i]);
        PropA[i] = GA(lambda[i], SigmaA[i]);
        PropK[i] = GK(PropR[i], PropA[i], SigmaK[i]);

        ssPropR[i] = SR(lambda[i], SigmaR[i]);
        ssPropA[i] = SA(lambda[i], SigmaA[i]);
        ssPropK[i] = SK(PropR[i], PropA[i], PropK[i]);
    }
}

//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(vector<Propagator> &SigmaR, vector<Propagator> &SigmaA, vector<Propagator> &SigmaK, vector<Propagator> &Sigma22, vector<Propagator> &PropR, vector<Propagator> &PropA, vector<Propagator> &PropK)
{
    for(int i=0; i<n_evolution; i++) {
        ostringstream self_energyR, self_energyA, self_energyK, self_energy22, propR, propA, propK;
        self_energyR << "self_energyR" << i << ".dat";
        self_energyA << "self_energyA" << i << ".dat";
        self_energyK << "self_energyK" << i << ".dat";
        self_energy22 << "self_energy22" << i << ".dat";

        propR << "propagatorR" << i << ".dat";
        propA << "propagatorA" << i << ".dat";
        propK << "propagatorK" << i << ".dat";

        ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_sigma22, my_file_propR, my_file_propA, my_file_propK;
        my_file_sigmaR.open(self_energyR.str());
        my_file_sigmaA.open(self_energyA.str());
        my_file_sigmaK.open(self_energyK.str());
        my_file_sigma22.open(self_energy22.str());

        my_file_propR.open(propR.str());
        my_file_propA.open(propA.str());
        my_file_propK.open(propK.str());


        for (int j = 0; j < omega.size(); j++) {
            my_file_sigmaR << omega[j] << " " << SigmaR[i][j].real() << " " << SigmaR[i][j].imag() << "\n";
            my_file_sigmaA << omega[j] << " " << SigmaA[i][j].real() << " " << SigmaA[i][j].imag() << "\n";
            my_file_sigmaK << omega[j] << " " << SigmaK[i][j].real() << " " << SigmaK[i][j].imag() << "\n";
            my_file_sigma22 << omega[j] << " " << Sigma22[i][j].real() << " " << Sigma22[i][j].imag() << "\n";

            my_file_propR << omega[j] << " " << PropR[i][j].real() << " " << PropR[i][j].imag() << "\n";
            my_file_propA << omega[j] << " " << PropA[i][j].real() << " " << PropA[i][j].imag() << "\n";
            my_file_propK << omega[j] << " " << PropK[i][j].real() << " " << PropK[i][j].imag() << "\n";
        }
        my_file_sigmaR.close();
        my_file_sigmaA.close();
        my_file_sigmaK.close();
        my_file_sigma22.close();

        my_file_propR.close();
        my_file_propA.close();
        my_file_propK.close();
    }
}