#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>

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
typedef matrix<comp> cmat;


const int n = 21;
const double omega_a = 0.0;
const double omega_b = 20.0;
const double domega = (omega_b-omega_a)/((double)(n-1));

const int grid = 100000;
const double x_a = -100.0;
const double x_b = 100.0;
const double dx = (x_b-x_a)/((double)(grid-1));

const comp epsilon = 0.; // NOLINT(cert-err58-cpp)
const comp Gamma = 1.; // NOLINT(cert-err58-cpp)

comp bubble(double omega, double nu);
comp integrate(cvec &Pi);
comp dot_product(vect &x, cvec &y);
void writeOutFile(cmat result, cmat result_lukas, vector<int> integ, cmat result_PAID, vect omega);

void setUpOmega();
void setUpFreq();
void setUpWeights();

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
            : vec(cvec(grid)), omega(0.0) {};

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

/* Initialize vectors for values, frequency grid, results (regular integrator), results (Lukas integrator), vector of ints to control integrator and
the matrix of values of the bubble (the integrand) for a given combination of omega and frequency. Every row corresponds to a fixed value of omega*/
vect omega(n);
vect freq(grid);
vect weights(grid);
IntegrandPAID integrand1;

class funct {
public:
    funct(){};

    comp operator() (double x){
        auto lindex = (int)((x-x_a)/dx);
        int uindex = lindex+1;
        double x1 = freq[lindex];
        double x2 = freq[uindex];
        comp y1 = integrand1[lindex];
        comp y2 = integrand1[uindex];

        return y1 + ((y2-y1)/dx)*(x-x1);
    }
};

int main()
{
    vector<int> cont (n);
    matrix<comp> result(1,n);
    matrix<comp> result_lukas(1,n);
    matrix<comp> result_PAID(1,n);
    matrix<comp> Pi(n,grid);
    Domain1D<comp> D(x_a, x_b);

    setUpOmega();
    setUpFreq();
    setUpWeights();



    for (int i=0; i<n; i++)     //Initialization of the matrix of values of the integrand
    {
        for (int j = 0; j<grid; j++)
        {
            Pi[i][j] = bubble(omega[i],freq[j]);
        }
    }


    //This loops iterates over the omega vector and for every iteration integrates the full function
    for (int i=0; i<n; i++)
    {
        cvec temp(grid);
        for (int j=0;j<grid;j++)
        {
            temp[j] = Pi[i][j];
        }
        double t0 = get_time();
        //Calculate the integral with regular method
        result[0][i] = integrate(temp);
        get_time(t0);

        double t1 = get_time();
        //Calculate the integral using Lukas' integrator
        Integrand temp1 (temp);
        cont[i] = intgk(result_lukas[0][i], x_a, x_b, 1e-7, 1e-2, 1e-9, temp1);
        get_time(t1);

        double t2 = get_time();
        //Calculate the integral using PAID
        integrand1.setVector(temp);
//        integrand1.setOmega(omega[i]);
        funct f;
        vector<PAIDInput> inputs;
        inputs.emplace_back(D, f, 1);

        PAID p(inputs);
        auto result_paid = p.solve();
        result_PAID[0][i] = result_paid[1];
        get_time(t2);
    }

    writeOutFile(result, result_lukas, cont, result_PAID, omega);

    return 0;
}

//Function calculating the value of the integrand in question
comp bubble(double omegax, double nu)
{
    return (1./(nu+omegax-epsilon+((comp)1.i)*Gamma))*(1./(nu-epsilon+((comp)1.i)*Gamma));
}

void writeOutFile(cmat result, cmat result_lukas, vector<int> integ, cmat result_PAID, vect omega)
{
    ofstream myfile;
    myfile.open("output.dat");
    for(int i=0; i<n; i++)
    {
        myfile << omega[i] << " " << result[0][i] << " " << result_lukas[0][i] << " " << integ[i] << " " << result_PAID[0][i] << "\n";
    }
    myfile.close();
}

//Integrates a matrix row-wise i.e. every row is integrated over the same grid and a vector of results is returned
comp integrate(cvec &Pi)
{
    comp resp;
    vect weights(grid);

    //Set up the weights
    for (int j=0;j<grid;j++)
    {
        weights[j] = 2 + 2 * (j % 2);
    }
    weights[0]=1.;
    weights[weights.size()-1]=1.;

    //Integrate using trapezoidal (Simpson's) rule
    resp = (dx/3.)*dot_product(weights,Pi);
    return resp;
}

//Self-explanatory
comp dot_product(vect &x, cvec &y){
    comp result =0.;
    for (int i=0; i<x.size(); i++)
        result += x[i]*y[i];
    return result;
}

void setUpOmega()
{
    for (int i=0; i<omega.size();i++)
        omega[i] = omega_a + i*domega;
}

void setUpFreq()
{
    for (int i=0; i<freq.size();i++)
        freq[i] = x_a + i*dx;
}

void setUpWeights()
{
    for (int i=0;i<weights.size();i++)
    {
        weights[i] = 2+2*(i%2);
    }
    weights[0]=1.;
    weights[weights.size()-1]=1.;
}