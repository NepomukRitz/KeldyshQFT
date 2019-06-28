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

using namespace std;
typedef complex<double> comp;
typedef vector<comp> cvec;
typedef vector<cvec> cmat;

//Number of integration grid points
const long int n_integration = 21;
//Number of evolution steps
const int n_evolution = 10;

//Frequency range on which it's integrated
const double x_a = -50.0;
const double x_b = 50.0;

//Values for the start and finish of the fRG flow
const double L_ini = 5.;
const double L_fin = 0.;

vector<double> lambda(n_evolution);     //Declaration of the evolution vector
vector<double> omega(n_integration), omega2(n_integration), omega3(n_integration);  //Declaration of the integration vectors
vector<double > weights(n_integration); //Declaration of the weights vector


const comp epsilon = 0.; // NOLINT(cert-err58-cpp)
const comp Gamma = 1.; // NOLINT(cert-err58-cpp)

class Propagator : public cvec
{
public:
    cvec vec;
    explicit Propagator (int size)
            :vec(cvec(size)) {}

    Propagator operator* (double dt){
        for (int i = 0; i<vec.size(); i++)
            vec[i] *= dt;
        return *this;
    }

    Propagator operator+ (Propagator& b)
    {
        for (int i=0; i<vec.size(); i++){
            vec[i] += b[i];
        }
        return *this;
    }

    comp& operator[] (int i)
    {
        return vec[i];
    }

//    Propagator operator= (Propagator rhs){
//        for(int i=0;i<vec.size();i++)
//        {
//            vec[i] = rhs[i];
//        }
//        return *this;
//    }

    int size(){
        return vec.size();
    }

    void setToZero(){
        for (int i=0; i<vec.size(); i++)
            vec[i] = 0.;
    }
};

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

//--------------------------------------------------------------------------------------------------------------------//
//Function predeclarations
void setUpLambda();
void setUpOmega();
void setUpWeights();
Propagator derivative_Sigma(Vertex &gamma, Propagator &ssProp);
Propagator derivativePropagator(double lambda, Propagator &Sigma);
Propagator L(Vertex &gamma, Propagator &ssProp);
Vertex derivative_gamma(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam);
Vertex bubble_dot_a(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam);
Vertex bubble_dot_p(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam);
Vertex bubble_dot_t(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam);
Propagator G(double lambda, Propagator &Sigma);
Propagator S(double lambda, Propagator &Sigma);
comp dot_product(cvec  &x, vector<double > &y);
void evolve(vector<Propagator> &Sigma, vector<Propagator> &Prop, vector<Propagator> &ssProp, vector<Vertex > &gamma);
void writeOutFile(vector<Propagator> &Sigma, vector<Propagator> &Prop, vector<Vertex > &gamma);

//--------------------------------------------------------------------------------------------------------------------//
//Main
int main() {

    setUpLambda();                          //Initialization of the evolution vector
    setUpOmega();                           //Initialization of the frequency (grid) vector(s)
    omega2 = omega;
    omega3 = omega;
    setUpWeights();                        // Initialization of the weights vector

    /*Vectors to initialize self_energy, propagator, single-scale propagator and vertex structures
     * Strictly auxiliary*/
    Propagator s_omega(omega.size());
    vector<comp> gamma_1(omega.size());
    vector<vector<comp> > gamma_2(omega.size(), gamma_1);
    Vertex gamma_3 (omega.size(), gamma_2);

    //Initialization of the self energy, the propagator, the single scale propagator and the vertex
    vector<Propagator> Sigma(lambda.size(), s_omega);
    vector<Propagator> Prop(lambda.size(), s_omega);
    vector<Propagator> ssProp(lambda.size(), s_omega);
    vector<Vertex > gamma(lambda.size(), gamma_3);

    //Initial conditions for lambda -> infty i.e. L_ini)
    for (int i=0;i<omega.size();i++)
    {
        Sigma[0][i] = 0.5;                  //Initial condition of the self-energy
        for (int j = 0; j<omega.size(); j++)
        {
            for (int k = 0; k<omega.size(); k++)
            {
                gamma[0][i][j][k] = 1.;         //Initial condition for the vertex
            }
        }
    }
    Prop[0] = G(L_ini,Sigma[0]);                    //Initial condition of the full propagator
    ssProp[0] = S(L_ini, Sigma[0]);                  //Initial condition of the single-scale propagator

    //Computes the fRG flow
    evolve(Sigma, Prop, ssProp, gamma);

    //Writes out the results in .dat files
    writeOutFile(Sigma, Prop, gamma);
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
    double domega = (x_b-x_a)/((double)(n_integration-1));
    for (int i=0; i<omega.size();i++)
        omega[i] = x_a + domega*i;

}

//Sets up Simpson's rule weights
void setUpWeights()
{
    for(int l=0;l<n_integration;l++)
        weights[l] = 2 + 2 * (l % 2);
    weights[0]=1.;
    weights[n_integration-1]=1.;
}


//Returns the single-scale propagator at a given scale Lambda
Propagator S(double l, Propagator &Sigma)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++)
        resp[i] = -((comp)0.5i)*G(l,Sigma)[i]*G(l,Sigma)[i];
    return resp;
}

//Returns the full propagator at a given scale Lambda
Propagator G(double l, Propagator &Sigma)
{
    Propagator resp(omega.size());
    for (int i=0; i<omega.size();i++) {

        resp[i] = 1. / (omega[i] - epsilon + ((comp) 0.5i) * (Gamma + l) - Sigma[i]);
    }
    return resp;
}

//Returns the value of the derivative of Sigma w.r.t to Lambda
Propagator derivative_Sigma(Vertex &gamma, Propagator &ssProp) {
    return L(gamma, ssProp);
}

//                                                                       .
//                                                                      / \
//Returns an effective propagator with a closed-loop self-energy      __\./___
Propagator L(Vertex &gamma, Propagator &ssProp)
{
    Propagator resp(ssProp.size());
    for(int i=0; i<resp.size(); i++) {
        vector<comp> y(n_integration);
        double delta_x = (x_b - x_a) / ((double) (n_integration));
        for (int l = 0; l < n_integration; l++) {
            y[l] = -gamma[i][l][l] * ssProp[l];
        }

        resp[i]= (delta_x / 3.) * dot_product(y, weights);
    }
    return resp;
}

//Self-explanatory
comp dot_product(cvec  &x, vector<double > &y)
{
    comp result = 0.;

    for (int i=0; i<y.size(); i++)
        result += x[i]*y[i];

    return result;
}

//Returns the value of the derivative of the full propagator G w.r.t to Lambda
Propagator derivativePropagator(double lamb, Propagator &Sigma)
{
    return S(lamb,Sigma); //+G*Sigma*G
}

//Returns the value of the derivative of the full vertex Gamma w.r.t to Lambda
Vertex derivative_gamma(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam)
{
    Vertex gamma_a = bubble_dot_a(gamma, Prop, dG, lam);
    Vertex gamma_p = bubble_dot_p(gamma, Prop, dG, lam);
    Vertex gamma_t = bubble_dot_t(gamma, Prop, dG, lam);

    return  gamma_a + gamma_p + gamma_t;
}

//Returns the contribution of a bubble diagramm in the a-channel
Vertex bubble_dot_a(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam)
{
    Vertex resp(gamma);
    vector<comp > y(n_integration);
    double delta_x = (x_b-x_a)/((double)(n_integration-1));
    Propagator zeroes (Prop);
    zeroes.setToZero();
    for(int i=0; i<omega.size(); i++) {
        for(int j=0; j<omega2.size();j++) {
            for(int k=0;k<omega3.size(); k++) {
                for (int l = 0; l < n_integration; l++) {
                    int prob = l+j-i;
                    if (prob > n_integration-1) {
                        y[l] = gamma[i][j][l]*(G(lam,zeroes)[l]*dG[n_integration-1])*gamma[l][n_integration-1][k];
                    }
                    else if (prob <= 0){
                        y[l] = gamma[i][j][l]*(G(lam,zeroes)[l]*dG[0])*gamma[l][0][k];
                    }
                    else {
                        y[l] = gamma[i][j][l]*(Prop[l]*dG[l+j-i])*gamma[l][l+j-i][k];
                    }
                }
                resp[i][j][k] = (delta_x / 3.) * dot_product(y, weights);       //No numerical prefactor in Ba
            }
        }
    }
    return resp;
}

//Returns the contribution of a bubble diagramm in the p-channel
Vertex bubble_dot_p(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam)
{
    Vertex resp(gamma);
    vector<comp > y(n_integration);
    double delta_x = (x_b-x_a)/((double)(n_integration-1));
    Propagator zeroes (Prop);
    zeroes.setToZero();
    for(int i=0; i<omega.size(); i++) {
        for(int j=0; j<omega2.size();j++) {
            for(int k=0;k<omega3.size(); k++) {
                for (int l = 0; l < n_integration; l++) {
                    int prob = j+k-l;
                    if (prob >= n_integration) {
                        y[l] =gamma[i][n_integration-1][l]*(G(lam, zeroes)[l]*dG[n_integration-1])*gamma[l][j][k];
                    }
                    else if(prob <= 0){
                        y[l] =gamma[i][0][l]*(G(lam, zeroes)[l]*dG[0])*gamma[l][j][k];
                    }
                    else {
                        y[l] = gamma[i][j+k-l][l]*(Prop[l]*dG[j+k-l])*gamma[l][j][k];
                    }
                }
                resp[i][j][k] = 0.5*(delta_x / 3.) * dot_product(y, weights);   //0.5 coming from the formulas for Bp
            }
        }
    }
    return resp;
}

//Returns the contribution of a bubble diagramm in the t-channel
Vertex bubble_dot_t(Vertex &gamma, Propagator &Prop, Propagator &dG, double lam)
{
    Vertex resp(gamma);
    vector<comp > y(n_integration);
    double delta_x = (x_b-x_a)/((double)(n_integration-1));
    Propagator zeroes (Prop);
    zeroes.setToZero();
    for(int i=0; i<omega.size(); i++) {
        for(int j=0; j<omega2.size();j++) {
            for(int k=0;k<omega3.size(); k++) {
                for (int l = 0; l < n_integration; l++) {
                    int prob = l+k-i;
                    if (prob >= n_integration) {
                        y[l] = gamma[i][l][k]*(G(lam, zeroes)[l]*dG[n_integration-1])*gamma[l][j][n_integration-1];
                    }
                    else if(prob<= 0) {
                        y[l] = gamma[i][l][k]*(Prop[l]*dG[0])*gamma[l][j][0];
                    }
                    else {
                        y[l] = gamma[i][l][k]*(Prop[l]*dG[l+k-i])*gamma[l][j][l+k-i];
                    }
                }
                resp[i][j][k] = -(delta_x / 3.) * dot_product(y, weights);      //- coming from the formulas for Bt
            }
        }
    }
    return resp;
}

//Computes the full fRG flow (here, one-loop fRG)
void evolve(vector<Propagator> &Sigma, vector<Propagator> &Prop, vector<Propagator> &ssProp, vector<Vertex > &gamma)
{
    double dt = (L_fin-L_ini)/((double)(n_evolution-1));
    for(int i=1; i<lambda.size(); i++) {
        Sigma[i] = Sigma[i-1] + derivative_Sigma(gamma[i - 1], ssProp[i - 1])*dt;
        Propagator dG = derivativePropagator(lambda[i - 1], Sigma[i - 1]);
        gamma[i] = gamma[i - 1] + derivative_gamma(gamma[i - 1], Prop[i - 1], dG, lambda[i - 1]) * dt;
        Prop[i] = G(lambda[i], Sigma[i]);
        ssProp[i] = S(lambda[i], Sigma[i]);
    }
}

//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(vector<Propagator> &Sigma, vector<Propagator> &Prop, vector<Vertex > &gamma)
{
    for(int i=0; i<n_evolution; i++) {
        ostringstream self_energy, prop, vertex;
        self_energy << "self_energy" << i << ".dat";
        prop << "propagator" << i << ".dat";
        vertex << "vertex" << i << ".dat";

        ofstream my_file_sigma, my_file_prop, my_file_vertex;
        my_file_sigma.open(self_energy.str());
        my_file_prop.open(prop.str());
        my_file_vertex.open(vertex.str());

        for (int j = 0; j < omega.size(); j++) {
            my_file_sigma << omega[j] << " " << Sigma[i][j].real() << " " << Sigma[i][j].imag() << "\n";
            my_file_prop << omega[j] << " " << Prop[i][j].real() << " " << Prop[i][j].imag() << "\n";
            my_file_vertex << omega[j] << " " << gamma[i][j][0][0].real() << " " << gamma[i][j][0][0].imag() << "\n";

        }
        my_file_sigma.close();
        my_file_prop.close();
        my_file_vertex.close();
    }
}
