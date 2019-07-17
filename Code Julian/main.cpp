#include <iostream>
#include<iomanip>
#include<cmath>
#include<fstream> 
#include<type_traits>
#include <string>
#include<tgmath.h>
#include "kagome.hpp"
#include <sstream>

using namespace std;
// using namespace boost::numeric::odeint;

//general contants:

const complex<double> zero(0.0,0.0);//zero written as complex number with real and imag part =0.
const complex<double> ci(0.0,1.0); //complex i

const double pi =3.14159265;

/********************************constants concerning the regularizer*************************/
const double sharp = 2;//governs the sharpness of the second cutoff


/********************************constants concerning the lattice*************************/

//The cutoff distance from the reference site d_c is then given by (nuc-1)/2 (in units of lattice spacings)
const double d_c = (nuc-1)/2;

/********************************constants concerning the frequency grid*************************/
vector<double> Lambdas;//grid for RG flow parameter lambda
const double Lambda_i = 10; //initial value of flow parameter
int nw = nw1;



//log grid:
const double wt = 4.5;//transition frequency. This freq is considered to be the last freq on the log grid. Its value should not be too high for good resolution
const double w0 = 0.05;//lowest freq on grid
const int nlog = 14;//number of freqs on log part . should be an even number
const double k = pow((wt/w0),2./(nlog-2));//exponent for log part
const double delw = wt - w0*pow(k,(nlog/2-2));//spacing between wt and next lower freq
const int nlin = nw - nlog;//compute number of points on linear parts of grid


//linear grid
//const double k =0.8;//some arbitrary factor that sets the density of the stored frequencies
//const double w0 = 0.5*k;
//const double wmax = nw/2*k;//value of maximal frequency that kan be reached for fixed k and nw.



const int reg = 2; //sets the regulator that is being used: 1: sharp cutoff, 2: smoothened cutoff
const int sym = 2;//turns on the symmetry relations 0: off, 1: generic operations, 2:physical symmetries




vector<double> ffreqs(nw);//fermionic matsubara freqs
vector<double> bfreqs(nw);//bosnonic matsubara freqs, NOTE: in the case T=0 (which is implemented here), these two grids are equivalent


state evaluate( double Lambda,state state_i){//this fucntion evaluates the RHS of the flow equation at iteration lambda

    state one_loop, two_loop,dstate_f;//result of this function
    dstate_f.Lambda = Lambda;


    dstate_f.selfenergy =  loop(Lambda, state_i.vertex,'s',state_i.selfenergy,state_i.selfenergy);//yields differentiated self energy. Note that the last argument is not used since this is always the single scale propagator


    cout << "finished self energy" << endl;
    cout << "computing bubbles.." << endl;

    //one loop:

    //s-channel
    cout << "first loop order:"<< endl;
    parvert<svert> sbubble_dif;//first argument specifies whcih of the vertices only contain complementary vertices
    sbubble_dif = sbubble(0,Lambda,state_i.vertex,state_i.vertex,'g','s',state_i.selfenergy,dstate_f.selfenergy)+ sbubble(0,Lambda,state_i.vertex,state_i.vertex,'s','g',state_i.selfenergy,dstate_f.selfenergy);
    one_loop.vertex.spinvertex.svertex = sbubble_dif.spinvertex;
    one_loop.vertex.densvertex.svertex = sbubble_dif.densvertex; //s-channel contribution to differentiated vertex
    cout << "finished s-channel" << endl;
    //    //  t-channel
    parvert<tvert> tbubble_dif;
    tbubble_dif = tbubble(0,Lambda,state_i.vertex,state_i.vertex,'g','s',state_i.selfenergy,dstate_f.selfenergy)+ tbubble(0,Lambda,state_i.vertex,state_i.vertex,'s','g',state_i.selfenergy,dstate_f.selfenergy);
    one_loop.vertex.spinvertex.tvertex = tbubble_dif.spinvertex;
    one_loop.vertex.densvertex.tvertex = tbubble_dif.densvertex; //t-channel contribution to differentiated vertex
    cout << "finished t-channel" << endl;
    ////     // u-channel
    parvert<uvert> ububble_dif;
    ububble_dif = ububble(0,Lambda,state_i.vertex,state_i.vertex,'g','s',state_i.selfenergy,dstate_f.selfenergy) + ububble(0,Lambda,state_i.vertex,state_i.vertex,'s','g',state_i.selfenergy,dstate_f.selfenergy);
    one_loop.vertex.spinvertex.uvertex = ububble_dif.spinvertex;
    one_loop.vertex.densvertex.uvertex = ububble_dif.densvertex; //u-channel contribution to differentiated vertex
    cout << "finished u-channel" << endl;



    //   //two loops:
    cout << "second loop order:"<< endl;
    //   //s-channel
//       sbubble_dif = sbubble(1,Lambda,one_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy)+ sbubble(2,Lambda,state_i.vertex,one_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
//       two_loop.vertex.spinvertex.svertex = sbubble_dif.spinvertex;
//       two_loop.vertex.densvertex.svertex = sbubble_dif.densvertex; //s-channel contribution to differentiated vertex
//       cout << "finished s-channel" << endl;
//       //    //  t-channel
//       tbubble_dif = tbubble(1,Lambda,one_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy)+ tbubble(2,Lambda,state_i.vertex,one_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
//       two_loop.vertex.spinvertex.tvertex = tbubble_dif.spinvertex;
//       two_loop.vertex.densvertex.tvertex = tbubble_dif.densvertex; //t-channel contribution to differentiated vertex
//       cout << "finished t-channel" << endl;
//       ////     // u-channel
//        ububble_dif = ububble(1,Lambda,one_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy) + ububble(2,Lambda,state_i.vertex,one_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
//       two_loop.vertex.spinvertex.uvertex = ububble_dif.spinvertex;
//       two_loop.vertex.densvertex.uvertex = ububble_dif.densvertex; //u-channel contribution to differentiated vertex
//       cout << "finished u-channel" << endl;

    cout << "finished bubbles.." << endl;
    cout << "Computing self energy corrections.." << endl;

    dstate_f = dstate_f+ one_loop + two_loop;

    cout << "Evaluation finished successfully. Ready for integration step.." << endl;
    return dstate_f;
}




state integrate(double Lambda_i,double Lambda_f, state& state_i){

    /*integration of vertex via Runge-Kutta-4*/
    double h = Lambda_f - Lambda_i; //negative step size since we are integrating downwards
    state k1 = h * evaluate(Lambda_i,state_i);
    state k2 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k1);
    state k3 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k2);
    state k4 = h * evaluate(Lambda_f,state_i + k3);//Note: Lambda_f = Lambda_i + h
    cout << "integration step completed.." << endl;
    return (state_i + 1./6 * (k1 + 2*k2 + 2*k3 + k4)); //returns state at Lambda_f
}



state integrate_euler(double Lambda_i,double Lambda_f, state& state_i){

    /*integration of vertex via Runge-Kutta-4*/
    double h = Lambda_f - Lambda_i; //negative step size since we are integrating downwards
    return (state_i + h*evaluate(Lambda_i,state_i));
} 








template<typename  T1,typename  T2>
struct lambda_params_bubble{

    T1& vert1;
    T2& vert2;
    char ptype1;
    char ptype2;
    self& selfen;
    self& diffselfen;

    int a; int b; int c;

    int d; int e; int f;

    double s;double w1; double w1p;

};




template<typename T1,typename T2>
double lambda_int_re_bubble(double Lambda, void * p){
    struct lambda_params_bubble<T1,T2> * params
            = static_cast< struct lambda_params_bubble<T1,T2> *>(p);
    complex<double> B;
    char p1 = (params ->ptype1);
    char p2 = (params ->ptype2);
    self& se = (params->selfen);
    self& dse = (params->diffselfen);

    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w1p = (params->w1p);

    return real(sbubble(Lambda,vert1,a,b,c,vert2,d,e,f,p1,p2,se,dse,s,w1,w1p));

}





inline bool exists_test (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}
int main(int argc, char *argv[]){

    stringstream stream;
    stream << std::fixed << std::setprecision(2) << argv[1];
    const float alpha = atof(argv[1]);
    const H5std_string	FILE_NAME("5_alpha" + stream.str() + "twoloopsm.h5");


    //**********************initialize frequency grid for ferm and bos case********************************
    // cout << "Initializing linear frequency grid with porp. factor: " << k << endl;

    setprecision(12);
    //set lin grid

//        for(int i=0; i<nw; i++){
//            ffreqs[i] = (i - nw/2 + 0.5) * k ;
//                bfreqs[i] = (i - nw/2 + 0.5) * k ;
//            };


    //set log grid:
    for (int i=0; i<nlog/2; i++){//set positive log part
        ffreqs[i+nw/2] = w0 * pow(k,i);
    };
    for (int i=0 ; i<nlin/2 ; i++){
        ffreqs[i + nw/2 + nlog/2 ] = wt + (i+1)* delw;//set positive linear part
    };
    for(int i=0; i<nw/2; i++){
        ffreqs[i] = - ffreqs[nw-i-1];//set entire negative part by copying positive part with reversed sign
    };

    for(int i=0; i<nw; i++){

        bfreqs[i] = ffreqs[i] ;
    };



    Lambdas.push_back(Lambda_i);
    double  Lambda = Lambda_i;

    while(Lambda > 2){
        Lambda *= 0.95;
        Lambdas.push_back(Lambda);
    };
    while(Lambda > 0.15){
        Lambda *= 0.99;
        Lambdas.push_back(Lambda);
    };






    /************************************************************* initialize bare vertices ********************************************************/

    state state_i;




    double bare = -sqrt(0.5)*cos(alpha*pi/180);
    double bare2= -sqrt(0.5)*sin(alpha*pi/180);
    state_i.vertex.spinvertex.irred.setvert(0,0,2,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(0,0,3,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(-1,0,3,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(0,-1,2,bare);//ititialize vertex with nearest neighbor hopping (J1=1)

    state_i.vertex.spinvertex.irred.setvert(0,-1,3,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(-1,0,2,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(-1,1,3,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
    state_i.vertex.spinvertex.irred.setvert(1,-1,2,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)

    /****PRINT ALL PARAMETERS TO OUTPUT FILE FOR LATER REFERENCE*******************/
    cout << "Writing to filename: " << FILE_NAME << endl;
    cout << "Total number of frequency points: " << nw << " with highest frequencies: " << ffreqs[nw-1] << setw(15) << ffreqs[(nw+nw2)/2-1] << setw(15) << ffreqs[(nw+nw3)/2-1] << endl;
    cout << "k=" << k << endl;
    cout << "Number of Lambda-iterations: " <<Lambdas.size() << "            Lowest Lambda value: " <<   Lambdas[Lambdas.size()-1] << endl;
    cout << "Used regulator: " << reg << endl;
    cout << "Used number of unit cells: " << nuc << endl;
    cout << "Used symmetry mode: " << sym << endl;
    cout << "wt: " << wt << "  w0: " << w0 << "  nlog: " << nlog << "   nw1: " << nw1 << "  nw2:" << nw2 << "  nw3: " << nw3 << endl;
    //cout << "  w0: " << w0 << "  k: " << k << "   nw1: " << nw1 << "  nw2:" << nw2 << "  nw3: " << nw3 << endl;
    cout << "Initial parameters: " << "J1=" << bare << "   ,J2=" << bare2 << endl;

    // state_i.verttbubble1(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + ex.spinvertex.irred.setvert(0,0,1,bare);//ititialize vertex with nearest neighbor hopping (J1=1)

    //   /*test for  lambda-derivative of bubbles**/
    double bubble_start;
    bubble_start = -0.5/(1-1/(pi*Lambda_i))+0.5;

    //double bubble_full_start_spin = -bare+bare/(1-(2*bare-3*bare*bare/(3.1415*Lambda_i))*(-1/(3.1415*Lambda_i)));
    //double bubble_full_start_dens = 3 * bare *(-1/(3.1415*Lambda_i)) *  bare/(1-(2*bare-3*bare*bare/(3.1415*Lambda_i))*(-1/(3.1415*Lambda_i)));


    //   state_i.vertex.spinvertex.uvertex.K1_setvert(0,0,2,nw/2,bubble_full_start_spin);
    //   state_i.vertex.spinvertex.uvertex.K1_setvert(0,0,3,nw/2,bubble_full_start_spin);

    //   state_i.vertex.densvertex.uvertex.K1_setvert(0,0,2,nw/2,bubble_full_start_dens);
    //   state_i.vertex.densvertex.uvertex.K1_setvert(0,0,3,nw/2,bubble_full_start_dens);
    //state_i.vertex.spinvertex.svertex.K1_setvert(-1,0,3,bfreqs[i],bubble_start);
    //state_i.vertex.spinvertex.svertex.K1_setvert(0,-1,2,bfreqs[i],bubble_start);

    //state_i.vertex.spinvertex.svertex.K1_setvert(0,-1,3,bfreqs[i],bubble_start);
    //state_i.vertex.spinvertex.uvertex.K1_setvert(-1,0,2,bfreqs[nw/2],bubble_start);
    //state_i.vertex.spinvertex.uvertex.K1_setvert(-1,1,3,bfreqs[nw/2],bubble_start);
    //state_i.vertex.spinvertex.svertex.K1_setvert(1,-1,2,bfreqs[i],bubble_start);


    ////    //differentiated and reintegrated:
    //    double result=0;
    //    gsl_function F;
    //    gsl_integration_workspace * w = gsl_integration_workspace_alloc (250);
    //    double result_re, error_re;
    //    struct lambda_params_bubble<irreducible,irreducible> lparams = {state_i.vertex.spinvertex.irred, state_i.vertex.spinvertex.irred,'g','s',state_i.selfenergy, dself0, 0,0,2,0,0,2 ,0,0,0};
    //    F.function = &lambda_int_re_bubble<irreducible,irreducible>;
    //    F.params = &lparams;
    //    gsl_integration_qag(&F,5,20, 1e-12, 1e-12,250,2,//adaptive integration in interval (5,20)
    //                         w, &result_re, &error_re);
    //    result += result_re;
    //    cout << "result with single scale propagator in propagator 1: " << result_re << endl;
    //    struct lambda_params_bubble<irreducible,irreducible> lparams2 = {state_i.vertex.spinvertex.irred, state_i.vertex.spinvertex.irred,'s','g',state_i.selfenergy, dself0, 0,0,2,0,0,2 ,0,0,0};
    //    F.params = &lparams2;
    //    gsl_integration_qag(&F,5,20, 1e-12, 1e-12,250,2,//adaptive integration in interval (5,20)
    //                        w, &result_re, &error_re);
    //cout <<  "result with single scale propagator in propagator 2: " <<result_re << endl;
    //    result += result_re;
    //    cout << "sum of reintegrated bubbles" << setw(30) << result << endl;
    //cout << "relative error" << setw(45) << (result-a)/a << endl;
    //    gsl_integration_workspace_free(w);


    int iterator = 0 ;
    int iter_start=0;


    string state_file = static_cast<string>(FILE_NAME);

    state_file.erase(state_file.end()-3, state_file.end());
    state_file += ".txt";




    state integrated = state_i;




    if(exists_test(state_file)==true){

        ifstream myfile (state_file);

        myfile >> iter_start >> iterator ;

        cout << "Program was restarted. Loading intermediate result from iteration " << iter_start <<"\n"<< "Continuing at iteration number: " <<iter_start+1 << endl;
        integrated = read_hdf(FILE_NAME,iterator,Lambdas.size()/3);

    }

    else{
        write_hdf(FILE_NAME,0,Lambdas.size()/3,state_i);
        ofstream myfile;
        myfile.open (state_file);
        myfile << 0  << '\n' << 0 << endl;

        myfile.close();
    };




    for(int i=iter_start+1; i<Lambdas.size();i++)
    {
        cout << "iteration number " << i << " with Lambda " << Lambdas[i] << endl;
        try{
            integrated = integrate(Lambdas[i-1],Lambdas[i],integrated);
            integrated.sus = suscept(Lambdas[i],integrated.vertex,integrated.selfenergy);
            if(i % 3 == 0){

                iterator += 1;
                add_hdf(FILE_NAME,iterator,Lambdas.size()/3, integrated);

                ofstream myfile;
                myfile.open (state_file, ios::out | ios::trunc );
                myfile << i <<"\n" << iterator;
                myfile.close();
            };


////                       integrated = integrate(Lambdas[i-1],Lambdas[i],integrated);
////                        integrated.sus = suscept(Lambdas[i],integrated.vertex,integrated.selfenergy);
////                        add_hdf(FILE_NAME,1,Lambdas.size()/3, integrated);

////                  double h = Lambdas[i] - Lambdas[i-1]; //negative step size since we are integrating downwards
////                state k1 = h * evaluate(Lambdas[i-1],integrated);
////             add_hdf(FILE_NAME,1,7,k1);
////               state k2 = h * evaluate(Lambdas[i-1]+h/2,integrated + 0.5 * k1);
////            add_hdf(FILE_NAME,2,7,k2);
////            state k3 = 0.5*k1;

////                state k3 = h * evaluate(Lambdas[i-1]+h/2,integrated + 0.5 * k2);
////               add_hdf(FILE_NAME,3,7,k3);
////             state k4 = h * evaluate(Lambdas[i],integrated + k3);//Note: Lambda_f = Lambda_i + h
////              add_hdf(FILE_NAME,4,7,k4);
////            ////   cout << "integration step completed.." << endl;
////            state final_state = (integrated + 1./6 * (k1 + 2*k2 + 2*k3 + k4));
////             add_hdf(FILE_NAME,5,7,final_state);
////            //// cout << final_state.vertex.spinvertex.irred.vval(0,0,2)<< endl;
////             final_state.sus = suscept(Lambdas[i],final_state.vertex,final_state.selfenergy);
////              add_hdf(FILE_NAME,6,7,final_state);

            //}
        }
        catch(...){cout << "exception raised. Program terminated" << endl; return 0;}
    };


    //for(int i=0; i<nw3; i++){
    //    double a = ffreqs[(nw-nw3)/2+i] ;
    //    double w = -a;
    //    cout << w << " " << fconv_n(w,nw3) << " " << nw3-1-fconv_n(ffreqs[(nw-nw3)/2+i],nw3) << " " << ffreqs[(nw-nw3)/2+fconv_n(w,nw3)] << endl;
    //};
    //cout << fconv_n(0.05,nw3) << endl;
    /*

    write_hdf(0,6,state_i);
    double h = Lambdas[1] - Lambda_i; //negative step size since we are integrating downwards
    state k1 = h * evaluate(Lambda_i,state_i);
    add_hdf(1,6, k1);


     state k2 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k1);
    add_hdf(2,6, k2);
    
       state k3 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k2);
    add_hdf(3,6, k3);
    state k4 = h * evaluate(Lambdas[1],state_i + k3);//Note: Lambda_f = Lambda_i + h
    add_hdf(4,6, k4);
    cout << "integration step completed.." << endl;
state state_f = state_i + 1./6 * (k1 + 2*k2 + 2*k3 + k4); //
add_hdf(5,6, state_f); */

    //for(int i=(nw-nw2)/2; i<(nw+nw2)/2-1; i++){
    //    cout << (nw+nw2)/2-1-i << " " << ffreqs[i] <<  " " << (nw+nw2)/2-1-i-fconv_n(-ffreqs[i],nw2) << " " <<endl;
    //};

    //    for(int i=0; i<nw1; i++){
    //        cout << i << "  " << ffreqs[i] << " " <<i- fconv_n(ffreqs[i],nw1) << " " <<  ffreqs[i]-ffreqs[fconv_n(ffreqs[i],nw1)] << endl;
    //    };
    //for(int i=0; i<50; i++){
    //    cout << i*0.8 << " " << fconv_n(i*0.8,nw1) << " " << ffreqs[fconv_n(i*0.8,nw1)] << " " << ffreqs[fconv_n(i*0.8,nw1)-1] << " " << ffreqs[fconv_n(i*0.8,nw1)+1] << endl;
    //};
    return 0;

}































