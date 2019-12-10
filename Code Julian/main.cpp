#include <iostream>
#include<iomanip>
#include<cmath>
#include<fstream> 
#include<type_traits>
#include <string>
//#include<tgmath.h>
#include "kagome_facil.hpp"
#if MPI==1
#include "bubbles_mpi.hpp"
#else
#include "bubbles.hpp"
#endif
#include <sstream>
#include<ctime>
using namespace std;
// using namespace boost::numeric::odeint;

//general contants:

const complex<double> zero(0.0,0.0);//zero written as complex number with real and imag part =0.
const complex<double> ci(0.0,1.0); //complex i

const double pi = 3.14159265359;

/********************************constants concerning the regularizer*************************/
const double sharp = 2;//governs the sharpness of the second cutoff


/********************************constants concerning the lattice*************************/

//The cutoff distance from the reference site d_c is then given by (nuc-1)/2 (in units of lattice spacings)
const double d_c = (nuc-1)/2;

/********************************constants concerning the frequency grid*************************/
vector<double> Lambdas;//grid for RG flow parameter lambda
const double Lambda_i = 10; //initial value of flow parameter
int nw = nw1;

#if grid==1
//linear grid
const double k =0.8;//some arbitrary factor that sets the density of the stored frequencies. Recommended: 0.8
const double w0 = 0.5*k;
const double wmax = nw/2*k;//value of maximal frequency that kan be reached for fixed k and nw.

#elif grid==2
//log grid:
const double wt = 6.5;//transition frequency. This freq is considered to be the last freq on the log grid. Its value should not be too high for good resolution
const double w0 = 0.05;//lowest freq on grid
const int nlog = 16;//number of freqs on log part . should be an even number
const double k = pow((wt/w0),2./(nlog-2));//exponent for log part
const double delw = wt - w0*pow(k,(nlog/2-2));//spacing between wt and next lower freq
const int nlin = nw - nlog;//compute number of points on linear parts of grid
#endif


vector<double> ffreqs(nw);//fermionic matsubara freqs
vector<double> bfreqs(nw);//bosnonic matsubara freqs, NOTE: in the case T=0 (which is implemented here), these two grids are equivalent


state evaluate( double Lambda,state state_i, int n_loops){//this fucntion evaluates the RHS of the flow equation at iteration lambda

    state  m_loop,dstate_f;//result of this function
 //   dstate_f.Lambda = Lambda;


    dstate_f.selfenergy =  loop(0,Lambda, state_i.vertex,'s',state_i.selfenergy,state_i.selfenergy);//yields differentiated self energy. Note that the last argument is not used since this is always the single scale propagator
    int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if(world_rank==0){

    cout << "finished self energy" << endl;
    cout << "computing bubbles.." << endl;};

    //one loop:

    //s-channel

      if(world_rank==0){
    cout << "1. loop order:"<< endl;};
      parvert<svert> sbubble_dif;//first argument specifies whcih of the vertices only contain complementary vertices
      sbubble_dif = sbubble(0,Lambda,state_i.vertex,state_i.vertex,'g','k',state_i.selfenergy,dstate_f.selfenergy);//+sbubble(0,Lambda,state_i.vertex,state_i.vertex,'k','g',state_i.selfenergy,dstate_f.selfenergy);
      dstate_f.vertex.spinvertex.svertex = sbubble_dif.spinvertex ;
      dstate_f.vertex.densvertex.svertex = sbubble_dif.densvertex ; //s-channel contribution to differentiated vertexs
    if(world_rank==0){
    cout << "finished s-channel" << endl;
   }

//    //    //  t-channel
    parvert<tvert> tbubble_dif;
    tbubble_dif = tbubble(0,Lambda,state_i.vertex,state_i.vertex,'g','k',state_i.selfenergy,dstate_f.selfenergy);//+ tbubble(0,Lambda,state_i.vertex,state_i.vertex,'k','g',state_i.selfenergy,dstate_f.selfenergy);
    dstate_f.vertex.spinvertex.tvertex = tbubble_dif.spinvertex ;
    dstate_f.vertex.densvertex.tvertex = tbubble_dif.densvertex ; //t-channel contribution to differentiated vertex
    if(world_rank==0){
    cout << "finished t-channel" << endl;};


    ////     // u-channel
  parvert<uvert> ububble_dif;
    ububble_dif = ububble(0,Lambda,state_i.vertex,state_i.vertex,'g','k',state_i.selfenergy,dstate_f.selfenergy);//+ububble(0,Lambda,state_i.vertex,state_i.vertex,'k','g',state_i.selfenergy,dstate_f.selfenergy);
    dstate_f.vertex.spinvertex.uvertex = ububble_dif.spinvertex;
    dstate_f.vertex.densvertex.uvertex = ububble_dif.densvertex; //u-channel contribution to differentiated vertex
    if(world_rank==0){
    cout << "finished u-channel" << endl;};

    //   //two loops:

     parvert<svert> sbubble_dif_left, sbubble_dif_right;
      parvert<tvert> tbubble_dif_left, tbubble_dif_right;
        parvert<uvert> ububble_dif_left, ububble_dif_right;
    if(n_loops >1){
    if(world_rank==0){
        cout << "2. loop order:"<< endl;};
    //   //s-channel
    sbubble_dif_left = sbubble(1,Lambda,dstate_f.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    sbubble_dif_right = sbubble(2,Lambda,state_i.vertex,dstate_f.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    m_loop.vertex.spinvertex.svertex = sbubble_dif_left.spinvertex + sbubble_dif_right.spinvertex;
    m_loop.vertex.densvertex.svertex = sbubble_dif_left.densvertex + sbubble_dif_right.densvertex; //s-channel contribution to differentiated vertex
   if(world_rank==0){ cout << "finished s-channel" << endl;};
    //    //  t-channel

    tbubble_dif_left = tbubble(1,Lambda,dstate_f.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    tbubble_dif_right =tbubble(2,Lambda,state_i.vertex,dstate_f.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    m_loop.vertex.spinvertex.tvertex = tbubble_dif_left.spinvertex + tbubble_dif_right.spinvertex ;
    m_loop.vertex.densvertex.tvertex = tbubble_dif_left.densvertex + tbubble_dif_right.densvertex; //t-channel contribution to differentiated vertex
    if(world_rank==0){ cout << "finished t-channel" << endl;};
    ////     // u-channel

    ububble_dif_left = ububble(1,Lambda,dstate_f.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    ububble_dif_right =  ububble(2,Lambda,state_i.vertex,dstate_f.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
    m_loop.vertex.spinvertex.uvertex = ububble_dif_left.spinvertex + ububble_dif_right.spinvertex;
    m_loop.vertex.densvertex.uvertex = ububble_dif_left.densvertex + ububble_dif_right.densvertex; //u-channel contribution to differentiated vertex
    if(world_rank==0){ cout << "finished u-channel" << endl;};
};
    dstate_f.vertex = dstate_f.vertex + m_loop.vertex;
    parvert<fullvert> center;
    if(n_loops>2){


parvert<fullvert> self_center;
    //multiloops:
    for(int l=3; l<n_loops+1; l++){
        if(world_rank==0){
            cout << l<<". loop order:"<< endl;};
        //s-channel:
        parvert<svert> scenter_buffer = sbubble(0,Lambda,state_i.vertex,sbubble_dif_left,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        center.spinvertex.svertex = scenter_buffer.spinvertex;
        center.densvertex.svertex = scenter_buffer.densvertex;
        sbubble_dif_left = sbubble(1,Lambda,m_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        sbubble_dif_right = sbubble(2,Lambda,state_i.vertex,m_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        m_loop.vertex.spinvertex.svertex = center.spinvertex.svertex + sbubble_dif_left.spinvertex + sbubble_dif_right.spinvertex;
        m_loop.vertex.densvertex.svertex = center.densvertex.svertex + sbubble_dif_left.densvertex + sbubble_dif_right.densvertex;

if(world_rank==0){ cout << "finished s-channel" << endl;};
        //t-channel:
        parvert<tvert> tcenter_buffer = tbubble(0,Lambda,state_i.vertex,tbubble_dif_left,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        center.spinvertex.tvertex = tcenter_buffer.spinvertex;
        center.densvertex.tvertex = tcenter_buffer.densvertex;
        tbubble_dif_left = tbubble(1,Lambda,m_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        tbubble_dif_right = tbubble(2,Lambda,state_i.vertex,m_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        m_loop.vertex.spinvertex.tvertex = center.spinvertex.tvertex + tbubble_dif_left.spinvertex + tbubble_dif_right.spinvertex;
        m_loop.vertex.densvertex.tvertex = center.densvertex.tvertex + tbubble_dif_left.densvertex + tbubble_dif_right.densvertex;
if(world_rank==0){ cout << "finished t-channel" << endl;};
        //u-channel:
        parvert<uvert> ucenter_buffer = ububble(0,Lambda,state_i.vertex,ububble_dif_left,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        center.spinvertex.uvertex = ucenter_buffer.spinvertex;
        center.densvertex.uvertex = ucenter_buffer.densvertex;
        ububble_dif_left = ububble(1,Lambda,m_loop.vertex,state_i.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        ububble_dif_right = ububble(2,Lambda,state_i.vertex,m_loop.vertex,'g','g',state_i.selfenergy,dstate_f.selfenergy);
        m_loop.vertex.spinvertex.uvertex = center.spinvertex.uvertex + ububble_dif_left.spinvertex + ububble_dif_right.spinvertex;
        m_loop.vertex.densvertex.uvertex = center.densvertex.uvertex + ububble_dif_left.densvertex + ububble_dif_right.densvertex;
        if(world_rank==0){ cout << "finished u-channel" << endl;};

        dstate_f.vertex = dstate_f.vertex + m_loop.vertex;
        self_center = self_center + center;
    };
    if(world_rank==0){
        cout << "finished bubbles.." << endl;
        cout << "Computing self energy corrections.." << endl;};
    self dself_tbar = loop(1,Lambda,self_center,'g',state_i.selfenergy,dstate_f.selfenergy);//the last argument specifies that only t-reduced contributions are read out from the vertex.
    if(world_rank==0){
        cout <<"fninished 1. self energy correction"<< endl;};
    self dself_t = loop(0,Lambda, state_i.vertex,'e',state_i.selfenergy,dself_tbar);
    if(world_rank==0){
        cout <<"fninished 2. self energy correction"<< endl;};
    dstate_f.selfenergy = dstate_f.selfenergy + dself_tbar + dself_t;
};



    if(world_rank==0){
    cout << "Evaluation finished successfully. Ready for integration step.." << endl;};

    return dstate_f;
}




state integrate(double Lambda_i,double Lambda_f, state& state_i,int n_loops){
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    /*integration of vertex via Runge-Kutta-4*/
    double h = Lambda_f - Lambda_i; //negative step size since we are integrating downwards
    state k1 = h * evaluate(Lambda_i,state_i,n_loops);
    state k2 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k1,n_loops);
    state k3 = h * evaluate(Lambda_i+h/2,state_i + 0.5 * k2,n_loops);
    state k4 = h * evaluate(Lambda_f,state_i + k3,n_loops);//Note: Lambda_f = Lambda_i + h
    if(world_rank==0){cout << "integration step completed.." << endl;};
    return (state_i + 1./6 * (k1 + 2*k2 + 2*k3 + k4)); //returns state at Lambda_f

}



state integrate_euler(double Lambda_i,double Lambda_f, state& state_i,int n_loops){

    /*integration of vertex via Runge-Kutta-4*/
    double h = Lambda_f - Lambda_i; //negative step size since we are integrating downwards
    return (state_i + h*evaluate(Lambda_i,state_i,n_loops));
} 


void rkck(state &state_i, state&dydx,state & result, double Lambda_i, double h, state& err,int n_loops){//Cash-Carp RK-step (cf. Numerical recipes, page 724)

    static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
                 b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3,b42 = -0.9, b43=1.2, b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,
                 b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
                 c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0, dc5 = -277.00/14336.0, dc1=c1-2825.0/27648.0,
                 dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc6=c6-0.25;

    state ak2,ak3,ak4,ak5,ak6,temporary;


    temporary = state_i + b21*h*dydx;
   ak2 = evaluate(Lambda_i+a2*h,temporary,n_loops);
   temporary=state_i+h*(b31*dydx+b32*ak2);
   ak3 = evaluate(Lambda_i+a3*h,temporary,n_loops);
   temporary=state_i+h*(b41*dydx+b42*ak2+b43*ak3);
   ak4 = evaluate(Lambda_i+a4*h,temporary,n_loops);
   temporary=state_i+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4);
   ak5 = evaluate(Lambda_i+a5*h,temporary,n_loops);
   temporary=state_i+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5);
   ak6 = evaluate(Lambda_i+a6*h,temporary,n_loops);
   result = state_i+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6);
   err = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);

}

void rkqs(state &state_i, state&dydx,state &yscal,double &Lambda_i, const double htry,double &hdid, double&hnext,const double eps, int n_loops){//quality-controlled RK-step (cf. Numerical recipes, page 723)
    const double SAFETY = 0.9, PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.e-3;
    double errmax,h,htemp,Lambda_new;
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    h= htry;
    state err,temporary;
    for(;;){
      if(world_rank==0){  cout << "Try stepsize " << h << endl;};
    rkck(state_i,dydx,temporary,Lambda_i,h,err,n_loops);



    errmax=max_err(err,yscal);
    errmax =errmax/eps;
     if(world_rank==0){   cout << "errmax: " <<errmax << endl;};
    if(errmax<=1.0){break;};
       if(world_rank==0){  cout << "Stepsize too big. Readjust.."<< endl;};
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h>=0.0? max({htemp,0.1*h},comp) : min({htemp,0.1*h},comp));
    Lambda_new = Lambda_i + h;
    if(Lambda_new==Lambda_i){cout << "Stepsize underflow in rkqs" << endl;};
    };
    if(errmax>ERRCON){hnext = SAFETY * h * pow(errmax,PGROW);}
    else{hnext = 1.2* h;};
    Lambda_i += (hdid=h);

    state_i=temporary;
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





//inline bool exists_test (const std::string& name) {
//    ifstream f(name.c_str());
//    return f.good();
//}
int main(int argc, char *argv[]){
clock_t begin = clock();
    MPI_Init(NULL, NULL);

//       // Get the number of processes
       int world_size;
       MPI_Comm_size(MPI_COMM_WORLD, &world_size);

       // Get the rank of the process
       int world_rank;
       MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

       // Get the name of the processor
       char processor_name[MPI_MAX_PROCESSOR_NAME];
       int name_len;
       MPI_Get_processor_name(processor_name, &name_len);
if(world_rank==0){
      cout << "Program executed on " << world_size << " processors." << endl;};
int n_loops = atoi(argv[2]);

    stringstream stream;
    stream << std::fixed << std::setprecision(2) << argv[1];
    const float alpha = atof(argv[1]);

#if temp==0
#if grid==1
    const H5std_string	FILE_NAME( to_string(nuc) + "sites_" + stream.str() + "degs_" + to_string(n_loops)+ "loop_" + to_string(nw3) + "nw3_lin.h5");
#elif grid==2
     const H5std_string	FILE_NAME(to_string(nuc) + "sites_" + stream.str() + "degs_" + to_string(n_loops)+ "loop_" + to_string(nw3) + "nw3_log.h5");
#endif
#elif temp==1
     const H5std_string	FILE_NAME(to_string(nuc) + "sites_" + stream.str() + "degs_" + to_string(n_loops)+ "loop_" + to_string(nw3) + "nw3_fin.h5");
#endif
    //**********************initialize frequency grid for ferm and bos case********************************
    // cout << "Initializing linear frequency grid with porp. factor: " << k << endl;

    setprecision(16);
#if temp==0
#if grid==1
    //set lin grid

        for(int i=0; i<nw; i++){
            ffreqs[i] = (i - nw/2 + 0.5) * k ;
                bfreqs[i] = (i - nw/2 + 0.5) * k ;
            };


#elif grid==2
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

#endif
#endif


    int iterator = 0 ;
     int iterator_sus = 0 ;
    int iter_start=0;
    int L_counts_sus = 240;//numver of saved steps for susceptibility
    int L_counts = 30;//number of saved steps for vertex

//    Lambdas.push_back(Lambda_i);
//    double  Lambda = Lambda_i;


//    while(Lambda > 4){
//        Lambda *= 0.9;
//        Lambdas.push_back(Lambda);

//    };
//    while(Lambda > 0.05){
//        Lambda *= 0.95;
//        Lambdas.push_back(Lambda);

//    };

//        while(Lambda > 0.8){
//            Lambda *= 0.99;
//            Lambdas.push_back(Lambda);
//        };
//    while(Lambda > 0.1){
//        Lambda *= 0.99;
//        Lambdas.push_back(Lambda);
//    };


    //count number of lambda iterations:
//    int L_counts=1;//first entry for initial state
//    int L_counts_sus =1;
//    for(int i=1; i<Lambdas.size();i++){
//        if(i%8==0){
//            L_counts += 1;
//        };
//        if(i%2){
//            L_counts_sus +=1;
//        };
//    };
    int site_counter=0;
 //count number of lattice sites
    if(world_rank==0){
    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= -(nuc_eff-1)/2; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){

                if( distance(a,b,c) <= d_c){
                    site_counter +=1;
                };};};};

    };

    double bare = 0;
    if(abs(-sqrt(0.5)*cos(alpha*pi/180))>1e-6){bare = -sqrt(0.5)*cos(alpha*pi/180);};
    double bare2= 0;
    if(abs(-sqrt(0.5)*sin(alpha*pi/180))>1e-6){bare2 = -sqrt(0.5)*sin(alpha*pi/180);};

    /************************************************************* initialize bare vertices ********************************************************/



    /****PRINT ALL PARAMETERS TO OUTPUT FILE FOR LATER REFERENCE*******************/
    if(world_rank==0){
    cout << "Writing to filename: " << FILE_NAME << endl;
    cout << "Total number of frequency points: " << nw << " with highest frequencies: " << ffreqs[nw-1] << setw(15) << ffreqs[(nw+nw2)/2-1] << setw(15) << ffreqs[(nw+nw3)/2-1] << endl;

  //  cout << "Number of Lambda-iterations: " <<Lambdas.size() << "    " << "Number of stored Lambda-iterations:  " << L_counts << "      Lowest Lambda value: " <<   Lambdas[Lambdas.size()-1] << endl;
    cout << "Used number of unit cells: " << nuc << " (nuc_eff =" << nuc_eff << "), which corresponds to a total of " << site_counter<<" lattice sites" <<  endl;
#if temp==0
#if grid==1
   cout <<"Used lin grid with settings:  k:" << k  << "  w0: " << w0 << "   nw1: " << nw1 << "  nw2:" << nw2 << "  nw3: " << nw3 << endl;
#elif grid==2

   cout <<"Used log+lin grid with settings:  wt: " << wt << "  w0: " << w0 << "  nlog: " << nlog << "  k:" << k <<endl;
#endif
#endif
#if reg==1
    cout << "Regulator: Heavyside Function" << endl;
#elif reg==2
    cout << "Regulator: Smoothened Step" << endl;
#endif
cout << "Numver of frequency points: nw1= " << nw1 << ", nw2 = " << nw2 << " nw3 = " << nw3 << endl;
#if temp==1
cout << "TEMPERATURE: " << T << endl;
cout << "Biggest frequencies: " << (nw1-1)*pi*T << " " <<(nw2-1)*pi*T << " " << (nw3-1)*pi*T << endl;
#endif
#if sym==0
     cout << "Symmetry mode: 0" << endl;
#elif sym==1
    cout << "Symmetry mode: 1" << endl;
#elif sym==2
    cout << "Symmetry mode: 2" << endl;
#endif
    cout << "Initial parameters: " << "J1=" << bare << "   ,J2=" << bare2 << endl;
    cout << "Number of loops: " << n_loops << endl;
};
    // state_i.verttbubble1(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + ex.spinvertex.irred.setvert(0,0,1,bare);//ititialize vertex with nearest neighbor hopping (J1=1)

    //   /*test for  lambda-derivative of bubbles**/
    double bubble_start;
    bubble_start = -0.5/(1-1/(pi*Lambda_i))+0.5;

double Lambda;
vector<double> Lambdas(L_counts_sus);
Lambdas[0]=Lambda_i;
   double h = -0.3;
    double hnext,hdid;
    string state_file = static_cast<string>(FILE_NAME);

    state_file.erase(state_file.end()-3, state_file.end());
    state_file += ".txt";




    state integrated;




    if(exists_test(state_file)==true){
for(int i=0; i<world_size; i++){
    if(world_rank==i){
        ifstream myfile (state_file);
        myfile >> iter_start >> iterator >> iterator_sus >> Lambda>>h;

     myfile.close();};
    MPI_Barrier(MPI_COMM_WORLD);
};

if(world_rank==0){
        cout << "Program was restarted. Loading intermediate result from iteration " << iter_start <<"\n"<< "Continuing at iteration number: " <<iter_start+1 << " with Lambda= "<<Lambda<< endl;
        integrated = read_hdf(FILE_NAME,iterator,L_counts,L_counts_sus,Lambdas);
};
 MPI_Bcast(&Lambdas[0],Lambdas.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);

 state_bcast(integrated);//boradcast read out data to all processes
    }



    else{
        state state_i;


        state_i.vertex.spinvertex.irred.setvert(0,0,2,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
        state_i.vertex.spinvertex.irred.setvert(0,0,3,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
       // state_i.vertex.spinvertex.irred.setvert(-1,0,3,bare);//ititialize vertex with nearest neighbor hopping (J1=1)
       // state_i.vertex.spinvertex.irred.setvert(0,-1,2,bare);//ititialize vertex with nearest neighbor hopping (J1=1)


     //   state_i.vertex.spinvertex.irred.setvert(0,-1,3,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
        state_i.vertex.spinvertex.irred.setvert(-1,0,2,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
        state_i.vertex.spinvertex.irred.setvert(-1,1,3,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
      //  state_i.vertex.spinvertex.irred.setvert(1,-1,2,bare2);//ititialize vertex with nearest neighbor hopping (J1=1)
        if(world_rank==0){
        write_hdf(FILE_NAME,Lambda_i,L_counts,L_counts_sus,state_i);
        ofstream myfile;
        myfile.open (state_file);
        myfile << 0  << '\n' << 0 << '\n' << 0  <<"\n" << Lambda_i <<"\n"<< h<<endl;//store starting values for i, iterator, and iterator_sus
        myfile.close();
        };
        integrated = state_i;
  Lambda= Lambda_i;

    };

//if(world_rank==0){ write_hdf(FILE_NAME,0,L_counts,state_i);}
    double eps=0.001;//relative error
    const int MAXSTP = 240;//maximal number of steps that is computed
    const double tiny = 1.0e-30;//to avoid division by zero;


    double Lambda_f = 0.05;

 const double hmin=1e-4;



    int nok=0,nbad=0;
for(int i=iter_start+1; i<MAXSTP; i++){
          if(world_rank==0){
        cout << "iteration number " << i << " with Lambda " << Lambda << endl;};
        try{


            state dydx = evaluate(Lambda,integrated,n_loops);

            state yscal= abs_sum_tiny(integrated,h*dydx, tiny);


            if((Lambda+h-Lambda_f)*(Lambda+h-Lambda_i)>0.0){h=Lambda_f-Lambda;};

            rkqs(integrated, dydx,yscal, Lambda,h,hdid, hnext,eps, n_loops);
            if(hdid==h){++nok;}else{++nbad;};

            if(abs(hnext)>hmin){h=hnext;}
            else{
            cout << "Step size too small "<< endl;
            break;
            };

            if(i%1==0 && iterator_sus < L_counts_sus -1){
                iterator_sus +=1;

            Susc susceptibility = suscept(Lambda,integrated.vertex,integrated.selfenergy);

 Lambdas[i]=Lambda;
 if (world_rank==0){
            add_hdf_sus(FILE_NAME,iterator_sus,Lambdas,L_counts_sus,susceptibility);//save susceptibility after every iteration

};
            };
 MPI_Barrier(MPI_COMM_WORLD);
            if(i % 8 == 0 &&iterator < L_counts -1){//save state only at every mth iteration

                iterator += 1;
                  if (world_rank==0){
                cout << "Storing state in Lambda-layer " << iterator << endl;
                add_hdf(FILE_NAME,iterator,L_counts, integrated);

                ofstream myfile;
                myfile.open (state_file, ios::out | ios::trunc );
               myfile << i <<"\n" << iterator <<"\n" << iterator_sus  <<"\n" <<  std::fixed << std::setprecision(20)<<Lambda<<"\n" <<h;
                myfile.close();

          };};





        }
        catch(...){cout << "exception raised. Program terminated" << endl; return 0;};
};



      MPI_Finalize();
      clock_t end = clock();
       double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       cout << "Time (sec) on processor " << world_rank << "/" << world_size <<" : " << elapsed_secs << endl;

    return 0;

}































