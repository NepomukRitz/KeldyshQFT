
#include <iostream>
#include<iomanip>
#include <complex>
#include<cmath>
#include <vector>
#include<fstream>
#include<type_traits>
#include <string>
#include<tgmath.h>
#include<cstdlib>
#include"kagome.hpp"
//using namespace std;
//using namespace H5;

/********************************constants concerning HDF5 data format*************************/
const int    NX = nw;                     // dataset dimensions
const int    NY = nw;
const int    NZ = nw;                       // dataset dimensions
const int    RANK_R = 8;
const int    RANK_K1 = 6;
const int    RANK_K2 = 7;
const int    RANK_K2b = 7;
const int    RANK_irreducible = 4;
const int    RANK_sus = 4;
const int    RANK_self = 2;
const int nw1 = 140;//put this to the same value as nw
const int nw2 = 100;
const int nw3 = 30;
const int nuc = 5; // numbwe of unit cells in x-direction: should be an odd number!

//define variable nuc_eff that allows to compute arbitrary lattice sizes (artifact of 
//if nuc<=5
const int nuc_eff = nuc;

////if nuc >5
//const int nuc_eff = nuc + 2;



////if sym==0
//const int nw3_q=nw3;
//const int nw3_w1 = nw3;
//const int nw3_w2 = nw3;
//const int nw2_q = nw2;
//const int nw2_w1 = nw2;
//const int nw1_q = nw1;
////if sym==1
//const int nw3_q=nw3/2;
//const int nw3_w1 = nw3/2;
//const int nw3_w2 = nw3;
//const int nw2_q = nw2;
//const int nw2_w1 = nw2/2;
//const int nw1_q = nw1/2;

//if sym=2
const int nw3_q=nw3/2;
const int nw3_w1 = nw3/2;
const int nw3_w2 = nw3;
const int nw2_q = nw2/2;
const int nw2_w1 = nw2/2;
const int nw1_q = nw1/2;

//small routine checking if a is greater than b. Needed for bubble functions
bool comp(int a, int b)
{
    return (a < b);
}


//function that write inital state into hdf5 format.  The second argument denotes the iteration number to which it is written. The third argument denotes the total number of saved iterations in the file
void write_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size,state& state_in){
    //Try block to detect exceptions raised by any of the calls inside it
    cout << "Starting to copy to buffer.." << endl;
    try
    {
        typedef struct complex_t{
            double spin_re;   /*real part */

            double dens_re;   /*real part */
            //   double dens_im;   /*imaginary part */
        }complex_t;

        typedef struct complex{
            double re;   /*real part */
            double im;   /*imaginary part */

        }complex;



        //buffer self energy
        complex selfenergy[nw];//irrdeucible vertex
        for(int i=0; i<nw; i++){
            selfenergy[i].re = real(state_in.selfenergy.sval(i));
            selfenergy[i].im = imag(state_in.selfenergy.sval(i));

        };



        //buffer irreducible vertex:

        complex_t irreducible_class[nuc_eff][nuc_eff][3];//irrdeucible vertex



        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if( distance(a,b,c) <= d_c){//only save sites in upper half of lattice

                        irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].spin_re = state_in.vertex.spinvertex.irred.vval(a,b,c);
                        irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].dens_re = state_in.vertex.densvertex.irred.vval(a,b,c);

                    };};};};





        //buffer all R-class-arrays:


        auto R_class = new complex_t[3][nuc_eff][nuc_eff][3][nw3_q][nw3_w1][nw3_w2];
#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){



                        for(int i=(nw3-nw3_q); i<nw3; i++){
                            for(int j=(nw3-nw3_w1) ; j<nw3; j++){
                               for(int k=(nw3-nw3_w2) ; k<nw3; k++){


                                    //rest class from s-channel
                                    R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.svertex.R_vval(a,b,c,i,j,k);
                                    R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.svertex.R_vval(a,b,c,i,j,k);


                                    //rest class from t-channel
                                    R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.tvertex.R_vval(a,b,c,i,j,k);
                                    R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.tvertex.R_vval(a,b,c,i,j,k);


                                    //rest class from u-channel
                                    R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.uvertex.R_vval(a,b,c,i,j,k);
                                    R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.uvertex.R_vval(a,b,c,i,j,k);



                                };};};};};};};

        //buffer all K1-class arrays:


        auto K1_class = new complex_t[3][nuc_eff][nuc_eff][3][nw1_q];

#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-11)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){

                        for(int i=(nw1-nw1_q) ; i<nw1; i++){


                            //  K1 class from s-channel
                            K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.svertex.K1_vval(a,b,c,i);
                            K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.svertex.K1_vval(a,b,c,i);

                            //K1 class from t-channel
                            K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.tvertex.K1_vval(a,b,c,i);
                            K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.tvertex.K1_vval(a,b,c,i);

                            //K1 class from u-channel
                            K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.uvertex.K1_vval(a,b,c,i);
                            K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.uvertex.K1_vval(a,b,c,i);

                        };};};};};

        //buffer all K2-class-arrays:


        auto K2_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){

                                                for(int i=(nw2-nw2_q) ; i<nw2; i++){
                            for(int j=(nw2-nw2_w1) ; j<nw2; j++){




                                //K2 class from s-channel
                                K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.svertex.K2_vval(a,b,c,i,j);
                                 K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re =state_in.vertex.densvertex.svertex.K2_vval(a,b,c,i,j);


                                //K2 class from t-channel
                                K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.tvertex.K2_vval(a,b,c,i,j);
 K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.tvertex.K2_vval(a,b,c,i,j);


                                //K2 class from u-channel
                                K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.uvertex.K2_vval(a,b,c,i,j);
                                          K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.uvertex.K2_vval(a,b,c,i,j);

                            };};};};};};



        //buffer all K2b-class-arrays:


        auto K2b_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){


                        for(int i=(nw2-nw2_q) ; i<nw2; i++){

                            for(int j=(nw2-nw2_w1) ; j<nw2; j++){


                                //rest class from s-channel
                                K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.svertex.K2b_vval(a,b,c,i,j);
                                      K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re =state_in.vertex.densvertex.svertex.K2b_vval(a,b,c,i,j);


                                //rest class from t-channel
                                K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.tvertex.K2b_vval(a,b,c,i,j);

                                K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.tvertex.K2b_vval(a,b,c,i,j);

                                //rest class from u-channel
                                K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.uvertex.K2b_vval(a,b,c,i,j);

                                K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.uvertex.K2b_vval(a,b,c,i,j);


                            };};};};};};





        //susceptibility
        double Suscept[nuc_eff][nuc_eff][3];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){


                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = state_in.sus.vval(a,b,c);

                    };};};};




        double ferm_freqs[nw];
        for(int i=0; i<nw; i++){
            ferm_freqs[i] = ffreqs[i];
        };

        double bos_freqs[nw];
        for(int i=0; i<nw; i++){
            bos_freqs[i] = ffreqs[i];
        };

        double Lambda_list[Lambdas.size()];
        for(int i=0; i<Lambdas.size(); i++){
            Lambda_list[i] = Lambdas[i];
        };



        cout << "Buffer ready. Preparing for saving into Hdf5 file.." << endl;
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        H5::Exception::dontPrint();

        // Create a new file using the default property lists.
        H5::H5File* file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);


        //Create the memory datatype
        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        //  mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
        //  mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);

        H5::CompType mtype2( sizeof(complex) );
        mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
        mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);

        complex_t fillvalue_vert;
        fillvalue_vert.spin_re = 0;
        fillvalue_vert.dens_re = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype1, &fillvalue_vert);

        complex fillvalue_self;
        fillvalue_self.re = 0;
        fillvalue_self.im = 0;

        H5::DSetCreatPropList plist_self;
        plist_self.setFillValue(mtype2, &fillvalue_self);

        // Create the dimension arrays for objects in file.
        hsize_t R_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2b_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims[]={static_cast<hsize_t>(Lambdas.size()),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};
        hsize_t freq_dims[]={static_cast<hsize_t>(nw)};                      // dataset dimensions for freqs grid (one time nw entries)
        hsize_t bos_freq_dims[]={static_cast<hsize_t>(nw)};
        hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nw)};
        hsize_t Lambda_dims[]={static_cast<hsize_t>(Lambdas.size())};

        // Create the dimension arrays for objects in buffer.
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t R_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};



        // Create the data space for the dataset in file.
        H5::DataSpace dataspacevertex_R(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K2(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
        H5::DataSpace dataspacevertex_K2b(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
        H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_sus(RANK_sus,sus_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacefreqs(1, freq_dims);
        H5::DataSpace dataspacefreqs_bos(1, bos_freq_dims);
        H5::DataSpace dataspaceself(RANK_self, self_dims);
        H5::DataSpace dataspacelambda(1, Lambda_dims);

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_R_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2b_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_sus_buffer(RANK_sus-1, sus_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);



        // Create the datasets in file:
        H5::DataSet* dataset_R;
        dataset_R = new H5::DataSet(file -> createDataSet(DATASET_R,mtype1, dataspacevertex_R,plist_vert));
        H5::DataSet* dataset_K1;
        dataset_K1 = new H5::DataSet(file -> createDataSet(DATASET_K1,mtype1, dataspacevertex_K1,plist_vert));
        H5::DataSet* dataset_K2;
        dataset_K2 = new H5::DataSet(file -> createDataSet(DATASET_K2,mtype1, dataspacevertex_K2,plist_vert));
        H5::DataSet* dataset_K2b;
        dataset_K2b = new H5::DataSet(file -> createDataSet(DATASET_K2b,mtype1, dataspacevertex_K2b,plist_vert));
        H5::DataSet* dataset_irred;
        dataset_irred = new H5::DataSet(file -> createDataSet(DATASET_irred,mtype1, dataspacevertex_irreducible,plist_vert));
        H5::DataSet* dataset_sus;
        dataset_sus = new H5::DataSet(file -> createDataSet(DATASET_sus,H5::PredType::NATIVE_DOUBLE, dataspacevertex_sus));
        H5::DataSet* dataset_ferm_freqs;
        dataset_ferm_freqs = new H5::DataSet(file -> createDataSet(FERM_FREQS_LIST,H5::PredType::NATIVE_DOUBLE, dataspacefreqs));
        H5::DataSet* dataset_bos_freqs;
        dataset_bos_freqs = new H5::DataSet(file -> createDataSet(BOS_FREQS_LIST,H5::PredType::NATIVE_DOUBLE, dataspacefreqs_bos));
        H5::DataSet* dataset_self;
        dataset_self = new H5::DataSet(file -> createDataSet(SELF_LIST,mtype2, dataspaceself,plist_self));
        H5::DataSet* dataset_lambda;
        dataset_lambda = new H5::DataSet(file -> createDataSet(LAMBDA_LIST,H5::PredType::NATIVE_DOUBLE, dataspacelambda));

        //Select hyperslabs:

        //R_class

        //Select hyperslab in the file where the data should be located
        hsize_t R_start[8];
        hsize_t R_stride[8];
        hsize_t R_count[8];
        hsize_t R_block[8];

        R_start[0] = Lambda_it;
        for(int i=1; i<8;i++){R_start[i] = 0;};
        for(int i=0; i<8;i++){
            R_stride[i] = 1;
            R_block[i] = 1;
        };
        R_count[0]= 1;
        R_count[1]= 3;
        R_count[2]= nuc_eff;
        R_count[3]= nuc_eff;
        R_count[4]= 3;
        R_count[5]= nw3_q;
        R_count[6]= nw3_w1;
        R_count[7]= nw3_w2;
        dataspacevertex_R.selectHyperslab(H5S_SELECT_SET, R_count,R_start,R_stride,R_block);
        //Select hyperslab in buffer
        hsize_t R_start_b[7];
        hsize_t R_stride_b[7];
        hsize_t R_count_b[7];
        hsize_t R_block_b[7];

        for(int i=0; i<7;i++){R_start_b[i] = 0;};
        for(int i=0; i<7;i++){
            R_stride_b[i] = 1;
            R_block_b[i] = 1;
        };
        R_count_b[0]= 3;
        R_count_b[1]= nuc_eff;
        R_count_b[2]= nuc_eff;
        R_count_b[3]= 3;
        R_count_b[4]= nw3_q;
        R_count_b[5]= nw3_w1;
        R_count_b[6]= nw3_w2;
        dataspacevertex_R_buffer.selectHyperslab(H5S_SELECT_SET, R_count_b,R_start_b,R_stride_b,R_block_b);




        //K1_class


        //Select hyperslab in the file where the data should be located
        hsize_t K1_start[6];
        hsize_t K1_stride[6];
        hsize_t K1_count[6];
        hsize_t K1_block[6];

        K1_start[0] = Lambda_it;
        for(int i=1; i<6;i++){K1_start[i] = 0;};
        for(int i=0; i<6;i++){
            K1_stride[i] = 1;
            K1_block[i] = 1;
        };
        K1_count[0]= 1;
        K1_count[1]= 3;
        K1_count[2]= nuc_eff;
        K1_count[3]= nuc_eff;
        K1_count[4]= 3;
        K1_count[5]= nw1_q;
        dataspacevertex_K1.selectHyperslab(H5S_SELECT_SET, K1_count,K1_start,K1_stride,K1_block);
        //Select hyperslab in buffer
        hsize_t K1_start_b[5];
        hsize_t K1_stride_b[5];
        hsize_t K1_count_b[5];
        hsize_t K1_block_b[5];

        for(int i=0; i<5;i++){
            K1_start_b[i] = 0;
            K1_stride_b[i] = 1;
            K1_block_b[i] = 1;
        };
        K1_count_b[0]= 3;
        K1_count_b[1]= nuc_eff;
        K1_count_b[2]= nuc_eff;
        K1_count_b[3]= 3;
        K1_count_b[4]= nw1_q;
        dataspacevertex_K1_buffer.selectHyperslab(H5S_SELECT_SET, K1_count_b,K1_start_b,K1_stride_b,K1_block_b);





        //K2_class

        //Select hyperslab in the file where the data should be located
        hsize_t K2_start[7];
        hsize_t K2_stride[7];
        hsize_t K2_count[7];
        hsize_t K2_block[7];

        K2_start[0] = Lambda_it;
        for(int i=1; i<7;i++){K2_start[i] = 0;};
        for(int i=0; i<7;i++){
            K2_stride[i] = 1;
            K2_block[i] = 1;
        };
        K2_count[0]= 1;
        K2_count[1]= 3;
        K2_count[2]= nuc_eff;
        K2_count[3]= nuc_eff;
        K2_count[4]= 3;
        K2_count[5]= nw2_q;
        K2_count[6]= nw2_w1;
        dataspacevertex_K2.selectHyperslab(H5S_SELECT_SET, K2_count,K2_start,K2_stride,K2_block);
        //Select hyperslab in buffer
        hsize_t K2_start_b[6];
        hsize_t K2_stride_b[6];
        hsize_t K2_count_b[6];
        hsize_t K2_block_b[6];

        for(int i=0; i<6;i++){K2_start_b[i] = 0;};
        for(int i=0; i<6;i++){
            K2_stride_b[i] = 1;
            K2_block_b[i] = 1;
        };
        K2_count_b[0]= 3;
        K2_count_b[1]= nuc_eff;
        K2_count_b[2]= nuc_eff;
        K2_count_b[3]= 3;
        K2_count_b[4]= nw2_q;
        K2_count_b[5]= nw2_w1;
        dataspacevertex_K2_buffer.selectHyperslab(H5S_SELECT_SET, K2_count_b,K2_start_b,K2_stride_b,K2_block_b);





        //K2b_class

        //Select hyperslab in the file where the data should be located
        hsize_t K2b_start[7];
        hsize_t K2b_stride[7];
        hsize_t K2b_count[7];
        hsize_t K2b_block[7];

        K2b_start[0] = Lambda_it;
        for(int i=1; i<7;i++){K2b_start[i] = 0;};
        for(int i=0; i<7;i++){
            K2b_stride[i] = 1;
            K2b_block[i] = 1;
        };
        K2b_count[0]= 1;
        K2b_count[1]= 3;
        K2b_count[2]= nuc_eff;
        K2b_count[3]= nuc_eff;;
        K2b_count[4]= 3;
        K2b_count[5]= nw2_q;
        K2b_count[6]= nw2_w1;
        dataspacevertex_K2b.selectHyperslab(H5S_SELECT_SET, K2b_count,K2b_start,K2b_stride,K2b_block);
        //Select hyperslab in buffer
        hsize_t K2b_start_b[6];
        hsize_t K2b_stride_b[6];
        hsize_t K2b_count_b[6];
        hsize_t K2b_block_b[6];

        for(int i=0; i<6;i++){K2b_start_b[i] = 0;};
        for(int i=0; i<6;i++){
            K2b_stride_b[i] = 1;
            K2b_block_b[i] = 1;
        };
        K2b_count_b[0]= 3;
        K2b_count_b[1]= nuc_eff;
        K2b_count_b[2]= nuc_eff;
        K2b_count_b[3]= 3;
        K2b_count_b[4]= nw2_q;
        K2b_count_b[5]= nw2_w1;
        dataspacevertex_K2b_buffer.selectHyperslab(H5S_SELECT_SET, K2b_count_b,K2b_start_b,K2b_stride_b,K2b_block_b);




        //irred_class


        //Select hyperslab in the file where the data should be located
        hsize_t irred_start[4];
        hsize_t irred_stride[4];
        hsize_t irred_count[4];
        hsize_t irred_block[4];

        irred_start[0] = Lambda_it;
        for(int i=1; i<4;i++){irred_start[i] = 0;};
        for(int i=0; i<4;i++){
            irred_stride[i] = 1;
            irred_block[i] = 1;
        };
        irred_count[0]= 1;
        irred_count[1]= nuc_eff;
        irred_count[2]= nuc_eff;
        irred_count[3]= 3;
        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, irred_count,irred_start,irred_stride,irred_block);
        //Select hyperslab in buffer
        hsize_t irred_start_b[3];
        hsize_t irred_stride_b[3];
        hsize_t irred_count_b[3];
        hsize_t irred_block_b[3];
        for(int i=0; i<3;i++){irred_start_b[i] = 0;};
        for(int i=0; i<3;i++){
            irred_stride_b[i] = 1;
            irred_block_b[i] = 1;
        };
        irred_count_b[0]= nuc_eff;
        irred_count_b[1]= nuc_eff;
        irred_count_b[2]= 3;
        dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, irred_count_b,irred_start_b,irred_stride_b,irred_block_b);



        //susceptibility

        //Select hyperslab in the file where the data should be located
        hsize_t sus_start[4];
        hsize_t sus_stride[4];
        hsize_t sus_count[4];
        hsize_t sus_block[4];

        sus_start[0] = Lambda_it;
        for(int i=1; i<4;i++){sus_start[i] = 0;};
        for(int i=0; i<4;i++){
            sus_stride[i] = 1;
            sus_block[i] = 1;
        };
        sus_count[0]= 1;
        sus_count[1]= nuc_eff;
        sus_count[2]= nuc_eff;
        sus_count[3]= 3;

        dataspacevertex_sus.selectHyperslab(H5S_SELECT_SET, sus_count,sus_start,sus_stride,sus_block);
        //Select hyperslab in buffer
        hsize_t sus_start_b[3];
        hsize_t sus_stride_b[3];
        hsize_t sus_count_b[3];
        hsize_t sus_block_b[3];

        for(int i=0; i<3;i++){
            sus_start_b[i] = 0;
            sus_stride_b[i] = 1;
            sus_block_b[i] = 1;
        };

        sus_count_b[0]= nuc_eff;
        sus_count_b[1]= nuc_eff;
        sus_count_b[2]= 3;

        dataspacevertex_sus_buffer.selectHyperslab(H5S_SELECT_SET, sus_count_b,sus_start_b,sus_stride_b,sus_block_b);





        //self energy
        //Select hyperslab in the file where the data should be located
        hsize_t self_start[2];
        hsize_t self_stride[2];
        hsize_t self_count[2];
        hsize_t self_block[2];

        self_start[0] = Lambda_it;
        for(int i=1; i<2;i++){self_start[i] = 0;};
        for(int i=0; i<2;i++){
            self_stride[i] = 1;
            self_block[i] = 1;
        };
        self_count[0]= 1;
        self_count[1]= nw;

        dataspaceself.selectHyperslab(H5S_SELECT_SET, self_count,self_start,self_stride,self_block);
        //Select hyperslab in buffer
        hsize_t self_start_b[1] ;
        hsize_t self_stride_b[1];
        hsize_t self_count_b[1] ;
        hsize_t self_block_b[1];
        self_start_b[0] = 0;
        self_stride_b[0] = 1;
        self_count_b[0] = nw;
        self_block_b[0] = 1;

        dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, self_count_b,self_start_b,self_stride_b,self_block_b);



        //write buffer data into file.


        dataset_R -> write( R_class, mtype1,dataspacevertex_R_buffer,dataspacevertex_R );
        dataset_K1 -> write( K1_class, mtype1 ,dataspacevertex_K1_buffer,dataspacevertex_K1 );
        dataset_K2 -> write( K2_class, mtype1 ,dataspacevertex_K2_buffer,dataspacevertex_K2);
        dataset_K2b -> write( K2b_class, mtype1 ,dataspacevertex_K2b_buffer,dataspacevertex_K2b);
        dataset_irred -> write( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);
        dataset_sus -> write( Suscept, H5::PredType::NATIVE_DOUBLE ,dataspacevertex_sus_buffer,dataspacevertex_sus );


        dataset_ferm_freqs -> write( ferm_freqs, H5::PredType::NATIVE_DOUBLE);
        dataset_bos_freqs -> write( bos_freqs, H5::PredType::NATIVE_DOUBLE);
        dataset_self -> write( selfenergy, mtype2, dataspaceself_buffer,dataspaceself);
        dataset_lambda -> write(Lambda_list, H5::PredType::NATIVE_DOUBLE );


        cout << "Successfully saved in hdf5 file: " << FILE_NAME << endl;
        //free R



        delete[] K1_class;
        delete[] R_class;
        delete[] K2_class;
        delete[] K2b_class;





        delete dataset_R;
        delete dataset_K1;
        delete dataset_K2;
        delete dataset_K2b;
        delete dataset_irred;
        delete dataset_sus;
        delete dataset_ferm_freqs;
        delete dataset_bos_freqs;
        delete dataset_self;
        delete dataset_lambda;


        delete file;

        
    }  // end of try block

    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printErrorStack();
        return ;
    }

    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
    {
        error.printErrorStack();
        return ;
    }

    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error)
    {
        error.printErrorStack();
        return ;
    }


}

//writes state of new iteration into en EXISTING Hdf5 file. The second argument denotes the iteration number to which it is written. The thrid arguemnt denots the total number of saved iterations in the file
void add_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size,state& state_in ){
    //Try block to detect exceptions raised by any of the calls inside it
    cout << "Starting to copy to buffer.." << endl;
    try
    {
        typedef struct complex_t{
            double spin_re;   /*real part */
            double dens_re;   /*real part */

        }complex_t;

        typedef struct complex{//complex data type for self energy
            double re;   /*real part */
            double im;   /*imaginary part */

        }complex;

        //buffer self energy
        complex selfenergy[nw];//irrdeucible vertex
        for(int i=0; i<nw; i++){
            selfenergy[i].re = real(state_in.selfenergy.sval(i));
            selfenergy[i].im = imag(state_in.selfenergy.sval(i));
        };


        //buffer irreducible vertex:

        complex_t irreducible_class[nuc_eff][nuc_eff][3];//irrdeucible vertex


        //split each component of vector into real and imaginary part and save result in array that is fed to HDF5

        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if( distance(a,b,c) <= d_c){//only save sites in upper half


                        irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].spin_re = state_in.vertex.spinvertex.irred.vval(a,b,c);
                        irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].dens_re = state_in.vertex.densvertex.irred.vval(a,b,c);
                                       };};};};








        //Note that the first index labels the channel with the correspondance ( 0 = s-channel, 1 = t-channel, 2 = u-channel)


        //buffer all R-class-arrays:






        auto R_class = new complex_t[3][nuc_eff][nuc_eff][3][nw3_q][nw3_w1][nw3_w2];
#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){

                        for(int i=(nw3-nw3_q); i<nw3; i++){
                                                      for(int j=(nw3-nw3_w1) ; j<nw3; j++){
                                                              for(int k=(nw3-nw3_w2) ; k<nw3; k++){


                                    //rest class from s-channel
                                    R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.svertex.R_vval(a,b,c,i,j,k);
                                   R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.svertex.R_vval(a,b,c,i,j,k);


                                    //rest class from t-channel
                                    R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.tvertex.R_vval(a,b,c,i,j,k);
                                    R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.tvertex.R_vval(a,b,c,i,j,k);


                                    //rest class from u-channel
                                    R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re = state_in.vertex.spinvertex.uvertex.R_vval(a,b,c,i,j,k);
                                    R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re = state_in.vertex.densvertex.uvertex.R_vval(a,b,c,i,j,k);


                                };};};};};};};

        //buffer all K1-class arrays:





        auto K1_class = new complex_t[3][nuc_eff][nuc_eff][3][nw1_q];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){


                        for(int i=(nw1-nw1_q) ; i<nw1; i++){






                            //  K1 class from s-channel
                            K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.svertex.K1_vval(a,b,c,i);
                                    K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.svertex.K1_vval(a,b,c,i);

                            //K1 class from t-channel
                            K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.tvertex.K1_vval(a,b,c,i);
                            K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.tvertex.K1_vval(a,b,c,i);

                            //K1 class from u-channel
                            K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re = state_in.vertex.spinvertex.uvertex.K1_vval(a,b,c,i);
                                         K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re = state_in.vertex.densvertex.uvertex.K1_vval(a,b,c,i);

                             };};};};};

        //buffer all K2-class-arrays:


        auto K2_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){

                        for(int i=(nw2-nw2_q) ; i<nw2; i++){
                                                   for(int j=(nw2-nw2_w1) ; j<nw2; j++){




                                //rest class from s-channel
                                K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.svertex.K2_vval(a,b,c,i,j);
                                                           K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re =state_in.vertex.densvertex.svertex.K2_vval(a,b,c,i,j);

                                                  //rest class from t-channel
                                K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.tvertex.K2_vval(a,b,c,i,j);
                                                                      K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.tvertex.K2_vval(a,b,c,i,j);

                                //rest class from u-channel
                                K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.uvertex.K2_vval(a,b,c,i,j);
                                 K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.uvertex.K2_vval(a,b,c,i,j);

                              };};};};};};



        //buffer all K2b-class-arrays:


        auto K2b_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){

                        for(int i=(nw2-nw2_q) ; i<nw2; i++){
                                               for(int j=(nw2-nw2_w1) ; j<nw2; j++){


                                //rest class from s-channel
                                K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.svertex.K2b_vval(a,b,c,i,j);
                                                     K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re =state_in.vertex.densvertex.svertex.K2b_vval(a,b,c,i,j);

                                //rest class from t-channel
                                K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.tvertex.K2b_vval(a,b,c,i,j);
                                K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.tvertex.K2b_vval(a,b,c,i,j);
                                   //rest class from u-channel
                                K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re = state_in.vertex.spinvertex.uvertex.K2b_vval(a,b,c,i,j);
                                     K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re = state_in.vertex.densvertex.uvertex.K2b_vval(a,b,c,i,j);

                            };};};};};};




        //susceptibility
        double Suscept[nuc_eff][nuc_eff][3];


#pragma omp parallel for
        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
                for (int c=1; c<4 ; c++){
                    if(distance(a,b,c) <= d_c){


                        //  Susceptibility in a-channel with bosonic transfer freq q
                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = state_in.sus.vval(a,b,c);

                    }
                    else{
                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = 0;
                    };};};};



        cout << "Buffer ready. Preparing for saving into Hdf5 file.." << endl;
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        H5::Exception::dontPrint();

        // Open an existing file and dataset.
        H5::H5File file(FILE_NAME, H5F_ACC_RDWR);
        H5::DataSet dataset_irred = file.openDataSet("irred");
        H5::DataSet dataset_R = file.openDataSet("R");
        H5::DataSet dataset_K1 = file.openDataSet("K1");
        H5::DataSet dataset_K2 = file.openDataSet("K2");
        H5::DataSet dataset_K2b = file.openDataSet("K2b");
        H5::DataSet dataset_sus = file.openDataSet("sus");
        H5::DataSet dataset_self = file.openDataSet("selflist");

        //Create the memory datatype
        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        //   mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
        //  mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);

        H5::CompType mtype2( sizeof(complex) );//for self energy
        mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
        mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);

        complex_t fillvalue_vert;
        fillvalue_vert.spin_re = 0;
        fillvalue_vert.dens_re = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype1, &fillvalue_vert);

        complex fillvalue_self;
        fillvalue_self.re = 0;
        fillvalue_self.im = 0;

        H5::DSetCreatPropList plist_self;
        plist_self.setFillValue(mtype2, &fillvalue_self);

        // Create the dimension arrays for objects in file.
        hsize_t R_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2b_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nw)};


        // Create the dimension arrays for objects in buffer.
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t R_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};



        // Create the data space for the dataset in file.
        H5::DataSpace dataspacevertex_R(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1(RANK_K1,K1_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2(RANK_K2,K2_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2b(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_sus(RANK_sus,sus_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspaceself(2, self_dims);

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_R_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2b_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_sus_buffer(RANK_sus-1, sus_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspaceself_buffer(1, self_dims_buffer);

        //Select hyperslabs:

        //R_class

        //Select hyperslab in the file where the data should be located
        hsize_t R_start[8];
        hsize_t R_stride[8];
        hsize_t R_count[8];
        hsize_t R_block[8];

        R_start[0] = Lambda_it;
        for(int i=1; i<8;i++){R_start[i] = 0;};
        for(int i=0; i<8;i++){
            R_stride[i] = 1;
            R_block[i] = 1;
        };
        R_count[0]= 1;
        R_count[1]= 3;//three different channels
        R_count[2]= nuc_eff;
        R_count[3]= nuc_eff;
        R_count[4]= 3;
        R_count[5]= nw3_q;
        R_count[6]= nw3_w1;
        R_count[7]= nw3_w2;
        dataspacevertex_R.selectHyperslab(H5S_SELECT_SET, R_count,R_start,R_stride,R_block);
        //Select hyperslab in buffer
        hsize_t R_start_b[7];
        hsize_t R_stride_b[7];
        hsize_t R_count_b[7];
        hsize_t R_block_b[7];

        for(int i=0; i<7;i++){R_start_b[i] = 0;};
        for(int i=0; i<7;i++){
            R_stride_b[i] = 1;
            R_block_b[i] = 1;
        };
        R_count_b[0]= 3;
        R_count_b[1]= nuc_eff;
        R_count_b[2]= nuc_eff;
        R_count_b[3]= 3;
        R_count_b[4]= nw3_q;
        R_count_b[5]= nw3_w1;
        R_count_b[6]= nw3_w2;
        dataspacevertex_R_buffer.selectHyperslab(H5S_SELECT_SET, R_count_b,R_start_b,R_stride_b,R_block_b);




        //K1_class


        //Select hyperslab in the file where the data should be located
        hsize_t K1_start[6];
        hsize_t K1_stride[6];
        hsize_t K1_count[6];
        hsize_t K1_block[6];

        K1_start[0] = Lambda_it;
        for(int i=1; i<6;i++){K1_start[i] = 0;};
        for(int i=0; i<6;i++){
            K1_stride[i] = 1;
            K1_block[i] = 1;
        };
        K1_count[0]= 1;
        K1_count[1]= 3;
        K1_count[2]= nuc_eff;
        K1_count[3]= nuc_eff;
        K1_count[4]= 3;
        K1_count[5]= nw1_q;
        dataspacevertex_K1.selectHyperslab(H5S_SELECT_SET, K1_count,K1_start,K1_stride,K1_block);
        //Select hyperslab in buffer
        hsize_t K1_start_b[5];
        hsize_t K1_stride_b[5];
        hsize_t K1_count_b[5];
        hsize_t K1_block_b[5];

        for(int i=0; i<5;i++){K1_start_b[i] = 0;};
        for(int i=0; i<5;i++){
            K1_stride_b[i] = 1;
            K1_block_b[i] = 1;
        };
        K1_count_b[0]= 3;
        K1_count_b[1]= nuc_eff;
        K1_count_b[2]= nuc_eff;
        K1_count_b[3]= 3;
        K1_count_b[4]= nw1_q;
        dataspacevertex_K1_buffer.selectHyperslab(H5S_SELECT_SET, K1_count_b,K1_start_b,K1_stride_b,K1_block_b);





        //K2_class

        //Select hyperslab in the file where the data should be located
        hsize_t K2_start[7];
        hsize_t K2_stride[7];
        hsize_t K2_count[7];
        hsize_t K2_block[7];

        K2_start[0] = Lambda_it;
        for(int i=1; i<7;i++){K2_start[i] = 0;};
        for(int i=0; i<7;i++){
            K2_stride[i] = 1;
            K2_block[i] = 1;
        };
        K2_count[0]= 1;
        K2_count[1]= 3;
        K2_count[2]= nuc_eff;
        K2_count[3]= nuc_eff;
        K2_count[4]= 3;
        K2_count[5]= nw2_q;
        K2_count[6]= nw2_w1;
        dataspacevertex_K2.selectHyperslab(H5S_SELECT_SET, K2_count,K2_start,K2_stride,K2_block);
        //Select hyperslab in buffer
        hsize_t K2_start_b[6];
        hsize_t K2_stride_b[6];
        hsize_t K2_count_b[6];
        hsize_t K2_block_b[6];

        for(int i=0; i<6;i++){K2_start_b[i] = 0;};
        for(int i=0; i<6;i++){
            K2_stride_b[i] = 1;
            K2_block_b[i] = 1;
        };
        K2_count_b[0]= 3;
        K2_count_b[1]= nuc_eff;
        K2_count_b[2]= nuc_eff;
        K2_count_b[3]= 3;
        K2_count_b[4]= nw2_q;
        K2_count_b[5]= nw2_w1;
        dataspacevertex_K2_buffer.selectHyperslab(H5S_SELECT_SET, K2_count_b,K2_start_b,K2_stride_b,K2_block_b);





        //K2b_class

        //Select hyperslab in the file where the data should be located
        hsize_t K2b_start[7];
        hsize_t K2b_stride[7];
        hsize_t K2b_count[7];
        hsize_t K2b_block[7];

        K2b_start[0] = Lambda_it;
        for(int i=1; i<7;i++){K2b_start[i] = 0;};
        for(int i=0; i<7;i++){
            K2b_stride[i] = 1;
            K2b_block[i] = 1;
        };
        K2b_count[0]= 1;
        K2b_count[1]= 3;
        K2b_count[2]= nuc_eff;
        K2b_count[3]= nuc_eff;
        K2b_count[4]= 3;
        K2b_count[5]= nw2_q;
        K2b_count[6]= nw2_w1;
        dataspacevertex_K2b.selectHyperslab(H5S_SELECT_SET, K2b_count,K2b_start,K2b_stride,K2b_block);
        //Select hyperslab in buffer
        hsize_t K2b_start_b[6];
        hsize_t K2b_stride_b[6];
        hsize_t K2b_count_b[6];
        hsize_t K2b_block_b[6];

        for(int i=0; i<6;i++){K2b_start_b[i] = 0;};
        for(int i=0; i<6;i++){
            K2b_stride_b[i] = 1;
            K2b_block_b[i] = 1;
        };
        K2b_count_b[0]= 3;
        K2b_count_b[1]= nuc_eff;
        K2b_count_b[2]= nuc_eff;
        K2b_count_b[3]= 3;
        K2b_count_b[4]= nw2_q;
        K2b_count_b[5]= nw2_w1;
        dataspacevertex_K2b_buffer.selectHyperslab(H5S_SELECT_SET, K2b_count_b,K2b_start_b,K2b_stride_b,K2b_block_b);





        //irred_class


        //Select hyperslab in the file where the data should be located
        hsize_t irred_start[4];
        hsize_t irred_stride[4];
        hsize_t irred_count[4];
        hsize_t irred_block[4];

        irred_start[0] = Lambda_it;
        for(int i=1; i<4;i++){irred_start[i] = 0;};
        for(int i=0; i<4;i++){
            irred_stride[i] = 1;
            irred_block[i] = 1;
        };
        irred_count[0]= 1;
        irred_count[1]= nuc_eff;
        irred_count[2]= nuc_eff;
        irred_count[3]= 3;
        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, irred_count,irred_start,irred_stride,irred_block);
        //Select hyperslab in buffer
        hsize_t irred_start_b[3];
        hsize_t irred_stride_b[3];
        hsize_t irred_count_b[3];
        hsize_t irred_block_b[3];
        for(int i=0; i<3;i++){irred_start_b[i] = 0;};
        for(int i=0; i<3;i++){
            irred_stride_b[i] = 1;
            irred_block_b[i] = 1;
        };
        irred_count_b[0]= nuc_eff;
        irred_count_b[1]= nuc_eff;
        irred_count_b[2]= 3;
        dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, irred_count_b,irred_start_b,irred_stride_b,irred_block_b);


        //susceptibility

        //Select hyperslab in the file where the data should be located
        hsize_t sus_start[4];
        hsize_t sus_stride[4];
        hsize_t sus_count[4];
        hsize_t sus_block[4];

        sus_start[0] = Lambda_it;
        for(int i=1; i<4;i++){sus_start[i] = 0;};
        for(int i=0; i<4;i++){
            sus_stride[i] = 1;
            sus_block[i] = 1;
        };
        sus_count[0]= 1;
        sus_count[1]= nuc_eff;
        sus_count[2]= nuc_eff;
        sus_count[3]= 3;

        dataspacevertex_sus.selectHyperslab(H5S_SELECT_SET, sus_count,sus_start,sus_stride,sus_block);
        //Select hyperslab in buffer
        hsize_t sus_start_b[3];
        hsize_t sus_stride_b[3];
        hsize_t sus_count_b[3];
        hsize_t sus_block_b[3];

        for(int i=0; i<3;i++){
            sus_start_b[i] = 0;
            sus_stride_b[i] = 1;
            sus_block_b[i] = 1;
        };

        sus_count_b[0]= nuc_eff;
        sus_count_b[1]= nuc_eff;
        sus_count_b[2]= 3;

        dataspacevertex_sus_buffer.selectHyperslab(H5S_SELECT_SET, sus_count_b,sus_start_b,sus_stride_b,sus_block_b);




        //self energy
        //Select hyperslab in the file where the data should be located
        hsize_t self_start[2];
        hsize_t self_stride[2];
        hsize_t self_count[2];
        hsize_t self_block[2];

        self_start[0] = Lambda_it;
        for(int i=1; i<2;i++){self_start[i] = 0;};
        for(int i=0; i<2;i++){
            self_stride[i] = 1;
            self_block[i] = 1;
        };
        self_count[0]= 1;
        self_count[1]= nw;

        dataspaceself.selectHyperslab(H5S_SELECT_SET, self_count,self_start,self_stride,self_block);
        //Select hyperslab in buffer
        hsize_t self_start_b[1] ;
        hsize_t self_stride_b[1];
        hsize_t self_count_b[1] ;
        hsize_t self_block_b[1];
        self_start_b[0] = 0;
        self_stride_b[0] = 1;
        self_count_b[0] = nw;
        self_block_b[0] = 1;

        dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, self_count_b,self_start_b,self_stride_b,self_block_b);



        //write buffer data into file.


        dataset_R.write( R_class, mtype1,dataspacevertex_R_buffer,dataspacevertex_R );
        dataset_K1.write( K1_class, mtype1 ,dataspacevertex_K1_buffer,dataspacevertex_K1 );
        dataset_K2.write( K2_class, mtype1 ,dataspacevertex_K2_buffer,dataspacevertex_K2);
        dataset_K2b.write( K2b_class, mtype1 ,dataspacevertex_K2b_buffer,dataspacevertex_K2b);
        dataset_irred.write( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);
        dataset_sus.write( Suscept, H5::PredType::NATIVE_DOUBLE ,dataspacevertex_sus_buffer,dataspacevertex_sus );
        dataset_self.write( selfenergy, mtype2, dataspaceself_buffer,dataspaceself);

        delete[] K1_class;
        delete[] R_class;
        delete[] K2_class;
        delete[] K2b_class;






        cout << "Successfully saved in hdf5 file: " << FILE_NAME << endl;
    }  // end of try block

    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printErrorStack();
        return ;
    }

    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
    {
        error.printErrorStack();
        return ;
    }

    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error)
    {
        error.printErrorStack();
        return ;
    }


}

//function to read out an exstiting Hdf5 file. Needed to resume computation if it has been interrupted during the flow
state read_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size){
    state result;


    typedef struct complex_t{
        double spin_re;   /*real part */
        double dens_re;   /*real part */
    }complex_t;

    typedef struct complex{//complex data type for self energy
        double re;   /*real part */
        double im;   /*imaginary part */

    }complex;

    //buffer self energy
    complex selfenergy[nw];//irrdeucible vertex
    complex_t irreducible_class[nuc_eff][nuc_eff][3];//irrdeucible vertex
    auto R_class = new complex_t[3][nuc_eff][nuc_eff][3][nw3_q][nw3_w1][nw3_w2];
    auto K1_class = new complex_t[3][nuc_eff][nuc_eff][3][nw1_q];
    auto K2_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];
    auto K2b_class = new complex_t[3][nuc_eff][nuc_eff][3][nw2_q][nw2_w1];
    double Suscept[nuc_eff][nuc_eff][3];





    H5::CompType mtype1( sizeof(complex_t) );
    mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);

    H5::CompType mtype2( sizeof(complex) );//for self energy
    mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
    mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);


    H5::H5File file(FILE_NAME, H5F_ACC_RDONLY );
    H5::DataSet dataset_irred = file.openDataSet("irred");
    H5::DataSet dataset_R = file.openDataSet("R");
    H5::DataSet dataset_K1 = file.openDataSet("K1");
    H5::DataSet dataset_K2 = file.openDataSet("K2");
    H5::DataSet dataset_K2b = file.openDataSet("K2b");
    H5::DataSet dataset_sus = file.openDataSet("sus");
    H5::DataSet dataset_self = file.openDataSet("selflist");


    H5::DataSpace dataspacevertex_R=dataset_R.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K1=dataset_K1.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K2=dataset_K2.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K2b=dataset_K2b.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_irreducible=dataset_irred.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_sus= dataset_sus.getSpace();//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspaceself=dataset_self.getSpace();

    hsize_t R_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t sus_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K2b_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nw)};





    // Create the dimension arrays for objects in buffer.
    hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t R_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw3_q),static_cast<hsize_t>(nw3_w1),static_cast<hsize_t>(nw3_w2)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t sus_dims_buffer[]= {static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw1_q)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(3),static_cast<hsize_t>(nuc_eff),static_cast<hsize_t>(nuc_eff),3,static_cast<hsize_t>(nw2_q),static_cast<hsize_t>(nw2_w1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};


    // Create the data space for buffer objects.
    H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_R_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K1_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K2_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_K2b_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspacevertex_sus_buffer(RANK_sus-1, sus_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    H5::DataSpace dataspaceself_buffer(1, self_dims_buffer);



    //Select hyperslabs:

    //R_class

    //Select hyperslab in the file
    hsize_t R_start[8];
    hsize_t R_stride[8];
    hsize_t R_count[8];
    hsize_t R_block[8];

    R_start[0] = Lambda_it;
    for(int i=1; i<8;i++){R_start[i] = 0;};
    for(int i=0; i<8;i++){
        R_stride[i] = 1;
        R_block[i] = 1;
    };
    R_count[0]= 1;
    R_count[1]= 3;//three different channels
    R_count[2]= nuc_eff;
    R_count[3]= nuc_eff;
    R_count[4]= 3;
    R_count[5]= nw3_q;
    R_count[6]= nw3_w1;
    R_count[7]= nw3_w2;
    dataspacevertex_R.selectHyperslab(H5S_SELECT_SET, R_count,R_start,R_stride,R_block);
    //Select hyperslab in buffer
    hsize_t R_start_b[7];
    hsize_t R_stride_b[7];
    hsize_t R_count_b[7];
    hsize_t R_block_b[7];

    for(int i=0; i<7;i++){R_start_b[i] = 0;};
    for(int i=0; i<7;i++){
        R_stride_b[i] = 1;
        R_block_b[i] = 1;
    };
    R_count_b[0]= 3;
    R_count_b[1]= nuc_eff;
    R_count_b[2]= nuc_eff;
    R_count_b[3]= 3;
    R_count_b[4]= nw3_q;
    R_count_b[5]= nw3_w1;
    R_count_b[6]= nw3_w2;
    dataspacevertex_R_buffer.selectHyperslab(H5S_SELECT_SET, R_count_b,R_start_b,R_stride_b,R_block_b);




    //K1_class


    //Select hyperslab in the file
    hsize_t K1_start[6];
    hsize_t K1_stride[6];
    hsize_t K1_count[6];
    hsize_t K1_block[6];

    K1_start[0] = Lambda_it;
    for(int i=1; i<6;i++){K1_start[i] = 0;};
    for(int i=0; i<6;i++){
        K1_stride[i] = 1;
        K1_block[i] = 1;
    };
    K1_count[0]= 1;
    K1_count[1]= 3;
    K1_count[2]= nuc_eff;
    K1_count[3]= nuc_eff;
    K1_count[4]= 3;
    K1_count[5]= nw1_q;
    dataspacevertex_K1.selectHyperslab(H5S_SELECT_SET, K1_count,K1_start,K1_stride,K1_block);
    //Select hyperslab in buffer
    hsize_t K1_start_b[5];
    hsize_t K1_stride_b[5];
    hsize_t K1_count_b[5];
    hsize_t K1_block_b[5];

    for(int i=0; i<5;i++){K1_start_b[i] = 0;};
    for(int i=0; i<5;i++){
        K1_stride_b[i] = 1;
        K1_block_b[i] = 1;
    };
    K1_count_b[0]= 3;
    K1_count_b[1]= nuc_eff;
    K1_count_b[2]= nuc_eff;
    K1_count_b[3]= 3;
    K1_count_b[4]= nw1_q;
    dataspacevertex_K1_buffer.selectHyperslab(H5S_SELECT_SET, K1_count_b,K1_start_b,K1_stride_b,K1_block_b);





    //K2_class

    //Select hyperslab in the file
    hsize_t K2_start[7];
    hsize_t K2_stride[7];
    hsize_t K2_count[7];
    hsize_t K2_block[7];

    K2_start[0] = Lambda_it;
    for(int i=1; i<7;i++){K2_start[i] = 0;};
    for(int i=0; i<7;i++){
        K2_stride[i] = 1;
        K2_block[i] = 1;
    };
    K2_count[0]= 1;
    K2_count[1]= 3;
    K2_count[2]= nuc_eff;
    K2_count[3]= nuc_eff;
    K2_count[4]= 3;
    K2_count[5]= nw2_q;
    K2_count[6]= nw2_w1;
    dataspacevertex_K2.selectHyperslab(H5S_SELECT_SET, K2_count,K2_start,K2_stride,K2_block);
    //Select hyperslab in buffer
    hsize_t K2_start_b[6];
    hsize_t K2_stride_b[6];
    hsize_t K2_count_b[6];
    hsize_t K2_block_b[6];

    for(int i=0; i<6;i++){K2_start_b[i] = 0;};
    for(int i=0; i<6;i++){
        K2_stride_b[i] = 1;
        K2_block_b[i] = 1;
    };
    K2_count_b[0]= 3;
    K2_count_b[1]= nuc_eff;
    K2_count_b[2]= nuc_eff;
    K2_count_b[3]= 3;
    K2_count_b[4]= nw2_q;
    K2_count_b[5]= nw2_w1;
    dataspacevertex_K2_buffer.selectHyperslab(H5S_SELECT_SET, K2_count_b,K2_start_b,K2_stride_b,K2_block_b);





    //K2b_class

    //Select hyperslab in the file
    hsize_t K2b_start[7];
    hsize_t K2b_stride[7];
    hsize_t K2b_count[7];
    hsize_t K2b_block[7];

    K2b_start[0] = Lambda_it;
    for(int i=1; i<7;i++){K2b_start[i] = 0;};
    for(int i=0; i<7;i++){
        K2b_stride[i] = 1;
        K2b_block[i] = 1;
    };
    K2b_count[0]= 1;
    K2b_count[1]= 3;
    K2b_count[2]= nuc_eff;
    K2b_count[3]= nuc_eff;
    K2b_count[4]= 3;
    K2b_count[5]= nw2_q;
    K2b_count[6]= nw2_w1;
    dataspacevertex_K2b.selectHyperslab(H5S_SELECT_SET, K2b_count,K2b_start,K2b_stride,K2b_block);
    //Select hyperslab in buffer
    hsize_t K2b_start_b[6];
    hsize_t K2b_stride_b[6];
    hsize_t K2b_count_b[6];
    hsize_t K2b_block_b[6];

    for(int i=0; i<6;i++){K2b_start_b[i] = 0;};
    for(int i=0; i<6;i++){
        K2b_stride_b[i] = 1;
        K2b_block_b[i] = 1;
    };
    K2b_count_b[0]= 3;
    K2b_count_b[1]= nuc_eff;
    K2b_count_b[2]= nuc_eff;
    K2b_count_b[3]= 3;
    K2b_count_b[4]= nw2_q;
    K2b_count_b[5]= nw2_w1;
    dataspacevertex_K2b_buffer.selectHyperslab(H5S_SELECT_SET, K2b_count_b,K2b_start_b,K2b_stride_b,K2b_block_b);





    //irred_class


    //Select hyperslab in the file
    hsize_t irred_start[4];
    hsize_t irred_stride[4];
    hsize_t irred_count[4];
    hsize_t irred_block[4];

    irred_start[0] = Lambda_it;
    for(int i=1; i<4;i++){irred_start[i] = 0;};
    for(int i=0; i<4;i++){
        irred_stride[i] = 1;
        irred_block[i] = 1;
    };
    irred_count[0]= 1;
    irred_count[1]= nuc_eff;
    irred_count[2]= nuc_eff;
    irred_count[3]= 3;
    dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, irred_count,irred_start,irred_stride,irred_block);
    //Select hyperslab in buffer
    hsize_t irred_start_b[3];
    hsize_t irred_stride_b[3];
    hsize_t irred_count_b[3];
    hsize_t irred_block_b[3];
    for(int i=0; i<3;i++){irred_start_b[i] = 0;};
    for(int i=0; i<3;i++){
        irred_stride_b[i] = 1;
        irred_block_b[i] = 1;
    };
    irred_count_b[0]= nuc_eff;
    irred_count_b[1]= nuc_eff;
    irred_count_b[2]= 3;
    dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, irred_count_b,irred_start_b,irred_stride_b,irred_block_b);


    //susceptibility

    //Select hyperslab in the file
    hsize_t sus_start[4];
    hsize_t sus_stride[4];
    hsize_t sus_count[4];
    hsize_t sus_block[4];

    sus_start[0] = Lambda_it;
    for(int i=1; i<4;i++){sus_start[i] = 0;};
    for(int i=0; i<4;i++){
        sus_stride[i] = 1;
        sus_block[i] = 1;
    };
    sus_count[0]= 1;
    sus_count[1]= nuc_eff;
    sus_count[2]= nuc_eff;
    sus_count[3]= 3;

    dataspacevertex_sus.selectHyperslab(H5S_SELECT_SET, sus_count,sus_start,sus_stride,sus_block);
    //Select hyperslab in buffer
    hsize_t sus_start_b[3];
    hsize_t sus_stride_b[3];
    hsize_t sus_count_b[3];
    hsize_t sus_block_b[3];

    for(int i=0; i<3;i++){
        sus_start_b[i] = 0;
        sus_stride_b[i] = 1;
        sus_block_b[i] = 1;
    };

    sus_count_b[0]= nuc_eff;
    sus_count_b[1]= nuc_eff;
    sus_count_b[2]= 3;

    dataspacevertex_sus_buffer.selectHyperslab(H5S_SELECT_SET, sus_count_b,sus_start_b,sus_stride_b,sus_block_b);




    //self energy
    //Select hyperslab in the file
    hsize_t self_start[2];
    hsize_t self_stride[2];
    hsize_t self_count[2];
    hsize_t self_block[2];

    self_start[0] = Lambda_it;
    for(int i=1; i<2;i++){self_start[i] = 0;};
    for(int i=0; i<2;i++){
        self_stride[i] = 1;
        self_block[i] = 1;
    };
    self_count[0]= 1;
    self_count[1]= nw;

    dataspaceself.selectHyperslab(H5S_SELECT_SET, self_count,self_start,self_stride,self_block);
    //Select hyperslab in buffer
    hsize_t self_start_b[1] ;
    hsize_t self_stride_b[1];
    hsize_t self_count_b[1] ;
    hsize_t self_block_b[1];
    self_start_b[0] = 0;
    self_stride_b[0] = 1;
    self_count_b[0] = nw;
    self_block_b[0] = 1;

    dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, self_count_b,self_start_b,self_stride_b,self_block_b);



    dataset_R.read( R_class, mtype1,dataspacevertex_R_buffer,dataspacevertex_R );
    dataset_K1.read( K1_class, mtype1 ,dataspacevertex_K1_buffer,dataspacevertex_K1 );
    dataset_K2.read( K2_class, mtype1 ,dataspacevertex_K2_buffer,dataspacevertex_K2);
    dataset_K2b.read( K2b_class, mtype1 ,dataspacevertex_K2b_buffer,dataspacevertex_K2b);
    dataset_irred.read( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);
    dataset_sus.read( Suscept, H5::PredType::NATIVE_DOUBLE ,dataspacevertex_sus_buffer,dataspacevertex_sus );
    dataset_self.read( selfenergy, mtype2, dataspaceself_buffer,dataspaceself);



    for(int i=0; i<nw; i++){
        result.selfenergy.setself(i,selfenergy[i].re + ci*selfenergy[i].im);
    };





    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half


                    result.vertex.spinvertex.irred.setvert(a,b,c,irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].spin_re);
                    result.vertex.densvertex.irred.setvert(a,b,c,irreducible_class[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1].dens_re);

                };};};};


    //Note that the first index labels the channel with the correspondance ( 0 = s-channel, 1 = t-channel, 2 = u-channel)


    //buffer all R-class-arrays:


#pragma omp parallel for
    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){

                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half


                    for(int i=(nw3-nw3_q); i<nw3; i++){
                        for(int j=(nw3-nw3_w1) ; j<nw3; j++){
                            for(int k=(nw3-nw3_w2) ; k<nw3; k++){


                                //rest class from s-channel
                                result.vertex.spinvertex.svertex.R_setvert(a,b,c,i,j,k, R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re );
                                                                result.vertex.densvertex.svertex.R_setvert(a,b,c,i,j,k, R_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re );

                                //rest class from t-channel
                                result.vertex.spinvertex.tvertex.R_setvert(a,b,c,i,j,k, R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re );
                               result.vertex.densvertex.tvertex.R_setvert(a,b,c,i,j,k, R_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re );


                                //rest class from u-channel
                                result.vertex.spinvertex.uvertex.R_setvert(a,b,c,i,j,k, R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].spin_re );
                                                                result.vertex.densvertex.uvertex.R_setvert(a,b,c,i,j,k, R_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)].dens_re );



                            };};};};};};};

    //buffer all K1-class arrays:






#pragma omp parallel for
    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

                    for(int i=(nw1-nw1_q) ; i<nw1; i++){





                        //  K1 class from s-channel
                        result.vertex.spinvertex.svertex.K1_setvert(a,b,c,i,K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re );
                        result.vertex.densvertex.svertex.K1_setvert(a,b,c,i,K1_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re );
                        //K1 class from t-channel
                        result.vertex.spinvertex.tvertex.K1_setvert(a,b,c,i,K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re );
                        result.vertex.densvertex.tvertex.K1_setvert(a,b,c,i,K1_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re );

                        //K1 class from u-channel
                        result.vertex.spinvertex.uvertex.K1_setvert(a,b,c,i,K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].spin_re );
                       result.vertex.densvertex.uvertex.K1_setvert(a,b,c,i,K1_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw1-nw1_q)].dens_re );
                    };};};};};

    //buffer all K2-class-arrays:



#pragma omp parallel for
    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

                    for(int i=(nw2-nw2_q) ; i<nw2; i++){
                                        for(int j=(nw2-nw2_w1) ; j<nw2; j++){




                            //rest class from s-channel
                            result.vertex.spinvertex.svertex.K2_setvert(a,b,c,i,j,K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re );
                            result.vertex.densvertex.svertex.K2_setvert(a,b,c,i,j,K2_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re );

                            //rest class from t-channel
                            result.vertex.spinvertex.tvertex.K2_setvert(a,b,c,i,j,K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re );
                            result.vertex.densvertex.tvertex.K2_setvert(a,b,c,i,j,K2_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re );

                            //rest class from u-channel
                            result.vertex.spinvertex.uvertex.K2_setvert(a,b,c,i,j,K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re );
                            result.vertex.densvertex.uvertex.K2_setvert(a,b,c,i,j,K2_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re );

                        };};};};};};



    //buffer all K2b-class-arrays:




#pragma omp parallel for
    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half



                    for(int i=(nw2-nw2_q) ; i<nw2; i++){
                        double q = bfreqs[(nw-nw2)/2+i];
                        for(int j=(nw2-nw2_w1) ; j<nw2; j++){
                            double w2 = ffreqs[(nw-nw2)/2+j];

                            //rest class from s-channel
                            result.vertex.spinvertex.svertex.K2b_setvert(a,b,c,i,j,K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re);
                            result.vertex.densvertex.svertex.K2b_setvert(a,b,c,i,j,K2b_class[0][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re);


                            //rest class from t-channel
                            result.vertex.spinvertex.tvertex.K2b_setvert(a,b,c,i,j,K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re);
                            result.vertex.densvertex.tvertex.K2b_setvert(a,b,c,i,j,K2b_class[1][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re);

                            //rest class from u-channel
                            result.vertex.spinvertex.uvertex.K2b_setvert(a,b,c,i,j,K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].spin_re);
                            result.vertex.densvertex.uvertex.K2b_setvert(a,b,c,i,j,K2b_class[2][a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)].dens_re);


                        };};};};};};







    //susceptibility



    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
            for (int c=1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) ){
                    if(distance(a,b,c) <= d_c){


                        //  Susceptibility in a-channel with bosonic transfer freq q
                        result.sus.write(a,b,c,Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] );

                    }
                    else {
                        result.sus.write(a,b,c,0);
                    };};};};};


    delete[] K1_class;
    delete[] R_class;
    delete[] K2_class;
    delete[] K2b_class;







    return result;



}



/****************************************** The following functions convert between physical frequencies and grid frequencies. Make sure to use right routine according to whether working on a lin or a log grid************************/


/****************************LINEAR GRID:***************/
//to convert on full frequency grid:



//int fconv(double w){//conversion  on linear grid
//    int i; //next lower lying lattice site on frequency vector
//    if(abs(w) < ffreqs[nw-1]){
//        i = w/k - 0.5 + nw/2;
//        if(i != nw-1 && w == ffreqs[i+1]){i+=1;}//avoid rounding errors:
//        else if(i != 0 && w == ffreqs[i-1]){i-=1;}
//        return i;}
//    else if(abs((w - ffreqs[nw-1]))<1e-12){return nw-1;}
//    else if(abs((w - ffreqs[0]))<1e-12){return 0;}
//    else{cout << "erroneous frequency conversion fconv:  max frequency is " << ffreqs[nw-1] << ". Given frequency: " <<  w << endl; return 0;};
//}



////to convert on reduced frequency grid:



//int fconv_n(double w, int n){//conversion  on linear grid of size n
//    int i; //next lower lying lattice site on frequency vector
//    if(fabs(w) < ffreqs[(nw+n)/2-1]){
//        i = w/k - 0.5 + n/2;
//        if(i != n-1 && w == ffreqs[(nw-n)/2+i+1]){i+=1;}//avoid rounding errors:
//        else if(i != 0 && w == ffreqs[(nw-n)/2+i-1]){i-=1;}
//        return i;}
//    else if(abs((w - ffreqs[(nw+n)/2-1]))<1e-12){return n-1;}
//    else if(abs((w - ffreqs[(nw-n)/2]))<1e-12){return 0;}
//    else{cout << setprecision (16) << endl;
//        cout << "erroneous frequency conversion fconv_n: n is " <<n << "  with the max frequency " << ffreqs[(nw+n)/2-1] << ". Given frequency: " <<  w << endl; return 0;};
//}


/****************************LOG GRID:***************/
//to convert on full frequency grid:

int fconv(double w){//conversion  on lcombination of linear and log grid. This function can only be called if it is ensured that w is in the range of the frequency grid and that abs(w) >= w0 where w0 is the smallest frequency saved.
    int i;



        if(abs(w) > wt){
            i = static_cast<int>((abs(w)-wt)/delw + nlog/2)-1 ;

        }
        else if(abs(w) <= wt){

            i = static_cast<int>(log(abs(w)/w0)/log(k));
        };

        if(w>0){
            i += nw/2;


        }
        else if(w<0){
            i = nw/2-1 - i;
        };


        if(i != nw-1 && abs(w - ffreqs[i+1])<1e-6){i+=1;}//avoid rounding errors:
        else if(i != 0 && abs( w- ffreqs[i-1])<1e-6){i-=1;};


    return i;

}



//to convert on reduced frequency grid:

int fconv_n(double w, int n){//conversion  on combination of linear and log grid. This function can only be called if it is ensured that w is in the range of the frequency grid and that abs(w) >= w0 where w0 is the smallest frequency saved.

    int i;


        if(abs(w) > wt){
            i = static_cast<int>((abs(w)-wt)/delw + nlog/2)-1 ;

        }
        else if(abs(w) <= wt){

            i = static_cast<int>(log(abs(w)/w0)/log(k));
        };

        if(w>0){
            i += n/2;
        }
        else if(w<0){
            i = n/2-1 - i;
        };


        if(i != n-1 && abs(w - ffreqs[(nw-n)/2+i+1]) <1e-6){i+=1;}//avoid rounding errors:
        else if(i != 0 && abs(w -ffreqs[(nw-n)/2+i-1])<1e-6){i-=1;};




    return i;


}





/***********************************************************************************************************************************************************************************/



/****************************FUNCTIONS CONCERNING KAGOME LATTICE***************/
site site_project(int a, int b, int c){//projects lattice site to its equivalent couterpart in upper half
    site result;
    if(b <0 || (b==0 && a < 0 && c != 2)){//only project if site lies in lower half

        if (c==1){
            result.a = -a;
            result.b = -b;
            result.c = 1;
        }
        else if(c==2){

            result.a = -a ;
            result.b = - b -1;
            result.c = 2;
        }
        else if(c==3){

            result.a = -a -1;
            result.b =  -b ;
            result.c = 3;
        };
    }
    else{
        result.a = a;
        result.b = b;
        result.c = c;
    };
    return result;
}

site site_switch(int a_raw, int b_raw, int c_raw){//switches two sites and projects onto upper plane

    int a,b,c;
    if (c_raw==1){
        a = - a_raw;
        b = - b_raw;
        c = 1;
    }
    else if(c_raw==2){

        a = - (- a_raw - b_raw );
        b = - ( 1 * a_raw);
        c = 3;
    }
    else if(c_raw==3){

        a = - ( b_raw );
        b = - ( -1 * a_raw - b_raw);
        c = 2;
    };

    return site_project(a,b,c);//project onto upper halfplane
}

double distance(int a, int b, int c){//determines the absolute distance between two points on the lattice
    double delta_1,delta_2;//shift due to different atoms in unit cell, given by index c
    if(c==1){
        delta_1 = 0;
        delta_2 = 0;}
    else if(c==2){
        delta_1 = 0.25 * 1;
        delta_2 = 0.25 * sqrt(3);}
    else if(c==3){
        delta_1 = 0.5;
        delta_2 = 0;}
    return (sqrt(pow( a+ 0.5*b + delta_1,2) + pow( 0.5 * sqrt(3) * b + delta_2,2)));

}



/*****************FUNCTIONS FOR THE S-VERTEX********************************************/

//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
double svert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (w1-w2+q)/2;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (-w1+w2-q)/2;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2;
            w2_s = w2-(q+w1)/2;};
        double value=0;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};


//        if(s > bfreqs[nw-1] && abs(s - bfreqs[nw-1]) < 1e-12){s = bfreqs[nw-1];}
//        else if(s<bfreqs[0] && abs(s - bfreqs[0]) < 1e-12){s = bfreqs[0];};


//        if(w1_s > ffreqs[nw-1] && abs(w1_s - ffreqs[nw-1]) < 1e-12){w1_s = ffreqs[nw-1];}
//        else if(w1_s<ffreqs[0] && abs(w1_s - ffreqs[0]) < 1e-12){w1_s = ffreqs[0];};


//        if(w2_s > ffreqs[nw-1] && abs(w2_s - ffreqs[nw-1]) < 1e-12){w2_s = ffreqs[nw-1];}
//        else if(w2_s<ffreqs[0] && abs(w2_s - ffreqs[0]) < 1e-12){w2_s = ffreqs[0];};

        value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations


        return value;}
    else{return 0;}

}
//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
double svert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;

        if(channel == 's'){
            s = q;
            w1_s = w1;
            w2_s = w2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (w1-w2+q)/2;}
        else if(channel == 'u'){

            s = w1+w2;
            w1_s = (w1-w2-q)/2;
            w2_s = (-w1+w2-q)/2;

        }
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            s = q+w1;
            w1_s = (q-w1)/2;
            w2_s = w2-(q+w1)/2;};
        double value=0;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};


//        if(s > bfreqs[nw-1] && abs(s - bfreqs[nw-1]) < 1e-12){s = bfreqs[nw-1];}
//        else if(s<bfreqs[0] && abs(s - bfreqs[0]) < 1e-12){s = bfreqs[0];};


//        if(w1_s > ffreqs[nw-1] && abs(w1_s - ffreqs[nw-1]) < 1e-12){w1_s = ffreqs[nw-1];}
//        else if(w1_s<ffreqs[0] && abs(w1_s - ffreqs[0]) < 1e-12){w1_s = ffreqs[0];};


//        if(w2_s > ffreqs[nw-1] && abs(w2_s - ffreqs[nw-1]) < 1e-12){w2_s = ffreqs[nw-1];}
//        else if(w2_s<ffreqs[0] && abs(w2_s - ffreqs[0]) < 1e-12){w2_s = ffreqs[0];};

        if(p==1){
            if(channel=='s'){
                if(f == 'R' || f == 'M'){value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2b_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'M'){value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        }

        else if(p==2){
            if(channel=='s'){
                if(f == 'R' || f == 'L'){value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K2_vvalsmooth(a,b,c,s,w2_s);}//if outer legs are conntected to differentbare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,s) + K2b_vvalsmooth(a,b,c,s,w1_s);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='t' || channel=='u'){
                if(f == 'R' || f== 'L'){value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s);};
            };
        };



        return value;}
    else{return 0;}

}
//overload with first extra firs two arguments: red_side = vertex number with only complementary channels (0,1,2), map = determines if channel mapping is turned on for this vertex (0=off/1=on)
double svert::vvalsmooth(int red_side,int map, int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){

    //THIS FUNCTION IS NEEDED ONLY WHEN BUBBLE FUNCTIONS ARE USED WITH VERTEX OF TYPE "PARVERT" INSTEAD OF "FULLVERT". LEADS TO PREVIOUS FUNCTION.
    return vvalsmooth(a, b, c, q, w1,w2, channel, p, f);
}
double svert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double s,w1_s,w2_s;
        s = q;
        w1_s = w1;
        w2_s = w2;

        if(abs(s) < bfreqs[nw/2]){ if (s >= 0) {s = bfreqs[nw/2];} else{s = bfreqs[nw/2-1];};};
        if(abs(w1_s) < ffreqs[nw/2]){if (w1_s >= 0) {w1_s = ffreqs[nw/2];} else{w1_s =  ffreqs[nw/2-1];};};
        if(abs(w2_s) < ffreqs[nw/2]){if (w2_s > 0) {w2_s =  ffreqs[nw/2];} else{w2_s =  ffreqs[nw/2-1];};};



        double value=0;
        value += R_vvalsmooth(a,b,c,s,w1_s,w2_s) + K1_vvalsmooth(a,b,c,s) + K2_vvalsmooth(a,b,c,s,w1_s) + K2b_vvalsmooth(a,b,c,s,w2_s)  ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}

}
void svert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){

            R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};
    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void svert::K1_setvert(int a, int b, int c,int i, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>=0 ){

            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};

    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
void svert::K2_setvert(int a, int b, int c, int i, int j, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value ;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void svert::K2b_setvert(int a, int b, int c, int i, int j,  double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){


            K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
double svert::R_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){

    
    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    double value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0){

            value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

        else if(sym == 1 ){

            if(i>=nw3/2){
                if(j>=nw3/2){
                    value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}
                else if(j<nw3/2){
                    site x = site_switch(a,b,c);
                    j=nw3-1-j;
                    k = nw3-1-k;
                    value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];
                };

            }
            else if(i<nw3/2){
                if(k>=nw3/2){
                    site x = site_switch(a,b,c);
                    i=nw3-1-i;

                    value += conj(R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);//note that j and k are interchanged
                }
                else if(k<nw3/2){
                    i=nw3-1-i;
                    j=nw3-1-j;
                    k = nw3-1-k;
                    value +=conj(R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);//note that j and k are interchanged
                }
            };
        }
        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){x = site_switch(x.a,x.b,x.c);i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};
            value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};


    };
    if(abs(value)<1e-100){value=0;};


    return value;
}
double svert::K1_vval(int a_raw, int b_raw, int c_raw,  int i){



    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i<nw1 && i>= 0 ){
        if (sym==0){
            
            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 || sym==2){
            if (i>= nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);
                i = nw1-1-i;
                value += conj(K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)]);};
        }

        ;};
    if(abs(value)<1e-100){value=0;};



    return value;
}
double svert::K2_vval(int a_raw, int b_raw, int c_raw,  int i, int j ){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j >= nw2/2){

                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
            else if (j< nw2/2){


                site x = site_switch(a,b,c);
                j = nw2 -1 - j;

                value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];

            };

        }
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i; x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j; x = site_switch(x.a,x.b,x.c);};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };
    };

    if(abs(value)<1e-100){value=0;};


    return value;
}
double svert::K2b_vval(int a_raw, int b_raw, int c_raw, int i, int j){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0){
        if(sym==0){
            
            value += K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            site x = site_switch(a,b,c);
            value += conj(K2_vval(x.a,x.b,x.c, nw2-1-i, j ));}
        else if(sym==2){
            site x = site_switch(a,b,c);
            value += K2_vval(x.a,x.b,x.c, i, nw2-1-j );
        };};
    if(abs(value)<1e-100){value=0;};

    return value;
}
double svert::R_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(s)+1e-6 < ffreqs[(nw3+nw)/2-1] && abs(w1)+1e-6<ffreqs[(nw3+nw)/2-1] && abs(w2)+1e-6 <ffreqs[(nw3+nw)/2-1]){//if frequency arguments are out of range, vertex vanishes
        if (sym==0){

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);
            int S_vert = fconv_n(s,nw3);
            int W1_vert = fconv_n(w1,nw3);
            int W2_vert = fconv_n(w2,nw3);


            value += ((R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                    ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){

            if(s > 0){
                if((w1)>0 ){
                    int S = fconv_n(s,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);

                    int S_vert = fconv_n(s,nw3)-(nw3-nw3_w1);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

                    value += ((R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                            ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }
                else if((w1)<0 ){
                    site x = site_switch(a,b,c);

                    int S = fconv_n(s,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);
                    int S_vert = fconv_n(s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;

                    value += ((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}

            else if (s < 0){
                if((w2)>0){

                    site x = site_switch(a,b,c);

                    int S = fconv_n(-s,nw3);
                    int W1 = fconv_n(w2,nw3);
                    int W2 = fconv_n(w1,nw3);
                    int S_vert = fconv_n(-s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;
                    s = -s;


                    value += conj((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }

                else if ((w2)<0){
                    int S = fconv_n(-s,nw3);
                    int W1 = fconv_n(-w2,nw3);
                    int W2 = fconv_n(-w1,nw3);

                    int S_vert = fconv_n(-s,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w2;
                    double w2_eff = -w1;
                    s = -s;


                    value += conj((R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

                };};

        }
        else if(sym==2){
            site x(a,b,c);

            double w1_eff = w1, w2_eff = w2, s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if((w1_eff) < 0){w1_eff = -w1_eff; w2_eff = -w2_eff; x = site_switch(x.a,x.b,x.c);};
            s = s_eff;

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int S_vert = fconv_n(s,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);

            value = ((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))
                    / ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));



        };}
    else{
        int i,j,k;
        if(abs(s)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(s,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value += R_vval(a,b,c,i,j,k);};};





    return value;
}
double svert::K1_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(s) + 1e-6 < ffreqs[(nw1+nw)/2-1]){
        if(sym==0){
            int S = fconv_n(s,nw1);
            int S_vert = fconv_n(s,nw1)-(nw1-nw1_q);
            value = ((K1[a+(nuc_eff-1)/2][b][c-1][S_vert]*(bfreqs[(nw-nw1)/2+S+1]-s)
                    +K1[a+(nuc_eff-1)/2][b][c-1][S_vert+1]*(-bfreqs[(nw-nw1)/2+S]+s))/(bfreqs[(nw-nw1)/2+S+1]-bfreqs[(nw-nw1)/2+S]));}
        else if(sym==1 || sym==2){


            if (s > 0){

                int S = fconv_n(s,nw1);
                int S_vert = fconv_n(s,nw1)-(nw1-nw1_q);
                value += ((K1[a+(nuc_eff-1)/2][b][c-1][S_vert]*(bfreqs[S+1]-s)
                        +K1[a+(nuc_eff-1)/2][b][c-1][S_vert+1]*(-bfreqs[S]+s))/(bfreqs[S+1]-bfreqs[S]));}
            else if (s<0){
                site x = site_switch(a,b,c);

                double s_eff = -s;
                int S = fconv_n(s_eff,nw1);
                int S_vert = fconv_n(s_eff,nw1)-(nw1-nw1_q);

                if(sym==1){

                    value = conj((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert]*(bfreqs[S+1]-s_eff)
                            +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));}
                else{
                    value = ((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert]*(bfreqs[S+1]-s_eff)
                            +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));};
            };};

    }
    else{
        int i;
        if(abs(s)<= bfreqs[nw1-1]){
            i= fconv_n(s,nw1);

            value += K1_vval(a,b,c,i);}

        ;};




    return value;

}
double svert::K2_vvalsmooth(int a_raw, int b_raw, int c_raw,  double s, double w1){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    double value=0;

    if(abs(s) + 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){
            int S = fconv_n(s,nw2);
            int W1 = fconv_n(w1,nw2);

            int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){
            if(w1 > 0){


                int S = fconv_n(s,nw2);
                int W1 = fconv_n(w1,nw2);


                int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
            else if(w1 < 0){

                site x = site_switch(a,b,c);

                int S = fconv_n(s,nw2);
                int W1 = fconv_n(-w1,nw2);

                int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(-w1,nw2)-(nw2-nw2_w1);
                double w1_eff = -w1 ;
                value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));};
        }
        else if(sym==2){

            site x(a,b,c);
            double w1_eff = w1;
            double s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};


            int S = fconv_n(s_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int S_vert = fconv_n(s_eff,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);

            value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][S_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

        };}
    else{
        if(abs(s)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            int i,j;

            i= fconv_n(s,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};
    };




    return value;
}
double svert::K2b_vvalsmooth(int a_raw, int b_raw, int c_raw,  double s, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    double value=0;
    if(abs(s)+ 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w2)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int S = fconv_n(s,nw2);
            int W2 = fconv_n(w2,nw2);

            int S_vert = fconv_n(s,nw2)-(nw2-nw2_q);
            int W2_vert = fconv_n(w2,nw2)-(nw2-nw2_w1);

            value += ((K2b[a+(nuc_eff-1)/2][b][c-1][S_vert][W2_vert]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W2_vert]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][S_vert][W2_vert+1]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][S_vert+1][W2_vert+1]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2])));}
        else if(sym==1){
            site x = site_switch(a,b,c);
            value += conj(K2_vvalsmooth(x.a,x.b,x.c, -s, w2 ));}
        else if(sym==2){
            site x = site_switch(a,b,c);
            value += K2_vvalsmooth(x.a,x.b,x.c, -s, w2 );
        };
         }
    else{
        if(abs(s)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            int i,j;
            i= fconv_n(s,nw2);
            j = fconv_n(w2,nw2);

            value += K2b_vval(a,b,c,i,j);};};


    return value;
}
//non-member functions
svert operator*(double alpha, const svert& vertex){
    svert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j -(nw2-nw2_w1); j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q) ; i++){
                    for(int j=sum_R-(nw3-nw3_w1)  ; j<nw3-(nw3-nw3_w1) ; j++){
                        for(int k=0 -(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2) ; k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1) ) <= abs(k-nw3/2+(nw3-nw3_w2) ))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;
}
svert operator*(const svert& vertex,double alpha){
    svert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;

}
svert operator+(const svert& vertex1,const svert& vertex2){
    svert vertex3;

    int sum_K1 = (sym==0?0:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);

#pragma omp for
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){

                //add K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    double value = vertex1.K1[a][b][c][i] +  vertex2.K1[a][b][c][i];
                    if(abs(value)>1e-16){
                        vertex3.K1[a][b][c][i] =  value;};};
                //add K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        double value = vertex1.K2[a][b][c][i][j] +  vertex2.K2[a][b][c][i][j];
                        if(abs(value)>1e-16){
                            vertex3.K2[a][b][c][i][j] = value;};
                        if(sym==0){
                            vertex3.K2b[a][b][c][i][j] = vertex1.K2b[a][b][c][i][j] + vertex2.K2b[a][b][c][i][j];};};};
                //add R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){

                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                double value = vertex1.R[a][b][c][i][j][k] + vertex2.R[a][b][c][i][j][k];
                                if(abs(value)>1e-16){
                                    vertex3.R[a][b][c][i][j][k] = value;};
                            };};};};};};};
    return vertex3;
}


/*****************FUNCTIONS FOR THE T-VERTEX********************************************/

double tvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2;
            w2_t = (-w1-w2+q)/2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
        else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2;
            w2_t =(w1+w2+q)/2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2;
            w2_t = w1+(q-w2)/2;

        };


        double value=0;

        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};



//        if(t > bfreqs[nw-1] && abs(t - bfreqs[nw-1]) < 1e-12){t = bfreqs[nw-1];}
//        else if(t<bfreqs[0] && abs(t - bfreqs[0]) < 1e-12){t = bfreqs[0];};


//        if(w1_t > ffreqs[nw-1] && abs(w1_t - ffreqs[nw-1]) < 1e-12){w1_t = ffreqs[nw-1];}
//        else if(w1_t<ffreqs[0] && abs(w1_t - ffreqs[0]) < 1e-12){w1_t = ffreqs[0];};


//        if(w2_t > ffreqs[nw-1] && abs(w2_t - ffreqs[nw-1]) < 1e-12){w2_t = ffreqs[nw-1];}
//        else if(w2_t<ffreqs[0] && abs(w2_t - ffreqs[0]) < 1e-12){w2_t = ffreqs[0];};


        value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t)+K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}


}
double tvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        if(channel == 's'){

            t = w2-w1;
            w1_t = (w1+w2+q)/2;
            w2_t = (-w1-w2+q)/2;
        }
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15

            t = q;
            w1_t = w1;
            w2_t = w2;}
        else if(channel == 'u'){

            t = w2-w1;
            w1_t = (w1+w2-q)/2;
            w2_t =(w1+w2+q)/2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)

            t = w2-q;
            w1_t = (q+w2)/2;
            w2_t = w1+(q-w2)/2;

        };


        double value=0;

        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};



//        if(t > bfreqs[nw-1] && abs(t - bfreqs[nw-1]) < 1e-12){t = bfreqs[nw-1];}
//        else if(t<bfreqs[0] && abs(t - bfreqs[0]) < 1e-12){t = bfreqs[0];};


//        if(w1_t > ffreqs[nw-1] && abs(w1_t - ffreqs[nw-1]) < 1e-12){w1_t = ffreqs[nw-1];}
//        else if(w1_t<ffreqs[0] && abs(w1_t - ffreqs[0]) < 1e-12){w1_t = ffreqs[0];};


//        if(w2_t > ffreqs[nw-1] && abs(w2_t - ffreqs[nw-1]) < 1e-12){w2_t = ffreqs[nw-1];}
//        else if(w2_t<ffreqs[0] && abs(w2_t - ffreqs[0]) < 1e-12){w2_t = ffreqs[0];};

        if(p==1){
            if(channel=='t'){
                if(f == 'R' || f == 'M'){value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2b_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='u'){
                if(f == 'R' || f== 'M'){value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);};
            };
        }

        else if(p==2){
            if(channel=='t'){
                if(f == 'R' || f == 'L'){value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K2_vvalsmooth(a,b,c,t,w2_t);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,t) + K2b_vvalsmooth(a,b,c,t,w1_t);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='u'){
                if(f == 'R' || f== 'L'){value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t);

                };
            };
        };



        return value;}
    else{return 0;}


}
//overload of previous function
double tvert::vvalsmooth(int red_side, int map,int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
    return vvalsmooth( a, b, c, q, w1,  w2,  channel,p,  f);
}
double tvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoother interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double t,w1_t,w2_t;
        t = q;
        w1_t = w1;
        w2_t = w2;
        double value=0;


        if(abs(t) < bfreqs[nw/2]){ if (t >= 0) {t = bfreqs[nw/2];} else{t = bfreqs[nw/2-1];};};
        if(abs(w1_t) < ffreqs[nw/2]){if (w1_t>= 0) {w1_t = ffreqs[nw/2];} else{w1_t =  ffreqs[nw/2-1];};};
        if(abs(w2_t) < ffreqs[nw/2]){if (w2_t > 0) {w2_t =  ffreqs[nw/2];} else{w2_t = ffreqs[nw/2-1];};};

        value += R_vvalsmooth(a,b,c,t,w1_t,w2_t) + K1_vvalsmooth(a,b,c,t) + K2_vvalsmooth(a,b,c,t,w1_t) + K2b_vvalsmooth(a,b,c,t,w2_t) ;//K2b is extracted from K2 by the symmetry relations
        return value;}
    else{return 0;}


}
void tvert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){

            R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};
    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void tvert::K1_setvert(int a, int b, int c,int i, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>=0 ){

            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};

    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
void tvert::K2_setvert(int a, int b, int c, int i, int j, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value ;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void tvert::K2b_setvert(int a, int b, int c, int i, int j,  double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0){


            K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
double tvert::R_vval(int a_raw, int b_raw,int c_raw,  int i, int j, int k){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0 ){

            value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}
        else if(sym == 1 ){

            if(i >= nw3/2){
                if(j >= nw3/2 ){

                    value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

                else if(i >= nw3/2  && j < nw3/2 ){

                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)]);};}


            else if(i < nw3/2){
                if(k >= nw3/2 ){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)];}//note that k and j are interchanged

                else if( k < nw3/2 ){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);};//note that k and j are interchanged
            };}

        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;x = site_switch(x.a,x.b,x.c);};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;};
            value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};
    };
    if(abs(value)<1e-100){value=0;};



    return value;
}
double tvert::K1_vval(int a_raw, int b_raw,int c_raw,  int i){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw1 && i >= 0){
        if(sym==0){

            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 ||sym==2){
            if (i>=nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);

                i = nw1-1-i;
                value += K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)];};
        };};
    if(abs(value)<1e-100){value=0;};



    return value;
}
double tvert::K2_vval(int a_raw, int b_raw,int c_raw,  int i, int j){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp

    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >=0){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j>=nw2/2){

                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}

            else if(j < nw2/2){


                j = nw2-1-j;
                value += conj(K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)]);
            };}
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;};
            if(j < nw2/2){j = nw2-1-j;};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };};
    if(abs(value)<1e-100){value=0;};


    return value;
}
double tvert::K2b_vval(int a_raw, int b_raw,int c_raw,  int i, int j){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >=0){
        if(sym==0){

            value += K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1 || sym==2){

            site x = site_switch(a,b,c);
            value += K2_vval(x.a,x.b,x.c, nw2-1-i, j );};};
    if(abs(value)<1e-100){value=0;};



    return value;
}
double tvert::R_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(t)+ 1e-6< bfreqs[(nw3+nw)/2-1] && abs(w1)+ 1e-6<ffreqs[(nw3+nw)/2-1]&& abs(w2)+ 1e-6<ffreqs[(nw3+nw)/2-1]){
        if (sym==0){


            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);

            int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

            value = ((R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                    ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){
            if(t > 0){
                if((w1)>0 ){

                    int T = fconv_n(t,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);

                    int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

                    value = ((R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                            ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }

                else if( (w1)<0){

                    int T = fconv_n(t,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);

                    int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;
                    value = conj((R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}


            else if (t < 0  ){
                if((w2)>0 ){
                    site x = site_switch(a,b,c);


                    int T = fconv_n(-t,nw3);
                    int W1 = fconv_n(w2,nw3);
                    int W2 = fconv_n(w1,nw3);

                    int T_vert = fconv_n(-t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;

                    t = -t;

                    value = (R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2]));
                }



                else if ((w2)<0 ){
                    site x = site_switch(a,b,c);


                    int T = fconv_n(-t,nw3);
                    int W1 = fconv_n(-w2,nw3);
                    int W2 = fconv_n(-w1,nw3);

                    int T_vert = fconv_n(-t,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1,nw3)-(nw3-nw3_w2);

                    double w1_eff = -w2;
                    double w2_eff = -w1;

                    t = -t;

                    value = conj((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };};}
        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, t_eff =t;
            if(t<0){t_eff = -t; };
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1; x = site_switch(x.a,x.b,x.c);};
            if((w1_eff)<0){w1_eff = -w1_eff; w2_eff = -w2_eff;};
            t = t_eff;


            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int T_vert = fconv_n(t,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);

            value = ((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                    ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

        };

    }
    else{
        int i,j,k;
        if(abs(t)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(t,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value += R_vval(a,b,c,i,j,k);};};
    if(abs(value)<1e-100){value=0;};



    return value;
}
double tvert::K1_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t){



    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(t)+ 1e-6< bfreqs[(nw1+nw)/2-1]){
        if(sym==0){

            int T = fconv_n(t,nw1);

            int T_vert = fconv_n(t,nw1)-(nw1-nw1_q);
            value = ((K1[a+(nuc_eff-1)/2][b][c-1][T_vert ]*(ffreqs[(nw-nw1)/2+T+1]-t)
                    +K1[a+(nuc_eff-1)/2][b][c-1][T_vert +1]*(-ffreqs[(nw-nw1)/2+T]+t))/((ffreqs[(nw-nw1)/2+T+1]-ffreqs[(nw-nw1)/2+T])));}
        else if(sym==1 || sym==2){
            if (t>0){
                int T = fconv_n(t,nw1);
                int T_vert = fconv_n(t,nw1)-(nw1-nw1_q);
                value = ((K1[a+(nuc_eff-1)/2][b][c-1][T_vert ]*(bfreqs[T+1]-t)
                        +K1[a+(nuc_eff-1)/2][b][c-1][T_vert +1]*(-bfreqs[T]+t))/((bfreqs[T+1]-bfreqs[T])));}
            else if (t<0){

                site x = site_switch(a,b,c);
                double t_eff = -t;

                int T = fconv_n(t_eff,nw1);
                int T_vert = fconv_n(t_eff,nw1)-(nw1-nw1_q);

                value = ((K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert ]*(bfreqs[T+1]-t_eff)
                        +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][T_vert +1]*(-bfreqs[T]+t_eff))/(bfreqs[T+1]-bfreqs[T]));
            };};}
    else{ int i;
        if(abs(t)<= bfreqs[nw1-1] ){
            i= fconv_n(t,nw1);

            value = K1_vval(a,b,c,i);};
    };



    return value;
}
double tvert::K2_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w1){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(t)+ 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int T = fconv_n(t,nw2);
            int W1 = fconv_n(w1,nw2);

            int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                    ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){
            if(w1 > 0){

                int T = fconv_n(t,nw2);
                int W1 = fconv_n(w1,nw2);

                int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                        ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
            else if (w1<0){

                double w1_eff = - w1;
                int T = fconv_n(t,nw2);
                int W1 =  fconv_n(w1_eff,nw2);

                int T_vert = fconv_n(t,nw2)-(nw2-nw2_q);
                int W1_vert =  fconv_n(w1_eff,nw2)-(nw2-nw2_w1);



                value += conj((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                        +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                        ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));};}
        else if(sym==2){
            double w1_eff = w1;
            double t_eff = t;
            if(t<0){t_eff = -t; };
            if(w1_eff<0){w1_eff = -w1_eff;};

            int T = fconv_n(t_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int T_vert = fconv_n(t_eff,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);
            value = ((K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                    +K2[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                    ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        };}
    else{ int i,j;
        if(abs(t)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(t,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};
    };
    if(abs(value)<1e-100){value=0;};



    return value;
}
double tvert::K2b_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(t)+ 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w2)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int T = fconv_n(T,nw2);
            int W2 = fconv_n(w2,nw2);

            int T_vert = fconv_n(T,nw2)-(nw2-nw2_q);
            int W2_vert = fconv_n(w2,nw2)-(nw2-nw2_w1);

            value += ((K2b[a+(nuc_eff-1)/2][b][c-1][T_vert][W2_vert]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W2_vert]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][T_vert][W2_vert+1]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][T_vert+1][W2_vert+1]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2])));}
        else if(sym==1 ||sym==2){
            site x = site_switch(a,b,c);

            value += K2_vvalsmooth(x.a,x.b,x.c, -t, w2 );};}
    else{int i,j;
        if(abs(t)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(t,nw2);
            j = fconv_n(w2,nw2);

            value += K2b_vval(a,b,c,i,j);};

    };



    return value;
}
//non-member functions
tvert operator*(double alpha, const tvert& vertex){
    tvert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;

}
tvert operator*(const tvert& vertex,double alpha){
    tvert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;
}
tvert operator+(const tvert& vertex1,const tvert& vertex2){
    tvert vertex3;
    int sum_K1 = (sym==0?0:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);

#pragma omp for
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){

                //add K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    double value = vertex1.K1[a][b][c][i] +  vertex2.K1[a][b][c][i];
                    if(abs(value)>1e-16){
                        vertex3.K1[a][b][c][i] =  value;};};
                //add K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        double value = vertex1.K2[a][b][c][i][j] +  vertex2.K2[a][b][c][i][j];
                        if(abs(value)>1e-16){
                            vertex3.K2[a][b][c][i][j] = value;};
                        if(sym==0){
                            vertex3.K2b[a][b][c][i][j] = vertex1.K2b[a][b][c][i][j] + vertex2.K2b[a][b][c][i][j];};};};
                //add R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){

                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                double value = vertex1.R[a][b][c][i][j][k] + vertex2.R[a][b][c][i][j][k];
                                if(abs(value)>1e-16){
                                    vertex3.R[a][b][c][i][j][k] = value;};
                            };};};};};};};
    return vertex3;
}


/*****************FUNCTIONS FOR THE U-VERTEX********************************************/

double uvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance

        double u,w1_u,w2_u;
        if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2;
            w2_u = (-w1+w2+q)/2;}
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2;
            w2_u = (w1+w2+q)/2;}
        else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2;
            w2_u = (w1+w2)/2;};

        double value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};



//        if(u > bfreqs[nw-1] && abs(u - bfreqs[nw-1]) < 1e-12){u = bfreqs[nw-1];}
//        else if(u<bfreqs[0] && abs(u - bfreqs[0]) < 1e-12){u = bfreqs[0];};


//        if(w1_u > ffreqs[nw-1] && abs(w1_u - ffreqs[nw-1]) < 1e-12){w1_u = ffreqs[nw-1];}
//        else if(w1_u<ffreqs[0] && abs(w1_u - ffreqs[0]) < 1e-12){w1_u = ffreqs[0];};


//        if(w2_u > ffreqs[nw-1] && abs(w2_u - ffreqs[nw-1]) < 1e-12){w2_u = ffreqs[nw-1];}
//        else if(w2_u<ffreqs[0] && abs(w2_u - ffreqs[0]) < 1e-12){w2_u = ffreqs[0];};


        value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u)  + K2b_vvalsmooth(a,b,c,u,w2_u)  ;//K2b is extracted from K2 by the symmetry relations

        return value;  }
    else{return 0;}



}
double uvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
    if(distance(a,b,c) <= d_c){//cutoff distance



        double u,w1_u,w2_u;
        if(channel == 's'){
            u = -w2-w1;
            w1_u = (w1-w2+q)/2;
            w2_u = (-w1+w2+q)/2;}
        else if(channel == 't'){//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
            u = w2-w1;
            w1_u = (w1+w2-q)/2;
            w2_u = (w1+w2+q)/2;}
        else if(channel == 'u'){
            u = q;
            w1_u = w1;
            w2_u = w2;}
        else if (channel == 'v'){//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
            u = w1-w2;
            w1_u = q + (w1-w2)/2;
            w2_u = (w1+w2)/2;};

        double value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};



//        if(u > bfreqs[nw-1] && abs(u - bfreqs[nw-1]) < 1e-12){u = bfreqs[nw-1];}
//        else if(u<bfreqs[0] && abs(u - bfreqs[0]) < 1e-12){u = bfreqs[0];};


//        if(w1_u > ffreqs[nw-1] && abs(w1_u - ffreqs[nw-1]) < 1e-12){w1_u = ffreqs[nw-1];}
//        else if(w1_u<ffreqs[0] && abs(w1_u - ffreqs[0]) < 1e-12){w1_u = ffreqs[0];};


//        if(w2_u > ffreqs[nw-1] && abs(w2_u - ffreqs[nw-1]) < 1e-12){w2_u = ffreqs[nw-1];}
//        else if(w2_u<ffreqs[0] && abs(w2_u - ffreqs[0]) < 1e-12){w2_u = ffreqs[0];};

        if(p==1){
            if(channel=='u'){
                if(f == 'R' || f == 'M'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
                else if(f == 'K' || f == 'L'){value += K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='t'){
                if(f == 'R' || f== 'M'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);};
            };
        }

        else if(p==2){
            if(channel=='u'){
                if(f == 'R' || f == 'L'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K2_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different bare vertex
                else if(f == 'K' || f == 'M'){value += K1_vvalsmooth(a,b,c,u) + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex

            }
            else if (channel=='s' || channel=='t'){
                if(f == 'R' || f== 'L'){value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u);
                };
            };
        };




        return value;  }
    else{return 0;}



}
//overload of previous function
double uvert::vvalsmooth(int red_side, int map,int a, int b, int c,  double q, double w1, double w2, char channel, int p, char f){
  return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
}
double uvert::vvalsmooth(int a, int b, int c,  double q, double w1, double w2){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
    if(distance(a,b,c) <= d_c){//cutoff distance


        double u,w1_u,w2_u;

        u = q;
        w1_u = w1;
        w2_u = w2;

        double value=0;

        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};

        value += R_vvalsmooth(a,b,c,u,w1_u,w2_u) + K1_vvalsmooth(a,b,c,u) + K2_vvalsmooth(a,b,c,u,w1_u) + K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
        return value;  }
    else{return 0;}



}
void uvert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0){


            R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void uvert::K1_setvert(int a, int b, int c, int i, double value){

    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>= 0 ){
            K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void uvert::K2_setvert(int a, int b, int c,int i, int j,double value){


    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0 ){

            K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void uvert::K2b_setvert(int a, int b, int c, int i, int j,  double value){


    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0 ){

            K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)] = value;};}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
double uvert::R_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){


    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        if(sym == 0 ){

            value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}


        else if(sym == 1 ){
            if(i>=nw3/2){
                if(j >=nw3/2){

                    value += R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)];}

                else if(j < nw3/2){
                    site x = site_switch(a,b,c);

                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][j-(nw3-nw3_w1)][k-(nw3-nw3_w2)]);};}

            else if (i<nw3/2 ){
                if(k >= nw3/2){
                    site x = site_switch(a,b,c);

                    i = nw3-1-i;
                    value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)];}//note that k and j are interchanged

                else if ( k < nw3/2){


                    i = nw3-1-i;
                    j = nw3-1-j;
                    k = nw3-1-k;
                    value += conj(R[a+(nuc_eff-1)/2][b][c-1][i-(nw3-nw3_q)][k-(nw3-nw3_w1)][j-(nw3-nw3_w2)]);};//note that k and j are interchanged

        };}
        else if(sym == 2 ){
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;x = site_switch(x.a,x.b,x.c);};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};
            value += R[x.a+(nuc_eff-1)/2][x.b][x.c-1][i_eff-(nw3-nw3_q)][j_eff-(nw3-nw3_w1)][k_eff-(nw3-nw3_w2)];};};
    if(abs(value)<1e-100){value=0;};



    return value;
}
double uvert::K1_vval(int a_raw, int b_raw, int c_raw,  int i){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw1 && i >= 0 ){
        if(sym==0){

            value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
        else if(sym==1 || sym==2){
            if (i>=nw1/2){

                value += K1[a+(nuc_eff-1)/2][b][c-1][i-(nw1-nw1_q)];}
            else if(i<nw1/2){

                site x = site_switch(a,b,c);
                i=nw1-1-i;
                value += K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw1-nw1_q)];};
        };};
    if(abs(value)<1e-100){value =0;};

    return value;
}
double uvert::K2_vval(int a_raw, int b_raw, int c_raw,  int i,int j){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0 ){
        if(sym==0){

            value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1){
            if(j >= nw2/2){


                value += K2[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
            else if (j<nw2/2){

                site x = site_switch(a,b,c);

                j = nw2-1-j;
                value += conj(K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)]);};
        }
        else if(sym==2){
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j;x = site_switch(x.a,x.b,x.c);};

            value += K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];
        };};
    if(abs(value)<1e-100){value =0;};



    return value;
}
double uvert::K2b_vval(int a_raw, int b_raw, int c_raw,   int i, int j){
    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(i < nw2 && i >= 0 && j < nw2 && j >= 0 ){
        if(sym==0){

            value += K2b[a+(nuc_eff-1)/2][b][c-1][i-(nw2-nw2_q)][j-(nw2-nw2_w1)];}
        else if(sym==1||sym==2){
            site x = site_switch(a,b,c);
            value += K2_vval(x.a,x.b,x.c, nw2-1-i, j );};};
    if(abs(value)<1e-100){value =0;};


    return value;
}
double uvert::R_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(u)+ 1e-6< ffreqs[(nw3+nw)/2-1] && abs(w1)+ 1e-6<ffreqs[(nw3+nw)/2-1]&& abs(w2)+ 1e-6<=ffreqs[(nw3+nw)/2-1]){
        if (sym==0){

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);


            int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);

            value += ((R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                    +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                    ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }
        else if(sym==1){
            if(u > 0 ){
                if( (w1)>0){

                    int U = fconv_n(u,nw3);
                    int W1 = fconv_n(w1,nw3);
                    int W2 = fconv_n(w2,nw3);


                    int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w2,nw3)-(nw3-nw3_w2);
                    value += ((R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                            ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));}
                else if((w1)<0){
                    site x = site_switch(a,b,c);

                    int U = fconv_n(u,nw3);
                    int W1 = fconv_n(-w1,nw3);
                    int W2 = fconv_n(-w2,nw3);

                    int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w1,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w2,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w1 ;
                    double w2_eff = -w2 ;
                    value += conj((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                };}

            else if (u < 0 ){
                if((w2)>0){
                    site x = site_switch(a,b,c);

                    int U = fconv_n(-u,nw3);
                    int W1 = fconv_n(w2 ,nw3);
                    int W2 = fconv_n(w1 ,nw3);


                    int U_vert = fconv_n(-u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(w2 ,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(w1 ,nw3)-(nw3-nw3_w2);
                    double w1_eff = w2;
                    double w2_eff = w1;
                    u = -u;

                    value += ((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
                }


                else if ((w2)<0){


                    int U = fconv_n(-u,nw3);
                    int W1 = fconv_n(-w2 ,nw3);
                    int W2 = fconv_n(-w1 ,nw3);

                    int U_vert = fconv_n(-u,nw3)-(nw3-nw3_q);
                    int W1_vert = fconv_n(-w2 ,nw3)-(nw3-nw3_w1);
                    int W2_vert = fconv_n(-w1 ,nw3)-(nw3-nw3_w2);
                    double w1_eff = -w2;
                    double w2_eff = -w1;
                    u = -u;

                    value += conj((R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                            +R[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                            ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));};

            };}

        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, u_eff =u;
            if(u<0){u_eff = -u; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if(w1_eff < 0){w1_eff = -w1_eff; w2_eff = -w2_eff;x = site_switch(x.a,x.b,x.c);};

            u = u_eff;

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int U_vert = fconv_n(u,nw3)-(nw3-nw3_q);
            int W1_vert = fconv_n(w1_eff,nw3)-(nw3-nw3_w1);
            int W2_vert = fconv_n(w2_eff,nw3)-(nw3-nw3_w2);
            value = ((R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1][W2_vert+1]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                    +R[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1][W2_vert+1]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                    ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        }


        ;}
    else{
        int i,j,k;
        if(abs(u)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(u,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value = R_vval(a,b,c,i,j,k);};
    };




    return value;
}
double uvert::K1_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(u)+ 1e-6<= ffreqs[(nw1+nw)/2-1]){
        if(sym==0){

            int U = fconv_n(u,nw1);

            int U_vert = fconv_n(u,nw1)-(nw1-nw1_q);

            value += ((K1[a+(nuc_eff-1)/2][b][c-1][U_vert ]*(bfreqs[(nw-nw1)/2+U+1]-u)
                    +K1[a+(nuc_eff-1)/2][b][c-1][U_vert +1]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));}
        else if(sym==1 ||sym==2){
            if (u>0){


                int U = fconv_n(u,nw1);

                int U_vert = fconv_n(u,nw1)-(nw1-nw1_q);

                value += ((K1[a+(nuc_eff-1)/2][b][c-1][U_vert]*(bfreqs[(nw-nw1)/2+U+1]-u)
                        +K1[a+(nuc_eff-1)/2][b][c-1][U_vert+1]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));}
            else if (u<0){

                site x = site_switch(a,b,c);
                int U = fconv_n(-u,nw1);

                int U_vert = fconv_n(-u,nw1)-(nw1-nw1_q);
                u = -u;

                value += (K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert]*(bfreqs[(nw-nw1)/2+U+1]-u)
                        +K1[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1]*(-bfreqs[(nw-nw1)/2+U]+u))/(bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U]);};};}
    else{
        int i;
        if(abs(u)<= bfreqs[nw1-1] ){
            i= fconv_n(u,nw1);

            value += K1_vval(a,b,c,i);};    };


    return value;
}
double uvert::K2_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;

    double value=0;
    if(abs(u)+ 1e-6<= bfreqs[(nw2+nw)/2-1] && abs(w1)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1,nw2);

            int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

            value += ((K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                    +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
        else if(sym==1){

            if(w1 > 0){

                int U = fconv_n(u,nw2);
                int W1 = fconv_n(w1,nw2);

                int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(w1,nw2)-(nw2-nw2_w1);

                value += ((K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                        +K2[a+(nuc_eff-1)/2][b][c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
            }
            else if (w1<0){

                site x = site_switch(a,b,c);

                int U = fconv_n(u,nw2);
                int W1 = fconv_n(-w1,nw2);

                int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
                int W1_vert = fconv_n(-w1,nw2)-(nw2-nw2_w1);
                double w1_eff = -w1 ;
                value += conj((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert ][W1_vert ]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert +1][W1_vert ]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert ][W1_vert +1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                        +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert +1][W1_vert +1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

            };
        }
        else if(sym==2){
            site x(a,b,c);
            double w1_eff = w1;
            if(u<0){u = -u;x = site_switch(x.a,x.b,x.c); };
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
            int W1_vert = fconv_n(w1_eff,nw2)-(nw2-nw2_w1);
            value += ((K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert][W1_vert+1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                    +K2[x.a+(nuc_eff-1)/2][x.b][x.c-1][U_vert+1][W1_vert+1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        };}
    else{ int i,j;
        if(abs(u)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(u,nw2);
            j = fconv_n(w1,nw2);

            value += K2_vval(a,b,c,i,j);};   };



    return value;
}
double uvert::K2b_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w2){

    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;
    if(abs(u)+ 1e-6< ffreqs[(nw2+nw)/2-1] && abs(w2)+ 1e-6< ffreqs[(nw2+nw)/2-1]){
        if(sym==0){

            int U = fconv_n(u,nw2);
            int W2 = fconv_n(w2,nw2);

            int U_vert = fconv_n(u,nw2)-(nw2-nw2_q);
            int W2_vert = fconv_n(w2,nw2)-(nw2-nw2_w1);

            value += ((K2b[a+(nuc_eff-1)/2][b][c-1][U_vert ][W2_vert ]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][U_vert +1][W2_vert ]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][U_vert ][W2_vert +1]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                    +K2b[a+(nuc_eff-1)/2][b][c-1][U_vert +1][W2_vert +1]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2])));}
        else if(sym==1 || sym==2){
            site x = site_switch(a,b,c);
            value += K2_vvalsmooth(x.a,x.b,x.c, -u, w2 );};}

    else{ int i,j;
        if(abs(u)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(u,nw2);
            j = fconv_n(w2,nw2);

            value += K2b_vval(a,b,c,i,j);};};



    return value;
}
//non-member functions
uvert operator*(double alpha, const uvert& vertex){
    uvert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;
}
uvert operator*(const uvert& vertex,double alpha){
    uvert vertex2;
    int sum_K1 = (sym==0?0:nw1/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);
#pragma omp for collapse(3)
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){


                //K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    vertex2.K1[a][b][c][i] = alpha * vertex.K1[a][b][c][i];};
                //K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        vertex2.K2[a][b][c][i][j] = alpha * vertex.K2[a][b][c][i][j];
                        if(sym==0){
                            vertex2.K2b[a][b][c][i][j] = alpha * vertex.K2b[a][b][c][i][j];};};};
                // R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){
                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                vertex2.R[a][b][c][i][j][k] = alpha * vertex.R[a][b][c][i][j][k];
                            };};};};};};};
    return vertex2;
}
uvert operator+(const uvert& vertex1,const uvert& vertex2){
    uvert vertex3;

    int sum_K1 = (sym==0?0:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = 0;}
    else if(sym==2){sum_K2_i = nw2/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = 0;}
    else if(sym==1 || sym==2){sum_K2_j = nw2/2;};

    int sum_R = (sym==0?0:nw3/2);

#pragma omp for
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){

                //add K1 contributions
                for(int i=sum_K1-(nw1-nw1_q) ; i<nw1-(nw1-nw1_q); i++){
                    double value = vertex1.K1[a][b][c][i] +  vertex2.K1[a][b][c][i];
                    if(abs(value)>1e-16){
                        vertex3.K1[a][b][c][i] =  value;};};
                //add K2 contributions
                for(int i=sum_K2_i-(nw2-nw2_q) ; i<nw2-(nw2-nw2_q); i++){
                    for(int j=sum_K2_j-(nw2-nw2_w1) ; j<nw2-(nw2-nw2_w1); j++){
                        double value = vertex1.K2[a][b][c][i][j] +  vertex2.K2[a][b][c][i][j];
                        if(abs(value)>1e-16){
                            vertex3.K2[a][b][c][i][j] = value;};
                        if(sym==0){
                            vertex3.K2b[a][b][c][i][j] = vertex1.K2b[a][b][c][i][j] + vertex2.K2b[a][b][c][i][j];};};};
                //add R contributions
                for(int i=sum_R-(nw3-nw3_q) ; i<nw3-(nw3-nw3_q); i++){
                    for(int j=sum_R -(nw3-nw3_w1); j<nw3-(nw3-nw3_w1); j++){
                        for(int k=0-(nw3-nw3_w2) ; k<nw3-(nw3-nw3_w2); k++){

                            if(sym !=2 || (abs(j-nw3/2+(nw3-nw3_w1)) <= abs(k-nw3/2+(nw3-nw3_w2)))){
                                double value = vertex1.R[a][b][c][i][j][k] + vertex2.R[a][b][c][i][j][k];
                                if(abs(value)>1e-16){
                                    vertex3.R[a][b][c][i][j][k] = value;};
                            };};};};};};};

    return vertex3;
}

/*****************FUNCTIONS FOR THE IRREDUCIBLE VERTEX********************************************/
double irreducible::vval(int a_raw, int b_raw, int c_raw){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
double irreducible::vvalsmooth(int a_raw, int b_raw, int c_raw){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        //cout << a << " " << b << " " << c << endl;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
double irreducible::vvalsmooth(int a_raw, int b_raw, int c_raw,double q, double w1, double w2, char channel, int par, char f){
    if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
        site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
        int a,b,c;
        a = project.a;
        b = project.b;
        c = project.c;
        return U_bare[a+(nuc_eff-1)/2][b][c-1];}
    else{return 0;};}
void irreducible::setvert(int a_raw, int b_raw, int c_raw,  double value){
    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    U_bare[a+(nuc_eff-1)/2][b][c-1] = value;}
//operators for irreducible vertex
irreducible operator*(double alpha, const irreducible& vertex){
    irreducible result;
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){
                result.U_bare[a][b][c] = alpha * vertex.U_bare[a][b][c];};};};
    return result;}
irreducible operator*(const irreducible& vertex,double alpha){
    irreducible result;
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){
                result.U_bare[a][b][c] = alpha * vertex.U_bare[a][b][c];};};};
    return result;}
irreducible operator+(const irreducible& vertex1,const irreducible& vertex2){
    irreducible result;
    for(int a=0; a<nuc_eff; a++){
        for(int b=0; b<(nuc_eff-1)/2+1; b++){
            for(int c=0; c<3; c++){
                result.U_bare[a][b][c] =  vertex1.U_bare[a][b][c] + vertex2.U_bare[a][b][c];};};};
    return result;}

/*****************MEMBER FUNCTIONS OF SUSCEPTIBILITY CLASS (WRITE/ADD_WRITE/READ)********************************************/

void Susc::write(int a ,int b,int c,double value){

    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) ){//only save sites in upper half
        if( distance(a,b,c) <= d_c){


            Sus[a+(nuc_eff-1)/2][b][c-1] = value;};
    }
    else{cout << "cannot write Susceptibility in lower half" << endl;};
}
void Susc::add_write(int a ,int b,int c,double value){

    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) ){//only save sites in upper half
        if( distance(a,b,c) <= d_c){
            Sus[a+(nuc_eff-1)/2][b][c-1] += value;};
    }
    else{cout << "cannot write Susceptibility in lower half" << endl;};
}
double Susc::vval(int a_raw, int b_raw, int c_raw){
    
    site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
    int a,b,c;
    a = project.a;
    b = project.b;
    c = project.c;
    double value=0;

    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        value = Sus[a+(nuc_eff-1)/2][b][c-1];};
    if(abs(value)<1e-100){value= 0;};

    return value;
}
/*************************************************************************************************************/








/*****************************************operators concerning parvert objetcs********************************************************/

parvert<svert> operator+(parvert<svert> vertex1,parvert<svert> vertex2){
    parvert<svert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<tvert> operator+(parvert<tvert> vertex1,parvert<tvert> vertex2){
    parvert<tvert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<uvert> operator+(parvert<uvert> vertex1,parvert<uvert> vertex2){
    parvert<uvert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<irreducible> operator+(parvert<irreducible> vertex1,parvert<irreducible> vertex2){
    parvert<irreducible> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<svert> operator+=(parvert<svert> vertex1,parvert<svert> vertex2){
    parvert<svert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<tvert> operator+=(parvert<tvert> vertex1,parvert<tvert> vertex2){
    parvert<tvert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<uvert> operator+=(parvert<uvert> vertex1,parvert<uvert> vertex2){
    parvert<uvert> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<irreducible> operator+=(parvert<irreducible> vertex1,parvert<irreducible> vertex2){
    parvert<irreducible> result;
    result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
    result.densvertex = vertex1.densvertex + vertex2.densvertex;
    return result;
}
parvert<svert> operator*(double alpha ,parvert<svert>& vertex){
    parvert<svert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<svert> operator*(parvert<svert>& vertex,double alpha){
    parvert<svert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<tvert> operator*(double alpha ,parvert<tvert>& vertex){
    parvert<tvert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<tvert> operator*(parvert<tvert>& vertex,double alpha){
    parvert<tvert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<uvert> operator*(double alpha ,parvert<uvert>& vertex){
    parvert<uvert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<uvert> operator*(parvert<uvert>& vertex,double alpha){
    parvert<uvert> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<irreducible> operator*(double alpha ,parvert<irreducible>& vertex){
    parvert<irreducible> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
parvert<irreducible> operator*(parvert<irreducible>& vertex,double alpha){
    parvert<irreducible> result;
    result.spinvertex = alpha * vertex.spinvertex;
    result.densvertex = alpha * vertex.densvertex;
    return result;
}
/*************************************************************************************************************/






/*****************************************FUNCTIONS FOR FULL VERTEX "FULLVERT"********************************************************/
//arguments are equivalent to those in the simple vertex functions

double fullvert::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel){
    double result = irred.vvalsmooth(a,b,c) + svertex.vvalsmooth(a,b,c,q,w1,w2,channel) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel);
    if(abs(result)>1e-16){
        return  result;}
    else{return 0;};
}
double fullvert::vvalsmooth(int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){
    double result=0;

    if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
        result += irred.vvalsmooth(a,b,c);

    }
    else if( p==2 && (f=='K' || f== 'M')){
        result += irred.vvalsmooth(a,b,c);

    };



    result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
    //if(p==2 && f=='L'){cout <<  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << " " << tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f)  << " " <<  uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) << endl;};
    if(abs(result)<1e-16){
        result  =0;};
    return result;
}
double fullvert::vvalsmooth(p, int a, int b, int c, double q, double w1, double w2, char channel, int p, char f){// red_side: if only complementary channel in one vertex, which one is reduced? (0/1/2), p: is this the left/upper (1) or the right/lower (2) vertex of the bubble?, f: diagrammatic class that is computed
    double result=0;

    if(red_side != p){
        if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
            result = irred.vvalsmooth(a,b,c);
        }
        else if( p==2 && (f=='K' || f== 'M')){
            result = irred.vvalsmooth(a,b,c);
        };
        result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
    }


    else if(red_side == p){

        if(map==0){
            if(channel == 's'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 'u'){  result =  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);};}
        else if(map==1){
            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);}
            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
            else if(channel == 'u'){  result =  uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) +  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) ;}
        };

    };

    if(abs(result)<1e-16){
        result  =0;};
    return result;
}
fullvert operator*(double alpha, const fullvert& vertex){
    fullvert result;
    result.irred = alpha * vertex.irred;
    result.svertex = alpha *vertex.svertex;
    result.tvertex = alpha * vertex.tvertex;
    result.uvertex = alpha * vertex.uvertex;
    return result;
}
fullvert operator+( const fullvert& vertex1, const fullvert& vertex2){
    fullvert result;
    result.irred = vertex1.irred + vertex2.irred ;
    result.svertex = vertex1.svertex + vertex2.svertex ;
    result.tvertex = vertex1.tvertex + vertex2.tvertex ;
    result.uvertex = vertex1.uvertex + vertex2.uvertex ;
    return result;
}








/*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
complex<double> self::sval(int i){
    if(i>=0 && i< nw){
        return selfenergy[i];}
    else{return zero;};
}
complex<double> self::svalsmooth(double w){//smoothly interpolates for values between discrete frequ values of mesh
    complex<double> value;
    if(abs(w) < ffreqs[nw-1]){
        int W = fconv(w);
            value += ((selfenergy[W]*(ffreqs[W+1]-w)+selfenergy[W+1]*(-ffreqs[W]+w))/(ffreqs[W+1]-ffreqs[W]));
        }

    else{
        return zero;};
    if(abs(real(value))<1e-100){value.real(0);};

    if(abs(imag(value))<1e-100){value.imag(0);};

    return value;
}
void self::setself(int i, complex<double> val){
    if(i>=0 && i< nw){

        selfenergy[i] = val ;};
}

//operators for self energy
self operator*(double alpha,const self& self1){//sum operator overloading
    self self2;
    for(int i=0; i<nw; i++){
        self2.selfenergy[i] = self1.selfenergy[i] * alpha;};
    return self2;
}
self operator*(const self& self1, double alpha){//sum operator overloading
    self self2;
    for(int i=0; i<nw; i++){
        self2.selfenergy[i] = self1.selfenergy[i] * alpha;};
    return self2;
}
self operator+(const self& self1, const self& self2){//sum operator overloading
    self self3;
    for(int i=0; i<nw; i++){
        complex<double> sum = self1.selfenergy[i] + self2.selfenergy[i] ;
        if(abs(sum)>1e-16){
        self3.selfenergy[i] = sum ;}
         else {self3.selfenergy[i] =zero;};
    };
    return self3;
}
self operator+=(const self& self1, const self& self2){//sum operator overloading
    self self3;
    for(int i=0; i<nw; i++){
        complex<double> sum = self1.selfenergy[i] + self2.selfenergy[i] ;
        if(abs(sum)>1e-16){
        self3.selfenergy[i] = sum ;}
         else {self3.selfenergy[i] =zero;};};
    return self3;
}



//operators containing state objetcs
state operator+(state state1, state state2){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred + state2.vertex.spinvertex.irred;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex + state2.vertex.spinvertex.svertex;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex + state2.vertex.spinvertex.tvertex;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex + state2.vertex.spinvertex.uvertex;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred + state2.vertex.densvertex.irred;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex + state2.vertex.densvertex.svertex;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex + state2.vertex.densvertex.tvertex;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex + state2.vertex.densvertex.uvertex;
    result.selfenergy = state1.selfenergy + state2.selfenergy;
    return result;
}
state operator*(double alpha, state state1){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex * alpha;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex * alpha;

    result.selfenergy = alpha * state1.selfenergy;
    return result;
}
state operator*(state state1, double alpha){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.svertex = state1.vertex.spinvertex.svertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;
    result.vertex.spinvertex.uvertex = state1.vertex.spinvertex.uvertex * alpha;
    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.svertex = state1.vertex.densvertex.svertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;
    result.vertex.densvertex.uvertex = state1.vertex.densvertex.uvertex * alpha;

    result.selfenergy = alpha * state1.selfenergy;
    return result;
}


/*****************************************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/
complex<double> propag(double Lambda, double w, self selfenergy,self diffselfenergy, char type){

    if(w!=0){
        complex<double> value;
        complex<double> iw(0.0,w);
        if(type == 'g'){//regular undifferentiated greensfunction
            complex<double> g0(0.0,0.0);
            if(reg==1){

                if(fabs(w) > Lambda){
                    g0 =1./iw;
                    value += 1./(1./g0-selfenergy.svalsmooth(w));

                }
                else if (fabs(w) == Lambda){
                    g0 = 1./iw;
                    value += 0.5/(1./g0-selfenergy.svalsmooth(w));
                    //this is the implementation of Morris Lemma, see SB.II, p. 9
                }
            }
            else if(reg==2){

                g0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw;

                value = 1./(1./g0-selfenergy.svalsmooth(w));
            };

        }
        else if(type == 's'){//single scale propagator

            if(reg==1){
                if(fabs(w) == Lambda){

                    value += (-1./(iw-selfenergy.svalsmooth(w)));}
                else{return zero;};
            }
            else if(reg==2){

                complex<double> G0 = (1.-exp(-pow(fabs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += pow((1. + G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * 1./iw;//see page 20 in SB II

            };
        }
        else if(type == 'k'){//katanin substitution
            if(reg==1){
                complex<double> S,G;
                if(fabs(w) == Lambda){ S = - 1./(iw-selfenergy.svalsmooth(w));};//single scale propagator
                if(fabs(w) > Lambda){ G = 1./(iw-selfenergy.svalsmooth(w));}
                else if(fabs(w) == Lambda) {G = 0.5/(iw-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9
                value = (S + G * G * diffselfenergy.svalsmooth(w));
            }
            else if(reg==2){
                complex<double> G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += pow((1. + G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * 1./iw
                        + G * G * diffselfenergy.svalsmooth(w);//see page 20 in SB II

            };

        }
        else if(type == 'e'){//only self energy extension of katanin substitution ( needed for sharp regulator)
            if(reg==1){
                complex<double> G;
                if(fabs(w) > Lambda){ G = 1./(iw-selfenergy.svalsmooth(w));}
                else if(fabs(w) == Lambda) {G = 0.5/(iw-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9
                value += (G * G * diffselfenergy.svalsmooth(w));
            }
            else if (reg==2){
                complex<double> G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * 1./iw; //bare greens function with smoothened cutoff at freq ffreqs[i]
                complex<double> G = 1./( 1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
                value += G * G * diffselfenergy.svalsmooth(w);
            };

        }
        else{cout << "could not resolve propagator type" << endl;};

        return value;
    }
    else{return 0;};}


/*****************************************FUNCTION TO COMPUTE SUSCEPTIBILITY********************************************************/
Susc suscept(double Lambda, parvert<fullvert>& vertex, self selfenergy){

    Susc result;
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:

    gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (1500);


    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);




    double one =1.0;
    fullvert bare;
    bare.irred.setvert(0,0,1,one);

    double first = real(0.5*ububble(0,0,0,w,Lambda, bare,0,0,1,bare,0,0,1,'g','g',selfenergy,selfenergy,bfreqs[nw/2],wlimit,wlimit,'K'));
    cout << "first susc-contr. computed.." << endl;

    //second and third diagram;

    fullvert second_third_left;

    fullvert spin_dens =  vertex.densvertex + (-1)*vertex.spinvertex;



    for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
        double w1 = ffreqs[j];
        second_third_left.uvertex.K2_setvert(0,0,1,nw2/2,j-(nw-nw2)/2, ububble(0,0,0,w,Lambda,bare,0,0,1, spin_dens,0,0,1,'g','g',selfenergy,selfenergy,bfreqs[nw/2],w1,wlimit,'L'));

    };



    double second_third = 0.5* ububble(0,0,0,w,Lambda,second_third_left,0,0,1,bare,0,0,1,'g','g',selfenergy,selfenergy,bfreqs[nw/2],wlimit,wlimit,'K');


    double firstthree = first+second_third;
    if(abs(firstthree)<1e-100){firstthree =0;};
    result.write(0,0,1,firstthree); //write diagrams 1-3

    cout << "second & third susc-contr. computed.." << endl;

    //fourth diagram:


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4; c++){



                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

                    tvert fourth_left;


                    fourth_left.K1_setvert(a,b,c,nw1/2, tbubble(0,0,1,w,Lambda,bare,0,0,1, vertex.spinvertex,a,b,c,'g','g',selfenergy,selfenergy,bfreqs[nw/2],wlimit,wlimit,'K'));
                    cout << fourth_left.K1_vval(a,b,c,nw1/2) << endl;


                    for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
                        double w1 = ffreqs[j];
                        fourth_left.K2_setvert(a,b,c,nw2/2,j-(nw-nw2)/2, tbubble(0,0,1,w,Lambda,bare,0,0,1, vertex.spinvertex,a,b,c,'g','g',selfenergy,selfenergy,bfreqs[nw/2],w1,wlimit,'L'));
                    };

                    result.add_write(a,b,c,(-1)*tbubble(0,0,0,w,Lambda, fourth_left,a,b,c,bare,0,0,1,'g','g',selfenergy,selfenergy,bfreqs[nw/2],wlimit,wlimit,'K'));//factor (-1) is implicitly contained in definition of t-bubble but since two t-bubble are computed here, they cancel and one must indeed include factor (-1) from combinatorics
                };};};};
    cout << "fourth susc-contr. computed.." << endl;

    return result;
}









