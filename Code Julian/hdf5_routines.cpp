#include"kagome_facil.hpp"
#include"bubbles_mpi.hpp"

#include <iostream>
#include<iomanip>
#include <complex>
#include<cmath>
#include <vector>
#include<fstream>
#include<type_traits>
#include <string>
//#include<tgmath.h>
#include<cstdlib>
/********************************constants concerning HDF5 data format*************************/
const int    NX = nw;                     // dataset dimensions
const int    NY = nw;
const int    NZ = nw;                       // dataset dimensions
const int    RANK_R = 2;
const int    RANK_K1 = 2;
const int    RANK_K2 = 2;
const int    RANK_K2b = 2;
const int    RANK_irreducible = 2;
const int    RANK_sus = 2;
const int    RANK_self = 2;


//function that write inital state into hdf5 format.  The second argument denotes the iteration number to which it is written. The third argument denotes the total number of saved iterations in the file
void write_hdf(const H5std_string FILE_NAME,double Lambda_i, long Lambda_size,int total_iterations, state& state_in){
    //Try block to detect exceptions raised by any of the calls inside it
    cout << "Starting to copy to buffer.." << endl;
    try
    {
        int Lambda_it=0;
        typedef struct complex_t{
            double spin_re;   /*real part */
            double dens_re;   /*real part */
        }complex_t;





        //buffer self energy
        auto selfenergy = new double[nw];//irrdeucible vertex
        for(int i=0; i<nw; i++){

            selfenergy[i] = state_in.selfenergy.sval(i);//saves the imaginary part of the purely imaginary self energy

        };




        //buffer irreducible vertex:

        auto irreducible_class = new complex_t[irred_dim1];//irrdeucible vertex


        for(int i=0;i<irred_dim1;i++){

            irreducible_class[i].spin_re = state_in.vertex.spinvertex.irred.acc(i);
            irreducible_class[i].dens_re = state_in.vertex.densvertex.irred.acc(i);
        };





        //buffer all R-class-arrays:


        auto R_class_s = new complex_t[K3_dim1];
        auto R_class_t = new complex_t[K3_dim1];
        auto R_class_u = new complex_t[K3_dim1];
        for(int i=0;i<K3_dim1;i++){
            R_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K3_acc(i);
            R_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K3_acc(i);

            R_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K3_acc(i);
            R_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K3_acc(i);

            R_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K3_acc(i);
            R_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K3_acc(i);
        };

        auto K1_class_s = new complex_t[K1_dim1];
        auto K1_class_t = new complex_t[K1_dim1];
        auto K1_class_u = new complex_t[K1_dim1];
        for(int i=0;i<K1_dim1;i++){
            K1_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K1_acc(i);
            K1_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K1_acc(i);

            K1_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K1_acc(i);
            K1_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K1_acc(i);

            K1_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K1_acc(i);
            K1_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K1_acc(i);
        };

        auto K2_class_s = new complex_t[K2_dim1];
        auto K2_class_t = new complex_t[K2_dim1];
        auto K2_class_u = new complex_t[K2_dim1];
        for(int i=0;i<K2_dim1;i++){
            K2_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K2_acc(i);
            K2_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K2_acc(i);

            K2_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K2_acc(i);
            K2_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K2_acc(i);

            K2_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K2_acc(i);
            K2_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K2_acc(i);
        };

#if sym==0
        auto K2b_class_s = new complex_t[K2b_dim1];
        auto K2b_class_t = new complex_t[K2b_dim1];
        auto K2b_class_u = new complex_t[K2b_dim1];
        for(int i=0;i<K2_dim1;i++){
            K2b_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K2b_acc(i);
            K2b_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K2b_acc(i);

            K2b_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K2b_acc(i);
            K2b_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K2b_acc(i);

            K2b_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K2b_acc(i);
            K2b_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K2b_acc(i);
        };
#endif

        //susceptibility
        double Suscept[nuc_eff*nuc_eff*3];
        //suscptibility has only zero-entries at beginning of flow

        //
        //        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        //            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
        //                for (int c=1; c<4 ; c++){
        //                    if(distance(a,b,c) <= d_c){


        //                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = state_in.sus.vval(a,b,c);

        //                    };};};};




        double ferm_freqs[nw];
        for(int i=0; i<nw; i++){
            ferm_freqs[i] = ffreqs[i];
        };

        double bos_freqs[nw];
        for(int i=0; i<nw; i++){
            bos_freqs[i] = ffreqs[i];
        };

        double Lambda_list[Lambda_size];

        Lambda_list[0] = Lambda_i;

        for(int i=1; i<Lambda_size; i++){
            Lambda_list[i]=0;};

        cout << "Buffer ready. Preparing for saving into Hdf5 file.." << endl;
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        //  H5::Exception::dontPrint();

        // Create a new file using the default property lists.
        H5::H5File* file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);


        //Create the memory datatype
        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        //  mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
        //  mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);



        complex_t fillvalue_vert;
        fillvalue_vert.spin_re = 0;
        fillvalue_vert.dens_re = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype1, &fillvalue_vert);





        // Create the dimension arrays for objects in file.
        hsize_t R_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K1_dim1)};   // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#if sym==0
        hsize_t K2b_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#endif
        hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims[]={static_cast<hsize_t>(total_iterations),static_cast<hsize_t>(nuc_eff*nuc_eff*3)};
        hsize_t freq_dims[]={static_cast<hsize_t>(nw)};                      // dataset dimensions for freqs grid (one time nw entries)
        hsize_t bos_freq_dims[]={static_cast<hsize_t>(nw)};
        hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nw)};
        hsize_t Lambda_dims[]={static_cast<hsize_t>(total_iterations)};

        // Create the dimension arrays for objects in buffer.
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t R_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t sus_dims_buffer[]= {static_cast<hsize_t>(nuc_eff*nuc_eff*3)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#if sym==0
        hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif
        hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};



        // Create the data space for the dataset in file.
        H5::DataSpace dataspacevertex_R_s(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_s(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K2_s(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_s(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
#endif
        H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)


        H5::DataSpace dataspacevertex_R_t(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K2_t(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_t(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
#endif


        H5::DataSpace dataspacevertex_R_u(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_u(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K2_u(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_u(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
#endif



        H5::DataSpace dataspacevertex_sus(RANK_sus,sus_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacefreqs(1, freq_dims);
        H5::DataSpace dataspacefreqs_bos(1, bos_freq_dims);
        H5::DataSpace dataspaceself(RANK_self, self_dims);
        H5::DataSpace dataspacelambda(1, Lambda_dims);

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataspacevertex_R_s_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_s_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_s_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_s_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


        H5::DataSpace dataspacevertex_R_t_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_t_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


        H5::DataSpace dataspacevertex_R_u_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_u_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_u_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_u_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


        H5::DataSpace dataspacevertex_sus_buffer(RANK_sus-1, sus_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);



        // Create the datasets in file:

        H5::DataSet* dataset_irred;
        dataset_irred = new H5::DataSet(file -> createDataSet(DATASET_irred,mtype1, dataspacevertex_irreducible,plist_vert));


        H5::DataSet* dataset_R_s;
        dataset_R_s = new H5::DataSet(file -> createDataSet(DATASET_R_s,mtype1, dataspacevertex_R_s,plist_vert));
        H5::DataSet* dataset_K1_s;
        dataset_K1_s = new H5::DataSet(file -> createDataSet(DATASET_K1_s,mtype1, dataspacevertex_K1_s,plist_vert));
        H5::DataSet* dataset_K2_s;
        dataset_K2_s = new H5::DataSet(file -> createDataSet(DATASET_K2_s,mtype1, dataspacevertex_K2_s,plist_vert));
#if sym==0
        H5::DataSet* dataset_K2b_s;
        dataset_K2b_s = new H5::DataSet(file -> createDataSet(DATASET_K2b_s,mtype1, dataspacevertex_K2b_s,plist_vert));
#endif


        H5::DataSet* dataset_R_t;
        dataset_R_t = new H5::DataSet(file -> createDataSet(DATASET_R_t,mtype1, dataspacevertex_R_t,plist_vert));
        H5::DataSet* dataset_K1_t;
        dataset_K1_t = new H5::DataSet(file -> createDataSet(DATASET_K1_t,mtype1, dataspacevertex_K1_t,plist_vert));
        H5::DataSet* dataset_K2_t;
        dataset_K2_t = new H5::DataSet(file -> createDataSet(DATASET_K2_t,mtype1, dataspacevertex_K2_t,plist_vert));
#if sym==0
        H5::DataSet* dataset_K2b_t;
        dataset_K2b_t = new H5::DataSet(file -> createDataSet(DATASET_K2b_t,mtype1, dataspacevertex_K2b_t,plist_vert));
#endif


        H5::DataSet* dataset_R_u;
        dataset_R_u = new H5::DataSet(file -> createDataSet(DATASET_R_u,mtype1, dataspacevertex_R_u,plist_vert));
        H5::DataSet* dataset_K1_u;
        dataset_K1_u = new H5::DataSet(file -> createDataSet(DATASET_K1_u,mtype1, dataspacevertex_K1_u,plist_vert));
        H5::DataSet* dataset_K2_u;
        dataset_K2_u = new H5::DataSet(file -> createDataSet(DATASET_K2_u,mtype1, dataspacevertex_K2_u,plist_vert));
#if sym==0
        H5::DataSet* dataset_K2b_u;
        dataset_K2b_u = new H5::DataSet(file -> createDataSet(DATASET_K2b_u,mtype1, dataspacevertex_K2b_u,plist_vert));
#endif


        H5::DataSet* dataset_sus;
        dataset_sus = new H5::DataSet(file -> createDataSet(DATASET_sus,H5::PredType::NATIVE_DOUBLE, dataspacevertex_sus));
        H5::DataSet* dataset_ferm_freqs;
        dataset_ferm_freqs = new H5::DataSet(file -> createDataSet(FERM_FREQS_LIST,H5::PredType::NATIVE_DOUBLE, dataspacefreqs));
        H5::DataSet* dataset_bos_freqs;
        dataset_bos_freqs = new H5::DataSet(file -> createDataSet(BOS_FREQS_LIST,H5::PredType::NATIVE_DOUBLE, dataspacefreqs_bos));
        H5::DataSet* dataset_self;
        dataset_self = new H5::DataSet(file -> createDataSet(SELF_LIST,H5::PredType::NATIVE_DOUBLE, dataspaceself));
        H5::DataSet* dataset_lambda;
        dataset_lambda = new H5::DataSet(file -> createDataSet(LAMBDA_LIST,H5::PredType::NATIVE_DOUBLE, dataspacelambda));

        //Select hyperslabs:

        //R_class

        //Select hyperslab in the file where the data should be located
        hsize_t start[2];
        hsize_t stride[2];
        hsize_t count[2];
        hsize_t block[2];

        start[0] = Lambda_it;
        start[1] = 0;
        for(int i=0; i<2;i++){
            stride[i] = 1;
            block[i] = 1;
        };
        count[0]= 1;


        count[1]= irred_dim1;

        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        count[1]= K3_dim1;

        dataspacevertex_R_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_R_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_R_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        count[1]= K1_dim1;

        dataspacevertex_K1_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        count[1]= K2_dim1;

        dataspacevertex_K2_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

#if sym==0
        dataspacevertex_K2b_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2b_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2b_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
#endif



        count[1]= nuc_eff*nuc_eff*3;

        dataspacevertex_sus.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        count[1]= nw1;
        dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);


        //Select hyperslab in buffer
//        hsize_t start_b[1];
//        hsize_t stride_b[1];
//        hsize_t count_b[1];
//        hsize_t block_b[1];

//        start_b[0] = 0;

//        stride_b[0] = 1;
//        block_b[0] = 1;



//        count_b[0]= irred_dim1;

    //    dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);


//        count_b[0]= K3_dim1;

//        dataspacevertex_R_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_R_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_R_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//        count_b[0]= K1_dim1;

//        dataspacevertex_K1_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K1_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K1_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//        count_b[0]= K2_dim1;

//        dataspacevertex_K2_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//#if sym==0
//        dataspacevertex_K2b_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2b_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2b_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//#endif

//        count_b[0]= nuc_eff*nuc_eff*3;
//        dataspacevertex_sus_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//        count_b[0]= nw1;

//        dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);




        //write buffer data into file.
        dataset_irred -> write(irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);

 cout << "passed irreeeed" << endl;
        dataset_R_s -> write( R_class_s, mtype1,dataspacevertex_R_s_buffer,dataspacevertex_R_s );
        dataset_K1_s -> write( K1_class_s, mtype1 ,dataspacevertex_K1_s_buffer,dataspacevertex_K1_s );
        dataset_K2_s -> write( K2_class_s, mtype1 ,dataspacevertex_K2_s_buffer,dataspacevertex_K2_s);
#if sym==0
        dataset_K2b_s -> write( K2b_class_s, mtype1 ,dataspacevertex_K2b_s_buffer,dataspacevertex_K2b_s);
#endif




        dataset_R_t -> write( R_class_t, mtype1,dataspacevertex_R_t_buffer,dataspacevertex_R_t );
        dataset_K1_t -> write( K1_class_t, mtype1 ,dataspacevertex_K1_t_buffer,dataspacevertex_K1_t );
        dataset_K2_t -> write( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#if sym==0
        dataset_K2b_t -> write( K2b_class_t, mtype1 ,dataspacevertex_K2b_t_buffer,dataspacevertex_K2b_t);
#endif



        dataset_R_u -> write( R_class_u, mtype1,dataspacevertex_R_u_buffer,dataspacevertex_R_u );
        dataset_K1_u -> write( K1_class_u, mtype1 ,dataspacevertex_K1_u_buffer,dataspacevertex_K1_u );
        dataset_K2_u -> write( K2_class_u, mtype1 ,dataspacevertex_K2_u_buffer,dataspacevertex_K2_u);
#if sym==0
        dataset_K2b_u -> write( K2b_class_u, mtype1 ,dataspacevertex_K2b_u_buffer,dataspacevertex_K2b_u);
#endif






        dataset_sus -> write( Suscept, H5::PredType::NATIVE_DOUBLE ,dataspacevertex_sus_buffer,dataspacevertex_sus );


        dataset_ferm_freqs -> write( ferm_freqs, H5::PredType::NATIVE_DOUBLE);
        dataset_bos_freqs -> write( bos_freqs, H5::PredType::NATIVE_DOUBLE);
        dataset_self -> write( selfenergy, H5::PredType::NATIVE_DOUBLE, dataspaceself_buffer,dataspaceself);
        dataset_lambda -> write(Lambda_list, H5::PredType::NATIVE_DOUBLE );


        cout << "Successfully saved in hdf5 file: " << FILE_NAME << endl;
        //free R



        delete[] K1_class_s;
        delete[] R_class_s;
        delete[] K2_class_s;
        delete[] K1_class_t;
        delete[] R_class_t;
        delete[] K2_class_t;
        delete[] K1_class_u;
        delete[] R_class_u;
        delete[] K2_class_u;
        delete[] selfenergy;
        delete[] irreducible_class;
#if sym==0
        delete[] K2b_class_s;
        delete[] K2b_class_t;
        delete[] K2b_class_u;

#endif



        //        delete dataset_R;
        //        delete dataset_K1;
        //        delete dataset_K2;
        //#if sym==0
        //        delete dataset_K2b;
        //#endif
        //        delete dataset_irred;
        //        delete dataset_sus;
        //        delete dataset_ferm_freqs;
        //        delete dataset_bos_freqs;
        //        delete dataset_self;
        //        delete dataset_lambda;



        //        delete file;
        dataset_bos_freqs->close();
        dataset_ferm_freqs->close();
        dataset_irred->close();

        dataset_K1_s->close();
        dataset_K2_s->close();

        dataset_K1_t->close();
        dataset_K2_t->close();

        dataset_K1_u->close();
        dataset_K2_u->close();
        dataset_lambda->close();
        dataset_R_s->close();
         dataset_R_t->close();
          dataset_R_u->close();
        dataset_self->close();
        dataset_sus->close();

#if sym==0
          dataset_K2b_s->close();
            dataset_K2b_t->close();
              dataset_K2b_u->close();
#endif
        file->close();
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
    if(Lambda_it < Lambda_size){
        try
        {
            typedef struct complex_t{
                double spin_re;   /*real part */
                double dens_re;   /*real part */

            }complex_t;

            //buffer self energy
            auto selfenergy = new double[nw];//irrdeucible vertex
            for(int i=0; i<nw; i++){
                selfenergy[i] = state_in.selfenergy.sval(i);
            };



            //buffer irreducible vertex:

            auto irreducible_class = new complex_t[irred_dim1];//irrdeucible vertex


            for(int i=0;i<irred_dim1;i++){

                irreducible_class[i].spin_re = state_in.vertex.spinvertex.irred.acc(i);
                irreducible_class[i].dens_re = state_in.vertex.densvertex.irred.acc(i);
            };



            cout << "passed irred" << endl;

            //buffer all R-class-arrays:


            auto R_class_s = new complex_t[K3_dim1];
            auto R_class_t = new complex_t[K3_dim1];
            auto R_class_u = new complex_t[K3_dim1];
            for(int i=0;i<K3_dim1;i++){
                R_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K3_acc(i);
                R_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K3_acc(i);

                R_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K3_acc(i);
                R_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K3_acc(i);

                R_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K3_acc(i);
                R_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K3_acc(i);
            };

            auto K1_class_s = new complex_t[K1_dim1];
            auto K1_class_t = new complex_t[K1_dim1];
            auto K1_class_u = new complex_t[K1_dim1];
            for(int i=0;i<K1_dim1;i++){
                K1_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K1_acc(i);
                K1_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K1_acc(i);

                K1_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K1_acc(i);
                K1_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K1_acc(i);

                K1_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K1_acc(i);
                K1_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K1_acc(i);
            };

            auto K2_class_s = new complex_t[K2_dim1];
            auto K2_class_t = new complex_t[K2_dim1];
            auto K2_class_u = new complex_t[K2_dim1];
            for(int i=0;i<K2_dim1;i++){
                K2_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K2_acc(i);
                K2_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K2_acc(i);

                K2_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K2_acc(i);
                K2_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K2_acc(i);

                K2_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K2_acc(i);
                K2_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K2_acc(i);
            };

    #if sym==0
            auto K2b_class_s = new complex_t[K2b_dim1];
            auto K2b_class_t = new complex_t[K2b_dim1];
            auto K2b_class_u = new complex_t[K2b_dim1];
            for(int i=0;i<K2_dim1;i++){
                K2b_class_s[i].spin_re=state_in.vertex.spinvertex.svertex.K2b_acc(i);
                K2b_class_s[i].dens_re=state_in.vertex.densvertex.svertex.K2b_acc(i);

                K2b_class_t[i].spin_re=state_in.vertex.spinvertex.tvertex.K2b_acc(i);
                K2b_class_t[i].dens_re=state_in.vertex.densvertex.tvertex.K2b_acc(i);

                K2b_class_u[i].spin_re=state_in.vertex.spinvertex.uvertex.K2b_acc(i);
                K2b_class_u[i].dens_re=state_in.vertex.densvertex.uvertex.K2b_acc(i);
            };
    #endif



            //        //susceptibility
            //        double Suscept[nuc_eff][nuc_eff][3];


            //
            //        for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
            //            for (int b=-(nuc_eff-1)/2; b<(nuc_eff-1)/2+1 ; b++){
            //                for (int c=1; c<4 ; c++){
            //                    if(distance(a,b,c) <= d_c){


            //                        //  Susceptibility in a-channel with bosonic transfer freq q
            //                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = state_in.sus.vval(a,b,c);

            //                    }
            //                    else{
            //                        Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] = 0;
            //                    };};};};



            cout << "Buffer ready. Preparing for saving into Hdf5 file.." << endl;
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            H5::Exception::dontPrint();

            // Open an existing file and dataset.
            //   H5::H5File file(FILE_NAME, H5F_ACC_RDWR);
            H5::H5File* file = 0;
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);

            H5::DataSet dataset_irred = file->openDataSet("irred");
            H5::DataSet dataset_R_s = file->openDataSet("R_s");
            H5::DataSet dataset_K1_s = file->openDataSet("K1_s");
            H5::DataSet dataset_K2_s = file->openDataSet("K2_s");
#if sym==0
            H5::DataSet dataset_K2b_s = file->openDataSet("K2b");
#endif


            H5::DataSet dataset_R_t = file->openDataSet("R_t");
            H5::DataSet dataset_K1_t = file->openDataSet("K1_t");
            H5::DataSet dataset_K2_t = file->openDataSet("K2_t");
#if sym==0
            H5::DataSet dataset_K2b_t = file->openDataSet("K2b_t");
#endif


            H5::DataSet dataset_R_u = file->openDataSet("R_u");
            H5::DataSet dataset_K1_u = file->openDataSet("K1_u");
            H5::DataSet dataset_K2_u = file->openDataSet("K2_u");
#if sym==0
            H5::DataSet dataset_K2b_u = file->openDataSet("K2b_u");
#endif
            // H5::DataSet dataset_sus = file->openDataSet("sus");
            H5::DataSet dataset_self = file->openDataSet("selflist");

            //Create the memory datatype
            H5::CompType mtype1( sizeof(complex_t) );
            mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
            //   mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
            mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
            //  mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);


            complex_t fillvalue_vert;
            fillvalue_vert.spin_re = 0;
            fillvalue_vert.dens_re = 0;
            H5::DSetCreatPropList plist_vert;
            plist_vert.setFillValue(mtype1, &fillvalue_vert);





            // Create the dimension arrays for objects in file.
            hsize_t R_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K1_dim1)};   // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    #if sym==0
            hsize_t K2b_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    #endif
            hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

            hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(nw)};

            // Create the dimension arrays for objects in buffer.
            hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t R_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
               hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
    #if sym==0
            hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

    #endif
            hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};



            // Create the data space for the dataset in file.

            H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)


            H5::DataSpace dataspacevertex_R_s(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_s(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
            H5::DataSpace dataspacevertex_K2_s(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_s(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
    #endif


            H5::DataSpace dataspacevertex_R_t(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_t(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
            H5::DataSpace dataspacevertex_K2_t(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_t(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
    #endif


            H5::DataSpace dataspacevertex_R_u(RANK_R, R_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_u(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
            H5::DataSpace dataspacevertex_K2_u(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_u(RANK_K2b,K2b_dims);//data space for vertex with three dimensions (two independent frequencies)
    #endif





            H5::DataSpace dataspaceself(RANK_self, self_dims);


            // Create the data space for buffer objects.
            H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

            H5::DataSpace dataspacevertex_R_s_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_s_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K2_s_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_s_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #endif


            H5::DataSpace dataspacevertex_R_t_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_t_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #endif


            H5::DataSpace dataspacevertex_R_u_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_u_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K2_u_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #if sym==0
            H5::DataSpace dataspacevertex_K2b_u_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
    #endif


               H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);


            //Select hyperslabs:

            //R_class

            //Select hyperslab in the file where the data should be located
            hsize_t start[2];
            hsize_t stride[2];
            hsize_t count[2];
            hsize_t block[2];

            start[0] = Lambda_it;
            start[1] = 0;
            for(int i=0; i<2;i++){
                stride[i] = 1;
                block[i] = 1;
            };
            count[0]= 1;


            count[1]= irred_dim1;

            dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);



            count[1]= K3_dim1;

            dataspacevertex_R_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_R_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_R_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            count[1]= K1_dim1;

            dataspacevertex_K1_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K1_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            count[1]= K2_dim1;

            dataspacevertex_K2_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

    #if sym==0
            dataspacevertex_K2b_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2b_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2b_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
    #endif




             count[1]= nw1;
            dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            //Select hyperslab in buffer
//            hsize_t start_b[1];
//            hsize_t stride_b[1];
//            hsize_t count_b[1];
//            hsize_t block_b[1];

//            start_b[0] = 0;

//            stride_b[0] = 1;
//            block_b[0] = 1;



//            count_b[0]= irred_dim1;

//            dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);



//            count_b[0]= K3_dim1;

//            dataspacevertex_R_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_R_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_R_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//            count_b[0]= K1_dim1;

//            dataspacevertex_K1_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K1_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K1_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//            count_b[0]= K2_dim1;

//            dataspacevertex_K2_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K2_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K2_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//    #if sym==0
//            dataspacevertex_K2b_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K2b_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//            dataspacevertex_K2b_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//    #endif


//            count_b[0]= nw1;

//            dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);






            //write buffer data into file.


            dataset_R_s.write( R_class_s, mtype1,dataspacevertex_R_s_buffer,dataspacevertex_R_s );
            dataset_K1_s.write( K1_class_s, mtype1 ,dataspacevertex_K1_s_buffer,dataspacevertex_K1_s );
            dataset_K2_s.write( K2_class_s, mtype1 ,dataspacevertex_K2_s_buffer,dataspacevertex_K2_s);
#if sym==0
            dataset_K2b_s.write( K2b_class_s, mtype1 ,dataspacevertex_K2b_s_buffer,dataspacevertex_K2b_s);
#endif


            dataset_R_t.write( R_class_t, mtype1,dataspacevertex_R_t_buffer,dataspacevertex_R_t );
            dataset_K1_t.write( K1_class_t, mtype1 ,dataspacevertex_K1_t_buffer,dataspacevertex_K1_t );
            dataset_K2_t.write( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#if sym==0
            dataset_K2b_t.write( K2b_class_t, mtype1 ,dataspacevertex_K2b_t_buffer,dataspacevertex_K2b_t);
#endif


            dataset_R_u.write( R_class_u, mtype1,dataspacevertex_R_u_buffer,dataspacevertex_R_u );
            dataset_K1_u.write( K1_class_u, mtype1 ,dataspacevertex_K1_u_buffer,dataspacevertex_K1_u );
            dataset_K2_u.write( K2_class_u, mtype1 ,dataspacevertex_K2_u_buffer,dataspacevertex_K2_u);
#if sym==0
            dataset_K2b_u.write( K2b_class_u, mtype1 ,dataspacevertex_K2b_u_buffer,dataspacevertex_K2b_u);
#endif

            dataset_irred.write( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);

              dataset_self.write( selfenergy, H5::PredType::NATIVE_DOUBLE, dataspaceself_buffer,dataspaceself);

            delete[] K1_class_s;
            delete[] R_class_s;
            delete[] K2_class_s;
            delete[] K1_class_t;
            delete[] R_class_t;
            delete[] K2_class_t;
            delete[] K1_class_u;
            delete[] R_class_u;
            delete[] K2_class_u;
            delete[] selfenergy;
            delete[] irreducible_class;
#if sym==0
            delete[] K2b_class_s;
             delete[] K2b_class_t;
             delete[] K2b_class_u;
#endif


            dataset_irred.close();
            dataset_R_s.close();
            dataset_K1_s.close();
            dataset_K2_s.close();
#if sym==0
            dataset_K2b_t.close();
#endif
            dataset_R_t.close();
            dataset_K1_t.close();
            dataset_K2_t.close();
#if sym==0
            dataset_K2b_t.close();
#endif
            dataset_R_u.close();
            dataset_K1_u.close();
            dataset_K2_u.close();
#if sym==0
            dataset_K2b_u.close();
#endif

            dataset_self.close();
            file->close();
            delete file;

            cout << "Successfully saved in hdf5 file: " << FILE_NAME << "  in Lambda-layer " << Lambda_it << endl;
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
    else{
        cout << "Cannot write to file " << FILE_NAME << " since Lambda layer is out of range." << endl;
    };
}

void add_hdf_sus(const H5std_string FILE_NAME,int iteration, vector<double>& Lambdas,long total_iterations,Susc& susceptibility){


    cout << "Starting to copy susceptibility to buffer.." << endl;
    if(iteration < total_iterations){
        try
        {

            //  H5::H5File file(FILE_NAME, H5F_ACC_RDWR);//open file
            H5::H5File* file = 0;
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);

            double Lambda_list[total_iterations];
            for(int i=0; i<total_iterations;i++){
                Lambda_list[i] = Lambdas[i];};

            H5::DataSet dataset_lambda = file->openDataSet("lambdas");

            dataset_lambda.write( Lambda_list, H5::PredType::NATIVE_DOUBLE  );//overwrite vector containing all values for lambda


            //susceptibility
            auto Suscept =new double[nuc_eff*nuc_eff*3];


int a=-(nuc_eff-1)/2;
int b=-(nuc_eff-1)/2;
int c=0;

for(int i=0;i<nuc_eff*nuc_eff*3;i++){
    c = c%3+1;

    Suscept[i] = susceptibility.vval(a,b,c);
    if((i+1)%3==0){b+=1;};
    if(b==(nuc_eff-1)/2+1){
        b= -(nuc_eff-1)/2;
        a += 1;
    };
};


            H5::DataSet dataset_sus = file->openDataSet("sus");
            hsize_t sus_dims[]= {static_cast<hsize_t>(total_iterations),static_cast<hsize_t>(nuc_eff*nuc_eff*3)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t sus_dims_buffer[]= {static_cast<hsize_t>(nuc_eff*nuc_eff*3)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

            H5::DataSpace dataspacevertex_sus(RANK_sus,sus_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_sus_buffer(RANK_sus-1, sus_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)






            //susceptibility

            //Select hyperslab in the file where the data should be located
            hsize_t sus_start[2];
            hsize_t sus_stride[2];
            hsize_t sus_count[2];
            hsize_t sus_block[2];

            sus_start[0] = iteration;
            sus_start[1] = 0;
            for(int i=0; i<2;i++){
                sus_stride[i] = 1;
                sus_block[i] = 1;
            };
            sus_count[0]= 1;
            sus_count[1]= nuc_eff*nuc_eff*3;

            dataspacevertex_sus.selectHyperslab(H5S_SELECT_SET, sus_count,sus_start,sus_stride,sus_block);
            //Select hyperslab in buffer
//            hsize_t sus_start_b[1];
//            hsize_t sus_stride_b[1];
//            hsize_t sus_count_b[1];
//            hsize_t sus_block_b[1];


//            sus_start_b[0] = 0;
//            sus_stride_b[0] = 1;
//            sus_block_b[0] = 1;


//            sus_count_b[0]= nuc_eff*nuc_eff*3;


 //           dataspacevertex_sus_buffer.selectHyperslab(H5S_SELECT_SET, sus_count_b,sus_start_b,sus_stride_b,sus_block_b);

            dataset_sus.write( Suscept, H5::PredType::NATIVE_DOUBLE ,dataspacevertex_sus_buffer,dataspacevertex_sus );

            delete[] Suscept;
            dataset_lambda.close();
            dataset_sus.close();
            file->close();
            delete file;



            cout << "Successfully saved susceptibility in hdf5 file: " << FILE_NAME << endl;
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
    else{
        cout << "Cannot write susceptibility to file " << FILE_NAME << " since iterator is out of range." << endl;
    };
}









//function to read out an exstiting Hdf5 file. Needed to resume computation if it has been interrupted during the flow
state read_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size,long total_iterations, vector<double> &Lambdas){
    state result;
    if(Lambda_it<Lambda_size){

cout << "OOOOOK" << endl;

        typedef struct complex_t{
            double spin_re;   /*real part */
            double dens_re;   /*real part */
        }complex_t;


        //buffer self energy

        auto irreducible_class = new complex_t[irred_dim1];//irrdeucible vertex
        auto R_class_s = new complex_t[K3_dim1];
        auto K1_class_s = new complex_t[K1_dim1];
        auto K2_class_s = new complex_t[K2_dim1];
        auto R_class_t = new complex_t[K3_dim1];
        auto K1_class_t = new complex_t[K1_dim1];
        auto K2_class_t = new complex_t[K2_dim1];
        auto R_class_u = new complex_t[K3_dim1];
        auto K1_class_u = new complex_t[K1_dim1];
        auto K2_class_u = new complex_t[K2_dim1];
        auto selfenergy = new double[nw];
        auto Lambdas_buff = new double[total_iterations];

#if sym==0
        auto K2b_class_s = new complex_t[K2_dim1];
          auto K2b_class_t = new complex_t[K2_dim1];
            auto K2b_class_u = new complex_t[K2_dim1];
#endif
        //  double Suscept[nuc_eff][nuc_eff][3];





        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);

        H5::H5File* file = 0;
        file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

        // H5::H5File file(FILE_NAME, H5F_ACC_RDONLY );


        H5::DataSet dataset_irred = file->openDataSet("irred");
        H5::DataSet dataset_R_s = file->openDataSet("R_s");
        H5::DataSet dataset_K1_s = file->openDataSet("K1_s");
        H5::DataSet dataset_K2_s = file->openDataSet("K2_s");
#if sym==0
        H5::DataSet dataset_K2b_s = file.openDataSet("K2b_s");
#endif
        H5::DataSet dataset_R_t = file->openDataSet("R_t");
        H5::DataSet dataset_K1_t = file->openDataSet("K1_t");
        H5::DataSet dataset_K2_t = file->openDataSet("K2_t");
#if sym==0
        H5::DataSet dataset_K2b_t = file.openDataSet("K2b_t");
#endif
        H5::DataSet dataset_R_u = file->openDataSet("R_u");
        H5::DataSet dataset_K1_u = file->openDataSet("K1_u");
        H5::DataSet dataset_K2_u = file->openDataSet("K2_u");
#if sym==0
        H5::DataSet dataset_K2b_u = file.openDataSet("K2b_u");
#endif

        H5::DataSet dataset_self = file->openDataSet("selflist");
        H5::DataSet dataset_lambdas = file->openDataSet("lambdas");




        H5::DataSpace dataspacevertex_R_s=dataset_R_s.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_s=dataset_K1_s.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_s=dataset_K2_s.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_s=dataset_K2b.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#endif

        H5::DataSpace dataspacevertex_R_t=dataset_R_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t=dataset_K1_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t=dataset_K2_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_t=dataset_K2b_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#endif

        H5::DataSpace dataspacevertex_R_u=dataset_R_u.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_u=dataset_K1_u.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_u=dataset_K2_u.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_u=dataset_K2b_u.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#endif
        H5::DataSpace dataspacevertex_irreducible=dataset_irred.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        //   H5::DataSpace dataspacevertex_sus= dataset_sus.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspaceself=dataset_self.getSpace();
        H5::DataSpace dataspacelambdas=dataset_lambdas.getSpace();







        // Create the dimension arrays for objects in buffer.
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t R_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#if sym==0
        hsize_t K2b_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif
        hsize_t self_dims_buffer[]={static_cast<hsize_t>(nw)};



        // Create the data space for buffer objects.
        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataspacevertex_R_s_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_s_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_s_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_s_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


        H5::DataSpace dataspacevertex_R_t_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_t_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


        H5::DataSpace dataspacevertex_R_u_buffer(RANK_R-1, R_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_u_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_u_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#if sym==0
        H5::DataSpace dataspacevertex_K2b_u_buffer(RANK_K2-1, K2b_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif


         H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);


        //Select hyperslabs:

        //R_class

        //Select hyperslab in the file where the data should be located
        hsize_t start[2];
        hsize_t stride[2];
        hsize_t count[2];
        hsize_t block[2];

        start[0] = Lambda_it;
        start[1] = 0;
        for(int i=0; i<2;i++){
            stride[i] = 1;
            block[i] = 1;
        };
        count[0]= 1;


        count[1]= irred_dim1;

        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);



        count[1]= K3_dim1;

        dataspacevertex_R_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_R_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_R_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        count[1]= K1_dim1;

        dataspacevertex_K1_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        count[1]= K2_dim1;

        dataspacevertex_K2_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

#if sym==0
        dataspacevertex_K2b_s.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2b_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2b_u.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
#endif



          count[1]= nw1;
        dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);


        //Select hyperslab in buffer
//        hsize_t start_b[1];
//        hsize_t stride_b[1];
//        hsize_t count_b[1];
//        hsize_t block_b[1];

//        start_b[0] = 0;

//        stride_b[0] = 1;
//        block_b[0] = 1;



//        count_b[0]= irred_dim1;

//        dataspacevertex_irreducible_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);



//        count_b[0]= K3_dim1;

//        dataspacevertex_R_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_R_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_R_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//        count_b[0]= K1_dim1;

//        dataspacevertex_K1_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K1_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K1_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//        count_b[0]= K2_dim1;

//        dataspacevertex_K2_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);

//#if sym==0
//        dataspacevertex_K2b_s_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2b_t_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//        dataspacevertex_K2b_u_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);
//#endif


//        count_b[0]= nw1;

//        dataspaceself_buffer.selectHyperslab(H5S_SELECT_SET, count_b,start_b,stride_b,block_b);



        dataset_R_s.read( R_class_s, mtype1,dataspacevertex_R_s_buffer,dataspacevertex_R_s );
        dataset_K1_s.read( K1_class_s, mtype1 ,dataspacevertex_K1_s_buffer,dataspacevertex_K1_s );
        dataset_K2_s.read( K2_class_s, mtype1 ,dataspacevertex_K2_s_buffer,dataspacevertex_K2_s);
#if sym==0
        dataset_K2b_s.read( K2b_class, mtype1 ,dataspacevertex_K2b_s_buffer,dataspacevertex_K2b_s);
#endif

        dataset_R_t.read( R_class_t, mtype1,dataspacevertex_R_t_buffer,dataspacevertex_R_t );
        dataset_K1_t.read( K1_class_t, mtype1 ,dataspacevertex_K1_t_buffer,dataspacevertex_K1_t );
        dataset_K2_t.read( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#if sym==0
        dataset_K2b_t.read( K2b_class, mtype1 ,dataspacevertex_K2b_t_buffer,dataspacevertex_K2b_t);
#endif

        dataset_R_u.read( R_class_u, mtype1,dataspacevertex_R_u_buffer,dataspacevertex_R_u );
        dataset_K1_u.read( K1_class_u, mtype1 ,dataspacevertex_K1_u_buffer,dataspacevertex_K1_u );
        dataset_K2_u.read( K2_class_u, mtype1 ,dataspacevertex_K2_u_buffer,dataspacevertex_K2_u);
#if sym==0
        dataset_K2b_u.read( K2b_class, mtype1 ,dataspacevertex_K2b_u_buffer,dataspacevertex_K2b_u);
#endif
        dataset_irred.read( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);

        dataset_self.read( selfenergy, H5::PredType::NATIVE_DOUBLE, dataspaceself_buffer,dataspaceself);
        dataset_lambdas.read(Lambdas_buff, H5::PredType::NATIVE_DOUBLE,dataspacelambdas);

        for(int i=0; i<total_iterations;i++){
            Lambdas[i] = Lambdas_buff[i];};

        for(int i=0; i<nw; i++){
            result.selfenergy.setself(i,selfenergy[i]);
        };



      for(int i=0; i<irred_dim1;i++){
            result.vertex.spinvertex.irred.direct_set(i,irreducible_class[i].spin_re);
            result.vertex.densvertex.irred.direct_set(i,irreducible_class[i].dens_re);
};



        //Note that the first index labels the channel with the correspondance ( 0 = s-channel, 1 = t-channel, 2 = u-channel)


        //buffer all R-class-arrays:



      for (int i=0; i<K3_dim1; i++){


          //rest class from s-channel
          result.vertex.spinvertex.svertex.R_direct_set(i,R_class_s[i].spin_re );
          result.vertex.densvertex.svertex.R_direct_set(i, R_class_s[i].dens_re );

          //rest class from t-channel
          result.vertex.spinvertex.tvertex.R_direct_set(i, R_class_t[i].spin_re );
          result.vertex.densvertex.tvertex.R_direct_set(i, R_class_t[i].dens_re );


          //rest class from u-channel
          result.vertex.spinvertex.uvertex.R_direct_set(i, R_class_t[i].spin_re );
          result.vertex.densvertex.uvertex.R_direct_set(i, R_class_t[i].dens_re );

      };






for (int i=0; i<K1_dim1; i++){

                                //  K1 class from s-channel
                                result.vertex.spinvertex.svertex.K1_direct_set(i,K1_class_s[i].spin_re );
                                result.vertex.densvertex.svertex.K1_direct_set(i,K1_class_s[i].dens_re );
                                //K1 class from t-channel
                                result.vertex.spinvertex.tvertex.K1_direct_set(i,K1_class_t[i].spin_re );
                                result.vertex.densvertex.tvertex.K1_direct_set(i,K1_class_t[i].dens_re );

                                //K1 class from u-channel
                                result.vertex.spinvertex.uvertex.K1_direct_set(i,K1_class_u[i].spin_re );
                                result.vertex.densvertex.uvertex.K1_direct_set(i,K1_class_u[i].dens_re );

};
        //buffer all K2-class-arrays:


for (int i=0; i<K2_dim1; i++){

                                //  K2 class from s-channel
                                result.vertex.spinvertex.svertex.K2_direct_set(i,K2_class_s[i].spin_re );
                                result.vertex.densvertex.svertex.K2_direct_set(i,K2_class_s[i].dens_re );
                                //K2 class from t-channel
                                result.vertex.spinvertex.tvertex.K2_direct_set(i,K2_class_t[i].spin_re );
                                result.vertex.densvertex.tvertex.K2_direct_set(i,K2_class_t[i].dens_re );

                                //K2 class from u-channel
                                result.vertex.spinvertex.uvertex.K2_direct_set(i,K2_class_u[i].spin_re );
                                result.vertex.densvertex.uvertex.K2_direct_set(i,K2_class_u[i].dens_re );
};

#if sym==0

        //buffer all K2b-class-arrays:
for (int i=0; i<K2_dim1; i++){

                                //  K2b class from s-channel
                                result.vertex.spinvertex.svertex.K2b_direct_set(i,K2b_class_s[i].spin_re );
                                result.vertex.densvertex.svertex.K2b_direct_set(i,K2b_class_s[i].dens_re );
                                //K2b class from t-channel
                                result.vertex.spinvertex.tvertex.K2b_direct_set(i,K2b_class_t[i].spin_re );
                                result.vertex.densvertex.tvertex.K2b_direct_set(i,K2b_class_t[i].dens_re );

                                //K2b class from u-channel
                                result.vertex.spinvertex.uvertex.K2b_direct_set(i,K2b_class_u[i].spin_re );
                                result.vertex.densvertex.uvertex.K2b_direct_set(i,K2b_class_u[i].dens_re );
};


#endif


        //susceptibility



        //    for (int a=-(nuc_eff-1)/2; a<(nuc_eff-1)/2+1 ; a++){
        //        for (int b=0; b<(nuc_eff-1)/2+1 ; b++){
        //            for (int c=1; c<4 ; c++){
        //                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) ){
        //                    if(distance(a,b,c) <= d_c){


        //                        //  Susceptibility in a-channel with bosonic transfer freq q
        //                        result.sus.write(a,b,c,Suscept[a+(nuc_eff-1)/2][b+(nuc_eff-1)/2][c-1] );

        //                    }
        //                    else {
        //                        result.sus.write(a,b,c,0);
        //                    };};};};};


        delete[] K1_class_s;
        delete[] R_class_s;
        delete[] K2_class_s;
delete[] K1_class_t;
delete[] R_class_t;
delete[] K2_class_t;
delete[] K1_class_u;
delete[] R_class_u;
delete[] K2_class_u;
        delete[] selfenergy;
        delete[] irreducible_class;

#if sym==0
        delete[] K2b_class_s;
 delete[] K2b_class_t;
 delete[] K2b_class_u;
#endif

        dataset_irred.close();
        dataset_K1_s.close();
        dataset_K2_s.close();
        dataset_K1_t.close();
        dataset_K2_t.close();
        dataset_K1_u.close();
        dataset_K2_u.close();
        dataset_lambdas.close();
        dataset_R_s.close();
         dataset_R_t.close();
          dataset_R_u.close();
        dataset_self.close();
        file->close();
        delete file;


        return result;

    }
    else{
        cout << "Cannot read from file " << FILE_NAME << " since Lambda layer out of range" << endl;
    };

}
