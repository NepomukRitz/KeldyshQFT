
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
#include "H5Cpp.h"

#include "parameters.h"
//#include "vertex.h"
//#include "state.h"
//#include "loop.h"
//#include "a_bubble.h"
//#include "p_bubble.h"
//#include "t_bubble.h"
//#include "propagator.h"
//#include "selfenergy.h"

/********************************constants concerning HDF5 data format*************************/
/*const int    NX = nw;                     // dataset dimensions
const int    NY = nw;
const int    NZ = nw;*/                       // dataset dimensions
const int    RANK_K3 = 2;
const int    RANK_K1 = 2;
const int    RANK_K2 = 2;
const int    RANK_K2b = 2;
const int    RANK_irreducible = 2;
const int    RANK_sus = 2;
const int    RANK_self = 2;



const H5std_string	DATASET_irred("irred");

const H5std_string	DATASET_K1_a("K1_a");
const H5std_string	DATASET_K1_p("K1_p");
const H5std_string	DATASET_K1_t("K1_t");

const H5std_string	DATASET_K2_a("K2_a");
const H5std_string	DATASET_K2_p("K2_p");
const H5std_string	DATASET_K2_t("K2_t");
const H5std_string	DATASET_K2b_a("K2b_a");
const H5std_string	DATASET_K2b_p("K2b_p");
const H5std_string	DATASET_K2b_t("K2b_t");

const H5std_string	DATASET_K3_a("K3_a");
const H5std_string	DATASET_K3_p("K3_p");
const H5std_string	DATASET_K3_t("K3_t");


const H5std_string	DATASET_sus("sus");
const H5std_string	FERM_FREQS_LIST("ferm_freqslist");
const H5std_string	BOS_FREQS_LIST("bos_freqslist");
const H5std_string	SELF_LIST("selflist");
const H5std_string	LAMBDA_LIST("lambdas");
const H5std_string  PARAM_LIST("parameters");
const H5std_string  MEMBER1( "spin_re" );
const H5std_string  MEMBER2( "spin_im" );
const H5std_string  MEMBER3( "dens_re" );
const H5std_string  MEMBER4( "dens_im" );
const H5std_string  MEMBER5( "re" );
const H5std_string  MEMBER6( "im" );
//function that write inital state into hdf5 format.  The second argument denotes the iteration number to which it is written. The third argument denotes the total number of saved iterations in the file

void write_hdf(const H5std_string FILE_NAME,double Lambda_i, long Lambda_size,State<complex<double>>& state_in){
    //Try block to detect exceptions raised by any of the calls inside it
    cout << "Starting to copy to buffer.." << endl;
    try
    {
        int Lambda_it=0;
        typedef struct complex_t{
            double spin_re;   /*real part */
            double spin_im;   /*imaginary part */
            double dens_re;   /*real part */
            double dens_im;   /*imaginary part */
        }complex_t;


        typedef struct complex{
            double re;   /*real part */
            double im;   /*imaginary part */

        }complex;


        //buffer self energy
        int self_dim =2*nSE;
        auto selfenergy = new complex[self_dim];//irrdeucible vertex
        for(int i=0; i<self_dim; i++){
            selfenergy[i].re = real(state_in.selfenergy.acc(i));
            selfenergy[i].im = imag(state_in.selfenergy.acc(i));
        };




        //buffer irreducible vertex:
        int irred_dim1 =16*n_in;
        auto irreducible_class = new complex_t[irred_dim1];//irrdeucible vertex


        for(int i=0;i<irred_dim1;i++){

            irreducible_class[i].spin_re = real(state_in.vertex.spinvertex.irred.acc(i));
            irreducible_class[i].dens_re = real(state_in.vertex.densvertex.irred.acc(i));
            irreducible_class[i].spin_im = imag(state_in.vertex.spinvertex.irred.acc(i));
            irreducible_class[i].dens_im = imag(state_in.vertex.densvertex.irred.acc(i));
        };



#if DIAG_CLASS >=1
        int K1_dim1=nK_K1 * nw1_wt * n_in;

        auto K1_class_p = new complex_t[K1_dim1];
        auto K1_class_t = new complex_t[K1_dim1];
        auto K1_class_a = new complex_t[K1_dim1];
        for(int i=0;i<K1_dim1;i++){
            K1_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K1_acc(i));
            K1_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K1_acc(i));
            K1_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K1_acc(i));
            K1_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K1_acc(i));

            K1_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K1_acc(i));
            K1_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K1_acc(i));
            K1_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K1_acc(i));
            K1_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K1_acc(i));

            K1_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K1_acc(i));
            K1_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K1_acc(i));
            K1_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K1_acc(i));
            K1_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K1_acc(i));
        }
#endif
#if DIAG_CLASS >=2
        int K2_dim1=nK_K2 * nw2_wt * nw2_nut * n_in;
        auto K2_class_p = new complex_t[K2_dim1];
        auto K2_class_t = new complex_t[K2_dim1];
        auto K2_class_a = new complex_t[K2_dim1];
        for(int i=0;i<K2_dim1;i++){
            K2_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K2_acc(i));
            K2_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K2_acc(i));
            K2_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K2_acc(i));
            K2_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K2_acc(i));

            K2_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K2_acc(i));
            K2_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K2_acc(i));
            K2_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K2_acc(i));
            K2_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K2_acc(i));

            K2_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K2_acc(i));
            K2_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K2_acc(i));
            K2_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K2_acc(i));
            K2_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K2_acc(i));
        }
#endif
#if DIAG_CLASS >=3
        //buffer all K3-class-arrays:
        int K3_dim1=nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in;

        auto K3_class_a = new complex_t[K3_dim1];
        auto K3_class_p = new complex_t[K3_dim1];
        auto K3_class_t = new complex_t[K3_dim1];
        for(int i=0;i<K3_dim1;i++){
            K3_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K3_acc(i));
            K3_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K3_acc(i));
            K3_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K3_acc(i));
            K3_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K3_acc(i));

            K3_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K3_acc(i));
            K3_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K3_acc(i));
            K3_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K3_acc(i));
            K3_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K3_acc(i));

            K3_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K3_acc(i));
            K3_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K3_acc(i));
            K3_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K3_acc(i));
            K3_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K3_acc(i));
        }
#endif



        double Lambda_list[Lambda_size];

        Lambda_list[0] = Lambda_i;

        for(int i=1; i<Lambda_size; i++){
            Lambda_list[i]=0;};

        cout << "Buffer ready. Preparing for saving into Hdf5 file..." << endl;
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        //  H5::Exception::dontPrint();

        // Create a new file using the default property lists.
        H5::H5File* file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);


        //Create the memory datatype
        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER3, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);


        H5::CompType mtype2( sizeof(complex) );
        mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
        mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);


        complex_t fillvalue_vert;
        fillvalue_vert.spin_re = 0;
        fillvalue_vert.dens_re = 0;
        fillvalue_vert.spin_im = 0;
        fillvalue_vert.dens_im = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype1, &fillvalue_vert);

        complex fillvalue_self;
        fillvalue_self.re = 0;
        fillvalue_self.im = 0;

        H5::DSetCreatPropList plist_self;
        plist_self.setFillValue(mtype2, &fillvalue_self);

        // Create the dimension arrays for objects in file and in buffer
        hsize_t Lambda_dims[]={static_cast<hsize_t>(Lambda_size)};
        hsize_t param_dims[]={static_cast<hsize_t>(param_size)};
        hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(self_dim)};
        hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        hsize_t self_dims_buffer[]={static_cast<hsize_t>(self_dim)};
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#if DIAG_CLASS >=1
        hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K1_dim1)};   // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif

#if DIAG_CLASS >=2
        hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif

#if DIAG_CLASS >=3
        hsize_t K3_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
        hsize_t K3_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif


        // Create the data space for the dataset in file and for buffer objects
        H5::DataSpace dataspacelambda(1, Lambda_dims);
        H5::DataSpace dataspaceparams(1, param_dims);
        H5::DataSpace dataspaceself(RANK_self, self_dims);
        H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        // Create the datasets in file:
        H5::DataSet* dataset_lambda;
        dataset_lambda = new H5::DataSet(file -> createDataSet(LAMBDA_LIST,H5::PredType::NATIVE_DOUBLE, dataspacelambda));

        H5::DataSet* dataset_params;
        dataset_params = new H5::DataSet(file -> createDataSet(PARAM_LIST,H5::PredType::NATIVE_DOUBLE, dataspaceparams));

        H5::DataSet* dataset_self;
        dataset_self = new H5::DataSet(file -> createDataSet(SELF_LIST,mtype2, dataspaceself));

        H5::DataSet* dataset_irred;
        dataset_irred = new H5::DataSet(file -> createDataSet(DATASET_irred,mtype1, dataspacevertex_irreducible,plist_vert));

#if DIAG_CLASS >=1
        H5::DataSpace dataspacevertex_K1_a(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K1_p(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataspacevertex_K1_t(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)

        H5::DataSpace dataspacevertex_K1_a_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_p_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        // Create the datasets in file:
        H5::DataSet* dataset_K1_a;
        dataset_K1_a = new H5::DataSet(file -> createDataSet(DATASET_K1_a,mtype1, dataspacevertex_K1_a,plist_vert));
        H5::DataSet* dataset_K1_p;
        dataset_K1_p = new H5::DataSet(file -> createDataSet(DATASET_K1_p,mtype1, dataspacevertex_K1_p,plist_vert));
        H5::DataSet* dataset_K1_t;
        dataset_K1_t = new H5::DataSet(file -> createDataSet(DATASET_K1_t,mtype1, dataspacevertex_K1_t,plist_vert));
#endif
#if DIAG_CLASS >=2
        H5::DataSpace dataspacevertex_K2_a(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
        H5::DataSpace dataspacevertex_K2_p(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
        H5::DataSpace dataspacevertex_K2_t(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)

        H5::DataSpace dataspacevertex_K2_a_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_p_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        // Create the datasets in file:
        H5::DataSet* dataset_K2_a;
        dataset_K2_a = new H5::DataSet(file -> createDataSet(DATASET_K2_a,mtype1, dataspacevertex_K2_a,plist_vert));
        H5::DataSet* dataset_K2_p;
        dataset_K2_p = new H5::DataSet(file -> createDataSet(DATASET_K2_p,mtype1, dataspacevertex_K2_p,plist_vert));
        H5::DataSet* dataset_K2_t;
        dataset_K2_t = new H5::DataSet(file -> createDataSet(DATASET_K2_t,mtype1, dataspacevertex_K2_t,plist_vert));
#endif

#if DIAG_CLASS >=3
        H5::DataSpace dataspacevertex_K3_a(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_p(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_t(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataspacevertex_K3_a_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_p_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_t_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

        // Create the datasets in file:
        H5::DataSet* dataset_K3_a;
        dataset_K3_a = new H5::DataSet(file -> createDataSet(DATASET_K3_a,mtype1, dataspacevertex_K3_a,plist_vert));
        H5::DataSet* dataset_K3_p;
        dataset_K3_p = new H5::DataSet(file -> createDataSet(DATASET_K3_p,mtype1, dataspacevertex_K3_p,plist_vert));
        H5::DataSet* dataset_K3_t;
        dataset_K3_t = new H5::DataSet(file -> createDataSet(DATASET_K3_t,mtype1, dataspacevertex_K3_t,plist_vert));
#endif




        //Select hyperslabs:

        //Select hyperslab in the file where the data should be located and after write buffer data into file.
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


        dataset_lambda -> write(Lambda_list, H5::PredType::NATIVE_DOUBLE );
        dataset_params -> write(parameter_list, H5::PredType::NATIVE_DOUBLE);

        count[1]= self_dim;
        dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataset_self -> write( selfenergy,mtype2, dataspaceself_buffer,dataspaceself);

        count[1]= irred_dim1;
        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataset_irred -> write(irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);

#if DIAG_CLASS >=1
        count[1]= K1_dim1;

        dataspacevertex_K1_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        dataset_K1_a -> write( K1_class_a, mtype1 ,dataspacevertex_K1_a_buffer,dataspacevertex_K1_a );
        dataset_K1_p -> write( K1_class_p, mtype1 ,dataspacevertex_K1_p_buffer,dataspacevertex_K1_p );
        dataset_K1_t -> write( K1_class_t, mtype1 ,dataspacevertex_K1_t_buffer,dataspacevertex_K1_t );
#endif

#if DIAG_CLASS >=2
        count[1]= K2_dim1;

        dataspacevertex_K2_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        dataset_K2_a -> write( K2_class_a, mtype1 ,dataspacevertex_K2_a_buffer,dataspacevertex_K2_a);
        dataset_K2_p -> write( K2_class_p, mtype1 ,dataspacevertex_K2_p_buffer,dataspacevertex_K2_p);
        dataset_K2_t -> write( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#endif

#if DIAG_CLASS >=3
        count[1]= K3_dim1;

        dataspacevertex_K3_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K3_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataspacevertex_K3_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        dataset_K3_a -> write( K3_class_a, mtype1,dataspacevertex_K3_a_buffer,dataspacevertex_K3_a );
        dataset_K3_p -> write( K3_class_p, mtype1,dataspacevertex_K3_p_buffer,dataspacevertex_K3_p );
        dataset_K3_t -> write( K3_class_t, mtype1,dataspacevertex_K3_t_buffer,dataspacevertex_K3_t );
#endif

        cout << "Successfully saved in hdf5 file: " << FILE_NAME << endl;


        //Terminate
        delete[] irreducible_class;
        delete[] selfenergy;

        dataset_lambda->close();
        dataset_params->close();
        dataset_irred->close();
        dataset_self->close();

#if DIAG_CLASS >=1
        delete[] K1_class_a;
        delete[] K1_class_p;
        delete[] K1_class_t;

        dataset_K1_a->close();
        dataset_K1_p->close();
        dataset_K1_t->close();
#endif

#if DIAG_CLASS >=2
        delete[] K2_class_a;
        delete[] K2_class_p;
        delete[] K2_class_t;

        dataset_K2_a->close();
        dataset_K2_p->close();
        dataset_K2_t->close();

#endif

#if DIAG_CLASS >=3
        delete[] K3_class_a;
        delete[] K3_class_p;
        delete[] K3_class_t;

        dataset_K3_a->close();
        dataset_K3_p->close();
        dataset_K3_t->close();
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

////writes state of new iteration into en EXISTING Hdf5 file. The second argument denotes the iteration number to which it is written. The thrid arguemnt denotes the total number of saved iterations in the file
void add_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size,State<complex<double>>& state_in, vector<double>& Lambdas ){
    //Try block to detect exceptions raised by any of the calls inside it
    cout << "Starting to copy to buffer.." << endl;
    if(Lambda_it < Lambda_size){
        try
        {
            typedef struct complex_t{
                double spin_re;   /*real part */
                double spin_im;   /*imaginary part */
                double dens_re;   /*real part */
                double dens_im;   /*imaginary part */
            }complex_t;


            typedef struct complex{
                double re;   /*real part */
                double im;   /*imaginary part */

            }complex;


           //buffer self energy
            int self_dim =2*nSE;
            auto selfenergy = new complex[self_dim];//irrdeucible vertex
            for(int i=0; i<self_dim; i++){
                selfenergy[i].re = real(state_in.selfenergy.acc(i));
                selfenergy[i].im = imag(state_in.selfenergy.acc(i));
            };




            //buffer irreducible vertex:
            int irred_dim1 =16*n_in;
            auto irreducible_class = new complex_t[irred_dim1];//irrdeucible vertex


            for(int i=0;i<irred_dim1;i++){

                irreducible_class[i].spin_re = real(state_in.vertex.spinvertex.irred.acc(i));
                irreducible_class[i].dens_re = real(state_in.vertex.densvertex.irred.acc(i));
                irreducible_class[i].spin_im = imag(state_in.vertex.spinvertex.irred.acc(i));
                irreducible_class[i].dens_im = imag(state_in.vertex.densvertex.irred.acc(i));
            }

#if DIAG_CLASS >=1
            int K1_dim1=nK_K1 * nw1_wt * n_in;

            auto K1_class_a = new complex_t[K1_dim1];
            auto K1_class_p = new complex_t[K1_dim1];
            auto K1_class_t = new complex_t[K1_dim1];
            for(int i=0;i<K1_dim1;i++){
                K1_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K1_acc(i));
                K1_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K1_acc(i));
                K1_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K1_acc(i));
                K1_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K1_acc(i));

                K1_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K1_acc(i));
                K1_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K1_acc(i));
                K1_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K1_acc(i));
                K1_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K1_acc(i));

                K1_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K1_acc(i));
                K1_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K1_acc(i));
                K1_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K1_acc(i));
                K1_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K1_acc(i));
            }

#endif

#if DIAG_CLASS >=2
            int K2_dim1=nK_K2 * nw2_wt * nw2_nut * n_in;

            auto K2_class_a = new complex_t[K2_dim1];
            auto K2_class_p = new complex_t[K2_dim1];
            auto K2_class_t = new complex_t[K2_dim1];
            for(int i=0;i<K2_dim1;i++){
                K2_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K2_acc(i));
                K2_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K2_acc(i));
                K2_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K2_acc(i));
                K2_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K2_acc(i));

                K2_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K2_acc(i));
                K2_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K2_acc(i));
                K2_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K2_acc(i));
                K2_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K2_acc(i));

                K2_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K2_acc(i));
                K2_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K2_acc(i));
                K2_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K2_acc(i));
                K2_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K2_acc(i));
            }
#endif

#if DIAG_CLASS >= 3
//buffer all K3-class-arrays:
            int K3_dim1=nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in;

            auto K3_class_a = new complex_t[K3_dim1];
            auto K3_class_p = new complex_t[K3_dim1];
            auto K3_class_t = new complex_t[K3_dim1];
            for(int i=0;i<K3_dim1;i++){
                K3_class_a[i].spin_re=real(state_in.vertex.spinvertex.avertex.K3_acc(i));
                K3_class_a[i].dens_re=real(state_in.vertex.densvertex.avertex.K3_acc(i));
                K3_class_a[i].spin_im=imag(state_in.vertex.spinvertex.avertex.K3_acc(i));
                K3_class_a[i].dens_im=imag(state_in.vertex.densvertex.avertex.K3_acc(i));

                K3_class_p[i].spin_re=real(state_in.vertex.spinvertex.pvertex.K3_acc(i));
                K3_class_p[i].dens_re=real(state_in.vertex.densvertex.pvertex.K3_acc(i));
                K3_class_p[i].spin_im=imag(state_in.vertex.spinvertex.pvertex.K3_acc(i));
                K3_class_p[i].dens_im=imag(state_in.vertex.densvertex.pvertex.K3_acc(i));

                K3_class_t[i].spin_re=real(state_in.vertex.spinvertex.tvertex.K3_acc(i));
                K3_class_t[i].dens_re=real(state_in.vertex.densvertex.tvertex.K3_acc(i));
                K3_class_t[i].spin_im=imag(state_in.vertex.spinvertex.tvertex.K3_acc(i));
                K3_class_t[i].dens_im=imag(state_in.vertex.densvertex.tvertex.K3_acc(i));
            }
#endif



            cout << "Buffer ready. Preparing for saving into Hdf5 file..." << endl;
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            H5::Exception::dontPrint();

            // Open an existing file and dataset.
            //   H5::H5File file(FILE_NAME, H5F_ACC_RDWR);
            H5::H5File* file = 0;
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);

            H5::DataSet dataset_self = file->openDataSet("selflist");
            H5::DataSet dataset_irred = file->openDataSet("irred");

#if DIAG_CLASS >=1
            H5::DataSet dataset_K1_a = file->openDataSet("K1_a");
            H5::DataSet dataset_K1_p = file->openDataSet("K1_p");
            H5::DataSet dataset_K1_t = file->openDataSet("K1_t");
#endif

#if DIAG_CLASS >=2
            H5::DataSet dataset_K2_a = file->openDataSet("K2_a");
            H5::DataSet dataset_K2_p = file->openDataSet("K2_p");
            H5::DataSet dataset_K2_t = file->openDataSet("K2_t");
#endif

#if DIAG_CLASS >=3
            H5::DataSet dataset_K3_a = file->openDataSet("K3_a");
            H5::DataSet dataset_K3_p = file->openDataSet("K3_p");
            H5::DataSet dataset_K3_t = file->openDataSet("K3_t");
#endif



            //Create memory datda type
            H5::CompType mtype1( sizeof(complex_t) );
            mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
            mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
            mtype1.insertMember( MEMBER3, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
            mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);


            H5::CompType mtype2( sizeof(complex) );
            mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
            mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);



            complex_t fillvalue_vert;
            fillvalue_vert.spin_re = 0;
            fillvalue_vert.dens_re = 0;
            fillvalue_vert.spin_im = 0;
            fillvalue_vert.dens_im = 0;
            H5::DSetCreatPropList plist_vert;
            plist_vert.setFillValue(mtype1, &fillvalue_vert);

            complex fillvalue_self;
            fillvalue_self.re = 0;
            fillvalue_self.im = 0;

            H5::DSetCreatPropList plist_self;
            plist_self.setFillValue(mtype2, &fillvalue_self);



            // Create the dimension arrays for objects in file and in buffer
            hsize_t self_dims[]={static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(self_dim)};
            hsize_t irreducible_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

            hsize_t self_dims_buffer[]={static_cast<hsize_t>(self_dim)};
            hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#if DIAG_CLASS >=1
            hsize_t K1_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K1_dim1)};   // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#endif
#if DIAG_CLASS >=2
            hsize_t K2_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K2_dim1)}; // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
#endif
#if DIAG_CLASS >=3
            hsize_t K3_dims[]= {static_cast<hsize_t>(Lambda_size),static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)
            hsize_t K3_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

#endif

            // Create the data space for the dataset in file and for the buffer objects
            H5::DataSpace dataspaceself(RANK_self, self_dims);
            H5::DataSpace dataspacevertex_irreducible(RANK_irreducible,irreducible_dims);//data space for vertex with three dimensions (three independent frequencies)

            H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);
            H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

#if DIAG_CLASS >=1
            H5::DataSpace dataspacevertex_K1_a(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
            H5::DataSpace dataspacevertex_K1_p(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)
            H5::DataSpace dataspacevertex_K1_t(RANK_K1,K1_dims);//data space for vertex with three dimensions (one independent frequencies)

            H5::DataSpace dataspacevertex_K1_a_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_p_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >=2
            H5::DataSpace dataspacevertex_K2_a(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
            H5::DataSpace dataspacevertex_K2_p(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)
            H5::DataSpace dataspacevertex_K2_t(RANK_K2,K2_dims);//data space for vertex with three dimensions (twoindependent frequencies)

            H5::DataSpace dataspacevertex_K2_a_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K2_p_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >=3
            H5::DataSpace dataspacevertex_K3_a(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K3_p(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K3_t(RANK_K3, K3_dims);//data space for vertex with three dimensions (three independent frequencies)

            H5::DataSpace dataspacevertex_K3_a_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K3_p_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataspacevertex_K3_t_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

   //         H5::DataSpace dataspacelambda(1, Lambda_dims);

            //Select hyperslab in the file where the data should be located and write buffer data into file
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


            count[1]= self_dim;
            dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataset_self.write( selfenergy, mtype2, dataspaceself_buffer,dataspaceself);

            count[1]= irred_dim1;
            dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataset_irred.write( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);

#if DIAG_CLASS >=1
            count[1]= K1_dim1;

            dataspacevertex_K1_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K1_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            dataset_K1_a.write( K1_class_a, mtype1 ,dataspacevertex_K1_a_buffer,dataspacevertex_K1_a );
            dataset_K1_p.write( K1_class_p, mtype1 ,dataspacevertex_K1_p_buffer,dataspacevertex_K1_p );
            dataset_K1_t.write( K1_class_t, mtype1 ,dataspacevertex_K1_t_buffer,dataspacevertex_K1_t );
#endif

#if DIAG_CLASS >=2
            count[1]= K2_dim1;

            dataspacevertex_K2_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            dataset_K2_a.write( K2_class_a, mtype1 ,dataspacevertex_K2_a_buffer,dataspacevertex_K2_a);
            dataset_K2_p.write( K2_class_p, mtype1 ,dataspacevertex_K2_p_buffer,dataspacevertex_K2_p);
            dataset_K2_t.write( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#endif

#if DIAG_CLASS >=3
            count[1]= K3_dim1;

            dataspacevertex_K3_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K3_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
            dataspacevertex_K3_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

            dataset_K3_a.write( K3_class_a, mtype1,dataspacevertex_K3_a_buffer,dataspacevertex_K3_a );
            dataset_K3_p.write( K3_class_p, mtype1,dataspacevertex_K3_p_buffer,dataspacevertex_K3_p );
            dataset_K3_t.write( K3_class_t, mtype1,dataspacevertex_K3_t_buffer,dataspacevertex_K3_t );

#endif

            double Lambda_list[Lambda_size];
            for(int i=0; i<Lambda_size;i++){
                Lambda_list[i] = Lambdas[i];
            }
            H5::DataSet dataset_lambda = file->openDataSet("lambdas");
            dataset_lambda.write( Lambda_list, H5::PredType::NATIVE_DOUBLE  );//overwrite vector containing all values for lambda



            delete[] irreducible_class;
            delete[] selfenergy;

            dataset_self.close();
            dataset_irred.close();

#if DIAG_CLASS >=1
            delete[] K1_class_a;
            delete[] K1_class_p;
            delete[] K1_class_t;

            dataset_K1_a.close();
            dataset_K1_p.close();
            dataset_K1_t.close();

#endif

#if DIAG_CLASS >=2
            delete[] K2_class_a;
            delete[] K2_class_p;
            delete[] K2_class_t;

            dataset_K2_a.close();
            dataset_K2_p.close();
            dataset_K2_t.close();
#endif

#if DIAG_CLASS >=3
            delete[] K3_class_a;
            delete[] K3_class_p;
            delete[] K3_class_t;

            dataset_K3_a.close();
            dataset_K3_p.close();
            dataset_K3_t.close();
#endif


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

//function to read out an exstiting Hdf5 file. Needed to resume computation if it has been interrupted during the flow
template<typename Q>
State<complex<double>> read_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size, vector<double> &Lambdas){
    State<complex<double>> result;
    if(Lambda_it<Lambda_size){

        result.Lambda = Lambdas[Lambda_it]; // TODO: or Lambdas[Lambda_it-1] ?

        H5::H5File* file = 0;
        file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

        typedef struct complex_t{
            double spin_re;   /*real part */
            double spin_im;   /*imaginary part */
            double dens_re;   /*real part */
            double dens_im;   /*imaginary part */
        }complex_t;


        typedef struct complex{
            double re;   /*real part */
            double im;   /*imaginary part */
        }complex;

        int irred_dim1 =16*n_in;
        int self_dim =2*nSE;

        auto Lambdas_buff = new double[Lambda_size];
        auto selfenergy = new complex[self_dim];
        auto irreducible_class = new complex_t[irred_dim1]; //irreducible vertex


#if DIAG_CLASS >=1
        int K1_dim1=nK_K1 * nw1_wt * n_in;
        auto K1_class_a = new complex_t[K1_dim1];
        auto K1_class_p = new complex_t[K1_dim1];
        auto K1_class_t = new complex_t[K1_dim1];

#endif
#if DIAG_CLASS >=2
        int K2_dim1=nK_K2 * nw2_wt * nw2_nut * n_in;

        auto K2_class_a = new complex_t[K2_dim1];
        auto K2_class_p = new complex_t[K2_dim1];
        auto K2_class_t = new complex_t[K2_dim1];

#endif
#if DIAG_CLASS >=3
        int K3_dim1=nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in;

        auto K3_class_a = new complex_t[K3_dim1];
        auto K3_class_p = new complex_t[K3_dim1];
        auto K3_class_t = new complex_t[K3_dim1];

#endif

        //Create memory datda type
        H5::CompType mtype1( sizeof(complex_t) );
        mtype1.insertMember( MEMBER1, HOFFSET(complex_t, spin_re),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER2, HOFFSET(complex_t, spin_im),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER3, HOFFSET(complex_t, dens_re),  H5::PredType::NATIVE_DOUBLE);
        mtype1.insertMember( MEMBER4, HOFFSET(complex_t, dens_im),  H5::PredType::NATIVE_DOUBLE);


        H5::CompType mtype2( sizeof(complex) );
        mtype2.insertMember( MEMBER5, HOFFSET(complex, re),  H5::PredType::NATIVE_DOUBLE);
        mtype2.insertMember( MEMBER6, HOFFSET(complex, im),  H5::PredType::NATIVE_DOUBLE);




        H5::DataSet dataset_lambdas = file->openDataSet("lambdas");
        H5::DataSet dataset_self = file->openDataSet("selflist");
        H5::DataSet dataset_irred = file->openDataSet("irred");


        H5::DataSpace dataspacelambdas = dataset_lambdas.getSpace();
        H5::DataSpace dataspaceself = dataset_self.getSpace();
        H5::DataSpace dataspacevertex_irreducible = dataset_irred.getSpace();//data space for vertex with three dimensions (three independent frequencies)


        // Create the dimension arrays for objects in buffer.
        hsize_t self_dims_buffer[]={static_cast<hsize_t>(self_dim)};
        hsize_t irreducible_dims_buffer[]= {static_cast<hsize_t>(irred_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        // Create the data space for buffer objects.
        H5::DataSpace dataspaceself_buffer(RANK_self-1, self_dims_buffer);
        H5::DataSpace dataspacevertex_irreducible_buffer(RANK_irreducible-1, irreducible_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

#ifdef SUSC
        H5::DataSet dataset_sus = file->openDataSet("sus");
        H5::DataSpace dataspacevertex_sus= dataset_sus.getSpace();//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >=1
        H5::DataSet dataset_K1_a = file->openDataSet("K1_a");
        H5::DataSet dataset_K1_p = file->openDataSet("K1_p");
        H5::DataSet dataset_K1_t = file->openDataSet("K1_t");

        H5::DataSpace dataspacevertex_K1_a=dataset_K1_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_p=dataset_K1_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t=dataset_K1_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)

        // Create the dimension arrays for objects in buffer.
        hsize_t K1_dims_buffer[]= {static_cast<hsize_t>(K1_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_K1_a_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_p_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K1_t_buffer(RANK_K1-1, K1_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >=2
        H5::DataSet dataset_K2_a = file->openDataSet("K2_a");
        H5::DataSet dataset_K2_p = file->openDataSet("K2_p");
        H5::DataSet dataset_K2_t = file->openDataSet("K2_t");

        H5::DataSpace dataspacevertex_K2_a=dataset_K2_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_p=dataset_K2_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t=dataset_K2_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)

        // Create the dimension arrays for objects in buffer.
        hsize_t K2_dims_buffer[]= {static_cast<hsize_t>(K2_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_K2_a_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_p_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K2_t_buffer(RANK_K2-1, K2_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >=3
        H5::DataSet dataset_K3_a = file->openDataSet("K3_a");
        H5::DataSet dataset_K3_p = file->openDataSet("K3_p");
        H5::DataSet dataset_K3_t = file->openDataSet("K3_t");

        H5::DataSpace dataspacevertex_K3_a=dataset_K3_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_p=dataset_K3_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_t=dataset_K3_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)

        // Create the dimension arrays for objects in buffer.
        hsize_t K3_dims_buffer[]= {static_cast<hsize_t>(K3_dim1)};  // dataset dimensions for vertex (it_nbr lambda iterations, nuc unit cells, 3 sity per unit cell, nw^3 freqs configurations)

        // Create the data space for buffer objects.
        H5::DataSpace dataspacevertex_K3_a_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_p_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataspacevertex_K3_t_buffer(RANK_K3-1, K3_dims_buffer);//data space for vertex with three dimensions (three independent frequencies)

#endif


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


        count[1]= self_dim;
        dataspaceself.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataset_self.read( selfenergy, mtype2, dataspaceself_buffer,dataspaceself);

        count[1]= irred_dim1;
        dataspacevertex_irreducible.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataset_irred.read( irreducible_class, mtype1 ,dataspacevertex_irreducible_buffer,dataspacevertex_irreducible);


#if DIAG_CLASS >=1
        count[1]= K1_dim1;
        dataspacevertex_K1_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K1_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K1_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataset_K1_a.read( K1_class_a, mtype1, dataspacevertex_K1_a_buffer, dataspacevertex_K1_a);
        dataset_K1_p.read( K1_class_p, mtype1, dataspacevertex_K1_p_buffer, dataspacevertex_K1_p);
        dataset_K1_t.read( K1_class_t, mtype1, dataspacevertex_K1_t_buffer, dataspacevertex_K1_t);
#endif
#if DIAG_CLASS >=2
        count[1]= K2_dim1;
        dataspacevertex_K2_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K2_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K2_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataset_K2_a.read( K2_class_a, mtype1 ,dataspacevertex_K2_a_buffer,dataspacevertex_K2_a);
        dataset_K2_p.read( K2_class_p, mtype1, dataspacevertex_K2_p_buffer,dataspacevertex_K2_p);
        dataset_K2_t.read( K2_class_t, mtype1 ,dataspacevertex_K2_t_buffer,dataspacevertex_K2_t);
#endif
#if DIAG_CLASS >=3
        count[1]= K3_dim1;
        dataspacevertex_K3_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K3_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataspacevertex_K3_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataset_K3_a.read( K3_class_a, mtype1, dataspacevertex_K3_a_buffer, dataspacevertex_K3_a);
        dataset_K3_p.read( K3_class_p, mtype1, dataspacevertex_K3_p_buffer, dataspacevertex_K3_p);
        dataset_K3_t.read( K3_class_t, mtype1, dataspacevertex_K3_t_buffer, dataspacevertex_K3_t);

#endif



        dataset_lambdas.read(Lambdas_buff, H5::PredType::NATIVE_DOUBLE,dataspacelambdas);

        for(int i=0; i<Lambda_size;i++){
            Lambdas[i] = Lambdas_buff[i];};



        Q val;//buffer value


        for(int i=0; i<self_dim; i++){
            val={selfenergy[i].re,selfenergy[i].im};
            result.selfenergy.direct_set(i,val);
        };



        for(int i=0; i<irred_dim1;i++){

            val={irreducible_class[i].spin_re,irreducible_class[i].spin_im};//assign complex number
            result.vertex.spinvertex.irred.direct_set(i,val);

            val={irreducible_class[i].dens_re,irreducible_class[i].dens_im};//assign complex number
            result.vertex.densvertex.irred.direct_set(i,val);
        };



        //Note that the first index labels the channel with the correspondance ( 0 = s-channel, 1 = t-channel, 2 = u-channel)


#if DIAG_CLASS >=1
        for (int i=0; i<K1_dim1; i++){
            val={K1_class_a[i].spin_re ,K1_class_a[i].spin_im};//assign complex number
            result.vertex.spinvertex.avertex.K1_direct_set(i,val);
            val={K1_class_a[i].dens_re,K1_class_a[i].dens_im};//assign complex number
            result.vertex.densvertex.avertex.K1_direct_set(i, val);

            val={K1_class_p[i].spin_re ,K1_class_p[i].spin_im};//assign complex number
            result.vertex.spinvertex.pvertex.K1_direct_set(i,val);
            val={K1_class_p[i].dens_re,K1_class_p[i].dens_im};//assign complex number
            result.vertex.densvertex.pvertex.K1_direct_set(i, val);

            val={K1_class_t[i].spin_re ,K1_class_t[i].spin_im};//assign complex number
            result.vertex.spinvertex.tvertex.K1_direct_set(i,val);
            val={K1_class_t[i].dens_re,K1_class_t[i].dens_im};//assign complex number
            result.vertex.densvertex.tvertex.K1_direct_set(i, val);
        };

        delete[] K1_class_a;
        delete[] K1_class_p;
        delete[] K1_class_t;

        dataset_K1_a.close();
        dataset_K1_p.close();
        dataset_K1_t.close();

#endif
#if DIAG_CLASS >=2
        for (int i=0; i<K2_dim1; i++){
            val={K2_class_a[i].spin_re ,K2_class_a[i].spin_im};//assign complex number
            result.vertex.spinvertex.avertex.K2_direct_set(i,val);
            val={K2_class_a[i].dens_re,K2_class_a[i].dens_im};//assign complex number
            result.vertex.densvertex.avertex.K2_direct_set(i, val);

            val={K2_class_p[i].spin_re ,K2_class_p[i].spin_im};//assign complex number
            result.vertex.spinvertex.pvertex.K2_direct_set(i,val);
            val={K2_class_p[i].dens_re,K2_class_p[i].dens_im};//assign complex number
            result.vertex.densvertex.pvertex.K2_direct_set(i, val);

            val={K2_class_t[i].spin_re ,K2_class_t[i].spin_im};//assign complex number
            result.vertex.spinvertex.tvertex.K2_direct_set(i,val);
            val={K2_class_t[i].dens_re,K2_class_t[i].dens_im};//assign complex number
            result.vertex.densvertex.tvertex.K2_direct_set(i, val);
        };

        delete[] K2_class_a;
        delete[] K2_class_p;
        delete[] K2_class_t;

        dataset_K2_a.close();
        dataset_K2_p.close();
        dataset_K2_t.close();
#endif
#if DIAG_CLASS >=3
        for (int i=0; i<K3_dim1; i++){
            val={K3_class_a[i].spin_re ,K3_class_a[i].spin_im};//assign complex number
            result.vertex.spinvertex.avertex.K3_direct_set(i,val);
            val={K3_class_a[i].dens_re,K3_class_a[i].dens_im};//assign complex number
            result.vertex.densvertex.avertex.K3_direct_set(i, val);

            val={K3_class_p[i].spin_re ,K3_class_p[i].spin_im};//assign complex number
            result.vertex.spinvertex.pvertex.K3_direct_set(i,val);
            val={K3_class_p[i].dens_re,K3_class_p[i].dens_im};//assign complex number
            result.vertex.densvertex.pvertex.K3_direct_set(i, val);

            val={K3_class_t[i].spin_re ,K3_class_t[i].spin_im};//assign complex number
            result.vertex.spinvertex.tvertex.K3_direct_set(i,val);
            val={K3_class_t[i].dens_re,K3_class_t[i].dens_im};//assign complex number
            result.vertex.densvertex.tvertex.K3_direct_set(i, val);


        };

        delete[] K3_class_a;
        delete[] K3_class_p;
        delete[] K3_class_t;

        dataset_K3_a.close();
        dataset_K3_p.close();
        dataset_K3_t.close();
#endif

        delete[] selfenergy;
        delete[] irreducible_class;
        dataset_lambdas.close();
        dataset_self.close();
        dataset_irred.close();

        file->close();
        delete file;

cout << "File '" << FILE_NAME<< "' successfully read out" << endl;
        return result;

    }
    else{
        cout << "Cannot read from file " << FILE_NAME << " since Lambda layer out of range" << endl;
    };

}

void test_hdf5(H5std_string FILE_NAME, int i, State<comp>& state) {
    // test hdf5: read files and compare to original file
    State<comp> out = read_hdf<comp>(FILE_NAME, i, nEVO, flow_grid);
    if (state.Lambda != out.Lambda) {
        cout << "Lambdas not equal" << endl;
    }
    for (int iK=0; iK<2; ++iK) {
        for (int iSE = 0; iSE < nSE; ++iSE) {
            if (state.selfenergy.sval(iK, iSE) != out.selfenergy.sval(iK, iSE)) {
                cout << "Self-energy not equal, " << iK << ", " << iSE << endl;
            }
        }
    }

    for (int iK=0; iK<nK_K1; ++iK) {
#if DIAG_CLASS >=1
        for (int iw1=0; iw1<nBOS; ++iw1) {
            if (state.vertex.densvertex.avertex.K1_vval(iK, iw1, 0) != out.vertex.densvertex.avertex.K1_vval(iK, iw1, 0)) {
                cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
            }
            if (state.vertex.densvertex.pvertex.K1_vval(iK, iw1, 0) != out.vertex.densvertex.pvertex.K1_vval(iK, iw1, 0)) {
                cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
            }
            if (state.vertex.densvertex.tvertex.K1_vval(iK, iw1, 0) != out.vertex.densvertex.tvertex.K1_vval(iK, iw1, 0)) {
                cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
            }
#if DIAG_CLASS >= 2
            for (int iw2=0; iw2<nFER; ++iw2) {
                if (state.vertex.densvertex.avertex.K2_vval(iK, iw1, iw2, 0) != out.vertex.densvertex.avertex.K2_vval(iK, iw1, iw2, 0)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                }
                if (state.vertex.densvertex.pvertex.K2_vval(iK, iw1, iw2, 0) != out.vertex.densvertex.pvertex.K2_vval(iK, iw1, iw2, 0)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                }
                if (state.vertex.densvertex.tvertex.K2_vval(iK, iw1, iw2, 0) != out.vertex.densvertex.tvertex.K2_vval(iK, iw1, iw2, 0)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                }
#if DIAG_CLASS == 3
                for (int iw3=0; iw3<nFER; ++iw3) {
                    if (state.vertex.densvertex.avertex.K3_vval(iK, iw1, iw2, iw3, 0) != out.vertex.densvertex.avertex.K3_vval(iK, iw1, iw2, iw3, 0)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                    }
                    if (state.vertex.densvertex.pvertex.K3_vval(iK, iw1, iw2, iw3, 0) != out.vertex.densvertex.pvertex.K3_vval(iK, iw1, iw2, iw3, 0)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                    }
                    if (state.vertex.densvertex.tvertex.K3_vval(iK, iw1, iw2, iw3, 0) != out.vertex.densvertex.tvertex.K3_vval(iK, iw1, iw2, iw3, 0)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                    }
                }
#endif
            }
#endif
        }
#endif
        if (state.vertex.densvertex.irred.bare[iK] != out.vertex.densvertex.irred.bare[iK]) {
            cout << "Vertex not equal, " << iK << endl;
        }
    }
}


