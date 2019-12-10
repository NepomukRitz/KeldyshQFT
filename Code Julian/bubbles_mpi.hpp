#ifndef BUBBLES_MPI_HPP
#define BUBBLES_MPI_HPP
//define bubble functions with parvert as argument and output:

#include"kagome_facil.hpp"

//for site-spin-freq parametrization from Reuther, define bubble function of vertices of the form (Gamma_spin *(Pauli)^2 + Gamma_dens * del^2). See Sb, p. 71
template<class T1, class T2> //last arguments defindes the object of type parvert<T> to which the result is written
parvert<svert> sbubble(int red_side, double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:

    //    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    //    int sum_K2_i;//concerns bos. freq.
    //    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    //    else if( sym== 2){sum_K2_i = nw/2;};

    //    int sum_K2_j;//concerns ferm. freq
    //    if(sym==0 ){sum_K2_j = (nw-nw2)/2;}
    //    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    //    int sum_R = (sym==0?(nw-nw3)/2:nw/2);




    vector<site> initials;
    site site_here1(0,0,2);
    initials.push_back(site_here1);
    site site_here2(0,0,3);
    initials.push_back(site_here2);
    site site_here3(-1,0,2);
    initials.push_back(site_here3);
    site site_here4(-1,1,3);
    initials.push_back(site_here4);



    int sum_K1,sum_K2_i,sum_K2_j,sum_R;
#if sym ==0
    sum_K1 = (nw-nw1)/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = (nw-nw2)/2;
    sum_R = (nw-nw3)/2;
#elif sym==1
    sum_K1 = nw/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#elif sym==2
    sum_K1 = nw/2;
    sum_K2_i = nw/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#endif

    parvert<svert> result;
    vector<site> upper_sites;

    int number=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};


    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);



    vector<double> K1_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));
    vector<double> K1_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));

    vector<double> K2_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));

#if sym==0
    vector<double> K2b_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2b_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));
#endif

    vector<double> K3_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));
    vector<double> K3_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));



    int iterator =0;


    for(int i=sum_K1 ; i< (nw+nw1)/2 ; i++){
        if((i-sum_K1) % world_size == world_rank){

            //site sum only over sites with initial condition nonzero

            for(int m=0; m<4; m++){



                int a = initials[m].a;
                int b = initials[m].b;
                int c = initials[m].c;

                double spinspin_K1;
                double spindens_K1;
                double densspin_K1;
                double densdens_K1;
                //note that the asymptotic behaviour is centered around s/2 (half the bosonic transfer freq in the s-channel)
             #if temp==0
                double s = bfreqs[i];
                gsl_integration_workspace * w
                        = gsl_integration_workspace_alloc (1500);
#elif temp==1

                        gsl_integration_workspace  * w = NULL;

                int s = (i-sum_K1)*2;
#endif
                spinspin_K1 = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                spindens_K1 = sbubble(red_side,0,0,w, Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');
                densspin_K1 = sbubble(red_side,0,0,w, Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');
                densdens_K1 = sbubble(red_side,0,0,w, Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,wlimit,'K');


                //set spin vertex. Combinatorial factors, see SBI
                //       result.spinvertex.K1_setvert(a,b,c,i,-4. * spinspin_K1 + 2. * densspin_K1 + 2. * spindens_K1);//combine vertices with and without tilde since u- and t-channel always appear as a sum

                //set density vertex
                //         result.densvertex.K1_setvert(a,b,c,i,6. * spinspin_K1 + 2. * densdens_K1);

                K1_spin_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = -4. * spinspin_K1 + 2. * densspin_K1 + 2. * spindens_K1;
                K1_dens_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] =  6. * spinspin_K1 + 2. * densdens_K1;



                // K2 and K2b
                //#pragma omp for
                //for(int m=0; m<number; m++){
                //    gsl_integration_workspace * w
                //            = gsl_integration_workspace_alloc (1500);

                //    int a = upper_sites[m].a;
                //    int b = upper_sites[m].b;
                //    int c = upper_sites[m].c;
                //        for(int i=sum_K2_i ; i<(nw+nw2)/2; i++)
#if temp==0
                gsl_integration_workspace_free(w);
           #endif
                if(i>= sum_K2_i && i <(nw+nw2)/2){
#pragma omp parallel
                    {
#if temp==0
                        gsl_integration_workspace * v
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * v = NULL;
#endif


#pragma omp for
                        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){

                            double spinspin_K2;
                            double spindens_K2;
                            double densspin_K2;
                            double densdens_K2;

#if temp==0
                            double s = bfreqs[i];
                            double  w1 = ffreqs[j];
#elif temp==1
 int s = (i-sum_K2_i)*2;
 int w1=(s%4==0? (j-sum_K2_j)*2+1:(j-sum_K2_j)*2);
#endif


                            spinspin_K2 = sbubble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                            spindens_K2 = sbubble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');
                            densspin_K2 = sbubble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');
                            densdens_K2 = sbubble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,wlimit,'L');                       //set spin vertex. Combinatorial factors, see SBI
                            //    result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-4. * spinspin_K2 + 2. * densspin_K2 + 2. * spindens_K2);
                            //   result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,6. * spinspin_K2 + 2. * densdens_K2);


                            K2_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] = -4. * spinspin_K2 + 2. * densspin_K2 + 2. * spindens_K2;
                            K2_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =  6. * spinspin_K2 + 2. * densdens_K2;


#if sym==0//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                            double bspinspin_K2;
                            double bspindens_K2;
                            double bdensspin_K2;
                            double bdensdens_K2;

                            bspinspin_K2 = sbubble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                            bspindens_K2 = sbubble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                            bdensspin_K2 = sbubble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                            bdensdens_K2 = sbubble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,wlimit,w1,'M');
                            //     result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-4. * bspinspin_K2 + 2. * bdensspin_K2 + 2. * bspindens_K2 );
                            //     result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,6 * bspinspin_K2 + 2. * bdensdens_K2);
                            K2b_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2-(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] = -4. * bspinspin_K2 + 2. * bdensspin_K2 + 2. * bspindens_K2;
                            K2b_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2-(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =  6. * bspinspin_K2 + 2. * bdensdens_K2;

#endif

                        };
                    #if temp==0
                        gsl_integration_workspace_free(v);
                    #endif
                    };
                };

            };
#pragma omp parallel
            {
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
#pragma omp for
                for(int m=0; m<number; m++){



                    int a = upper_sites[m].a;
                    int b = upper_sites[m].b;
                    int c = upper_sites[m].c;
                    //R (rest function):

                    if(i>= sum_R && i <(nw+nw3)/2){
                        //   for(int i=sum_R ; i<(nw+nw3)/2; i++){
                        for(int j=sum_R ; j<(nw+nw3)/2; j++){
                            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble

#if sym==2
                                if((abs(j-nw/2) <= abs(k-nw/2))){
                                    double spinspin_K3;
                                    double spindens_K3;
                                    double densspin_K3;
                                    double densdens_K3;

#if temp==0
                                    double s = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
#elif temp==1
                                    int s = (i-sum_R)*2;
                                    int w1,w2=0;
                                    if(s%4==0){
                                        w1=(j-sum_R)*2+1;
                                        w2=(k-sum_R)*2+1;}
                                    else{
                                        w1=(j-sum_R)*2;
                                        w2=(k-sum_R)*2;
                                    };
#endif
                                    spinspin_K3 = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                                    spindens_K3 = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
                                    densspin_K3 = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
                                    densdens_K3 = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');

                                    //set spin vertex. Combinatorial factors, see SBI

                                    //   result.spinvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2  ,-4. * spinspin_K3 + 2. * densspin_K3 + 2. * spindens_K3);
                                    //set density vertex

                                    //        result.densvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2 ,6. * spinspin_K3 + 2. * densdens_K3  );
                                    K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = -4. * spinspin_K3 + 2. * densspin_K3 + 2. * spindens_K3;
                                    K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] =  6. * spinspin_K3 + 2. * densdens_K3;

                                };

#else
                                double spinspin_K3;
                                double spindens_K3;
                                double densspin_K3;
                                double densdens_K3;

                                double s = bfreqs[i];
                                double w1 = ffreqs[j];
                                double w2 = ffreqs[k];

                                spinspin_K3 = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');//these four grids contain all frequency values for fixed (a,b,c) and Lambda
                                spindens_K3 = sbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
                                densspin_K3 = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');
                                densdens_K3 = sbubble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,s,w1,w2,'R');

                                //set spin vertex. Combinatorial factors, see SBI

                                //   result.spinvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2  ,-4. * spinspin_K3 + 2. * densspin_K3 + 2. * spindens_K3);
                                //set density vertex

                                //        result.densvertex.R_setvert(a,b,c,i-(nw1-nw3)/2,j-(nw1-nw3)/2,k-(nw1-nw3)/2 ,6. * spinspin_K3 + 2. * densdens_K3  );

                                K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = -4. * spinspin_K3 + 2. * densspin_K3 + 2. * spindens_K3;
                                K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] =  6. * spinspin_K3 + 2. * densdens_K3;


#endif
                            };};};



                };
#if temp==0
                gsl_integration_workspace_free(w);
            #endif
            };
            iterator += 1;
        };
    };


    // //collect data and write it to result:
    vector<double> K1_spin_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));
    vector<double> K1_dens_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));

    vector<double> K2_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));

#if sym==0
    vector<double> K2b_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2b_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
#endif

    vector<double> K3_spin_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));
    vector<double> K3_dens_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));


    ////   //Send all information to all nodes
    MPI_Allgather(&K1_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_spin_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K1_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_dens_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgather(&K2_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

#if sym==0
    MPI_Allgather(&K2b_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2b_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
#endif

    MPI_Allgather(&K3_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K3_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    vector<double> final_spin_K1;
    vector<double> final_dens_K1;

    vector<double> final_spin_K2;
    vector<double> final_dens_K2;
#if sym==0
    vector<double> final_spin_K2b;
    vector<double> final_dens_K2b;
#endif
    vector<double> final_spin_K3;
    vector<double> final_dens_K3;

    ////reorder the vector such that it can directly be copied to the class objects
    ////K1:
    for(int i=0; i<(nw1_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw1_q%world_size > j){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);};};

    //K2:
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);};};

    //K2b
#if sym==0
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);};};
#endif
    //K3
    for(int i=0; i<(nw3_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw3_q%world_size > j){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);};};

    ////write into class object of result

    result.spinvertex.K1_copy(final_spin_K1);
    result.densvertex.K1_copy(final_dens_K1);

    result.spinvertex.K2_copy(final_spin_K2);
    result.densvertex.K2_copy(final_dens_K2);

#if sym==0
    result.spinvertex.K2b_copy(final_spin_K2b);
    result.densvertex.K2b_copy(final_dens_K2b);
#endif

    result.spinvertex.K3_copy(final_spin_K3);
    result.densvertex.K3_copy(final_dens_K3);


    return result;
}

//****************three types of t-bubbles in full parametrization ***********************
//NOTE: the last argument is a parvert<T> vertex to which is result is written in the ZEROTH component
template<class T1, class T2>//this is the RPA term, compare to SB1, page 73.
parvert<tvert> tbubble1(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;
    //mixed combinations like spindens or densspin do not contribute since they cancel out (see Sb1, p. 66)

#pragma omp parallel for
    //  K1:
    for(int i=sum_K1 ; i< (nw+nw1)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){

                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspinval, densdensval;
                        double t = bfreqs[i] ;
                        for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.


                            if(f==1){
                                for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                spinspinval += tbubble(red_side, 0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                            };};};};}

                            else if(f==2){
                                for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        int ap = a-d;
                                        int bp = b-e;//translate such that j lies in first Wigner Seitz cell

                                        int app = -ap-bp;
                                        int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                            };};};};}
                            else if(f==3){
                                for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                    for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                        int ap = a - d;
                                        int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                        int app = bp;
                                        int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                        if(distance(d,e,f) <= d_c){//cutoff distance
                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                                densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                            };

                                        };};};};
                        };
                        result.spinvertex.K1_setvert(a,b,c,i,2.*spinspinval);
                        result.densvertex.K1_setvert(a,b,c,i,2.*densdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                        spinspinval = 0;//reset values for next iteration step with different frequencies
                        densdensval = 0;
                    };

                };};};
#if temp==0
        gsl_integration_workspace_free(w);
#endif
    };


    //   K2 and K2b:


#pragma omp parallel for

    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){

                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspinval, densdensval;
                            double t = bfreqs[i];
                            double  w1 = ffreqs[j];
                            for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                if(f==1){
                                    for(int d =-(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};};};}

                                else if(f==2){
                                    for(int d =-(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a- d;
                                            int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                            int app = -ap-bp;
                                            int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};};};}
                                else if(f==3){
                                    for(int d = -(nuc_eff-1)/2; d <(nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a - d;
                                            int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                            int app = bp;
                                            int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                    spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                    densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                                };};

                                        };};};
                            };

                            result.spinvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*spinspinval );
                            result.densvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*densdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                            spinspinval = 0;//reset values for next iteration step with different frequencies
                            densdensval =0;


                            //K2b (interchange k and j in bubble calculations):


#if sym==0//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                            double bspinspinval, bdensdensval;
                            for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                if(f==1){
                                    for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                    bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                };};};};}

                                else if(f==2){
                                    for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a- d;
                                            int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                            int app = -ap-bp;
                                            int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                    bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                };};};};}
                                else if(f==3){
                                    for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                        for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                            int ap = a - d;
                                            int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                            int app = bp;
                                            int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)
                                            if(distance(d,e,f) <= d_c){//cutoff distance
                                                if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                    bspinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                    bdensdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                                };};

                                        };};};
                            };
                            result.spinvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bspinspinval);//
                            result.densvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bdensdensval);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                            bspinspinval = 0;//reset values for next iteration step with different frequencies
                            bdensdensval = 0;


#endif
                        };};

                };};};
#if temp==0
        gsl_integration_workspace_free(w);
    #endif
    };

    // R (rest function):

#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w= NULL;
#endif
        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){

                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspinval, densdensval;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];

                                    double w2 = ffreqs[k];

                                    for(int f=1; f<4 ; f++){//iterates through different atoms within each unit cell for dummy site j.
                                        if(f==1){
                                            for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};};};}

                                        else if(f==2){
                                            for(int d =-(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e=-(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    int ap = a- d;
                                                    int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                                    int app = -ap-bp;
                                                    int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                    int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)
                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};};};}
                                        else if(f==3){
                                            for(int d = -(nuc_eff-1)/2; d < (nuc_eff-1)/2+1;d++){//iterates through different unit cells for dummy site j.
                                                for(int e= -(nuc_eff-1)/2; e<(nuc_eff-1)/2+1; e++){
                                                    int ap = a - d;
                                                    int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                                    int app = bp;
                                                    int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                                    int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                                                    if(distance(d,e,f) <= d_c){//cutoff distance
                                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                            spinspinval += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                            densdensval += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                        };};

                                                };};};
                                    };
                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspinval );
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*densdensval  );

                                    spinspinval = 0;//reset values for next iteration step with different frequencies
                                    densdensval = 0;


                                };};};};};};};
#if temp==0
        gsl_integration_workspace_free(w);
    #endif
    };





    return result;
}

template<class T1, class T2>
parvert<tvert> tbubble2(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2,char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;



#pragma omp parallel for
    //K1:
    for(int i=sum_K1; i< (nw+nw1)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){


                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspin, spindens, densspin, densdens;
                        double t = bfreqs[i];
                        spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                        spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                        result.spinvertex.K1_setvert(a,b,c,i,spinspin - densspin);
                        result.densvertex.K1_setvert(a,b,c,i-(nw-nw1)/2,-3.*spindens -1.*densdens);
                    };
                };};};
#if temp==0
        gsl_integration_workspace_free(w);
    #endif
    };

#pragma omp parallel for
    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif

        for(int j=sum_K2_j; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){


                        //K2 and K2b:

                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspin, spindens, densspin, densdens;
                            double t = bfreqs[i];
                            double w1 = ffreqs[j] ;
                            spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                            spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                            result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*spinspin -1.*densspin );
                            result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*spindens -1.*densdens);

#if sym==0//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                            double bspinspin, bspindens, bdensspin, bdensdens;
                            bspinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                            bspindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');

                            result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*bspinspin -1.*bdensspin );
                            result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*bspindens -1.*bdensdens);
#endif
                        };};
                };};};
#if temp==0
        gsl_integration_workspace_free(w);
    #endif
    };

#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in t-bubble
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){
                                //rest function:

                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspin, spindens, densspin, densdens;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
                                    spinspin = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,1.*spinspin -1.*densspin);
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,-3.*spindens -1.*densdens);
                                };};};

                    };};};};
      #if temp==0
        gsl_integration_workspace_free(w);
    #endif
    };

    return result;
}

template<class T1, class T2>
parvert<tvert> tbubble3(int red_side, double Lambda, parvert<T1>& vert1, parvert<T2>& vert2,char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    int sum_K2_i;
    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    else if(sym==2){sum_K2_i = nw/2;};

    int sum_K2_j;
    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    parvert<tvert> result;



#pragma omp parallel for

    //K1:
    for(int i=sum_K1; i< (nw+nw1)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif

        for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                for(int c= 1; c<4; c++){


                    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                        double spinspin, spindens, densspin, densdens;
                        double t = bfreqs[i];
                        spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                        spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                        densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                        result.spinvertex.K1_setvert(a,b,c,i,spinspin - spindens);
                        result.densvertex.K1_setvert(a,b,c,i,-3.*densspin -1.*densdens);
                    };
                };};};
#if temp==0
        gsl_integration_workspace_free(w);
#endif
    };
    //K2 and K2b:

#pragma omp parallel for
    for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                    for(int c= 1; c<4; c++){



                        if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                            double spinspin, spindens, densspin, densdens;
                            double t = bfreqs[i];
                            double w1 = ffreqs[j] ;
                            spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                            spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

                            result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*spinspin -1.*spindens );
                            result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*densspin -1.*densdens);

#if sym==0//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                            double bspinspin, bspindens, bdensspin, bdensdens;
                            bspinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                            bspindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');

                            result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,1.*bspinspin -1.*bspindens);
                            result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,-3.*bdensspin -1.*bdensdens );
#endif
                        };};
                };};};
#if temp==0
        gsl_integration_workspace_free(w);
#endif
    };

    //R:


#pragma omp parallel for
    for(int i=sum_R ; i<(nw+nw3)/2; i++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif

        for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble
                if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
                    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
                        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
                            for(int c= 1; c<4; c++){



                                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                                    double spinspin, spindens, densspin, densdens;
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
                                    spinspin = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                    result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,1.*spinspin -1.*spindens );
                                    result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,-3.*densspin -1.*densdens );

                                };};};};};};};
#if temp==0
        gsl_integration_workspace_free(w);

#endif
    };

    return result;
}

//template<class T1, class T2>
//parvert<tvert> tbubble(double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
//    return(tbubble1(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + tbubble2(Lambda, vert1, vert2,p1,p2,selfenergy, diffselfenergy) + tbubble3(Lambda, vert1, vert2,  p1,p2,selfenergy, diffselfenergy) );
//}





template<class T1, class T2>//this is the RPA term, compare to SB1, page 73.
parvert<tvert> tbubble(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){
    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:

    vector<site> all_sites,upper_sites;

    int number_all=0, number_upper=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= -(nuc_eff-1)/2; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){

                if( distance(a,b,c) <= d_c){
                    number_all +=1;
                    site site_here(a,b,c);
                    all_sites.push_back(site_here);};};};};


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){

                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number_upper +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};




    //    int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

    //    int sum_K2_i;
    //    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    //    else if(sym==2){sum_K2_i = nw/2;};

    //    int sum_K2_j;
    //    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    //    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    //    int sum_R = (sym==0?(nw-nw3)/2:nw/2);


    int sum_K1,sum_K2_i,sum_K2_j,sum_R;
#if sym ==0
    sum_K1 = (nw-nw1)/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = (nw-nw2)/2;
    sum_R = (nw-nw3)/2;
#elif sym==1
    sum_K1 = nw/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#elif sym==2
    sum_K1 = nw/2;
    sum_K2_i = nw/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#endif


    parvert<tvert> result;
    //mixed combinations like spindens or densspin do not contribute since they cancel out (see Sb1, p. 66)

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);



    vector<double> K1_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));
    vector<double> K1_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));

    vector<double> K2_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));

#if sym==0
    vector<double> K2b_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2b_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));
#endif

    vector<double> K3_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));
    vector<double> K3_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));



    int iterator =0;


    for(int i=sum_K1 ; i< (nw+nw1)/2 ; i++){
        if((i-sum_K1) % world_size == world_rank){

#pragma omp parallel
            {

#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
#pragma omp for
                for(int m=0; m<number_upper; m++){



                    int a = upper_sites[m].a;
                    int b = upper_sites[m].b;
                    int c = upper_sites[m].c;


                    //K1 class

#if temp==0
   double t = bfreqs[i];

#elif temp==1
   int t = (i-sum_K1)*2;
#endif
                    //first diagrammatic type (RPA)
                    double spinspinval_K1=0, densdensval_K1=0;
                    double spinspin_2_K1=0,spindens_2_K1=0,densspin_2_K1=0,densdens_2_K1=0;double spinspin_3_K1=0,spindens_3_K1=0,densspin_3_K1=0,densdens_3_K1=0;

                    for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                        int d = all_sites[s].a;
                        int e = all_sites[s].b;
                        int f = all_sites[s].c;

                        //    if(f==1 && distance(a-d,b-e,c) <= d_c){//cutoff distance
                        //stronger condition that knows about initial conditions
                        if(f==1 && distance(a-d,b-e,c) <= d_c){//cutoff distance

                            spinspinval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                            //        densdensval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                        }

                        else if(f==2){

                            int ap = a-d;
                            int bp = b-e;//translate such that j lies in first Wigner Seitz cell

                            int app = -ap-bp;
                            int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                            int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)


                            if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                spinspinval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                //       densdensval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                            };}
                        else if(f==3){

                            int ap = a - d;
                            int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                            int app = bp;
                            int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                            int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                            if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                spinspinval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                                //       densdensval_K1 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                            };

                        };
                    };
if(world_rank==0){
cout << spinspinval_K1 << endl;
};

                    //second diagrammatic type

                    //no contribution to K1-class since they have bare vertices at the endes which are always zero for site combination (0,0,1) from initial condition

                    //                if(distance(a,b,c) <= distance(-1,0,2)){//only contributions when initial condition nonzero
                    //                spinspin_2_K1 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the fisrt vertex is taken with the distance vector 0.(interaction of a site with itself)
                    //                spindens_2_K1 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    //                densspin_2_K1 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    //                densdens_2_K1 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');

                   //third diagrammatic type

                    //                spinspin_3_K1 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                    //                spindens_3_K1 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    //                densspin_3_K1 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');
                    //                densdens_3_K1 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,wlimit,'K');


                    //                //    result.spinvertex.K1_setvert(a,b,c,i,2.*spinspinval+spinspin_2-densspin_2+spinspin_3-spindens_3);
                    //                //    result.densvertex.K1_setvert(a,b,c,i-(nw-nw1)/2,2.*densdensval-3*spindens_2-densdens_2-3*densspin_3-densdens_3);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                    //                };

                    K1_spin_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = 2.*spinspinval_K1+spinspin_2_K1-densspin_2_K1+spinspin_3_K1-spindens_3_K1;
                    K1_dens_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = 2.*densdensval_K1-3*spindens_2_K1-densdens_2_K1-3*densspin_3_K1-densdens_3_K1;



                    //K2-class:
                    if( i>=sum_K2_i && i<(nw+nw2)/2){
                        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){

                            double spinspinval_K2=0, densdensval_K2=0;
                            double bspinspinval_K2, bdensdensval_K2;
                            double spinspin_2_K2=0,spindens_2_K2=0,densspin_2_K2=0,densdens_2_K2=0;double spinspin_3_K2=0,spindens_3_K2=0,densspin_3_K2=0,densdens_3_K2=0;
                            double bspinspin_2_K2=0,bspindens_2_K2=0,bdensspin_2_K2=0,bdensdens_2_K2=0;double bspinspin_3_K2=0,bspindens_3_K2=0,bdensspin_3_K2=0,bdensdens_3_K2=0;

#if temp==0
                            double t = bfreqs[i];
                            double  w1 = ffreqs[j];
#elif temp==1
 int t = (i-sum_K2_i)*2;
 int w1=(t%4==0? (j-sum_K2_j)*2+1:(j-sum_K2_j)*2);
#endif

                            //first diagrammatic contribution (RPA):
                            for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                                int d = all_sites[s].a;
                                int e = all_sites[s].b;
                                int f = all_sites[s].c;
                                if(f==1 && distance(a-d,b-e,c) <= d_c){//cutoff distance
                                    spinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                    //                     densdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
#if sym==0
                                    bspinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                    //                      bdensdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
#endif
                                }

                                else if(f==2){

                                    int ap = a- d;
                                    int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                    int app = -ap-bp;
                                    int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                    int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                                    if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                        spinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                        //                             densdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
#if sym==0
                                        bspinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                        //                            bdensdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
#endif
                                    };}
                                else if(f==3){

                                    int ap = a - d;
                                    int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                    int app = bp;
                                    int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                    int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)

                                    if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                        spinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                                        //                            densdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
#if sym==0
                                        bspinspinval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                                        //                            bdensdensval_K2 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
#endif
                                    };};



                            };

                            //second diagrammatic contribution

                            //    spinspin_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                            //    spindens_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            //    densspin_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            //    densdens_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

#if sym==0
                            bspinspin_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                            //               bspindens_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            bdensspin_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            //               bdensdens_2_K2 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
#endif

                            //third diagrammatic contribution
                            spinspin_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                            spindens_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            //                 densspin_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');
                            //                 densdens_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,wlimit,'L');

#if sym==0
                            //    bspinspin_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                            //    bspindens_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            //    bdensspin_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
                            //    bdensdens_3_K2 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,wlimit,w1,'M');
#endif

#if sym==0
                            //   result.spinvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bspinspinval_K2+bspinspin_2_K2-bdensspin_2_K2+1.*bspinspin_3_K2 -1.*bspindens_3_K2);//
                            //   result.densvertex.K2b_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*bdensdensval_K2-3*bspindens_2_K2 - bdensdens_2_K2- 3.*bdensspin_3_K2 -1.*bdensdens_3_K2);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed
                            K2b_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =2.*bspinspinval_K2 + bspinspin_2_K2 - bdensspin_2_K2 + 1.*bspinspin_3_K2 - 1.*bspindens_3_K2;
                            K2b_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =  2.*bdensdensval_K2 - 3.*bspindens_2_K2 - bdensdens_2_K2 - 3.*bdensspin_3_K2 - 1.*bdensdens_3_K2;

#endif
                            //  result.spinvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*spinspinval_K2 + spinspin_2_K2 - densspin_2_K2 + 1.*spinspin_3_K2 - 1.*spindens_3_K2 );
                            //  result.densvertex.K2_setvert(a,b,c,i-(nw-nw2)/2,j-(nw-nw2)/2,2.*densdensval_K2 - 3.*spindens_2_K2 - densdens_2_K2 - 3.*densspin_3_K2 - 1.*densdens_3_K2);//set vertex for one specific set of frequencies where the spacial sum of the RPA term has been performed

                            K2_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =2.*spinspinval_K2 + spinspin_2_K2 - densspin_2_K2 + 1.*spinspin_3_K2 - 1.*spindens_3_K2;
                            K2_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =  2.*densdensval_K2 - 3.*spindens_2_K2 - densdens_2_K2 - 3.*densspin_3_K2 - 1.*densdens_3_K2;

                        };};





                    // R (rest function):

                    if( i>=sum_R && i<(nw+nw3)/2){
                        for(int j=sum_R ; j<(nw+nw3)/2; j++){
                            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){
#if sym==2
                                if((abs(j-nw/2) <= abs(k-nw/2))){

#if temp==0
                                    double t = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
#elif temp==1
                                    int t = (i-sum_R)*2;
                                    int w1,w2=0;
                                    if(t%4==0){
                                        w1=(j-sum_R)*2+1;
                                        w2=(k-sum_R)*2+1;}
                                    else{
                                        w1=(j-sum_R)*2;
                                        w2=(k-sum_R)*2;
                                    };
#endif
                                    double spinspinval_K3=0, densdensval_K3=0;
                                    double spinspin_2_K3=0,spindens_2_K3=0,densspin_2_K3=0,densdens_2_K3=0;double spinspin_3_K3=0,spindens_3_K3=0,densspin_3_K3=0,densdens_3_K3=0;
                                    //first diagrammatic contribution (RPA):

                                    for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                                        int d = all_sites[s].a;
                                        int e = all_sites[s].b;
                                        int f = all_sites[s].c;
                                        if(f==1){
                                            if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                                spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                            };}

                                        else if(f==2){

                                            int ap = a- d;
                                            int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                            int app = -ap-bp;
                                            int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                                spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                            };}
                                        else if(f==3){

                                            int ap = a - d;
                                            int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                            int app = bp;
                                            int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                            int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)


                                            if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                                spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                                densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


                                            };};
                                    };

                                    //second diagrammatic conmtribution
                                    spinspin_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                    //third diagrammatic contribution
                                    spinspin_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                                    spindens_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densspin_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                    densdens_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


                                    // result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspinval + spinspin_2 - densspin_2 + 1.*spinspin_3 - 1.*spindens_3 );
                                    //result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*densdensval - 3*spindens_2 - densdens_2 - 3.*densspin_3 - 1.*densdens_3 );

                                    K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*spinspinval_K3 + spinspin_2_K3 - densspin_2_K3 + 1.*spinspin_3_K3 - 1.*spindens_3_K3 ;
                                    K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*densdensval_K3 - 3*spindens_2_K3 - densdens_2_K3 - 3.*densspin_3_K3 - 1.*densdens_3_K3;



                                };
#else
                                double t = bfreqs[i];
                                double w1 = ffreqs[j];
                                double w2 = ffreqs[k];


                                double spinspinval_K3=0, densdensval_K3=0;
                                double spinspin_2_K3=0,spindens_2_K3=0,densspin_2_K3=0,densdens_2_K3=0;double spinspin_3_K3=0,spindens_3_K3=0,densspin_3_K3=0,densdens_3_K3=0;
                                //first diagrammatic contribution (RPA):

                                for(int s=0; s<number_all ; s++){//iterates through different atoms within each unit cell for dummy site j.
                                    int d = all_sites[s].a;
                                    int e = all_sites[s].b;
                                    int f = all_sites[s].c;
                                    if(f==1){
                                        if(distance(a-d,b-e,c) <= d_c){//cutoff distance
                                            spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,a-d,b-e,c,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                            densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,a-d,b-e,c,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                        };}

                                    else if(f==2){

                                        int ap = a- d;
                                        int bp = b- e;//translate such that j lies in first Wigner Seitz cell

                                        int app = -ap-bp;
                                        int bpp = ap;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 240 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c+1) % 3 + 1;//each such rotation about 240 degr. (clockwise) moves two steps forward (or equivalently one step backw.) in cylic ordering of (1,2,3)

                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance

                                            spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                            densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                        };}
                                    else if(f==3){

                                        int ap = a - d;
                                        int bp = b - e;//translate such that j lies in first Wigner Seitz cell
                                        int app = bp;
                                        int bpp = -ap-bp;//rotate coordinates such that j is the first site in the three-atom basis, i.e. the vector connecting site j and site (d,e,f) can be desribed in terms of our original coordinate system (rotated about 120 degrees clockwise)-> see page 73 in SB1 for derivation of this coordinate trafo
                                        int cp = (c % 3)+1;//each such rotation about 120 degr.(clockwise) moves on step forward in cylic ordering of (1,2,3)


                                        if(distance(app,bpp,cp) <= d_c){//cutoff distance
                                            spinspinval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.spinvertex,app,bpp,cp,vert2.spinvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                            densdensval_K3 += tbubble(red_side,0,0,w,Lambda, vert1.densvertex,app,bpp,cp,vert2.densvertex,d,e,f,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


                                        };};
                                };

                                //second diagrammatic conmtribution
                                spinspin_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the second vertex is taken with the distance vector 0.(interaction of a site with itself)
                                spindens_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.spinvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                densspin_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                densdens_2_K3 = tbubble(red_side,1,0,w,Lambda, vert1.densvertex,0,0,1,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');

                                //third diagrammatic contribution
                                spinspin_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');//Note that the first vertex is taken with the distance vector 0.(interaction of a site with itself)
                                spindens_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                densspin_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');
                                densdens_3_K3 = tbubble(red_side,0,1,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,0,0,1,p1,p2,selfenergy, diffselfenergy,t,w1,w2,'R');


                                // result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspinval + spinspin_2 - densspin_2 + 1.*spinspin_3 - 1.*spindens_3 );
                                //result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*densdensval - 3*spindens_2 - densdens_2 - 3.*densspin_3 - 1.*densdens_3 );

                                K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*spinspinval_K3 + spinspin_2_K3 - densspin_2_K3 + 1.*spinspin_3_K3 - 1.*spindens_3_K3 ;
                                K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*densdensval_K3 - 3*spindens_2_K3 - densdens_2_K3 - 3.*densspin_3_K3 - 1.*densdens_3_K3;

#endif

                            };};};


                };
#if temp==0
                gsl_integration_workspace_free(w);

            #endif
            };

            iterator += 1;
        };
    };

    // //collect data and write it to result:
    vector<double> K1_spin_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));
    vector<double> K1_dens_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));

    vector<double> K2_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));

#if sym==0
    vector<double> K2b_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2b_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
#endif

    vector<double> K3_spin_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));
    vector<double> K3_dens_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));


    ////   //Send all information to all nodes
    MPI_Allgather(&K1_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_spin_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K1_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_dens_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgather(&K2_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

#if sym==0
    MPI_Allgather(&K2b_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2b_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
#endif

    MPI_Allgather(&K3_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K3_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    vector<double> final_spin_K1;
    vector<double> final_dens_K1;

    vector<double> final_spin_K2;
    vector<double> final_dens_K2;
#if sym==0
    vector<double> final_spin_K2b;
    vector<double> final_dens_K2b;
#endif
    vector<double> final_spin_K3;
    vector<double> final_dens_K3;


    ////reorder the vector such that it can directly be copied to the class objects
    ////K1:
    for(int i=0; i<(nw1_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw1_q%world_size > j){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);};};

    //K2:
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);};};

    //K2b
#if sym==0
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);};};
#endif
    //K3
    for(int i=0; i<(nw3_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw3_q%world_size > j){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);};};
    ////write into class object of result

    result.spinvertex.K1_copy(final_spin_K1);
    result.densvertex.K1_copy(final_dens_K1);

    result.spinvertex.K2_copy(final_spin_K2);
    result.densvertex.K2_copy(final_dens_K2);

#if sym==0
    result.spinvertex.K2b_copy(final_spin_K2b);
    result.densvertex.K2b_copy(final_dens_K2b);
#endif

    result.spinvertex.K3_copy(final_spin_K3);
}








//*******************************************************************************************
template<class T1, class T2>
parvert<uvert> ububble(int red_side,double Lambda, parvert<T1>& vert1, parvert<T2>& vert2, char p1, char p2, self selfenergy, self diffselfenergy){


    //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
    //    int sum_K1 = (sym==0?0:nw/2);

    //    int sum_K2_i;
    //    if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
    //    else if(sym==2){sum_K2_i = nw/2;};

    //    int sum_K2_j;
    //    if(sym==0){sum_K2_j = (nw-nw2)/2;}
    //    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

    //    int sum_R = (sym==0?(nw-nw3)/2:nw/2);

    int sum_K1,sum_K2_i,sum_K2_j,sum_R;
#if sym ==0
    sum_K1 = (nw-nw1)/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = (nw-nw2)/2;
    sum_R = (nw-nw3)/2;
#elif sym==1
    sum_K1 = nw/2;
    sum_K2_i = (nw-nw2)/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#elif sym==2
    sum_K1 = nw/2;
    sum_K2_i = nw/2;
    sum_K2_j = nw/2;
    sum_R = nw/2;
#endif



    vector<site> initials;
    site site_here1(0,0,2);
    initials.push_back(site_here1);
    site site_here2(0,0,3);
    initials.push_back(site_here2);
    site site_here3(-1,0,2);
    initials.push_back(site_here3);
    site site_here4(-1,1,3);
    initials.push_back(site_here4);



    parvert<uvert> result;
    vector<site> upper_sites;

    int number=0;


    for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
        for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4 ; c++){
                if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
                    number +=1;
                    site site_here(a,b,c);
                    upper_sites.push_back(site_here);};};};};


    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);



    vector<double> K1_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));
    vector<double> K1_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1));

    vector<double> K2_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));

#if sym==0
    vector<double> K2b_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * (nw2_q/world_size+1));
    vector<double> K2b_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *(nw2_q/world_size+1));
#endif

    vector<double> K3_spin_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));
    vector<double> K3_dens_buffer(3*(nuc_eff+1)/2*nuc_eff * nw3_w2 *nw3_w1 *(nw3_q/world_size+1));



    int iterator =0;


    for(int i=sum_K1 ; i< (nw+nw1)/2 ; i++){
        if((i-sum_K1) % world_size == world_rank){




            for(int m=0; m<4; m++){
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif

                int a = initials[m].a;
                int b = initials[m].b;
                int c = initials[m].c;

                //K1:



                double spinspin_K1, spindens_K1, densspin_K1, densdens_K1;
#if temp==0
   double u = bfreqs[i];

#elif temp==1
   int u = (i-sum_K1)*2;
#endif
                //compute all four possible vertex combinations in s-bubble
                spinspin_K1 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
                spindens_K1 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
                densspin_K1 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
                densdens_K1 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,wlimit,'K');
                //  result.spinvertex.K1_setvert(a,b,c,i,2.*spinspin_K1 + densspin_K1 + spindens_K1);
                //  result.densvertex.K1_setvert(a,b,c,i, 3.* spinspin_K1 + densdens_K1);


                K1_spin_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = 2.*spinspin_K1 + densspin_K1 + spindens_K1;
                K1_dens_buffer[iterator * K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = 3.* spinspin_K1 + densdens_K1;


                //K2 and K2b:


#if temp==0
                gsl_integration_workspace_free(w);
#endif
                if(i>=sum_K2_i && i<(nw+nw2)/2){
#pragma omp parallel
                    {
#if temp==0
                        gsl_integration_workspace * v
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * v= NULL;
#endif
#pragma omp for
                        for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){


                            double spinspin_K2, spindens_K2, densspin_K2, densdens_K2;
#if temp==0
                            double u = bfreqs[i];
                            double  w1 = ffreqs[j];
#elif temp==1
                            int u = (i-sum_K2_i)*2;
                            int w1=(u%4==0? (j-sum_K2_j)*2+1:(j-sum_K2_j)*2);
#endif
                            spinspin_K2 = ububble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                            spindens_K2 = ububble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                            densspin_K2 = ububble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                            densdens_K2 = ububble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,wlimit,'L');
                            // result.spinvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,2.*spinspin_K2  +  densspin_K2 + spindens_K2);
                            // result.densvertex.K2_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2, 3.* spinspin_K2 + densdens_K2 );

                            K2_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] = 2.*spinspin_K2  +  densspin_K2 + spindens_K2;
                            K2_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2 -(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] = 3.* spinspin_K2 + densdens_K2;


#if sym==0//if symmetry relations are used, the vertex K2b does not need to be computed explicitely
                            double bspinspin_K2, bspindens_K2, bdensspin_K2, bdensdens_K2;

                            bspinspin_K2 = ububble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                            bspindens_K2 = ububble(red_side,0,0,v,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                            bdensspin_K2 = ububble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                            bdensdens_K2 = ububble(red_side,0,0,v,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,wlimit,w1,'M');
                            //    result.spinvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2,2.*bspinspin_K2 +  bdensspin_K2 + bspindens_K2 );
                            //    result.densvertex.K2b_setvert(a,b,c,i-(nw1-nw2)/2,j-(nw1-nw2)/2, 3.* bspinspin_K2 + bdensdens_K2 );

                            K2b_spin_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2-(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] = 2.*bspinspin_K2 +  bdensspin_K2 + bspindens_K2 ;
                            K2b_dens_buffer[(iterator-(sum_K2_i - sum_K1)) * K2_dim2+ (j-(nw1-nw2)/2-(nw2-nw2_w1)) * K2_dim3 + (a+(nuc_eff-1)/2) * K2_dim4 + b * K2_dim5 + (c-1)] =  3.* bspinspin_K2 + bdensdens_K2;

#endif

                        };
#if temp==0
                        gsl_integration_workspace_free(v);
           #endif
                    };
                };
            };
#pragma omp parallel
            {
#if temp==0
                        gsl_integration_workspace * w
                                = gsl_integration_workspace_alloc (1500);
#elif temp==1
                        gsl_integration_workspace  * w = NULL;
#endif
#pragma omp for
                for(int m=0; m<number; m++){



                    int a = upper_sites[m].a;
                    int b = upper_sites[m].b;
                    int c = upper_sites[m].c;
                    //rest class

                    if(i>=sum_R && i<(nw+nw3)/2){
                        for(int j=sum_R ; j<(nw+nw3)/2; j++){
                            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){    //compute all four possible vertex combinations in s-bubble
#if sym==2
                                if((abs(j-nw/2) <= abs(k-nw/2))){


                                    double spinspin_K3, spindens_K3, densspin_K3, densdens_K3;
#if temp==0
                                    double u = bfreqs[i];
                                    double w1 = ffreqs[j];
                                    double w2 = ffreqs[k];
#elif temp==1
                                    int u = (i-sum_R)*2;
                                    int w1,w2=0;
                                    if(u%4==0){
                                        w1=(j-sum_R)*2+1;
                                        w2=(k-sum_R)*2+1;}
                                    else{
                                        w1=(j-sum_R)*2;
                                        w2=(k-sum_R)*2;
                                    };
#endif
                                    spinspin_K3 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                    spindens_K3 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                    densspin_K3 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                    densdens_K3 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                    //  result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspin_K3 +  densspin_K3 + spindens_K3 );
                                    // result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2, 3.* spinspin_K3 + densdens_K3 );

                                    K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*spinspin_K3 +  densspin_K3 + spindens_K3 ;
                                    K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 3.* spinspin_K3 + densdens_K3 ;
                                };

#else
                                double spinspin_K3, spindens_K3, densspin_K3, densdens_K3;
                                double u = bfreqs[i];
                                double w1 = ffreqs[j];
                                double w2 = ffreqs[k];
                                spinspin_K3 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                spindens_K3 = ububble(red_side,0,0,w,Lambda, vert1.spinvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                densspin_K3 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.spinvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                densdens_K3 = ububble(red_side,0,0,w,Lambda, vert1.densvertex,a,b,c,vert2.densvertex,a,b,c,p1,p2,selfenergy, diffselfenergy,u,w1,w2,'R');
                                //  result.spinvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2,2.*spinspin_K3 +  densspin_K3 + spindens_K3 );
                                // result.densvertex.R_setvert(a,b,c,i-(nw-nw3)/2,j-(nw-nw3)/2,k-(nw-nw3)/2, 3.* spinspin_K3 + densdens_K3 );

                                K3_spin_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 2.*spinspin_K3 +  densspin_K3 + spindens_K3 ;
                                K3_dens_buffer[(iterator-(sum_R - sum_K1)) * K3_dim2+ (j-(nw1-nw3)/2-(nw3-nw3_w1)) * K3_dim3 +(k-(nw1-nw3)/2-(nw3-nw3_w2))*K3_dim4 + (a+(nuc_eff-1)/2) * K3_dim5 + b * K3_dim6 + (c-1)] = 3.* spinspin_K3 + densdens_K3 ;

#endif
                            };};};




                };
#if temp==0
                gsl_integration_workspace_free(w);
            #endif
            };

            iterator += 1;
        };
    };


    // //collect data and write it to result:
    vector<double> K1_spin_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));
    vector<double> K1_dens_results(3*(nuc_eff+1)/2*nuc_eff * ( (nw1_q-(nw1_q%world_size))+world_size));

    vector<double> K2_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 * ( (nw2_q-(nw2_q%world_size))+world_size));

#if sym==0
    vector<double> K2b_spin_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
    vector<double> K2b_dens_results(3*(nuc_eff+1)/2*nuc_eff * nw2_w1 *( (nw2_q-(nw2_q%world_size))+world_size));
#endif

    vector<double> K3_spin_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));
    vector<double> K3_dens_results(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2 * nw3_w1 *( (nw3_q-(nw3_q%world_size))+world_size));


    ////   //Send all information to all nodes
    MPI_Allgather(&K1_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_spin_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K1_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, &K1_dens_results[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff * (nw1_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgather(&K2_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

#if sym==0
    MPI_Allgather(&K2b_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K2b_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, &K2b_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw2_w1 *(nw2_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
#endif

    MPI_Allgather(&K3_spin_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_spin_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&K3_dens_buffer[0], static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, &K3_dens_results[0],  static_cast<int>(3*(nuc_eff+1)/2*nuc_eff *  nw3_w2*  nw3_w1 *(nw3_q/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);

    vector<double> final_spin_K1;
    vector<double> final_dens_K1;

    vector<double> final_spin_K2;
    vector<double> final_dens_K2;
#if sym==0
    vector<double> final_spin_K2b;
    vector<double> final_dens_K2b;
#endif
    vector<double> final_spin_K3;
    vector<double> final_dens_K3;


    ////reorder the vector such that it can directly be copied to the class objects
    ////K1:
    for(int i=0; i<(nw1_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+i) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw1_q%world_size > j){
            vector<double>::const_iterator first_spin = K1_spin_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K1_dim2;

            vector<double>::const_iterator first_dens = K1_dens_results.begin()+ ((j*(nw1_q/world_size+1))+nw1_q/world_size) *K1_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K1_dim2;

            final_spin_K1.insert(final_spin_K1.end(), first_spin, last_spin);
            final_dens_K1.insert(final_dens_K1.end(),first_dens, last_dens);};};

    //K2:
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2.insert(final_spin_K2.end(), first_spin, last_spin);
            final_dens_K2.insert(final_dens_K2.end(),first_dens, last_dens);};};

    //K2b
#if sym==0
    for(int i=0; i<(nw2_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+i) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw2_q%world_size > j){
            vector<double>::const_iterator first_spin = K2b_spin_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K2_dim2;

            vector<double>::const_iterator first_dens = K2b_dens_results.begin()+ ((j*(nw2_q/world_size+1))+nw2_q/world_size) *K2_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K2_dim2;

            final_spin_K2b.insert(final_spin_K2b.end(), first_spin, last_spin);
            final_dens_K2b.insert(final_dens_K2b.end(),first_dens, last_dens);};};
#endif
    //K3
    for(int i=0; i<(nw3_q/world_size); i++){
        for(int j=0; j<world_size; j++){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+i) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);
        };
    };
    for(int j=0; j<world_size; j++){
        if(nw3_q%world_size > j){
            vector<double>::const_iterator first_spin = K3_spin_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_spin = first_spin + K3_dim2;

            vector<double>::const_iterator first_dens = K3_dens_results.begin()+ ((j*(nw3_q/world_size+1))+nw3_q/world_size) *K3_dim2 ;
            vector<double>::const_iterator last_dens = first_dens + K3_dim2;

            final_spin_K3.insert(final_spin_K3.end(), first_spin, last_spin);
            final_dens_K3.insert(final_dens_K3.end(),first_dens, last_dens);};};
    ////write into class object of result

    result.spinvertex.K1_copy(final_spin_K1);
    result.densvertex.K1_copy(final_dens_K1);

    result.spinvertex.K2_copy(final_spin_K2);
    result.densvertex.K2_copy(final_dens_K2);

#if sym==0
    result.spinvertex.K2b_copy(final_spin_K2b);
    result.densvertex.K2b_copy(final_dens_K2b);
#endif

    result.spinvertex.K3_copy(final_spin_K3);
    result.densvertex.K3_copy(final_dens_K3);



    return result;
}

#endif // BUBBLES_MPI_HPP

