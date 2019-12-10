
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
#include"kagome_facil.hpp"
#include"bubbles_mpi.hpp"
//using namespace std;
//using namespace H5;


const int nw1 = 80;//put this to the same value as nw. normally set to 240
const int nw2 = 60;//normally set to 160
const int nw3 = 4;//normally set to 34
const int nuc = 5; // number of unit cells in x-direction: should be an odd number!

const double T=1;

//define variable nuc_eff that allows to compute arbitrary lattice sizes (artifact of lattice parametrization)

const int nuc_eff = (nuc<=7?nuc:(nuc+2));



const int K3_dim6=3;
const int K3_dim5=(nuc_eff+1)/2 * K3_dim6;
const int K3_dim4=nuc_eff  * K3_dim5;
const int K3_dim3=nw3_w2*K3_dim4;
const int K3_dim2=nw3_w1 * K3_dim3;
const int K3_dim1=nw3_q * K3_dim2;


//K1:


const int K1_dim4=3 ;
const int K1_dim3=(nuc_eff+1)/2 * K1_dim4 ;
const int K1_dim2 = nuc_eff * K1_dim3;
const int K1_dim1=nw1_q*K1_dim2;


//K2:

const int K2_dim5=3 ;
const int K2_dim4=(nuc_eff+1)/2 * K2_dim5 ;
const int K2_dim3 = nuc_eff * K2_dim4;
const int K2_dim2=nw2_w1*K2_dim3;
const int K2_dim1=nw2_q * K2_dim2;

//irred:
const int irred_dim3=3;
const int irred_dim2=(nuc_eff+1)/2 * irred_dim3;
const int irred_dim1=nuc_eff*irred_dim2;

#if sym==0
const int nw3_q=nw3;
const int nw3_w1 = nw3;
const int nw3_w2 = nw3;
const int nw2_q = nw2;
const int nw2_w1 = nw2;
const int nw1_q = nw1;
#elif sym==1
const int nw3_q=nw3/2;
const int nw3_w1 = nw3/2;
const int nw3_w2 = nw3;
const int nw2_q = nw2;
const int nw2_w1 = nw2/2;
const int nw1_q = nw1/2;

#elif sym==2
const int nw3_q=nw3/2;
const int nw3_w1 = nw3/2;
const int nw3_w2 = nw3;
const int nw2_q = nw2/2;
const int nw2_w1 = nw2/2;
const int nw1_q = nw1/2;
#endif






/****************************************** The following functions convert between physical frequencies and grid frequencies. Make sure to use right routine according to whether working on a lin or a log grid************************/
#if temp==0
#if grid==1
/****************************LINEAR GRID:***************/
//to convert on full frequency grid:



int fconv(double w){//conversion  on linear grid
    int i; //next lower lying lattice site on frequency vector

    i = w/k - 0.5 + nw/2;
    //        if(i != nw-1 && abs(w- ffreqs[i+1])<1e-12){i+=1;}//avoid rounding errors:
    //        else if(i != 0 && abs(w -ffreqs[i-1])<1e-12){i-=1;}
    return i;
}



//to convert on reduced frequency grid:



int fconv_n(double w, int n){//conversion  on linear grid of size n
    int i; //next lower lying lattice site on frequency vector

    i = w/k - 0.5 + n/2;
    //        if(i != n-1 && abs(w - ffreqs[(nw-n)/2+i+1])<1e-12){i+=1;}//avoid rounding errors:
    //        else if(i != 0 && (w - ffreqs[(nw-n)/2+i-1])<1e-12){i-=1;}
    return i;

}

#elif grid == 2
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


    //    if(i != nw-1 && abs(w - ffreqs[i+1])<1e-12){i+=1;}//avoid rounding errors:
    //    else if(i != 0 && abs( w- ffreqs[i-1])<1e-12){i-=1;};


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


    //    if(i != n-1 && abs(w - ffreqs[(nw-n)/2+i+1]) <1e-12){i+=1;}//avoid rounding errors:
    //    else if(i != 0 && abs(w -ffreqs[(nw-n)/2+i-1])<1e-12){i-=1;};




    return i;


}
#endif

#elif temp==1
int fconv(int w){
    return static_cast<int>(w/2.+(nw-1)/2.);
}

//to convert on reduced frequency grid:

int fconv_n(int w, int n){
  return static_cast<int>(w/2.+(n-1)/2.);

}
#endif


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


void svert::R_direct_set(int i, double value){
    if(i>=0 && i<K3_dim1){
        K3[i] = value;
    };
}
void svert::K1_direct_set(int i, double value){
    if(i>=0 && i<K1_dim1){
        K1[i] = value;
    };
}
void svert::K2_direct_set(int i, double value){
    if(i>=0 && i<K2_dim1){
        K2[i] = value;
    };
}

void svert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

        if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0 /*&&abs(value) > 1e-16*/ ){

                      K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ] = value;
        };
    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
void svert::K1_setvert(int a, int b, int c,int i, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw1 && i>=0 /*&&abs(value) > 1e-16*/ ){

            K1[(i-(nw1-nw1_q))* K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = value;
        };

    }
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
void svert::K1_copy(vector<double> vec){
    if(vec.size()== K1_dim1){
        for(int i=0; i<K1_dim1; i++){
            K1[i] = vec[i];};
    }
    else{cout << "error: K1-Vector cannot be copied due to mismatching dimensions."<< endl;};
}

void svert::K2_copy(vector<double> vec){
    if(vec.size()== K2_dim1){
        for(int i=0; i<K2_dim1; i++){
                   K2[i] = vec[i];
        };
    }
    else{cout << "error: K2-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
#if sym==0
void svert::K2b_copy(vector<double> vec){
    if(vec.size()== K2_dim1){
        for(int i=0; i<K2_dim1; i++){
            Kb2[i] = vec[i];};
    }
    else{cout << "error: K2b-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
#endif
void svert::K3_copy(vector<double> vec){
    if(vec.size()== K3_dim1){
        for(int i=0; i<K3_dim1; i++){
                K3[i] = vec[i];
        };
    }
    else{cout << "error: K3-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
void svert::K2_setvert(int a, int b, int c, int i, int j, double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
        if(i< nw2 && i>= 0 && j< nw2 && j>= 0 /*&&abs(value) > 1e-16*/){

               K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] = value ;
        };}
    else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
#if sym==0
void svert::K2b_setvert(int a, int b, int c, int i, int j,  double value){
    if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c{//only save sites in upper half
            if(i< nw2 && i>= 0 && j< nw2 && j>= 0 /*&&abs(value) > 1e-16*/){


            K2b[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] = value;};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}

            void svert::K2b_direct_set(int i, double value){
                if(i>=0 && i<K2_dim1){
                    K2b[i] = value;
                };
            }
        #endif
            double svert::R_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){


            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;

            double value=0;
            if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        #if sym==0

            value = K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ];


        #elif sym == 1
            if(i>=nw3/2){
            if(j>=nw3/2){
            value = K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ];}
            else if(j<nw3/2){
            site x = site_switch(a,b,c);
            j=nw3-1-j;
            k = nw3-1-k;
            value = K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ];
};

}
            else if(i<nw3/2){
            if(k>=nw3/2){
            site x = site_switch(a,b,c);
            i=nw3-1-i;

            value = conj( K3[ (i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]);//note that j and k are interchanged
}
            else if(k<nw3/2){
            i=nw3-1-i;
            j=nw3-1-j;
            k = nw3-1-k;
            value =conj(K3[ (i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]);//note that j and k are interchanged
}
};
        #elif sym == 2
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;

            if(i<nw3/2){x = site_switch(x.a,x.b,x.c);i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};            

            value = K3[ (i_eff-(nw3-nw3_q))*K3_dim2 + (j_eff-(nw3-nw3_w1)) * K3_dim3 +(k_eff-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ];
        #endif


};



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
        #if sym==0

            value = K1[(i-(nw1-nw1_q))* K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] ;
        #elif sym==1 || sym==2


            if (i>= nw1/2){

            value = K1[(i-(nw1-nw1_q))* K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] ;}
            else if(i<nw1/2){

            site x = site_switch(a,b,c);
            i = nw1-1-i;
            value =  K1[(i-(nw1-nw1_q))* K1_dim2 + (x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1)] ;};
        #endif
};


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
        #if sym==0

            value = K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)];
        #elif sym==1
            if(j >= nw2/2){

            value = K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)];}
            else if (j< nw2/2){


            site x = site_switch(a,b,c);
            j = nw2 -1 - j;

            value = K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)];

};

        #elif sym==2
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i; x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j; x = site_switch(x.a,x.b,x.c);};

            value = K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)];
        #endif
};




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
        #if sym==0

            value = K2b[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)];
        #elif sym==1
            site x = site_switch(a,b,c);
            value = conj(K2_vval(x.a,x.b,x.c, nw2-1-i, j ));
        #elif sym==2

            site x = site_switch(a,b,c);
            value = K2_vval(x.a,x.b,x.c, i, nw2-1-j );

        #endif
};


            return value;
}
            double svert::R_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s, double w1, double w2){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(s)+1e-12 <= ffreqs[(nw3+nw)/2-1] && abs(w1)+1e-12 <=ffreqs[(nw3+nw)/2-1] && abs(w2)+1e-12 <=ffreqs[(nw3+nw)/2-1]){//if frequency arguments are out of range, vertex vanishes
        #if sym==0

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);
            int S_vert = S;
            int W1_vert = W1;
            int W2_vert = W2;


            value = ((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(ffreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-ffreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(ffreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-ffreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+S+1]-ffreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        #elif sym==1

            if(s > 0){
            if((w1)>0 ){
            int S = fconv_n(s,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);

            int S_vert = S-(nw3-nw3_w1);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);

            value = ((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}
            else if((w1)<0 ){
            site x = site_switch(a,b,c);

            int S = fconv_n(s,nw3);
            int W1 = fconv_n(-w1,nw3);
            int W2 = fconv_n(-w2,nw3);
            int S_vert = S-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = -w1 ;
            double w2_eff = -w2 ;

            value = ((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
};}

            else if (s < 0){
            if((w2)>0){

            site x = site_switch(a,b,c);

            int S = fconv_n(-s,nw3);
            int W1 = fconv_n(w2,nw3);
            int W2 = fconv_n(w1,nw3);
            int S_vert =S-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = w2;
            double w2_eff = w1;
            s = -s;


            value = conj((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}

            else if ((w2)<0){
            int S = fconv_n(-s,nw3);
            int W1 = fconv_n(-w2,nw3);
            int W2 = fconv_n(-w1,nw3);

            int S_vert = S-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = -w2;
            double w2_eff = -w1;
            s = -s;

            value = conj((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

};};

        #elif sym==2
            site x(a,b,c);

            double w1_eff = w1, w2_eff = w2, s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if((w1_eff) < 0){w1_eff = -w1_eff; w2_eff = -w2_eff; x = site_switch(x.a,x.b,x.c);};


            int S = fconv_n(s_eff,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int S_vert = S-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);

            value = ((K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(bfreqs[(nw-nw3)/2+S+1]-s_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(S_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ]*(-bfreqs[(nw-nw3)/2+S]+s_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((bfreqs[(nw-nw3)/2+S+1]-bfreqs[(nw-nw3)/2+S])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));



        #endif
}
            else{
            int i,j,k;
            if(abs(s)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(s,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value = R_vval(a,b,c,i,j,k);};};



            return value;
}
            double svert::K1_vvalsmooth(int a_raw, int b_raw, int c_raw,   double s){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(s)  < ffreqs[(nw1+nw)/2-1]){
        #if sym==0
            int S = fconv_n(s,nw1);
            int S_vert = S-(nw1-nw1_q);

            value = ((K1[(S_vert)*K1_dim2 + (a+(nuc_eff-1)/2)*K1_dim3 + b*K1_dim4 +(c-1) ]*(bfreqs[(nw-nw1)/2+S+1]-s)
                      +K1[(S_vert+1)*K1_dim2 + (a+(nuc_eff-1)/2)*K1_dim3 + b*K1_dim4 +(c-1) ]*(-bfreqs[(nw-nw1)/2+S]+s))/(bfreqs[(nw-nw1)/2+S+1]-bfreqs[(nw-nw1)/2+S]));
        #elif sym==1

            if (s > 0){

            int S = fconv_n(s,nw1);
            int S_vert = S-(nw1-nw1_q);
            value = ((K1[(S_vert)*K1_dim2+(a+(nuc_eff-1)/2)*K1_dim3 +b*K1_dim4+(c-1) ]*(bfreqs[S+1]-s)
                      +K1[(S_vert+1)*K1_dim2+(a+(nuc_eff-1)/2)*K1_dim3 +b*K1_dim4+(c-1) ]*(-bfreqs[S]+s))/(bfreqs[S+1]-bfreqs[S]));}
            else if (s<0){
            site x = site_switch(a,b,c);

            double s_eff = -s;
            int S = fconv_n(s_eff,nw1);
            int S_vert = S-(nw1-nw1_q);


            value = conj((K1[(S_vert)*K1_dim2 + (x.a+(nuc_eff-1)/2)*K1_dim3 + x.b*K1_dim4 + (x.c-1) ]*(bfreqs[S+1]-s_eff)
                          +K1[(S_vert+1)*K1_dim2 + (x.a+(nuc_eff-1)/2)*K1_dim3 + x.b*K1_dim4 + (x.c-1) ]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));
};

        #elif sym==2

            if (s > 0){

            int S = fconv_n(s,nw1);
            int S_vert = S-(nw1-nw1_q);
            value = ((K1[(S_vert)*K1_dim2+(a+(nuc_eff-1)/2)*K1_dim3 +b*K1_dim4+(c-1) ]*(bfreqs[S+1]-s)
                      +K1[(S_vert+1)*K1_dim2+(a+(nuc_eff-1)/2)*K1_dim3 +b*K1_dim4+(c-1) ]*(-bfreqs[S]+s))/(bfreqs[S+1]-bfreqs[S]));}
            else if (s<0){
            site x = site_switch(a,b,c);

            double s_eff = -s;
            int S = fconv_n(s_eff,nw1);
            int S_vert = S-(nw1-nw1_q);


            value = ((K1[(S_vert)*K1_dim2 + (x.a+(nuc_eff-1)/2)*K1_dim3 + x.b*K1_dim4 + (x.c-1) ]*(bfreqs[S+1]-s_eff)
                      +K1[(S_vert+1)*K1_dim2 + (x.a+(nuc_eff-1)/2)*K1_dim3 + x.b*K1_dim4 + (x.c-1) ]*(-bfreqs[S]+s_eff))/(bfreqs[S+1]-bfreqs[S]));

};
        #endif

}
            else{
            int i;
            if(abs(s)<= bfreqs[nw1-1]){
            i= fconv_n(s,nw1);

            value = K1_vval(a,b,c,i);
}

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

            if(abs(s) < ffreqs[(nw2+nw)/2-1] && abs(w1) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0
            int S = fconv_n(s,nw2);
            int W1 = fconv_n(w1,nw2);

            int S_vert =S-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(S_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(S_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        #elif sym==1
            if(w1 > 0){


            int S = fconv_n(s,nw2);
            int W1 = fconv_n(w1,nw2);


            int S_vert = S-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(S_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(S_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
}
            else if(w1 < 0){

            site x = site_switch(a,b,c);

            int S = fconv_n(s,nw2);
            int W1 = fconv_n(-w1,nw2);

            int S_vert = S-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);
            double w1_eff = -w1 ;
            value = ((K2[(S_vert)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(S_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
};
        #elif sym==2

            site x(a,b,c);
            double w1_eff = w1;
            double s_eff = s;
            if(s<0){s_eff = -s; x = site_switch(x.a,x.b,x.c);};
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};


            int S = fconv_n(s_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int S_vert = S-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(S_vert)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(S_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                      +K2[(S_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+S]+s_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

        #endif
}
            else{
            if(abs(s)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            int i,j;

            i= fconv_n(s,nw2);
            j = fconv_n(w1,nw2);

            value = K2_vval(a,b,c,i,j);};
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
            if(abs(s) < ffreqs[(nw2+nw)/2-1] && abs(w2) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0

            int S = fconv_n(s,nw2);
            int W2 = fconv_n(w2,nw2);

            int S_vert =S-(nw2-nw2_q);
            int W2_vert = W2-(nw2-nw2_w1);

            value = (K2b[(S_vert)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim3 + b * K2_dim4 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                     +K2b[(S_vert+1)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim3 + b * K2_dim4 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                     +K2b[(S_vert)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim3 + b * K2_dim4 + (c-1)]*(ffreqs[(nw-nw2)/2+S+1]-s)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                     +K2b[(S_vert+1)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim3 + b * K2_dim4 + (c-1)]*(-ffreqs[(nw-nw2)/2+S]+s)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+S+1]-ffreqs[(nw-nw2)/2+S])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2]));
        #elif sym==1
            site x = site_switch(a,b,c);
            value = conj(K2_vvalsmooth(x.a,x.b,x.c, -s, w2 ));
        #elif sym==2
            site x = site_switch(a,b,c);
            value = K2_vvalsmooth(x.a,x.b,x.c, -s, w2 );
        #endif
}
            else{
            if(abs(s)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            int i,j;
            i= fconv_n(s,nw2);
            j = fconv_n(w2,nw2);

            value = K2b_vval(a,b,c,i,j);};};



            return value;
}
            //non-member functions
            svert operator*(double alpha, const svert& vertex){
            svert vertex2;


            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha* vertex.K1[i];
};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif

            // R contributions

            for(int i=0 ;i<K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}
            svert operator*(const svert& vertex,double alpha){
            svert vertex2;


            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha * vertex.K1[i];

};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif


            // R contributions

            for(int i=0 ;i< K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}
            svert operator+(const svert& vertex1,const svert& vertex2){
            svert vertex3;




            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i] = vertex1.K1[i] +  vertex2.K1[i];

};

            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i] = vertex1.K2[i] +  vertex2.K2[i];

        #if sym==0
            vertex3.K2b[i] = vertex1.K2b[i] +  vertex2.K2b[i];

        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i] = vertex1.K3[i] +  vertex2.K3[i];

};

            return vertex3;
}

            svert abs_sum_tiny(const svert& vertex1,const svert& vertex2, double tiny){
            svert vertex3;


            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i] = abs(vertex1.K1[i]) +  abs(vertex2.K1[i])+tiny;

};

            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i] = abs(vertex1.K2[i]) +  abs(vertex2.K2[i])+tiny;

        #if sym==0
            vertex3.K2b[i] = abs(vertex1.K2b[i]) +  abs(vertex2.K2b[i])+tiny;


        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i] = abs(vertex1.K3[i]) +  abs(vertex2.K3[i])+tiny;

};
            return vertex3;
}


            /*****************FUNCTIONS FOR THE T-VERTEX********************************************/


            void tvert::R_direct_set(int i, double value){
                if(i>=0 && i<K3_dim1){
                    K3[i] = value;
                };
            }
            void tvert::K1_direct_set(int i, double value){
                if(i>=0 && i<K1_dim1){
                    K1[i] = value;
                };
            }
            void tvert::K2_direct_set(int i, double value){
                if(i>=0 && i<K2_dim1){
                    K2[i] = value;
                };
            }
            void tvert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0 /*&&abs(value) > 1e-16*/){

                    K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]  = value;
};
}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
            void tvert::K1_setvert(int a, int b, int c,int i, double value){
            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw1 && i>=0  /*&&abs(value) > 1e-16*/){

            K1[(i-(nw1-nw1_q))*K1_dim2 + (a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)] = value;
};

}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;};
}
            void tvert::K2_setvert(int a, int b, int c, int i, int j, double value){
            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw2 && i>= 0 && j< nw2 && j>= 0 /*&&abs(value) > 1e-16*/){

            K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]  = value ;
};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
        #if sym==0
            void tvert::K2b_setvert(int a, int b, int c, int i, int j,  double value){
            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw2 && i>= 0 && j< nw2 && j>= 0 /*&&abs(value) > 1e-16*/){


            K2b[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]  = value;};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
            void tvert::K2b_direct_set(int i, double value){
                if(i>=0 && i<K2_dim1){
                    K2b[i] = value;
                };
            }
        #endif


            void tvert::K1_copy(vector<double> vec){
            if(vec.size()== K1_dim1){
            for(int i=0; i<K1_dim1; i++){
            K1[i] = vec[i];};
}
            else{cout << "error: K1-Vector cannot be copied due to mismatching dimensions."<< endl;};
}

            void tvert::K2_copy(vector<double> vec){
            if(vec.size()== K2_dim1){
            for(int i=0; i<K2_dim1; i++){
                   K2[i] = vec[i];
};
}
            else{cout << "error: K2-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
        #if sym==0
            void tvert::K2b_copy(vector<double> vec){
            if(vec.size()== K2_dim1){
            for(int i=0; i<K2_dim1; i++){
            K2[i] = vec[i];};
}
            else{cout << "error: K2b-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
        #endif
            void tvert::K3_copy(vector<double> vec){
            if(vec.size()== K3_dim1){
            for(int i=0; i<K3_dim1; i++){
                 K3[i] = vec[i];
};
}
            else{cout << "error: K3-Vector cannot be copied due to mismatching dimensions."<< endl;};
}



            double tvert::R_vval(int a_raw, int b_raw,int c_raw,  int i, int j, int k){


            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        #if sym==0

            value =  K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] ;
        #elif sym==1

            if(i >= nw3/2){
            if(j >= nw3/2 ){

            value =  K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] ;}

            else if(i >= nw3/2  && j < nw3/2 ){

            j = nw3-1-j;
            k = nw3-1-k;
            value = conj( K3[ (i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] );};}


            else if(i < nw3/2){
            if(k >= nw3/2 ){
            site x = site_switch(a,b,c);

            i = nw3-1-i;
            value =  K3[ (i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] ;}//note that k and j are interchanged

            else if( k < nw3/2 ){
            site x = site_switch(a,b,c);

            i = nw3-1-i;
            j = nw3-1-j;
            k = nw3-1-k;
            value = conj(K3[ (i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] );};//note that k and j are interchanged
};
        #elif sym==2
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;x = site_switch(x.a,x.b,x.c);};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;};
            value =  K3[ (i_eff-(nw3-nw3_q))*K3_dim2 + (j_eff-(nw3-nw3_w1)) * K3_dim3 +(k_eff-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1) ];
        #endif
};




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
        #if sym==0

            value =  K1[(i-(nw1-nw1_q))*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)  ] ;
        #elif sym==1 || sym==2
            if (i>=nw1/2){

            value =  K1[(i-(nw1-nw1_q))*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)  ] ;}
            else if(i<nw1/2){

            site x = site_switch(a,b,c);

            i = nw1-1-i;
            value =   K1[(i-(nw1-nw1_q))*K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1)  ] ;};
        #endif
};




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
        #if sym==0

            value = K2[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;
        #elif sym==1
            if(j>=nw2/2){

            value = K2[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;}

            else if(j < nw2/2){


            j = nw2-1-j;
            value = conj(K2[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] );
};
        #elif sym==2
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;};
            if(j < nw2/2){j = nw2-1-j;};

            value = K2[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)] ;
        #endif
};


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
        #if sym==0

            value = K2b[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;
        #elif sym==1 || sym==2

            site x = site_switch(a,b,c);
            value = K2_vval(x.a,x.b,x.c, nw2-1-i, j );
        #endif
};




            return value;
}
            double tvert::R_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w1, double w2){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(t) < bfreqs[(nw3+nw)/2-1] && abs(w1)<ffreqs[(nw3+nw)/2-1]&& abs(w2) <ffreqs[(nw3+nw)/2-1]){
        #if sym==0


            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert =W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);




            value = ((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        #elif sym==1
            if(t > 0){
            if((w1)>0 ){

            int T = fconv_n(t,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert =W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);

            value = ((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}

            else if( (w1)<0){

            int T = fconv_n(t,nw3);
            int W1 = fconv_n(-w1,nw3);
            int W2 = fconv_n(-w2,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = -w1 ;
            double w2_eff = -w2 ;
            value = conj((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
};}


            else if (t < 0  ){
            if((w2)>0 ){
            site x = site_switch(a,b,c);


            int T = fconv_n(-t,nw3);
            int W1 = fconv_n(w2,nw3);
            int W2 = fconv_n(w1,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = w2;
            double w2_eff = w1;

            t = -t;

            value = ((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}



            else if ((w2)<0 ){
            site x = site_switch(a,b,c);


            int T = fconv_n(-t,nw3);
            int W1 = fconv_n(-w2,nw3);
            int W2 = fconv_n(-w1,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);

            double w1_eff = -w2;
            double w2_eff = -w1;

            t = -t;

            value = conj((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                          +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
};};
        #elif sym ==2
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, t_eff =t;
            if(t<0){t_eff = -t; };
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1; x = site_switch(x.a,x.b,x.c);};
            if((w1_eff)<0){w1_eff = -w1_eff; w2_eff = -w2_eff;};


            int T = fconv_n(t_eff,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int T_vert = T-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);

            value = ((K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(ffreqs[(nw-nw3)/2+T+1]-t_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                      +K3[(T_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] *(-ffreqs[(nw-nw3)/2+T]+t_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+T+1]-ffreqs[(nw-nw3)/2+T])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));

        #endif

}
            else{
            int i,j,k;
            if(abs(t)<= bfreqs[(nw1+nw3)/2-1] && abs(w1)<= ffreqs[(nw1+nw3)/2-1]&& abs(w2)<= ffreqs[(nw1+nw3)/2-1]){
            i= fconv_n(t,nw3);
            j = fconv_n(w1,nw3);
            k = fconv_n(w2,nw3);

            value = R_vval(a,b,c,i,j,k);};};




            return value;
}
            double tvert::K1_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t){



            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(t) < bfreqs[(nw1+nw)/2-1]){
        #if sym==0

            int T = fconv_n(t,nw1);

            int T_vert = T-(nw1-nw1_q);
            value = ((K1[ (T_vert)* K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1) ]*(ffreqs[(nw-nw1)/2+T+1]-t)
                      +K1[(T_vert+1)* K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1) ]*(-ffreqs[(nw-nw1)/2+T]+t))/((ffreqs[(nw-nw1)/2+T+1]-ffreqs[(nw-nw1)/2+T])));
        #else
            if (t>0){
            int T = fconv_n(t,nw1);
            int T_vert = T-(nw1-nw1_q);
            value = ((K1[ (T_vert)* K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1) ]*(bfreqs[T+1]-t)
                      +K1[ (T_vert+1)* K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1) ]*(-bfreqs[T]+t))/((bfreqs[T+1]-bfreqs[T])));}
            else if (t<0){

            site x = site_switch(a,b,c);
            double t_eff = -t;

            int T = fconv_n(t_eff,nw1);
            int T_vert = T-(nw1-nw1_q);

            value = ((K1[ (T_vert)* K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1) ]*(bfreqs[T+1]-t_eff)
                      +K1[ (T_vert+1)* K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1) ]*(-bfreqs[T]+t_eff))/(bfreqs[T+1]-bfreqs[T]));
};
        #endif
}
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
            if(abs(t) < ffreqs[(nw2+nw)/2-1] && abs(w1) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0
            int T = fconv_n(t,nw2);
            int W1 = fconv_n(w1,nw2);

            int T_vert =T-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(T_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(T_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                     ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        #elif sym==1
            if(w1 > 0){

            int T = fconv_n(t,nw2);
            int W1 = fconv_n(w1,nw2);

            int T_vert = T-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(T_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(T_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1))/
                     ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));}
            else if (w1<0){

            double w1_eff = - w1;
            int T = fconv_n(t,nw2);
            int W1 =  fconv_n(w1_eff,nw2);

            int T_vert =T-(nw2-nw2_q);
            int W1_vert =  W1-(nw2-nw2_w1);



            value = conj((K2[(T_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                          +K2[(T_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                          +K2[(T_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                          +K2[(T_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                         ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
};
        #elif sym==2
            double w1_eff = w1;
            double t_eff = t;
            if(t<0){t_eff = -t; };
            if(w1_eff<0){w1_eff = -w1_eff;};

            int T = fconv_n(t_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int T_vert = T-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);
            value = ((K2[(T_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(T_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                      +K2[(T_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/
                     ((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        #endif
}
            else{ int i,j;
            if(abs(t)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(t,nw2);
            j = fconv_n(w1,nw2);

            value = K2_vval(a,b,c,i,j);};
};




            return value;
}
            double tvert::K2b_vvalsmooth(int a_raw, int b_raw,int c_raw,  double t, double w2){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(t) < ffreqs[(nw2+nw)/2-1] && abs(w2) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0
            int T = fconv_n(T,nw2);
            int W2 = fconv_n(w2,nw2);

            int T_vert =T-(nw2-nw2_q);
            int W2_vert = W2-(nw2-nw2_w1);

            value = ((K2b[(T_vert)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                      +K2b[(T_vert+1)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                      +K2b[(T_vert)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+T+1]-t)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                      +K2b[(T_vert+1)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+T]+t)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+T+1]-ffreqs[(nw-nw2)/2+T])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2])));
        #elif sym==1 ||sym==2
            site x = site_switch(a,b,c);

            value = K2_vvalsmooth(x.a,x.b,x.c, -t, w2 );
        #endif
}
            else{int i,j;
            if(abs(t)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(t,nw2);
            j = fconv_n(w2,nw2);

            value = K2b_vval(a,b,c,i,j);};

};


            return value;
}
            //non-member functions
            tvert operator*(double alpha, const tvert& vertex){
            tvert vertex2;

            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha* vertex.K1[i];
};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif


            // R contributions

            for(int i=0 ;i< K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}






            tvert operator*(const tvert& vertex,double alpha){
            tvert vertex2;

            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha* vertex.K1[i];
};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif


            // R contributions

            for(int i=0 ;i< K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}

            tvert operator+(const tvert& vertex1,const tvert& vertex2){
            tvert vertex3;


            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i]= vertex1.K1[i] +  vertex2.K1[i];};


            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i] = vertex1.K2[i] +  vertex2.K2[i];

        #if sym==0
            vertex3.K2b[i] = vertex1.K2b[i] +  vertex2.K2b[i];

        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i]  = vertex1.K3[i] +  vertex2.K3[i];

};

            return vertex3;
}


            tvert abs_sum_tiny(const tvert& vertex1,const tvert& vertex2,double tiny){
            tvert vertex3;


            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i] = abs(vertex1.K1[i]) +  abs(vertex2.K1[i])+tiny;

};

            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i]  = abs(vertex1.K2[i]) +  abs(vertex2.K2[i])+tiny;

        #if sym==0
            vertex3.K2b[i] = abs(vertex1.K2b[i]) +  abs(vertex2.K2b[i])+tiny;

        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i] = abs(vertex1.K3[i]) +  abs(vertex2.K3[i])+tiny;
};

            return vertex3;
}


            /*****************FUNCTIONS FOR THE U-VERTEX********************************************/

            void uvert::R_direct_set(int i, double value){
                if(i>=0 && i<K3_dim1){
                    K3[i] = value;
                };
            }
            void uvert::K1_direct_set(int i, double value){
                if(i>=0 && i<K1_dim1){
                    K1[i] = value;
                };
            }
            void uvert::K2_direct_set(int i, double value){
                if(i>=0 && i<K2_dim1){
                    K2[i] = value;
                };
            }
            void uvert::R_setvert(int a, int b, int c, int i, int j, int k, double value){
            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw3 && i>= 0 && j< nw3 && j>= 0 &&  k< nw3 && k>= 0 /*&&abs(value) > 1e-16*/){


                 K3[(i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] = value;
};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
            void uvert::K1_setvert(int a, int b, int c, int i, double value){

            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half
            if(i< nw1 && i>= 0  /*&&abs(value) > 1e-16*/){

            K1[(i-(nw1-nw1_q))*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)]  = value;
};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
            void uvert::K2_setvert(int a, int b, int c,int i, int j,double value){


            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c ){//only save sites in upper half
            if(i< nw2 && i>= 0 && j< nw2 && j>= 0  /*&&abs(value) > 1e-16*/){

                   K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] = value;
};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}
        #if sym==0
            void uvert::K2b_setvert(int a, int b, int c, int i, int j,  double value){


            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c ){//only save sites in upper half
            if(i< nw2 && i>= 0 && j< nw2 && j>= 0  /*&&abs(value) > 1e-16*/){

            K2b[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1) ] = value;};}
            else{cout << "error: Site (" << a<< " , " <<b<<" , " <<c << ") cannot be written since it lies in the lower half plane.."<< endl;}
}

            void uvert::K2b_direct_set(int i, double value){
                if(i>=0 && i<K2_dim1){
                    K2b[i] = value;
                };
            }
        #endif


            void uvert::K1_copy(vector<double> vec){
            if(vec.size()== K1_dim1){
            for(int i=0; i<K1_dim1; i++){
            K1[i] = vec[i];};
}
            else{cout << "error: K1-Vector cannot be copied due to mismatching dimensions."<< endl;};
}

            void uvert::K2_copy(vector<double> vec){
            if(vec.size()== K2_dim1){
            for(int i=0; i<K2_dim1; i++){
                   K2[i] = vec[i];
};
}
            else{cout << "error: K2-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
        #if sym==0
            void uvert::K2b_copy(vector<double> vec){
            if(vec.size()== K2_dim1){
            for(int i=0; i<K2_dim1; i++){
            K2b[i] = vec[i];};
}
            else{cout << "error: K2b-Vector cannot be copied due to mismatching dimensions."<< endl;};
}
        #endif
            void uvert::K3_copy(vector<double> vec){
            if(vec.size()== K3_dim1){
            for(int i=0; i<K3_dim1; i++){
             K3[i] = vec[i];
};
}
            else{cout << "error: K3-Vector cannot be copied due to mismatching dimensions."<< endl;};
}




            double uvert::R_vval(int a_raw, int b_raw, int c_raw,  int i, int j, int k){


            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(i < nw3 && i >= 0 && j < nw3 && j >= 0 && k < nw3 && k >= 0){
        #if sym==0
            value = K3[(i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] ;
        #elif sym==1
            if(i>=nw3/2){
            if(j >=nw3/2){

            value = K3[(i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)] ;}

            else if(j < nw3/2){
            site x = site_switch(a,b,c);

            j = nw3-1-j;
            k = nw3-1-k;
            value = conj(K3[(i-(nw3-nw3_q))*K3_dim2 + (j-(nw3-nw3_w1)) * K3_dim3 +(k-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]);};}

            else if (i<nw3/2 ){
            if(k >= nw3/2){
            site x = site_switch(a,b,c);

            i = nw3-1-i;
            value = K3[(i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] ;}//note that k and j are interchanged

            else if ( k < nw3/2){


            i = nw3-1-i;
            j = nw3-1-j;
            k = nw3-1-k;
            value = conj(K3[(i-(nw3-nw3_q))*K3_dim2 + (k-(nw3-nw3_w1)) * K3_dim3 +(j-(nw3-nw3_w2))*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]  );};//note that k and j are interchanged

};
        #elif sym==2
            site x(a,b,c);
            int i_eff = i, j_eff = j, k_eff = k;
            if(i<nw3/2){i_eff = nw3-1-i;x = site_switch(x.a,x.b,x.c);};
            if(abs(j-nw3/2) > abs(k-nw3/2)){j_eff = k; k_eff = j;};
            if(j_eff-nw3/2 < 0){j_eff = nw3-1-j_eff; k_eff = nw3-1-k_eff;x = site_switch(x.a,x.b,x.c);};
            value = K3[ (i_eff-(nw3-nw3_q))*K3_dim2 + (j_eff-(nw3-nw3_w1)) * K3_dim3 +(k_eff-(nw3-nw3_w2))*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)] ;
        #endif
};




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
        #if sym==0

            value = K1[(i-(nw1-nw1_q))*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)];
        #elif sym==1 || sym==2
            if (i>=nw1/2){

            value =K1[(i-(nw1-nw1_q))*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)];}
            else if(i<nw1/2){

            site x = site_switch(a,b,c);
            i=nw1-1-i;
            value =K1[(i-(nw1-nw1_q))*K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1)];};
        #endif
};


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
        #if sym==0
            value =  K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;
        #elif sym==1
            if(j >= nw2/2){


            value =  K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;}
            else if (j<nw2/2){

            site x = site_switch(a,b,c);

            j = nw2-1-j;
            value = conj( K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]);};
        #elif sym==2
            site x(a,b,c);
            if(i < nw2/2){i = nw2-1-i;x = site_switch(x.a,x.b,x.c);};
            if(j < nw2/2){j = nw2-1-j;x = site_switch(x.a,x.b,x.c);};

            value =  K2[(i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)] ;
        #endif
};


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
        #if sym==0

            value =  K2b[ (i-(nw2-nw2_q))*K2_dim2 + (j-(nw2-nw2_w1))*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)] ;
        #elif sym==1 || sym==2
            site x = site_switch(a,b,c);
            value = K2_vval(x.a,x.b,x.c, nw2-1-i, j );
        #endif
};

            return value;
}
            double uvert::R_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1, double w2){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(u) < ffreqs[(nw3+nw)/2-1] && abs(w1) <ffreqs[(nw3+nw)/2-1]&& abs(w2)<ffreqs[(nw3+nw)/2-1]){
        #if sym==0

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);


            int U_vert =U-(nw3-nw3_q);
            int W1_vert =W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);



            value = (( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        #elif sym==1
            if(u > 0 ){
            if( (w1)>0){

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(w1,nw3);
            int W2 = fconv_n(w2,nw3);


            int U_vert = U-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert =W2-(nw3-nw3_w2);

            value = (( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(ffreqs[(nw-nw3)/2+W2+1]-w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1)*(-ffreqs[(nw-nw3)/2+W2]+w2))/
                     ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}
            else if((w1)<0){
            site x = site_switch(a,b,c);

            int U = fconv_n(u,nw3);
            int W1 = fconv_n(-w1,nw3);
            int W2 = fconv_n(-w2,nw3);

            int U_vert = U-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = -w1 ;
            double w2_eff = -w2 ;
            value = conj(( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
};}

            else if (u < 0 ){
            if((w2)>0){
            site x = site_switch(a,b,c);

            int U = fconv_n(-u,nw3);
            int W1 = fconv_n(w2 ,nw3);
            int W2 = fconv_n(w1 ,nw3);


            int U_vert = U-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert =W2-(nw3-nw3_w2);
            double w1_eff = w2;
            double w2_eff = w1;
            u = -u;

            value = (( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
}


            else if ((w2)<0){


            int U = fconv_n(-u,nw3);
            int W1 = fconv_n(-w2 ,nw3);
            int W2 = fconv_n(-w1 ,nw3);

            int U_vert = U-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert = W2-(nw3-nw3_w2);
            double w1_eff = -w2;
            double w2_eff = -w1;
            u = -u;

            value = conj(( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                           +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(a+(nuc_eff-1)/2)*K3_dim5 + b * K3_dim6 + (c-1)]*(-ffreqs[(nw-nw3)/2+U]+u)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                         ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
};

};
        #elif sym==2
            site x(a,b,c);
            double w1_eff = w1, w2_eff = w2, u_eff =u;
            if(u<0){u_eff = -u; x = site_switch(x.a,x.b,x.c);};
            if(abs(w1_eff) > abs(w2_eff)){w1_eff = w2; w2_eff = w1;};
            if(w1_eff < 0){w1_eff = -w1_eff; w2_eff = -w2_eff;x = site_switch(x.a,x.b,x.c);};



            int U = fconv_n(u_eff,nw3);
            int W1 = fconv_n(w1_eff,nw3);
            int W2 = fconv_n(w2_eff,nw3);

            int U_vert = U-(nw3-nw3_q);
            int W1_vert = W1-(nw3-nw3_w1);
            int W2_vert =W2-(nw3-nw3_w2);
            value = (( K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(ffreqs[(nw-nw3)/2+W2+1]-w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u_eff)*(ffreqs[(nw-nw3)/2+W1+1]-w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(ffreqs[(nw-nw3)/2+U+1]-u_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff)
                       +K3[(U_vert+1)*K3_dim2 + (W1_vert+1) * K3_dim3 +(W2_vert+1)*K3_dim4+(x.a+(nuc_eff-1)/2)*K3_dim5 + x.b * K3_dim6 + (x.c-1)]*(-ffreqs[(nw-nw3)/2+U]+u_eff)*(-ffreqs[(nw-nw3)/2+W1]+w1_eff)*(-ffreqs[(nw-nw3)/2+W2]+w2_eff))/
                     ((ffreqs[(nw-nw3)/2+U+1]-ffreqs[(nw-nw3)/2+U])*(ffreqs[(nw-nw3)/2+W1+1]-ffreqs[(nw-nw3)/2+W1])*(ffreqs[(nw-nw3)/2+W2+1]-ffreqs[(nw-nw3)/2+W2])));
        #endif
}
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
            if(abs(u) < ffreqs[(nw1+nw)/2-1]){
        #if sym==0
            int U = fconv_n(u,nw1);

            int U_vert = U-(nw1-nw1_q);


            value = (( K1[(U_vert)*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)]*(bfreqs[(nw-nw1)/2+U+1]-u)
                       + K1[(U_vert+1)*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));
        #elif sym==1 || sym==2
            if (u>0){


            int U = fconv_n(u,nw1);

            int U_vert = U-(nw1-nw1_q);

            value = (( K1[(U_vert)*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)]*(bfreqs[(nw-nw1)/2+U+1]-u)
                       + K1[(U_vert+1)*K1_dim2+(a+(nuc_eff-1)/2) * K1_dim3 + b * K1_dim4 + (c-1)]*(-bfreqs[(nw-nw1)/2+U]+u))/((bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U])));}
            else if (u<0){
            double u_eff = -u;

            site x = site_switch(a,b,c);
            int U = fconv_n(u_eff,nw1);

            int U_vert = U-(nw1-nw1_q);


            value = ( K1[(U_vert)*K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1)]*(bfreqs[(nw-nw1)/2+U+1]-u_eff)
                      + K1[(U_vert+1)*K1_dim2+(x.a+(nuc_eff-1)/2) * K1_dim3 + x.b * K1_dim4 + (x.c-1)]*(-bfreqs[(nw-nw1)/2+U]+u_eff))/(bfreqs[(nw-nw1)/2+U+1]-bfreqs[(nw-nw1)/2+U]);
};
        #endif
}
            else{
            int i;
            if(abs(u)<= bfreqs[nw1-1] ){
            i= fconv_n(u,nw1);

            value = K1_vval(a,b,c,i);};    };


            return value;
}
            double uvert::K2_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w1){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;

            double value=0;
            if(abs(u) < bfreqs[(nw2+nw)/2-1] && abs(w1) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1,nw2);

            int U_vert = U-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(U_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(U_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        #elif sym==1

            if(w1 > 0){

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(w1,nw2);

            int U_vert = U-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);

            value = ((K2[(U_vert)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1)
                      +K2[(U_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
}
            else if (w1<0){

            site x = site_switch(a,b,c);

            int U = fconv_n(u,nw2);
            int W1 = fconv_n(-w1,nw2);

            int U_vert =U-(nw2-nw2_q);
            int W1_vert = W1-(nw2-nw2_w1);
            double w1_eff = -w1 ;
            value = conj((K2[(U_vert)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                          +K2[(U_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                          +K2[(U_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                          +K2[(U_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));

};
        #elif sym==2
            site x(a,b,c);
            double w1_eff = w1;
            double u_eff =u;
            if(u<0){u_eff = -u;x = site_switch(x.a,x.b,x.c); };
            if(w1_eff<0){w1_eff = -w1_eff; x = site_switch(x.a,x.b,x.c);};

            int U = fconv_n(u_eff,nw2);
            int W1 = fconv_n(w1_eff,nw2);

            int U_vert = U-(nw2-nw2_q);
            int W1_vert =W1-(nw2-nw2_w1);
            value = ((K2[(U_vert)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+U]+u_eff)*(ffreqs[(nw-nw2)/2+W1+1]-w1_eff)
                      +K2[(U_vert)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff)
                      +K2[(U_vert+1)*K2_dim2 + (W1_vert+1)*K2_dim3+(x.a+(nuc_eff-1)/2)*K2_dim4 + x.b * K2_dim5 + (x.c-1)]*(-ffreqs[(nw-nw2)/2+U]+u_eff)*(-ffreqs[(nw-nw2)/2+W1]+w1_eff))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W1+1]-ffreqs[(nw-nw2)/2+W1])));
        #endif
}
            else{ int i,j;
            if(abs(u)<= bfreqs[(nw1+nw2)/2-1] && abs(w1)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(u,nw2);
            j = fconv_n(w1,nw2);

            value = K2_vval(a,b,c,i,j);};   };


            return value;
}
            double uvert::K2b_vvalsmooth(int a_raw, int b_raw, int c_raw,  double u, double w2){

            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            double value=0;
            if(abs(u) < ffreqs[(nw2+nw)/2-1] && abs(w2) < ffreqs[(nw2+nw)/2-1]){
        #if sym==0

            int U = fconv_n(u,nw2);
            int W2 = fconv_n(w2,nw2);

            int U_vert = U-(nw2-nw2_q);
            int W2_vert = W2-(nw2-nw2_w1);

            value = ((K2b[(U_vert)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                      +K2b[(U_vert+1)*K2_dim2 + (W2_vert)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(ffreqs[(nw-nw2)/2+W2+1]-w2)
                      +K2b[(U_vert)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(ffreqs[(nw-nw2)/2+U+1]-u)*(-ffreqs[(nw-nw2)/2+W2]+w2)
                      +K2b[(U_vert+1)*K2_dim2 + (W2_vert+1)*K2_dim3+(a+(nuc_eff-1)/2)*K2_dim4 + b * K2_dim5 + (c-1)]*(-ffreqs[(nw-nw2)/2+U]+u)*(-ffreqs[(nw-nw2)/2+W2]+w2))/((ffreqs[(nw-nw2)/2+U+1]-ffreqs[(nw-nw2)/2+U])*(ffreqs[(nw-nw2)/2+W2+1]-ffreqs[(nw-nw2)/2+W2])));
        #elif sym==1 || sym==2
            site x = site_switch(a,b,c);
            value = K2_vvalsmooth(x.a,x.b,x.c, -u, w2 );
        #endif
}

            else{ int i,j;
            if(abs(u)<= bfreqs[(nw1+nw2)/2-1] && abs(w2)<= ffreqs[(nw1+nw2)/2-1]){
            i= fconv_n(u,nw2);
            j = fconv_n(w2,nw2);

            value = K2b_vval(a,b,c,i,j);};};



            return value;
}
            //non-member functions
            uvert operator*(double alpha, const uvert& vertex){
            uvert vertex2;

            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha* vertex.K1[i];
};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif


            // R contributions

            for(int i=0 ;i< K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}


            uvert operator*(const uvert& vertex,double alpha){
            uvert vertex2;

            //K1 contributions

            for( int i=0 ; i<K1_dim1 ; i++){
            vertex2.K1[i] = alpha* vertex.K1[i];
};
            //K2 contributions

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2[i] = alpha * vertex.K2[i];
};

        #if sym==0

            for(int i=0;i<K2_dim1; i++){

            vertex2.K2b[i] = alpha * vertex.K2b[i];
};
        #endif


            // R contributions

            for(int i=0 ;i< K3_dim1; i++){

            vertex2.K3[i] = alpha * vertex.K3[i];};

            return vertex2;

}

            uvert operator+(const uvert& vertex1,const uvert& vertex2){
            uvert vertex3;



            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i] = vertex1.K1[i] +  vertex2.K1[i];
};

            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i] = vertex1.K2[i] +  vertex2.K2[i];

        #if sym==0
            vertex3.K2b[i]  = vertex1.K2b[i] +  vertex2.K2b[i];

        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i]= vertex1.K3[i] +  vertex2.K3[i];

};

            return vertex3;
}



            uvert abs_sum_tiny(const uvert& vertex1,const uvert& vertex2, double tiny){
            uvert vertex3;



            for(int i=0; i<K1_dim1; i++){
            vertex3.K1[i] = abs(vertex1.K1[i]) + abs( vertex2.K1[i])+tiny;

};

            //add K2 contributions

            for(int i=0; i<K2_dim1; i++){
            vertex3.K2[i]  = abs(vertex1.K2[i]) +  abs(vertex2.K2[i])+tiny;

        #if sym==0
            vertex3.K2b[i] =abs(vertex1.K2b[i]) +  abs(vertex2.K2b[i])+tiny;

        #endif
};
            //add R contributions

            for(int i=0; i<K3_dim1; i++){
            vertex3.K3[i] = abs(vertex1.K3[i]) +  abs(vertex2.K3[i])+tiny;
};
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
            return U_bare[(a+(nuc_eff-1)/2)*irred_dim2 +b*irred_dim3+(c-1)];}
            else{return 0;};}
            double irreducible::vvalsmooth(int a_raw, int b_raw, int c_raw){
            if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            //cout << a << " " << b << " " << c << endl;

            return U_bare[(a+(nuc_eff-1)/2)*irred_dim2 +b*irred_dim3+(c-1)];}

            else{return 0;};}
            double irreducible::vvalsmooth(int a_raw, int b_raw, int c_raw,double q, double w1, double w2, char channel, int par, char f){
            if(distance(a_raw,b_raw,c_raw) <= d_c){//cutoff distance
            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;
            return U_bare[(a+(nuc_eff-1)/2)*irred_dim2 +b*irred_dim3+(c-1)];}
            else{return 0;};}

            void irreducible::setvert(int a_raw, int b_raw, int c_raw,  double value){
            if(distance(a_raw,b_raw,c_raw)<=d_c ){
            site project = site_project(a_raw,b_raw,c_raw);//project onto upper hp if site xs located in lower hp
            int a,b,c;
            a = project.a;
            b = project.b;
            c = project.c;

            U_bare[(a+(nuc_eff-1)/2)*irred_dim2 +b*irred_dim3+(c-1)] = value;
}
            else{
            cout << "Cannot write irred: Out of range" << endl;
};
}

            void irreducible::direct_set(int i,  double value){
            if(i>=0&&i<irred_dim1 ){

            U_bare[i] = value;
}
            else{
            cout << "Cannot write irred: Out of range" << endl;
};
}


            void irreducible::copy(vector<double> vec){
            if(vec.size()==irred_dim1){
            U_bare = vec;
};
}


            //operators for irreducible vertex
            irreducible operator*(double alpha, const irreducible& vertex){
            irreducible result;
            for(int i=0; i<irred_dim1; i++){
            result.U_bare[i] = alpha * vertex.U_bare[i];};
            return result;}
            irreducible operator*(const irreducible& vertex,double alpha){
            irreducible result;
            for(int i=0; i<irred_dim1; i++){
            result.U_bare[i] = alpha * vertex.U_bare[i];};
            return result;}
            irreducible operator+(const irreducible& vertex1,const irreducible& vertex2){
            irreducible result;
            for(int i=0; i<irred_dim1; i++){
            result.U_bare[i] =  vertex1.U_bare[i] + vertex2.U_bare[i];};
            return result;}

            irreducible abs_sum_tiny(const irreducible& vertex1,const irreducible& vertex2, double tiny){
            irreducible result;
            for(int i=0; i<irred_dim1; i++){
            result.U_bare[i] =  abs(vertex1.U_bare[i]) + abs(vertex2.U_bare[i])+ tiny;};
            return result;}

            /*****************MEMBER FUNCTIONS OF SUSCEPTIBILITY CLASS (WRITE/ADD_WRITE/READ)********************************************/

            void Susc::write(int a ,int b,int c,double value){

            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){


            Sus[(a+(nuc_eff-1)/2)*irred_dim2+b*irred_dim3+(c-1)] = value;
}
         else{cout << "cannot write Susceptibility in lower half" << endl;};
}
            void Susc::add_write(int a ,int b,int c,double value){

           if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){
            Sus[(a+(nuc_eff-1)/2)*irred_dim2+b*irred_dim3+(c-1)] += value;}

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
            value = Sus[(a+(nuc_eff-1)/2)*irred_dim2+b*irred_dim3+(c-1)];};


            return value;
}

            double Susc::acc(int i){
            if(i>=0 && i<irred_dim1){
            return Sus[i];
}
            else{cout <<"Error: Trying to access susceptibility value beyond range" << endl;
return 0;};
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

            parvert<fullvert> operator+(parvert<fullvert> vertex1,parvert<fullvert> vertex2){
            parvert<fullvert> result;
            result.spinvertex = vertex1.spinvertex + vertex2.spinvertex;
            result.densvertex = vertex1.densvertex + vertex2.densvertex;
            return result;
}

            parvert<svert> abs_sum_tiny(parvert<svert> vertex1,parvert<svert> vertex2,double tiny){
            parvert<svert> result;
            result.spinvertex = abs_sum_tiny(vertex1.spinvertex, vertex2.spinvertex, tiny);
            result.densvertex = abs_sum_tiny(vertex1.densvertex, vertex2.densvertex, tiny);
            return result;
}
            parvert<tvert> abs_sum_tiny(parvert<tvert> vertex1,parvert<tvert> vertex2,double tiny){
            parvert<tvert> result;
            result.spinvertex = abs_sum_tiny(vertex1.spinvertex, vertex2.spinvertex, tiny);
            result.densvertex = abs_sum_tiny(vertex1.densvertex, vertex2.densvertex, tiny);
            return result;
}
            parvert<uvert> abs_sum_tiny(parvert<uvert> vertex1,parvert<uvert> vertex2,double tiny){
            parvert<uvert> result;
            result.spinvertex = abs_sum_tiny(vertex1.spinvertex, vertex2.spinvertex, tiny);
            result.densvertex = abs_sum_tiny(vertex1.densvertex, vertex2.densvertex, tiny);
            return result;
}
            parvert<irreducible> abs_sum_tiny(parvert<irreducible> vertex1,parvert<irreducible> vertex2,double tiny){
            parvert<irreducible> result;
            result.spinvertex = abs_sum_tiny(vertex1.spinvertex, vertex2.spinvertex, tiny);
            result.densvertex = abs_sum_tiny(vertex1.densvertex, vertex2.densvertex, tiny);
            return result;
}

            parvert<fullvert>abs_sum_tiny(parvert<fullvert> vertex1,parvert<fullvert> vertex2,double tiny){
            parvert<fullvert> result;
            result.spinvertex = abs_sum_tiny(vertex1.spinvertex, vertex2.spinvertex, tiny);
            result.densvertex = abs_sum_tiny(vertex1.densvertex, vertex2.densvertex, tiny);
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
template<typename T0>
            double fullvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2, char channel){

            double result = irred.vvalsmooth(a,b,c) + svertex.vvalsmooth(a,b,c,q,w1,w2,channel) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel);

            return  result;

}

template<typename T0>
            double fullvert::vvalsmooth(int a, int b, int c, T0 q, T0 w1, T0 w2, char channel, int p, char f){



            double result=0;

            if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
            result += irred.vvalsmooth(a,b,c);

}
            else if( p==2 && (f=='K' || f== 'M')){
            result += irred.vvalsmooth(a,b,c);

};


            result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);

            return result;
}
template<typename T0>
            double fullvert::vvalsmooth(int red_side, int map, int a, int b, int c, T0 q, T0 w1, T0 w2, char channel, int p, char f){// red_side: if only complementary channel in one vertex, which one is reduced? (0/1/2), p: is this the left/upper (1) or the right/lower (2) vertex of the bubble?, f: diagrammatic class that is computed

            double result=0;
            //if vertex is mapped, also frequency arguments must be adjusted (1' and 2' are interchanged): for the p-channel, this results in: w2 -> -w2, for a- and t-channel, this just interchanges the channel arguements (t<->a)

            if(red_side != p){
            if( p==1 && (f=='K' || f== 'L')){//only yield irred part if both legs are connected to the same bare vertex
            result = irred.vvalsmooth(a,b,c);
}
            else if( p==2 && (f=='K' || f== 'M')){
            result = irred.vvalsmooth(a,b,c);
};
            if(map==0){
            result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,channel,p,f);
}
            else if(map==1){
            if(channel =='s'){
            //if channels are mapped,  the arguments must be adjusted appropriately
            result +=  svertex.vvalsmooth(a,b,c,q,w1,-w2,'s',p,f) + tvertex.vvalsmooth(a,b,c,q,w1,-w2,'s',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,-w2,'s',p,f);
}
            else if(channel =='t'){
            result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f);
}
            else if(channel =='u'){
            result +=  svertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f);
};
};
}

            else if(red_side == p){

            if(map==0){
            if(channel == 's'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,'s',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,'s',p,f);}
            else if(channel == 't'){  result =  svertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f);}
            else if(channel == 'u'){  result =  svertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f) + tvertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f);};}
            else if(map==1){

            if(channel == 's' ){  result =  tvertex.vvalsmooth(a,b,c,q,w1,-w2,'s',p,f) + uvertex.vvalsmooth(a,b,c,q,w1,-w2,'s',p,f);}
            else if(channel == 't'){  result =  tvertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f) +  svertex.vvalsmooth(a,b,c,q,w1,w2,'u',p,f) ;}
            else if(channel == 'u'){  result =  uvertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f) +  svertex.vvalsmooth(a,b,c,q,w1,w2,'t',p,f) ;}
};

};


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
            fullvert abs_sum_tiny( const fullvert& vertex1, const fullvert& vertex2, double tiny){
            fullvert result;
            result.irred = abs_sum_tiny(vertex1.irred, vertex2.irred, tiny);
            result.svertex =  abs_sum_tiny(vertex1.svertex , vertex2.svertex, tiny);
            result.tvertex =  abs_sum_tiny(vertex1.tvertex , vertex2.tvertex, tiny);
            result.uvertex =  abs_sum_tiny(vertex1.uvertex , vertex2.uvertex , tiny);
            return result;
}







            /*****************************************FUNCTIONS FOR SELF ENERGY********************************************************/
            double self::sval(int i){
            double value = 0;
            if(i>=0 && i< nw){
            value = selfenergy[i];};
            return value;
}
            void self::self_copy(vector<double> vec){
            if(vec.size()==nw1){
            for(int i=0; i<nw1; i++){
            selfenergy[i] = vec[i];};
}
            else{cout << "error: Self energy cannot be copied due to mismatching dimensions."<< endl;};
}



            void self::setself(int i, double val){
            if(i>=0 && i< nw && abs(val) > 1e-16){

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
            double sum = self1.selfenergy[i] + self2.selfenergy[i] ;

            self3.selfenergy[i] = sum ;

};
            return self3;
}

            self abs_sum_tiny(const self& self1, const self& self2, double tiny){//sum operator overloading
            self self3;
            for(int i=0; i<nw; i++){
            self3.selfenergy[i] = abs(self1.selfenergy[i]) + abs(self2.selfenergy[i])+tiny ;
};
            return self3;
}
            self operator+=(const self& self1, const self& self2){//sum operator overloading
            self self3;
            for(int i=0; i<nw; i++){
            double sum = self1.selfenergy[i] + self2.selfenergy[i] ;

            self3.selfenergy[i] = sum ;};
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

            state abs_sum_tiny(state state1, state state2, double tiny){
            state result;
            result.vertex.spinvertex.irred = abs_sum_tiny(state1.vertex.spinvertex.irred ,state2.vertex.spinvertex.irred, tiny);
            result.vertex.spinvertex.svertex = abs_sum_tiny(state1.vertex.spinvertex.svertex , state2.vertex.spinvertex.svertex, tiny);
            result.vertex.spinvertex.tvertex = abs_sum_tiny(state1.vertex.spinvertex.tvertex , state2.vertex.spinvertex.tvertex, tiny);
            result.vertex.spinvertex.uvertex = abs_sum_tiny(state1.vertex.spinvertex.uvertex , state2.vertex.spinvertex.uvertex, tiny);
            result.vertex.densvertex.irred = abs_sum_tiny(state1.vertex.densvertex.irred , state2.vertex.densvertex.irred, tiny);
            result.vertex.densvertex.svertex = abs_sum_tiny(state1.vertex.densvertex.svertex, state2.vertex.densvertex.svertex, tiny);
            result.vertex.densvertex.tvertex = abs_sum_tiny(state1.vertex.densvertex.tvertex , state2.vertex.densvertex.tvertex, tiny);
            result.vertex.densvertex.uvertex = abs_sum_tiny(state1.vertex.densvertex.uvertex , state2.vertex.densvertex.uvertex, tiny);
            result.selfenergy =abs_sum_tiny( state1.selfenergy , state2.selfenergy, tiny);
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


            double max_err(state& err, state& scale){
            int sum_K1 = (sym==0?(nw-nw1)/2:nw/2);

            int sum_K2_i;
            if(sym==0 || sym==1){sum_K2_i = (nw-nw2)/2;}
            else if(sym==2){sum_K2_i = nw/2;};

            int sum_K2_j;
            if(sym==0){sum_K2_j = (nw-nw2)/2;}
            else if(sym==1 || sym==2){sum_K2_j = nw/2;};

            int sum_R = (sym==0?(nw-nw3)/2:nw/2);

            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            double max_error1, max_error2,max_error3,max_error4=0.0;

            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= 0; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4; c++){

            if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half


            for(int i=sum_K1 ; i< (nw+nw1)/2; i++){
            //spin Vertex:
            double s_rel_err_spin=(abs(scale.vertex.spinvertex.svertex.K1_vval(a,b,c,i))>1e-6? err.vertex.spinvertex.svertex.K1_vval(a,b,c,i)/scale.vertex.spinvertex.svertex.K1_vval(a,b,c,i) : 0) ;
            double t_rel_err_spin=(abs(scale.vertex.spinvertex.tvertex.K1_vval(a,b,c,i))>1e-6? err.vertex.spinvertex.tvertex.K1_vval(a,b,c,i)/scale.vertex.spinvertex.tvertex.K1_vval(a,b,c,i) : 0) ;
            double u_rel_err_spin=(abs(scale.vertex.spinvertex.uvertex.K1_vval(a,b,c,i))>1e-6? err.vertex.spinvertex.uvertex.K1_vval(a,b,c,i)/scale.vertex.spinvertex.uvertex.K1_vval(a,b,c,i) : 0) ;
            //dens Vertex:
            double s_rel_err_dens=(abs(scale.vertex.densvertex.svertex.K1_vval(a,b,c,i))>1e-6? err.vertex.densvertex.svertex.K1_vval(a,b,c,i)/scale.vertex.densvertex.svertex.K1_vval(a,b,c,i) : 0) ;
            double t_rel_err_dens=(abs(scale.vertex.densvertex.tvertex.K1_vval(a,b,c,i))>1e-6? err.vertex.densvertex.tvertex.K1_vval(a,b,c,i)/scale.vertex.densvertex.tvertex.K1_vval(a,b,c,i) : 0) ;
            double u_rel_err_dens=(abs(scale.vertex.densvertex.uvertex.K1_vval(a,b,c,i))>1e-6? err.vertex.densvertex.uvertex.K1_vval(a,b,c,i)/scale.vertex.densvertex.uvertex.K1_vval(a,b,c,i) : 0) ;

            max_error1 = max({abs(max_error1),abs(s_rel_err_spin),abs(t_rel_err_spin),abs(u_rel_err_spin),abs(s_rel_err_dens),abs(t_rel_err_dens),abs(u_rel_err_dens)},comp);


};


            for(int i=sum_K2_i ; i<(nw+nw2)/2; i++){
            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
            //spin Vertex:
            //spin Vertex:
            double s_rel_err_spin=(abs(scale.vertex.spinvertex.svertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.svertex.K2_vval(a,b,c,i,j)/scale.vertex.spinvertex.svertex.K2_vval(a,b,c,i,j) : 0) ;
            double t_rel_err_spin=(abs(scale.vertex.spinvertex.tvertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.tvertex.K2_vval(a,b,c,i,j)/scale.vertex.spinvertex.tvertex.K2_vval(a,b,c,i,j) : 0) ;
            double u_rel_err_spin=(abs(scale.vertex.spinvertex.uvertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.uvertex.K2_vval(a,b,c,i,j)/scale.vertex.spinvertex.uvertex.K2_vval(a,b,c,i,j) : 0) ;
            //dens Vertex
            double s_rel_err_dens=(abs(scale.vertex.densvertex.svertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.svertex.K2_vval(a,b,c,i,j)/scale.vertex.densvertex.svertex.K2_vval(a,b,c,i,j) : 0) ;
            double t_rel_err_dens=(abs(scale.vertex.densvertex.tvertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.tvertex.K2_vval(a,b,c,i,j)/scale.vertex.densvertex.tvertex.K2_vval(a,b,c,i,j) : 0) ;
            double u_rel_err_dens=(abs(scale.vertex.densvertex.uvertex.K2_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.uvertex.K2_vval(a,b,c,i,j)/scale.vertex.densvertex.uvertex.K2_vval(a,b,c,i,j) : 0) ;

            max_error2 = max({abs(max_error2),abs(s_rel_err_spin),abs(t_rel_err_spin),abs(u_rel_err_spin),abs(s_rel_err_dens),abs(t_rel_err_dens),abs(u_rel_err_dens)},comp);


        #if sym==0
            //spin Vertex:
            double s_rel_err_spin=(abs(scale.vertex.spinvertex.svertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.svertex.K2b_vval(a,b,c,i,j)/scale.vertex.spinvertex.svertex.K2b_vval(a,b,c,i,j) : 0) ;
            double t_rel_err_spin=(abs(scale.vertex.spinvertex.tvertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.tvertex.K2b_vval(a,b,c,i,j)/scale.vertex.spinvertex.tvertex.K2b_vval(a,b,c,i,j) : 0) ;
            double u_rel_err_spin=(abs(scale.vertex.spinvertex.uvertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.spinvertex.uvertex.K2b_vval(a,b,c,i,j)/scale.vertex.spinvertex.uvertex.K2b_vval(a,b,c,i,j) : 0) ;
            //dens Vertex
            double s_rel_err_dens=(abs(scale.vertex.densvertex.svertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.svertex.K2b_vval(a,b,c,i,j)/scale.vertex.densvertex.svertex.K2b_vval(a,b,c,i,j) : 0) ;
            double t_rel_err_dens=(abs(scale.vertex.densvertex.tvertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.tvertex.K2b_vval(a,b,c,i,j)/scale.vertex.densvertex.tvertex.K2b_vval(a,b,c,i,j) : 0) ;
            double u_rel_err_dens=(abs(scale.vertex.densvertex.uvertex.K2b_vval(a,b,c,i,j))>1e-6? err.vertex.densvertex.uvertex.K2b_vval(a,b,c,i,j)/scale.vertex.densvertex.uvertex.K2b_vval(a,b,c,i,j) : 0) ;
            max_error = max({abs(max_error),abs(s_rel_err_spin),abs(t_rel_err_spin),abs(u_rel_err_spin),abs(s_rel_err_dens),abs(t_rel_err_dens),abs(u_rel_err_dens)},comp);
        #endif


};};



            for(int i=sum_R ; i<(nw+nw3)/2; i++){
            for(int j=sum_R ; j<(nw+nw3)/2; j++){
            for(int k=(nw-nw3)/2 ; k<(nw+nw3)/2; k++){
            if(sym !=2 || (abs(j-nw/2) <= abs(k-nw/2))){
            double s_rel_err_spin=(abs(scale.vertex.spinvertex.svertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.spinvertex.svertex.R_vval(a,b,c,i,j,k)/scale.vertex.spinvertex.svertex.R_vval(a,b,c,i,j,k) : 0) ;
            double t_rel_err_spin=(abs(scale.vertex.spinvertex.tvertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.spinvertex.tvertex.R_vval(a,b,c,i,j,k)/scale.vertex.spinvertex.tvertex.R_vval(a,b,c,i,j,k) : 0) ;
            double u_rel_err_spin=(abs(scale.vertex.spinvertex.uvertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.spinvertex.uvertex.R_vval(a,b,c,i,j,k)/scale.vertex.spinvertex.uvertex.R_vval(a,b,c,i,j,k) : 0) ;
            //dens Vertex
            double s_rel_err_dens=(abs(scale.vertex.densvertex.svertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.densvertex.svertex.R_vval(a,b,c,i,j,k)/scale.vertex.densvertex.svertex.R_vval(a,b,c,i,j,k) : 0) ;
            double t_rel_err_dens=(abs(scale.vertex.densvertex.tvertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.densvertex.tvertex.R_vval(a,b,c,i,j,k)/scale.vertex.densvertex.tvertex.R_vval(a,b,c,i,j,k) : 0) ;
            double u_rel_err_dens=(abs(scale.vertex.densvertex.uvertex.R_vval(a,b,c,i,j,k))>1e-6? err.vertex.densvertex.uvertex.R_vval(a,b,c,i,j,k)/scale.vertex.densvertex.uvertex.R_vval(a,b,c,i,j,k) : 0) ;
            max_error3 = max({abs(max_error3),abs(s_rel_err_spin),abs(t_rel_err_spin),abs(u_rel_err_spin),abs(s_rel_err_dens),abs(t_rel_err_dens),abs(u_rel_err_dens)},comp);


};};};};

            //close loop over sites
};};};};

            for(int i=0; i<nw;i++){
            double self = (abs(scale.selfenergy.acc(i))>1e-6? err.selfenergy.acc(i)/scale.selfenergy.acc(i) :0);

            max_error4 = max({abs(max_error4),abs(self)},comp);
};
            double max_error = max({max_error1,max_error2,max_error3,max_error4},comp);
            if(world_rank==0 && max({max_error1,max_error2,max_error3,max_error4},comp)==max_error1){
            cout << "K1 biggest" << endl;
}
            else if(world_rank==0 && max({max_error1,max_error2,max_error3,max_error4},comp)==max_error2){
            cout << "K2 biggest" << endl;
}
            else if(world_rank==0 && max({max_error1,max_error2,max_error3,max_error4},comp)==max_error3){
            cout << "R biggest" << endl;
}
            else if(world_rank==0 && max({max_error1,max_error2,max_error3,max_error4},comp)==max_error4){
            cout << "selfbiggest" << endl;
};
            if(world_rank==0){
            cout << "Maximal error detected" << max_error << endl;
};
            return max_error;

}

            /*****************************************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/
            //ATTENTION: The propagator is purely imaginary so the double-value is the imaginary part
template<typename T0>
            double propag(double Lambda, T0 w, self selfenergy,self diffselfenergy, char type){
#if temp==0
            if(w!=0){
            double value=0;

            if(type == 'g'){//regular undifferentiated greensfunction
            double g0=0;
        #if reg==1
            if(abs(w)-1e-12 > Lambda){
            g0 = - 1./w;//minus sign since we compute the imaginary part
            value += - 1./(-1./g0-selfenergy.svalsmooth(w));//minus sign structure since we compute the imaginary part

}
            else if (abs(abs(w)-abs(Lambda)) <1e-12){
            g0 = -1./w;//minus sign since we compute the imaginary part
            value += -0.5/(-1./g0-selfenergy.svalsmooth(w));
            //this is the implementation of Morris Lemma, see SB.II, p. 9
}
        #elif reg==2
            g0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * (-1.)/w;//minus sign since we compute the imaginary part
            value = -1./(-1./g0-selfenergy.svalsmooth(w));
        #elif reg==3
            g0 = - 1./(w*(1. - exp(-pow(abs(w)/Lambda,sharp))) + Lambda * sgn(w) * exp(-pow(abs(w)/Lambda,sharp)));
            value = -1./(-1./g0-selfenergy.svalsmooth(w));
        #endif

}
            else if(type == 's'){//single scale propagator

        #if reg==1
            if(abs(w) == Lambda){

            value += -(-1./(w-selfenergy.svalsmooth(w)));}//one minus sign is from the differentiation of the theta-function, one is from the fact the we compute the imaginary poart of the propagator
            else{return 0;};
        #elif reg==2
            double G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * (-1.)/w; //bare greens function with smoothened cutoff at freq ffreqs[i]. minus sign since we compute the imaginary part
            double G = -1./( -1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
            value += pow((1. - G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * (-1.)/w;//see page 20 in SB II. In the first bracket, there is a minues sign since we multiply two purely imaginary quantities. Minues sign in last factor due to fact that we compute imag part of the propagator
        #elif reg==3
            double G0 = - 1./(w*(1. - exp(-pow(abs(w)/Lambda,sharp))) + Lambda * sgn(w) * exp(-pow(abs(w)/Lambda,sharp)));
            double G = -1./(-1./G0-selfenergy.svalsmooth(w));//this is the imaginary part of the propagator
            value = pow(G,2) * (-2*pow(w,3)/pow(Lambda,3)+sgn(w)+2* pow(w,2)/pow(Lambda,2)) * exp(-pow(abs(w)/Lambda,sharp));// (see SB3, p. 3).. in the overall sign, the (i)^2 that stems from twice the imaginary part of G is taken into account
        #endif
}
            else if(type == 'k'){//katanin substitution
        #if reg==1
            double S=0,G=0;
            if(abs(w) == Lambda){ S = -(- 1.)/(w-selfenergy.svalsmooth(w));//single scale propagator. //one minus sign is from the differentiation of the theate-function, one is from the fact the we compute the imaginary poart of the propagator
            G = -0.5/(w-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9. //minus sign since we compute the imaginary part
            else if(abs(w) > Lambda){ G = -1./(w-selfenergy.svalsmooth(w));};//minus sign since we compute the imaginary part

            value = (S - G * G * diffselfenergy.svalsmooth(w));//minus sign since we mulitply two purely imaginary greensfunctions G*G where G denotes the imaginary part of G
        #elif reg==2
            double G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * (-1.)/w; //bare greens function with smoothened cutoff at freq ffreqs[i]//minus sign since we compute the imaginary part
            double G = -1./( -1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
            value += pow((1. - G * selfenergy.svalsmooth(w)),2) * (-sharp/Lambda) * pow(abs(w)/Lambda,sharp) * exp(-pow(abs(w)/Lambda,sharp)) * (-1.)/w//In the first bracket, there is a minues sign since we multiply two purely imaginary quantities. Minues sign in last factor due to fact that we compute imag part of the propagator
            - G * G * diffselfenergy.svalsmooth(w);//see page 20 in SB II. //minus sign since we mulitply two purely imaginary greensfunctions G*G where G denotes the imaginary part of G
        #elif reg==3
            double G0 = - 1./(w*(1. - exp(-pow(abs(w)/Lambda,sharp))) + Lambda * sgn(w) * exp(-pow(abs(w)/Lambda,sharp)));
            double G = -1./(-1./G0-selfenergy.svalsmooth(w));//this is the imaginary part of the propagator
            value = pow(G,2) * (-2*pow(w,3)/pow(Lambda,3)+sgn(w)+2* pow(w,2)/pow(Lambda,2)) * exp(-pow(abs(w)/Lambda,sharp)) - G * G * diffselfenergy.svalsmooth(w);// (see SB3, p. 3).. in the overall sign, the (i)^2 that stems from twice the imaginary part of G is taken into account
        #endif

}
            else if(type == 'e'){//only self energy extension of katanin substitution ( needed for sharp regulator)
        #if reg==1
            double G=0;
            if(abs(w)-1e-12 > Lambda){ G = -1./(-w-selfenergy.svalsmooth(w));}
            else if(abs(abs(w)-abs(Lambda)) <1e-12) {G = -0.5/(-w-selfenergy.svalsmooth(w));}//this is the implementation of Morris Lemma, see SB.II, p. 9
            value = -(G * G * diffselfenergy.svalsmooth(w)); //there is a minus sign since we multiply two purely imaginary quantities.
        #elif reg==2
            double G0 = (1.-exp(-pow(abs(w)/Lambda,sharp))) * (-1.)/w; //bare greens function with smoothened cutoff at freq ffreqs[i]
            double G = -1./( -1./G0 - selfenergy.svalsmooth(w)); //full greens function at ffreqs[i]
            value = - G * G * diffselfenergy.svalsmooth(w);
        #elif reg==3
            double G0 = - 1./(w*(1. - exp(-pow(abs(w)/Lambda,sharp))) + Lambda * sgn(w) * exp(-pow(abs(w)/Lambda,sharp)));
            double G = -1./(-1./G0-selfenergy.svalsmooth(w));//this is the imaginary part of the propagator
            value = - G * G * diffselfenergy.svalsmooth(w);// (see SB3, p. 3).. in the overall sign, the (i)^2 that stems from twice the imaginary part of G is taken into account

        #endif


}
            else{cout << "could not resolve propagator type" << endl;};

            return value;
}
            else{return 0;};

        #elif temp==1

            double value=0;
            if((w-1)%2==0){
            if(type =='g'){
            double g0 = - 1./(w*pi*T);//minus sign since we compute the imaginary part
            if(abs(Lambda/(pi*T) - abs(w) )<= 1.){
            value += - 1./(-1./g0-selfenergy.svalsmooth(w)) * (abs(w*pi*T)-(Lambda-pi*T))/(2*pi*T);//minus sign structure since we compute the imaginary part
}
            else if(abs(w*pi*T)>Lambda){
            value += - 1./(-1./g0-selfenergy.svalsmooth(w));//minus sign structure since we compute the imaginary part
};
}
       else if (type=='s'){
            if(abs(Lambda/(pi*T) - abs(w) )<= 1.){
            value += -(-1./(w*pi*T-selfenergy.svalsmooth(w))) * sgn(w)/(2*pi*T);};//one minus sign is from the differentiation of the theta-function, one is from the fact the we compute the imaginary poart of the propagator

}
            else if (type=='k'){
            double G,S=0;
             double g0 = - 1./(w*pi*T);//minus sign since we compute the imaginary part
                 if(abs(Lambda/(pi*T) - abs(w) )<= 1.){
           G = - 1./(-1./g0-selfenergy.svalsmooth(w)) * (abs(w*pi*T)-(Lambda-pi*T))/(2*pi*T);}
            else if(abs(w*pi*T)>Lambda){
            G= - 1./(-1./g0-selfenergy.svalsmooth(w));
};

            if(abs(Lambda/(pi*T) - w )<= 1.){
            S += -(-1./(w*pi*T-selfenergy.svalsmooth(w))) * sgn(w)/(2*pi*T);};//one minus sign is from the differentiation of the theta-function, one is from the fact the we compute the imaginary poart of the propagator
             value = S - G * G * diffselfenergy.svalsmooth(w);
     }
            else if (type=='e'){
            double G=0;
             double g0 = - 1./(w*pi*T);//minus sign since we compute the imaginary part
                 if(abs(Lambda/(pi*T) - abs(w) )<= 1.){
           G = - 1./(-1./g0-selfenergy.svalsmooth(w)) * (abs(w*pi*T)-(Lambda-pi*T))/(2*pi*T);}
            else if(abs(w*pi*T)>Lambda){
            G= - 1./(-1./g0-selfenergy.svalsmooth(w));
};

                  value =  - G * G * diffselfenergy.svalsmooth(w);
     };
        }
    else{cout << "PROPAG ERROR " << w   << " " << (w-1)%2<< endl;
};
            return value;
        #endif
}






            /*****************************************FUNCTION TO COMPUTE SUSCEPTIBILITY********************************************************/
            Susc suscept(double Lambda, parvert<fullvert>& vertex, self selfenergy){
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            Susc result;

            //range of frequency integrations can be reduced by using symmetry relations. The lower frequency bounds for the vertices are:
#if temp==0
            gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (1500);
        #elif temp==1
            gsl_integration_workspace * w =NULL;
        #endif
            if(world_rank==0){
            cout << "SUS_WORLDRANK"<< world_size << endl;};
            int sum_K2_j;
            //    if(sym==0){sum_K2_j = (nw-nw2)/2;}
            //    else if(sym==1 || sym==2){sum_K2_j = nw/2;};

        #if sym==0
            sum_K2_j = (nw-nw2)/2;
        #elif sym==1 || sym==2
            sum_K2_j = nw/2;
        #endif


            double one =1.0;
            fullvert bare;
            bare.irred.setvert(0,0,1,one);

        #if temp==0
            double u = bfreqs[nw/2];
        #elif temp==1
            int u = 0;
        #endif
            double first = real(0.5*ububble(0,0,0,w,Lambda, bare,0,0,1,bare,0,0,1,'g','g',selfenergy,selfenergy,u,wlimit,wlimit,'K'));
            if(world_rank==0){
            cout << "first susc-contr. computed.." << endl;
};
            //second and third diagram;

            fullvert second_third_left;

            fullvert spin_dens =  vertex.densvertex + (-1)*vertex.spinvertex;



            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
        #if temp==0

             double w1 = ffreqs[j];
        #elif temp==1

            int w1 = -nw1+1+j*2;
        #endif

            second_third_left.uvertex.K2_setvert(0,0,1,nw2/2,j-(nw-nw2)/2, ububble(0,0,0,w,Lambda,bare,0,0,1, spin_dens,0,0,1,'g','g',selfenergy,selfenergy,u,w1,wlimit,'L'));

};



            double second_third = 0.5* ububble(0,0,0,w,Lambda,second_third_left,0,0,1,bare,0,0,1,'g','g',selfenergy,selfenergy,u,wlimit,wlimit,'K');


            double firstthree = first+second_third;
            if(abs(firstthree)<1e-100){firstthree =0;};
            result.write(0,0,1,firstthree); //write diagrams 1-3
            if(world_rank==0){
            cout << "second & third susc-contr. computed.." << endl;
};
            //fourth diagram:


            for(int a= -(nuc_eff-1)/2; a<(nuc_eff-1)/2+1; a++){
            for(int b= -(nuc_eff-1)/2; b<(nuc_eff-1)/2+1; b++){
            for(int c= 1; c<4; c++){



         if((b > 0 || (b==0 && a>=0) || (b==0 && a<0 && c==2)) && distance(a,b,c) <= d_c){//only save sites in upper half

            tvert fourth_left;
        #if temp==0
            double t = bfreqs[nw/2];
        #elif temp==1
            int t = 0;
        #endif

            fourth_left.K1_setvert(a,b,c,nw1/2, tbubble(0,0,0,w,Lambda,bare,0,0,1, vertex.spinvertex,a,b,c,'g','g',selfenergy,selfenergy,t,wlimit,wlimit,'K'));
            if(world_rank==0){
            cout << fourth_left.K1_vval(a,b,c,nw1/2) << endl;
};

            for(int j=sum_K2_j ; j<(nw+nw2)/2; j++){
        #if temp==0

             double w1 = ffreqs[j];
        #elif temp==1

            int w1 = -nw1+1+j*2;
        #endif
            fourth_left.K2_setvert(a,b,c,nw2/2,j-(nw-nw2)/2, tbubble(0,0,0,w,Lambda,bare,0,0,1, vertex.spinvertex,a,b,c,'g','g',selfenergy,selfenergy,t,w1,wlimit,'L'));
};

            result.add_write(a,b,c,(-1)*tbubble(0,0,0,w,Lambda, fourth_left,a,b,c,bare,0,0,1,'g','g',selfenergy,selfenergy,t,wlimit,wlimit,'K'));//factor (-1) is implicitly contained in definition of t-bubble but since two t-bubble are computed here, they cancel and one must indeed include factor (-1) from combinatorics
};};};};
            if(world_rank==0){
            cout << "fourth susc-contr. computed.." << endl;
};
        #if temp==0
                    gsl_integration_workspace_free(w);
        #endif

            return result;
}






            void state_bcast(state & data){

            vector<double> self_rcv(nw1);

            vector<double> K1_s_rcv_spin(K1_dim1);
            vector<double> K1_t_rcv_spin(K1_dim1);
            vector<double> K1_u_rcv_spin(K1_dim1);

            vector<double> K2_s_rcv_spin(K2_dim1);
            vector<double> K2_t_rcv_spin(K2_dim1);
            vector<double> K2_u_rcv_spin(K2_dim1);
        #if sym ==0
            vector<double> K2b_s_rcv_spin(K2_dim1);
            vector<double> K2b_t_rcv_spin(K2_dim1);
            vector<double> K2b_u_rcv_spin(K2_dim1);
        #endif
            vector<double> K3_s_rcv_spin(K3_dim1);
            vector<double> K3_t_rcv_spin(K3_dim1);
            vector<double> K3_u_rcv_spin(K3_dim1);

            vector<double> irred_rcv_spin(irred_dim1);

            vector<double> K1_s_rcv_dens(K1_dim1);
            vector<double> K1_t_rcv_dens(K1_dim1);
            vector<double> K1_u_rcv_dens(K1_dim1);

            vector<double> K2_s_rcv_dens(K2_dim1);
            vector<double> K2_t_rcv_dens(K2_dim1);
            vector<double> K2_u_rcv_dens(K2_dim1);
        #if sym ==0
            vector<double> K2b_s_rcv_dens(K2_dim1);
            vector<double> K2b_t_rcv_dens(K2_dim1);
            vector<double> K2b_u_rcv_dens(K2_dim1);
        #endif
            vector<double> K3_s_rcv_dens(K3_dim1);
            vector<double> K3_t_rcv_dens(K3_dim1);
            vector<double> K3_u_rcv_dens(K3_dim1);

            vector<double> irred_rcv_dens(irred_dim1);
            // Get the rank of the process
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            if(world_rank==0){

            self_rcv = data.selfenergy.allout();

            K1_s_rcv_spin = data.vertex.spinvertex.svertex.allout_K1();
            K1_s_rcv_dens = data.vertex.densvertex.svertex.allout_K1();
            K1_t_rcv_spin = data.vertex.spinvertex.tvertex.allout_K1();
            K1_t_rcv_dens = data.vertex.densvertex.tvertex.allout_K1();
            K1_u_rcv_spin= data.vertex.spinvertex.uvertex.allout_K1();
            K1_u_rcv_dens= data.vertex.densvertex.uvertex.allout_K1();

            K2_s_rcv_spin = data.vertex.spinvertex.svertex.allout_K2();
            K2_s_rcv_dens= data.vertex.densvertex.svertex.allout_K2();
            K2_t_rcv_spin = data.vertex.spinvertex.tvertex.allout_K2();
            K2_t_rcv_dens = data.vertex.densvertex.tvertex.allout_K2();
            K2_u_rcv_spin = data.vertex.spinvertex.uvertex.allout_K2();
            K2_u_rcv_dens =data.vertex.densvertex.uvertex.allout_K2();
        #if sym ==0
            K2b_s_rcv_spin = data.vertex.spinvertex.svertex.allout_K2b();
            K2b_s_rcv_dens= data.vertex.densvertex.svertex.allout_K2b();
            K2b_t_rcv_spin = data.vertex.spinvertex.tvertex.allout_K2b();
            K2b_t_rcv_dens = data.vertex.densvertex.tvertex.allout_K2b();
            K2b_u_rcv_spin = data.vertex.spinvertex.uvertex.allout_K2b();
            K2b_u_rcv_dens =data.vertex.densvertex.uvertex.allout_K2b();
        #endif
            K3_s_rcv_spin = data.vertex.spinvertex.svertex.allout_K3();
            K3_s_rcv_dens= data.vertex.densvertex.svertex.allout_K3();
            K3_t_rcv_spin = data.vertex.spinvertex.tvertex.allout_K3();
            K3_t_rcv_dens = data.vertex.densvertex.tvertex.allout_K3();
            K3_u_rcv_spin = data.vertex.spinvertex.uvertex.allout_K3();
            K3_u_rcv_dens =data.vertex.densvertex.uvertex.allout_K3();

            irred_rcv_spin =  data.vertex.spinvertex.irred.allout();
            irred_rcv_dens =  data.vertex.densvertex.irred.allout();

            cout << "ok1" << endl;
};

            MPI_Bcast(&self_rcv[0],nw1,MPI_DOUBLE,0,MPI_COMM_WORLD);

            MPI_Bcast(&K1_s_rcv_spin[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K1_s_rcv_dens[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K1_t_rcv_spin[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K1_t_rcv_dens[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K1_u_rcv_spin[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K1_u_rcv_dens[0],K1_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);

            MPI_Bcast(&K2_s_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2_s_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2_t_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2_t_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2_u_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2_u_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);

        #if sym==0
            MPI_Bcast(&K2b_s_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2b_s_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2b_t_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2b_t_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2b_u_rcv_spin[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K2b_u_rcv_dens[0],K2_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);


        #endif

            MPI_Bcast(&K3_s_rcv_spin[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K3_s_rcv_dens[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K3_t_rcv_spin[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K3_t_rcv_dens[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K3_u_rcv_spin[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&K3_u_rcv_dens[0],K3_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);


            MPI_Bcast(&irred_rcv_spin[0],irred_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(&irred_rcv_dens[0],irred_dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            cout << "ok2" << endl;
            if(world_rank!=0){
            data.selfenergy.self_copy(self_rcv);

            data.vertex.spinvertex.svertex.K1_copy(K1_s_rcv_spin);
            data.vertex.densvertex.svertex.K1_copy(K1_s_rcv_dens);
            data.vertex.spinvertex.tvertex.K1_copy(K1_t_rcv_spin);
            data.vertex.densvertex.tvertex.K1_copy(K1_t_rcv_dens);
            data.vertex.spinvertex.uvertex.K1_copy(K1_u_rcv_spin);
            data.vertex.densvertex.uvertex.K1_copy(K1_u_rcv_dens);

            data.vertex.spinvertex.svertex.K2_copy(K2_s_rcv_spin);
            data.vertex.densvertex.svertex.K2_copy(K2_s_rcv_dens);
            data.vertex.spinvertex.tvertex.K2_copy(K2_t_rcv_spin);
            data.vertex.densvertex.tvertex.K2_copy(K2_t_rcv_dens);
            data.vertex.spinvertex.uvertex.K2_copy(K2_u_rcv_spin);
            data.vertex.densvertex.uvertex.K2_copy(K2_u_rcv_dens);
        #if sym==0
            data.vertex.spinvertex.svertex.K2b_copy(K2b_s_rcv_spin);
            data.vertex.densvertex.svertex.K2b_copy(K2b_s_rcv_dens);
            data.vertex.spinvertex.tvertex.K2b_copy(K2b_t_rcv_spin);
            data.vertex.densvertex.tvertex.K2b_copy(K2b_t_rcv_dens);
            data.vertex.spinvertex.uvertex.K2b_copy(K2b_u_rcv_spin);
            data.vertex.densvertex.uvertex.K2b_copy(K2b_u_rcv_dens);
        #endif

            data.vertex.spinvertex.svertex.K3_copy(K3_s_rcv_spin);
            data.vertex.densvertex.svertex.K3_copy(K3_s_rcv_dens);
            data.vertex.spinvertex.tvertex.K3_copy(K3_t_rcv_spin);
            data.vertex.densvertex.tvertex.K3_copy(K3_t_rcv_dens);
            data.vertex.spinvertex.uvertex.K3_copy(K3_u_rcv_spin);
            data.vertex.densvertex.uvertex.K3_copy(K3_u_rcv_dens);

            data.vertex.spinvertex.irred.copy(irred_rcv_spin);
            data.vertex.densvertex.irred.copy(irred_rcv_dens);
};



            if(world_rank==0){
            cout << "Data distributed to all processes.." << endl;
};
}


