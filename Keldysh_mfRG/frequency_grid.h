//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_H
#define KELDYSH_MFRG_FREQUENCY_GRID_H

using namespace std;

//TODO: revise this
//TODO: find an elegant way to set the grid (maybe via preprocessor macro)

/**************************** LINEAR GRID ***************/
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

/**************************** LOG GRID ******************/
//to convert on full frequency grid:

/*

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

*/

#endif //KELDYSH_MFRG_FREQUENCY_GRID_H
