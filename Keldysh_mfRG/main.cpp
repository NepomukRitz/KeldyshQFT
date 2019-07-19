#include <iostream>
#include <vector>
#include <complex>
#include "vertex.h"
#include "util.h"

using namespace std;

typedef complex<double> comp; // TODO: redefine comp() function in Julian's kagome.cpp


//class internal_structure : public vector<vector<vector<double> > > {
//public:
//    internal_structure() : vector<vector<vector<double> > > () {};
//    internal_structure(int input[3]) : vector<vector<vector<double> > > (input[0], vector<vector<double> > (input[1], vector<double> (input[2]))) {};
//};





int main() {
    for(int i=0; i<1; ++i) {
        std::cout << "Hello, World!" << std::endl;
    }

    int nuc_eff = 5;
    int input[3] = {nuc_eff, (nuc_eff+1)/2, 3};


    int d_in = 3;
    int K_in[] = {100,100,100};
    int i_in[] = {3,71,25};

    double t0 = get_time();


    for (int k=0; k<1000000; ++k) {
        int index = 0;
        for (int i = 0; i < d_in; ++i) {
            int mult = 1;
            for (int j = i + 1; j < d_in; ++j) {
                mult *= K_in[j];
            }
            index += mult * i_in[i];
        }
    }
    get_time(t0, "us");

    double t1 = get_time();

    for (int k=0; k<1000000; ++k) {
        int index2 = i_in[0] * K_in[1] * K_in[2] + i_in[1] * K_in[2] + i_in[2];
    }

    get_time(t1, "us");

    //cout << index << " " << index2 << endl;


    parvert<fullvert<comp> > test;


    return 0;
}