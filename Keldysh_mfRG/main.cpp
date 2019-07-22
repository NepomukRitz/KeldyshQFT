#include <iostream>
#include <vector>
#include <complex>
#include "vertex.h"
#include "util.h"
#include "parameters.h"

using namespace std;

typedef complex<double> comp; // TODO: redefine comp() function in Julian's kagome.cpp


//class internal_structure : public vector<vector<vector<double> > > {
//public:
//    internal_structure() : vector<vector<vector<double> > > () {};
//    internal_structure(int input[3]) : vector<vector<vector<double> > > (input[0], vector<vector<double> > (input[1], vector<double> (input[2]))) {};
//};


int main() {


    //cout << index << " " << index2 << endl;


    parvert<fullvert<comp> > test;


    return 0;
}