#include <iostream>
#include <vector>
#include <complex>
#include "vertex.h"

using namespace std;

typedef complex<double> comp; // TODO: redefine comp() function in Julian's kagome.cpp

class internal_structure : public vector<vector<vector<double> > > {
public:
    internal_structure() : vector<vector<vector<double> > > () {};
    internal_structure(int input[3]) : vector<vector<vector<double> > > (input[0], vector<vector<double> > (input[1], vector<double> (input[2]))) {};
};

int main() {
    std::cout << "Hello, World!" << std::endl;

    int nuc_eff = 5;
    int input[3] = {nuc_eff, (nuc_eff+1)/2, 3};

    parvert<fullvert<internal_structure> > test (input);




    return 0;
}