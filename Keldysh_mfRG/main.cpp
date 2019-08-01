#include <cstdlib>
#include <iostream>
#include <complex>

#include "vertex.h"

using namespace std;

typedef complex<double> comp; // TODO: redefine comp() function in Julian's kagome.cpp


int main() {


    Vertex<fullvert<comp> > test;
    test.spinvertex.avertex.K1_setvert(0, 1, 2, 3);

    Vertex<avert<comp> > test1;
    Vertex<avert<comp> > test2;
    Vertex<avert<comp> > test3;

    test3 = test1 + test2;


    return 0;
}