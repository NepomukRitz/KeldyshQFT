#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace H5;

const H5std_string FILE_NAME ("test.h5");

int main() {
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
}
