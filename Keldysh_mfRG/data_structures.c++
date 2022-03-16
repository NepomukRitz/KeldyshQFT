#include "data_structures.hpp"

auto isfinite(comp z) -> bool {
    return std::isfinite(real(z)) and std::isfinite(imag(z));
};


std::ostream& operator << (std::ostream& out, K_class k) {
    if (k == k1 or k == k2) {out << "K" << static_cast<int>(k+1);}
    else if(k == k3) {out << "K" << static_cast<int>(k);}
    else {out << "K" << static_cast<int>(k) << "'";}
    return out;
}

std::ostream& operator << (std::ostream& out, vertexType symmtype) {
    if (symmtype == symmetric_full) {out << "symmetric_full";}
    else {out << "non-symmetric_full" ;}
    return out;
}

