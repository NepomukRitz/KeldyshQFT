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

std::ostream& operator << (std::ostream& out, symmetryType symmtype) {
    if (symmtype == symmetric) {out << "symmetric";}
    else {out << "non-symmetric" ;}
    return out;
}

