#include "data_structures.hpp"

auto isfinite(comp z) -> bool {
    return std::isfinite(real(z)) and std::isfinite(imag(z));
};


std::ostream& operator << (std::ostream& out, K_class k) {
    if (k == k1 or k == k2) {out << "K" << static_cast<int>(k+1);}
    else if(k == k3) {out << "K" << static_cast<int>(k);}
    else if(k == k2b) {out << "K" << static_cast<int>(k) << "'";}
    else  {out << "selfenergy";}
    return out;
}

std::ostream& operator << (std::ostream& out, vertexType symmtype) {
    if (symmtype == symmetric_full) {out << "symmetric_full";}
    else if (symmtype == symmetric_r_irred) {out << "symmetric_r_irred";}
    else if (symmtype == non_symmetric_diffleft) {out << "non_symmetric_diffleft";}
    else if (symmtype == non_symmetric_diffright) {out << "non_symmetric_diffright";}
    else {assert(false);}
    return out;
}

