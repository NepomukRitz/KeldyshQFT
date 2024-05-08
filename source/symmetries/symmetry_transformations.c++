#include "symmetry_transformations.hpp"

void switch_channel(IndicesSymmetryTransformations& indices) {
    switch (indices.channel_rvert) {
        case 'a':
            indices.channel_rvert = 't';
            break;
        case 't':
            indices.channel_rvert = 'a';
            break;
        default:;
    }
    if (INTERPOL2D_FOR_K3) {
        switch (indices.channel_bubble) {
            case 'a':
                indices.channel_bubble = 't';
                break;
            case 't':
                indices.channel_bubble = 'a';
                break;
            default: ;
        }
    }
}

