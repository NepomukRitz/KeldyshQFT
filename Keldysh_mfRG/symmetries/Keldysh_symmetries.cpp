# include "Keldysh_symmetries.hpp"

auto convertToIndepIndex(int iK) -> int
{
    if(isInList(iK, non_zero_Keldysh_K3)) {
        if (iK == 0 || iK == 1)                                             //  iK(indep)   Real index
            return iK;                                                      //  0           0
        else if (iK == 3)                                                   //  1           1
            return 2;                                                       //  2           3
        else                                                                //  3           5
            return iK - 2;                                                  //  4           6
    } else{                                                                 //  5           7
        print("convertToIndepIndex is throwing this error");
        return -1;
    }
}

auto alphas(int index) -> std::vector<int> {
    std::vector<int> alphas (4);
    int base = 8;
    for (int i = 0; i < 4; i++) {
        alphas[i] = index / base;
        index -= alphas[i]*base;
        base /=2;
        alphas[i]++;
    }
    //alphas[0] = (index % 16)/8 + 1;
    //alphas[1] = (index % 8)/4 + 1;
    //alphas[2] = (index % 4)/2 + 1;
    //alphas[3] = (index % 2) + 1;
    return alphas;
}

auto indices_sum(int i0, int i2, const char channel) -> std::vector<int> {
    std::vector<int> indices ({0, 0});           // Return std::vector, already correct for Matsubara
    if (KELDYSH){
        std::vector<int> alphas_i0 = alphas(i0);   // Calculate the alphas of each input. Refer to these alphas as (1'2'|12)
        std::vector<int> alphas_i2 = alphas(i2);   // Calculate the alphas of each input. Refer to these alphas as (34|3'4')

        //Distribute the alphas of indices i0 and i2 into i1 and i3
        switch (channel) {
            case 'a':
                indices[0] = 8*(alphas_i0[0]-1) + 4*(alphas_i2[3]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i0[3]-1);  // i1 = (1'4'|32)
                indices[1] = 8*(alphas_i2[2]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i2[1]-1);  // i3 = (3'2'|14)
                break;
            case 'p':
                indices[0] = 8*(alphas_i0[0]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i2[1]-1);  // i1 = (1'2'|34)
                indices[1] = 8*(alphas_i2[2]-1) + 4*(alphas_i2[3]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i0[3]-1);  // i3 = (3'4'|12)
                break;
            case 't':
                indices[0] = 8*(alphas_i2[3]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i0[3]-1);  // i1 = (4'2'|32)
                indices[1] = 8*(alphas_i0[0]-1) + 4*(alphas_i2[2]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i2[1]-1);  // i3 = (1'3'|14)
                break;
            default:;
        }
    }
    return indices;
}

