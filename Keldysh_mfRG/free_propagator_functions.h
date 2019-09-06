//
// Created by Sa.Aguirre on 9/4/19.
//

#ifndef KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H
#define KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H

#include "data_structures.h"
#include "parameters.h"

/*******PROPAGATOR FUNCTIONS***********/
comp gR(double Lambda, double omega)
{
#if REG==1
    return 1./(omega-epsilon);
#elif REG==2
    return 1./(omega-epsilon + (comp)0.5i*(GAMMA_REG+Lambda));
#endif
}

comp gA(double Lambda, double omega)
{
#if REG==1
    return 1./(omega-epsilon);
#elif REG==2
    return 1./(omega-epsilon - (comp)0.5i*(GAMMA_REG+Lambda));
#endif
}


#if REG==2
//Self-explanatory
double Fermi_distribution(double omega)
{
    return 1./(exp((omega-mu)/T)+1.);
}

comp gK(double Lambda, double omega)
{
    return (1.-2.*Fermi_distribution(omega))*(gR(Lambda,omega) - gA(Lambda,omega));
}


comp sR(double Lambda, double omega)
{
    return -(comp)0.5i*pow(gR(Lambda,omega),2.);
}
comp sA(double Lambda, double omega)
{
    return +(comp)0.5i*pow(gA(Lambda,omega),2.);
}
comp sK(double Lambda, double omega)
{
    return (1.-2.*Fermi_distribution(omega))*(sR(Lambda,omega) - sA(Lambda,omega));
}

#endif

#endif //KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H