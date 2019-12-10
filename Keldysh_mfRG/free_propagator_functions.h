//
// Created by Sa.Aguirre on 9/4/19.
//

#ifndef KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H
#define KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H

#include "data_structures.h"
#include "parameters.h"

/*******PROPAGATOR FUNCTIONS***********/
auto gR(double Lambda, double omega) -> comp
{
#if REG==1
    return 1./(omega-epsilon);
#elif REG==2
    return 1./(omega-epsilon + ((comp)0.5i*(GAMMA_REG+Lambda)));
#endif
}

auto gA(double Lambda, double omega) -> comp
{
#if REG==1
    return 1./(omega-epsilon);
#elif REG==2
    return 1./(omega-epsilon - ((comp)0.5i*(GAMMA_REG+Lambda)));
#endif
}

//Self-explanatory
auto Fermi_distribution(double omega) -> double
{
    return 1./(exp((omega-mu)/T)+1.);
}

auto gK(double Lambda, double omega) -> comp
{
    return (1.-2.*Fermi_distribution(omega))*(gR(Lambda,omega) - gA(Lambda,omega));
}

#if REG==2

auto sR(double Lambda, double omega) -> comp
{
    return -((comp)0.5i)*gR(Lambda,omega)*gR(Lambda,omega);
}
auto sA(double Lambda, double omega) -> comp
{
    return ((comp)0.5i)*gA(Lambda,omega)*gA(Lambda,omega);
}
auto sK(double Lambda, double omega) -> comp
{
    return (1.-2.*Fermi_distribution(omega))*(sR(Lambda,omega) - sA(Lambda,omega));
}

#endif

#endif //KELDYSH_MFRG_FREE_PROPAGATOR_FUNCTIONS_H
