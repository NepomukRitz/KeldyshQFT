/*
 * compute physical objects such as Fermi function, free propagator
 */

#ifndef PHYSICS_FUNCTIONS_H
#define PHYSICS_FUNCTIONS_H

// Fermi-Dirac distribution function at temperature T and chem. potential mu, as a function of frequency omega
rvec fermiFunction(rvec& omega, double T, double mu) {
  rvec nF (omega.size());
  for (int i_omega=0; i_omega<omega.size(); ++i_omega) {
    nF(i_omega) = 1./(exp((omega(i_omega) - mu)/T) + 1.);
  }
  return nF;
}

// Fermi-Dirac distribution function at temperature T and chem. potential mu, for fixed frequency omega
double fermiFunction(double omega, double T, double mu) {
  return 1./(exp((omega - mu)/T) + 1.);
}

// reservoir-dressed non-interacting impurity Green's function as a function of frequency,
// for given impurity energies eps_sigma, temperature T, chem. pot. mu, hybridization Gamma, flow parameter Lambda
// TODO: adjust Keldysh structure to fermionic convention of Kamenev
V2P get_G0(rvec& omega, rvec& eps_sigma, double T, double mu, double Gamma, double Lambda) {
  int N_omega = omega.size();
  V2P G0 (N_omega);

  for (int i_sigma=0; i_sigma<2; ++i_sigma) {
    for (int i_omega=0; i_omega<N_omega; ++i_omega) {
      G0(i_sigma, 1, 0, i_omega) = 1./(omega(i_omega) - eps_sigma(i_sigma) + I*(Gamma + Lambda)/2.); // retarded G0
      G0(i_sigma, 0, 1, i_omega) = conj(G0(i_sigma, 1, 0, i_omega));                                 // advanced G0
      G0(i_sigma, 1, 1, i_omega) = (1 - 2*fermiFunction(omega(i_omega), T, mu))
                                    * (G0(i_sigma, 1, 0, i_omega) - G0(i_sigma, 0, 1, i_omega));
    }
  }
  return G0;
}

// return self-energy with initial value Sigma^R = Sigma^A = U/2
// TODO: adjust Keldysh structure to fermionic convention of Kamenev
V2P initialize_SE(int N_omega, double U) {
  V2P SE (N_omega);
  for (int i_sigma=0; i_sigma<2; ++i_sigma) {
    for (int i_omega=0; i_omega<N_omega; ++i_omega) {
      SE(i_sigma, 0, 1, i_omega) = U / 2;  // retarded self-energy
      SE(i_sigma, 1, 0, i_omega) = U / 2;  // advanced self-energy
    }
  }
  return SE;
}

// compute the impurity Green's function using the non-interacting one and the self-energy
V2P get_G(V2P& G0, V2P& SE) {
  V2P temp (G0.size());
  temp = (G0.inv() + SE.inv()).inv();
  return temp;
}

// identity matrix in Keldysh space
V2P identity(int N_omega) {
  V2P temp (N_omega);
  for (int i_sigma=0; i_sigma<2; ++i_sigma) {
    for (int i_omega=0; i_omega<N_omega; ++i_omega) {
      temp(i_sigma, 0, 0, i_omega) = 1;
      temp(i_sigma, 1, 1, i_omega) = 1;
    }
  }
  return temp;
}

// compute single-scale propagator
// TODO: need \dot{G0}, or more compact formula as derived by Severin Jakobs
V2P get_S(V2P& G0, V2P& SE) {
  int N_omega = G0.size();
  V2P G = get_G(G0, SE);
  return (identity(N_omega) + G * SE) * G0 * (SE * G + identity(N_omega));
}


#endif // PHYSICS_FUNCTIONS_H
