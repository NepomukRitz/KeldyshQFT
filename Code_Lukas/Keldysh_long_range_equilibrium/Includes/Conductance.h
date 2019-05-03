#ifndef CONDUCTANCE_07042017  
#define CONDUCTANCE_07042017 

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "approxxpp.h"
#include "basic.h"

#include "Physics.h"
#include "Numerics.h"

using namespace std;

//Conductance fuer gebrochene Paritaetssymmetrie:

class Conductance_eq_nonint{
 public:
 Physics phy;
 int Nges;
 int number_of_freq;
 matrix<double> frequencies;
 matrix<complex<double> > hamiltonian;
 double c_up;
 double c_down;
 matrix<double> conductance_integrand_up;
 Conductance_eq_nonint(Physics phy_in, int Nges, int number_of_freq, matrix<complex<double> > &hamiltonian_in);
 void save(char *filename);
};

Conductance_eq_nonint::Conductance_eq_nonint(Physics phy_in, int Nges_in, int number_of_freq_in, matrix<complex<double> > &hamiltonian_in): phy(phy_in), Nges(Nges_in), number_of_freq(number_of_freq_in), frequencies(number_of_freq), hamiltonian(hamiltonian_in), conductance_integrand_up(number_of_freq){
 double delta_freq=4./(number_of_freq -1.0);
 for(int fi=0; fi<number_of_freq; ++fi){
  frequencies(fi)=-2 +delta_freq*fi;
 }

 matrix<double> der(number_of_freq);
 
 //Erzeuge eine triviale Distributionsfunktion
 matrix<double> tmp_x(3);
 matrix<matrix<complex<double> > > tmp_y(3);
 matrix<complex<double> > tmp_y_sites(Nges+2);
 complex<double> I(0.0,1.0);
 tmp_y_sites=0.0*I;
 tmp_x(0)=0.0;
 tmp_x(1)=1.0;
 tmp_x(2)=2.0;
 tmp_y(0)=tmp_y_sites;
 tmp_y(1)=tmp_y_sites;
 tmp_y(2)=tmp_y_sites;
 linear_ipol_bin<matrix<complex<double> > > distribution(tmp_x,tmp_y);
 
 //Erzeuge triviale Keldyshkomponente
 matrix<complex<double> > EKel(Nges, Nges);
 EKel=0.0*I;
 double taul=1.0;
 double Lambda=1e-8;

 if(phy.T!=0.0){
  matrix<double> cu(number_of_freq);
  matrix<double> cd(number_of_freq);
  for(int fi=0; fi<number_of_freq; ++fi){
   matrix<matrix<complex<double> > > Gu=green_and_single_various_distribution_functions(frequencies(fi), hamiltonian, EKel, phy.h, distribution, taul, Lambda);
   der(fi) = 1./phy.T/(1.+exp((frequencies(fi)-phy.mu)/phy.T))/(1.+exp(-(frequencies(fi)-phy.mu)/phy.T));
   double gu;
   if(frequencies(fi) + phy.h/2. <-2. || frequencies(fi) + phy.h/2. > 2.){
    gu = .0;
   }
   else{
    gu=1/(2.)*sqrt(4.-(frequencies(fi)+phy.h/2.)*(frequencies(fi)+phy.h/2.));
   }

   cu(fi) = (2.*abs(Gu(0)(Nges-1,0))*gu);
   cu(fi) = cu(fi)*cu(fi);
   
   matrix<matrix<complex<double> > > Gd=green_and_single_various_distribution_functions(frequencies(fi), hamiltonian, EKel, -phy.h, distribution, taul, Lambda);
   der(fi) = 1./phy.T/(1.+exp((frequencies(fi)-phy.mu)/phy.T))/(1.+exp(-(frequencies(fi)-phy.mu)/phy.T));
   double gd;
   if(frequencies(fi)- phy.h/2. <-2. || frequencies(fi) - phy.h/2. > 2.){
    gd = .0;
   }
   else{
    gd=1/(2.)*sqrt(4.-(frequencies(fi)-phy.h/2.)*(frequencies(fi)-phy.h/2.));
   }

   cd(fi) = (2.*abs(Gd(0)(Nges-1,0))*gd);
   cd(fi) = cd(fi)*cd(fi);
  }

  //Integration
  c_up=0.0;
  c_down=0.0;
  double norm=0.0;
  for(int fi=0; fi<number_of_freq-1; ++fi){
   conductance_integrand_up(fi)=.5*(der(fi)*cu(fi) +der(fi+1)*cu(fi+1) );
   c_up   += .5*(der(fi)*cu(fi) +der(fi+1)*cu(fi+1) )*(frequencies(fi+1)-frequencies(fi));
   c_down += .5*(der(fi)*cd(fi) +der(fi+1)*cd(fi+1) )*(frequencies(fi+1)-frequencies(fi));
   norm   += .5*(der(fi)        +der(fi+1)          )*(frequencies(fi+1)-frequencies(fi));
  }
  //c_up=c_up/norm; /* Das ist nur fuer kleine Temperaturen valide*/
  //c_down=c_down/norm;
  c_up=c_up; /* Das ist nur fuer kleine Temperaturen valide*/
  c_down=c_down;
 }
 else{
   matrix<matrix<complex<double> > > Gu=green_and_single_various_distribution_functions(phy.mu, hamiltonian, EKel, phy.h, distribution, taul, Lambda);
   Gu(0).save("Conductance.mat","Gu");
   double gu;
   if(phy.mu+phy.h/2. <-2. || phy.mu+phy.h/2. >2.){
    gu = .0;
   }
   else{
    gu=1/(2.)*sqrt(4.-(phy.mu+phy.h/2.)*(phy.mu+phy.h/2.));
   }

   c_up = (2.*abs(Gu(0)(Nges-1,0))*gu);
   c_up = c_up*c_up;
   
   matrix<matrix<complex<double> > > Gd=green_and_single_various_distribution_functions(phy.mu, hamiltonian, EKel, -phy.h, distribution, taul, Lambda);
   double gd;
   if(phy.mu-phy.h/2. <-2. || phy.mu-phy.h/2. >2.){
    gd = .0;
   }
   else{
    gd=1/(2.)*sqrt(4.-(phy.mu-phy.h/2.)*(phy.mu-phy.h/2.));
   }

   c_down = (2.*abs(Gd(0)(Nges-1,0))*gd);
   c_down = c_down*c_down;
 }
}



#endif
