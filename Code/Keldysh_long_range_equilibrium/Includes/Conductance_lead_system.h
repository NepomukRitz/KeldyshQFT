#ifndef CONDUCTANCE_LEAD_SYSTEM_28072017
#define CONDUCTANCE_LEAD_SYSTEM_28072017


#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Vertex.h"
#include "Barevertex.h"

#define ACCURACY_VERTEX_CONTRIBUTION 1e-5
#define ACCURACY_CONDUCTANCE_INTEGRAND 1e-6
#define NUMBER_OF_FREQUENCIES_FOR_VERTEX_CONTRIBUTION 101

template <int mode> class Integrand_lead_system_vertex_contribution_zero_mag{
	public:
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Vertex<mode> &gamma;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_lead_system_vertex_contribution_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode> Integrand_lead_system_vertex_contribution_zero_mag<mode>::Integrand_lead_system_vertex_contribution_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), gamma(gamma_in){} 

template<int mode> matrix<complex<double> > Integrand_lead_system_vertex_contribution_zero_mag<mode>::operator()(double internal){
	matrix<complex<double> > ret(num.Nges,num.Nges);
	double intline = sub.resu_concatenated(internal);
	double sum = intline + external_freq;
	double diff =intline - external_freq;
	double hyb = sqrt(4.- pow(intline,2));
 	double prefactor_p =(1./tanh((sum/2.-phy.mu)/phy.T)-1.+2.*fermi(intline,phy.mu,phy.T))/M_PI;
 	double prefactor_x =(-1./tanh(diff/(2.*phy.T))-1.+2.*fermi(intline,phy.mu,phy.T))/M_PI;
	cout<<"intline="<<intline<<endl;
	cout<<"sum="<<sum<<endl;
	cout<<"diff="<<diff<<endl;
	cout<<"prefactor_p="<<prefactor_p<<endl;
	cout<<"prefactor_x="<<prefactor_x<<endl;
	syma<complex<double> > Gu = pre.iGu(internal);
	//dynamic central contribution:	
	matrix<complex<double> > Pud = Trafo(gamma.aPud_central_ipol(sum));
	matrix<complex<double> > Pdd = Trafo(gamma.aPdd_central_ipol(sum));
	matrix<complex<double> > Xud = Trafo(gamma.aXud_central_ipol(diff));
	matrix<complex<double> > Ddd = Trafo(gamma.aDdd_central_ipol(-diff));
	for(int j=0; j<num.Nges; ++j){
		for(int i=0; i<num.Nges; ++i){
			ret(j,i) = (imag(Pud(j,i)) + imag(Pdd(j,i)))*conj(Gu(num.Nges-1,i))*Gu(num.Nges-1,j)*hyb*prefactor_p;
			ret(j,i)+= (imag(Xud(i,j)) - imag(Ddd(i,j)))*conj(Gu(num.Nges-1,j))*Gu(num.Nges-1,i)*hyb*prefactor_x;
		}
	}
	cout<<"sub.weight_concatenated(internal)="<<sub.weight_concatenated(internal)<<endl;
	return sub.weight_concatenated(internal)*ret;
	//matrix<complex<double> > rueck(num.Nges,num.Nges);
	//rueck=(complex<double>) prefactor_p;
	//rueck=(complex<double>) 1.0;
	//return sub.weight_concatenated(internal)*rueck;
}

template<int mode> matrix<double> Integrand_lead_system_vertex_contribution_zero_mag<mode>::select(matrix<complex<double> > &M){
	int eff= 2.*M.dim_r*M.dim_c;
	int eff_half= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<num.Nges;++j){
		for(int i=0; i<num.Nges;++i,++z){
			n(z) = M(j,i).imag();
			n(z+eff_half) = M(j,i).real();
		}
	}
	return n;
}


template <int mode> class Lead_system_vertex_contribution_zero_mag{
	public:
		static const double eps = 1e-16;
		static const double accuracy=ACCURACY_VERTEX_CONTRIBUTION; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		static const double Lambda = 1e-8;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Vertex<mode> gamma;
		Stops<mode> stops_obj;
		Lead_system_vertex_contribution_zero_mag(Physics &phy_in, 
		                                         Numerics &num_in, 
		                                         Precomputation_zeromag<mode> &pre_in, 
		                                         Substitution<mode> &sub_in, 
		                                         Vertex<mode> &gamma_in);
		matrix<complex<double> > operator()(double external_freq);
};

template<int mode> Lead_system_vertex_contribution_zero_mag<mode>::Lead_system_vertex_contribution_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), gamma(gamma_in), stops_obj(phy, sub, Lambda){} 

template<int mode> matrix<complex<double> > Lead_system_vertex_contribution_zero_mag<mode>::operator()(double external_freq){
	
	Integrand_lead_system_vertex_contribution_zero_mag<mode> vertex_int(external_freq, phy, num, pre, sub, gamma);
	matrix<double> additional_stops(9); 
	additional_stops(0) = sub.subst_concatenated(max(-2.,min(external_freq,2.)));
	additional_stops(1) = sub.subst_concatenated(max(-2.,min(external_freq-10.*phy.T,2.)));
	additional_stops(2) = sub.subst_concatenated(max(-2.,min(external_freq+10.*phy.T,2.)));
	additional_stops(3) = sub.subst_concatenated(max(-2.,min(2.*phy.mu-external_freq,2.)));
	additional_stops(4) = sub.subst_concatenated(max(-2.,min(2.*phy.mu-external_freq-10.*phy.T,2.)));
	additional_stops(5) = sub.subst_concatenated(max(-2.,min(2.*phy.mu-external_freq+10.*phy.T,2.)));
	additional_stops(6) = sub.subst_concatenated(max(-2.,min(phy.mu,2.)));
	additional_stops(7) = sub.subst_concatenated(max(-2.,min(phy.mu-5.*phy.T,2.)));
	additional_stops(8) = sub.subst_concatenated(max(-2.,min(phy.mu+5.*phy.T,2.)));
	matrix<complex<double> > erg(num.Nges, num.Nges);
	erg = (complex<double>) .0;
	additional_stops(0)=0.0;
	matrix<double> stops = stops_obj.Conductance_stops(additional_stops);
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,vertex_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}


template <int mode> class Integrand_conductance_lead_system_trivial_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_conductance_lead_system_trivial_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode> Integrand_conductance_lead_system_trivial_zero_mag<mode>::Integrand_conductance_lead_system_trivial_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in){} 

template<int mode> matrix<double> Integrand_conductance_lead_system_trivial_zero_mag<mode>::operator()(double internal){
	//Trivial Contribution:	
	matrix<double> ret(1);
	double intline = sub.resu_concatenated(internal);
	double dfermi = -1./phy.T*fermi(intline,phy.mu,phy.T)*fermi(-intline+2.*phy.mu,phy.mu,phy.T); 
	double hyb = sqrt(4.- pow(intline,2));
	syma<complex<double> > G = pre.iGu(internal);
	double g_trivial = hyb*abs(G(num.Nges-1,0));
	g_trivial = g_trivial*g_trivial;
	ret(0) = -g_trivial*dfermi*sub.weight_concatenated(internal);
	cout<<"dfermi="<<dfermi<<endl;

	return ret;
}

template<int mode> matrix<double> Integrand_conductance_lead_system_trivial_zero_mag<mode>::select(matrix<double> &M){
	return M;
}

template <int mode> class Integrand_conductance_lead_system_vertex_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		linear_ipol_bin<matrix<complex<double> > > &vertex_cont;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_conductance_lead_system_vertex_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, linear_ipol_bin<matrix<complex<double> > > &vertex_cont_in);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode> Integrand_conductance_lead_system_vertex_zero_mag<mode>::Integrand_conductance_lead_system_vertex_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, linear_ipol_bin<matrix<complex<double> > > &vertex_cont_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), vertex_cont(vertex_cont_in){} 

template<int mode> matrix<complex<double> > Integrand_conductance_lead_system_vertex_zero_mag<mode>::operator()(double internal){
	//Vertex Contribution:	
	matrix<complex<double> > ret(1);
	double intline = sub.resu_concatenated(internal);
	double dfermi = -1./phy.T*fermi(intline,phy.mu,phy.T)*fermi(-intline+2.*phy.mu,phy.mu,phy.T); 
	double hyb = sqrt(4.- pow(intline,2));
	syma<complex<double> > G = pre.iGu(internal);
	matrix<complex<double> > Vertex_contribution = vertex_cont(intline);
	matrix<complex<double> > tmp = G.conj()*Vertex_contribution*G;  
	complex<double> g_vertex = tmp(0,0);
	ret(0) = -g_vertex*dfermi*hyb*sub.weight_concatenated(internal);

	
	return ret;
}

template<int mode> matrix<double> Integrand_conductance_lead_system_vertex_zero_mag<mode>::select(matrix<complex<double> > &M){
	matrix<double> n(2);
	n(0) = real(M(0));
	n(1) = imag(M(0));
	return n;
}


template <int mode> class Conductance_lead_system_zero_mag{
	public:
		static const double eps = 1e-7;
		static const double Lambda = 1e-8;
		static const double accuracy=ACCURACY_CONDUCTANCE_INTEGRAND; //Baue hier eventuell noch die dynamische Genauigkeit ein!

		Physics &phy;
		Numerics &num;
		Vertex<mode> &gamma; 
		Substitution<mode> &sub;
		Precomputation_zeromag<mode> &pre;
		Stops<mode> stops_obj;
		Conductance_lead_system_zero_mag(Physics &phy_in, Numerics &num_in, Vertex<mode> &gamma_in, Substitution<mode> &sub_in, Precomputation_zeromag<mode> &pre_in);
		complex<double> compute_conductance();
};

template <int mode> Conductance_lead_system_zero_mag<mode>::Conductance_lead_system_zero_mag(Physics &phy_in, Numerics &num_in, Vertex<mode> &gamma_in, Substitution<mode> &sub_in, Precomputation_zeromag<mode> &pre_in): phy(phy_in), num(num_in), gamma(gamma_in), sub(sub_in), pre(pre_in), stops_obj(phy, sub, Lambda){}

template <int mode> complex<double> Conductance_lead_system_zero_mag<mode>::compute_conductance(){
	pre.precompute_non_ps(Lambda, sub, gamma.ERetu_ipol_subst); 
	
	//Compute Vertex contribution:
	int number_freq_vertex_cont=NUMBER_OF_FREQUENCIES_FOR_VERTEX_CONTRIBUTION;
	double lower_bound=-2.;
	double upper_bound=+2.;
	matrix<double> freq_vertex_cont(number_freq_vertex_cont);
	double delta_freq = (upper_bound - lower_bound)/(number_freq_vertex_cont-1);
	for(int i=0; i<number_freq_vertex_cont; ++i){
		freq_vertex_cont(i) = lower_bound +i*delta_freq; 
	}
	Lead_system_vertex_contribution_zero_mag<mode> Vertex_contribution(phy, num, pre, sub, gamma);
	matrix<matrix<complex<double> > > vertex_contribution(number_freq_vertex_cont);
	for(int i=0; i<number_freq_vertex_cont; ++i){
		cout<<"i="<<i<<endl;
		cout<<"freq_vertex_cont="<<freq_vertex_cont(i)<<endl;
		vertex_contribution(i) = Vertex_contribution(freq_vertex_cont(i)); 
	}
	linear_ipol_bin<matrix<complex<double> > > vertex_cont(freq_vertex_cont, vertex_contribution); 
	cout<<"Vertex contribution computed"<<endl;

	//Compute total conductance:
	Integrand_conductance_lead_system_trivial_zero_mag<mode> cond_triv(phy, num, pre, sub);
	Integrand_conductance_lead_system_vertex_zero_mag<mode> cond_vertex(phy, num, pre, sub, vertex_cont);
	
	matrix<double> additional_stops(7); 
	matrix<double> erg_triv(1);
	matrix<complex<double> > erg_vertex(1);
	erg_triv =  .0;
	erg_vertex = (complex<double>) .0;

	additional_stops(0)=sub.subst_concatenated(max(-2.,min(phy.mu,2.)));
	additional_stops(1)=sub.subst_concatenated(max(-2.,min(phy.mu-15.*phy.T,2.)));
	additional_stops(2)=sub.subst_concatenated(max(-2.,min(phy.mu-10.*phy.T,2.)));
	additional_stops(3)=sub.subst_concatenated(max(-2.,min(phy.mu-5.*phy.T,2.)));
	additional_stops(4)=sub.subst_concatenated(max(-2.,min(phy.mu+5.*phy.T,2.)));
	additional_stops(5)=sub.subst_concatenated(max(-2.,min(phy.mu+10.*phy.T,2.)));
	additional_stops(6)=sub.subst_concatenated(max(-2.,min(phy.mu+15.*phy.T,2.)));
	matrix<double> stops = stops_obj.Conductance_stops(additional_stops);
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg_triv,stops(i),stops(i+1),accuracy,1e-4,1e-14,cond_triv);     //params: result, start, stop, tolerance, initial step, minimum step, function
			intgk(erg_vertex,stops(i),stops(i+1),accuracy,1e-4,1e-14,cond_vertex);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	cout<<"erg_triv="<<erg_triv(0)<<endl;
	cout<<"erg_vertex="<<erg_vertex(0)<<endl;
	return erg_triv(0) + erg_vertex(0);
	
}




#endif

