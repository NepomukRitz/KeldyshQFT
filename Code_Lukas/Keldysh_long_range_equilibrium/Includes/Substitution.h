#ifndef SUBSTITUTION_15032017
#define SUBSTITUTION_15032017

using namespace std;

template <int mode> class Substitution{
	public:
		double Lambda;
		Substitution(double Lambda_in);
		double subst_concatenated (double omega);
		double resu_concatenated (double y);
		double weight_concatenated (double y);
		//rpa_mode:
		static const double alpha = 1.1; //Use the alpha dependent substitution in the RPA case.
		static const double beta =1.0; //This slows the substitution down! do not use the beta dependent substitution for production jobs!
};
 
 
template<int mode> Substitution<mode>::Substitution(double Lambda_in): Lambda(Lambda_in){};

//line to finite line segments
template <>
double Substitution<0>::subst_concatenated (double omega)
{
	if (omega<-6.)
	{
		double oshifted = (omega +6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted)) -6.;
	}
	else if (omega==-6.)
	{
		return -6.;
	}
	else if (omega<-2.)
	{
		return -2.*sqrt(-omega-2.)-2.;
	}
	else if (omega==-2.)
	{
		return -2.;
	}
	else if (omega<2.)
	{
		double x = omega;
		return (sqrt(2.+x)-sqrt(2.-x));
	}
	else if (omega==2.)
	{
		return 2.;
	}
	else if (omega<6.)
	{
		return 2.*sqrt(omega-2.)+2.;
	}
	else if (omega==6.)
	{
		return 6.;
	}
	else
	{
		double oshifted = (omega-6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted)) +6.;
	}
}

//finite line segments to line
template <>
double Substitution<0>::resu_concatenated (double y)
{
	if (y<-6.)
	{
		double yshifted = (y+6.);
		return (-2.*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)-6.;
	}
	else if (y==-6.)
	{
		return -6.;
	}
	else if (y<-2.)
	{
		return -2.-(y+2.)*(y+2.)/4.;
	}
	else if (y==-2.)
	{
		return -2.;
	}
	else if (y<2.)
	{
		if (y==.0)
			return .0;
		else
			return y*sqrt(4./y/y-(y*y-4.)*(y*y-4.)/4./y/y);
	}
	else if (y==2.)
	{
		return 2.;
	}
	else if (y<6.)
	{
		return 2.+(y-2.)*(y-2.)/4.;
	}
	else if (y==6.)
	{
		return 6.;
	}
	else
	{
		double yshifted = (y-6.);
		return (-2.*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)+6.;
	}
}

//measure on finite line segments
template <>
double Substitution<0>::weight_concatenated (double y)
{
	if (y<-6.)
	{
		double yshifted = (y+6.);
		yshifted *= yshifted;
		return 2.*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)*(1.+Lambda);
	}
	else if (y==-6.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-6. This should not happen" << endl;
		return 1e20;
	}
	else if (y<-2.)
	{
		return -1.-y/2.;
	}
	else if (y==-2.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-2. This should not happen" << endl;
		return 1e20;
	}
	else if (y<2.)
	{
		double x = Substitution<0>::resu_concatenated(y);
		return 2./((1./sqrt(2.+x)+1./sqrt(2.-x)));
	}
	else if (y==2.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=2. This should not happen" << endl;
		return 1e20;
	}
	else if (y<6.)
	{
		return -1.+y/2.;
	}
	else if (y==6.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=6. This should not happen" << endl;
		return 1e20;
	}
	else
	{
		double yshifted = (y-6.);
		yshifted *= yshifted;
		return 2.*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)*(1.+Lambda);
	}
}

//Substitution for rpa mode:

//line to finite line segments
template <>
double Substitution<1>::subst_concatenated (double omega)
{
	if (omega<-6.)
	{
	 	return -7. + pow(6./(6.-(omega+6.)*beta),alpha-1.);
	}
	else if (omega==-6.)
	{
		return -6.;
	}
	else if (omega<-2.)
	{
		return -2.*sqrt(-omega-2.)-2.;
	}
	else if (omega==-2.)
	{
		return -2.;
	}
	else if (omega<2.)
	{
		double x = omega;
		return (sqrt(2.+x)-sqrt(2.-x));
	}
	else if (omega==2.)
	{
		return 2.;
	}
	else if (omega<6.)
	{
		return 2.*sqrt(omega-2.)+2.;
	}
	else if (omega==6.)
	{
		return 6.;
	}
	else
	{
	 	return 7. - pow(6./(6.+(omega-6.)*beta),alpha-1.);
	}
}

//finite line segments to line
template <>
double Substitution<1>::resu_concatenated (double y)
{
	if (y<-6.)
	{
	 	return 6./beta*(1. - 1./pow(y+7.,1./(alpha-1.)))-6.;
	}
	else if (y==-6.)
	{
		return -6.;
	}
	else if (y<-2.)
	{
		return -2.-(y+2.)*(y+2.)/4.;
	}
	else if (y==-2.)
	{
		return -2.;
	}
	else if (y<2.)
	{
		if (y==.0)
			return .0;
		else
			return y*sqrt(4./y/y-(y*y-4.)*(y*y-4.)/4./y/y);
	}
	else if (y==2.)
	{
		return 2.;
	}
	else if (y<6.)
	{
		return 2.+(y-2.)*(y-2.)/4.;
	}
	else if (y==6.)
	{
		return 6.;
	}
	else
	{
	 	return 6./beta*(1./pow(7.-y,1./(alpha-1.))-1.) + 6.;
	}
}

//measure on finite line segments
template <>
double Substitution<1>::weight_concatenated (double y)
{
	if (y<-6.)
	{
	 	return 6./(beta*(alpha - 1.)*pow(7.+y,alpha/(alpha-1.)));
	}
	else if (y==-6.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-6. This should not happen" << endl;
		return 1e20;
	}
	else if (y<-2.)
	{
		return -1.-y/2.;
	}
	else if (y==-2.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-2. This should not happen" << endl;
		return 1e20;
	}
	else if (y<2.)
	{
		double x = Substitution<1>::resu_concatenated(y);
		return 2./((1./sqrt(2.+x)+1./sqrt(2.-x)));
	}
	else if (y==2.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=2. This should not happen" << endl;
		return 1e20;
	}
	else if (y<6.)
	{
		return -1.+y/2.;
	}
	else if (y==6.)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=6. This should not happen" << endl;
		return 1e20;
	}
	else
	{
	 	return 6./(beta*(alpha - 1.)*pow(7.-y,alpha/(alpha-1.)));
	}
}



#endif
