#ifndef SUBSTITUTION_FLOW_28072017
#define SUBSTITUTION_FLOW_28072017
 
using namespace std;

class Substitution_flow{
	public:
		Substitution_flow();
		double subst (double Lambda);
		double resu  (double x);
		double weight(double x);
};

Substitution_flow::Substitution_flow(){};

double Substitution_flow::subst(double Lambda){
	return log(Lambda/(1.+Lambda));
}
 
double Substitution_flow::resu(double x){
	return exp(x)/(1.-exp(x));
}

double Substitution_flow::weight(double x){
	return exp(x)/(1.-exp(x))/(1.-exp(x));
}
 


#endif
