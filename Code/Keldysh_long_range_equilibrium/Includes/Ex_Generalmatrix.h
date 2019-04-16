#ifndef Ex_GENERALMATRIX_30102018
#define Ex_GENERALMATRIX_30102018

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "Numerics.h"
#include "Ex_functions.h"

using namespace std;

class Ex_Generalmatrix{
	public:
		matrix<matrix<syma<complex<double> > > > self_str;
		matrix< matrix < matrix < matrix < complex < double > > > > > dynamic_str;
		matrix < matrix < matrix < double > > > static_str;
		Ex_Generalmatrix(){};
		Ex_Generalmatrix(Numerics &num_in);
		void resize(Numerics &num_in);
		void initialize(double d);
		void initialize_random(int D);
		void save(string filename, string variable);
		void load(string filename, string variable);
		Ex_Generalmatrix operator+ (const Ex_Generalmatrix &gm2);
		double errnorm(double atol,double rtol, Ex_Generalmatrix &y1, Ex_Generalmatrix &y2); 
};

Ex_Generalmatrix::Ex_Generalmatrix(Numerics &num){
	resize(num);
}

void Ex_Generalmatrix::resize(Numerics &num){
	self_str.resize(2);
	dynamic_str.resize(7);
	static_str.resize(7);
	resize_str(self_str(0),num.Nff,num.N);
	resize_str(self_str(1),num.Nff,num.N);
	resize_str(dynamic_str(0),num.Lp_structure,num.N);
	resize_str(dynamic_str(1),num.Lp_structure,num.N);
	resize_str(dynamic_str(2),num.Lp_structure,num.N);
	resize_str(dynamic_str(3),num.Lx_structure,num.N);
	resize_str(dynamic_str(4),num.Lx_structure,num.N);
	resize_str(dynamic_str(5),num.Lx_structure,num.N);
	resize_str(dynamic_str(6),num.Lx_structure,num.N);
	for(int i=0; i<7; ++i){
		resize_str(static_str(i),num.L,num.N);
	}
}

void Ex_Generalmatrix::initialize(double d){
	init(self_str, (complex<double>) d);
	init(dynamic_str, (complex<double>) d);
	init(static_str, d);
}

void Ex_Generalmatrix::initialize_random(int D){
	init_random(self_str, D);
	init_random(dynamic_str, D);
	init_random(static_str, D);
}

void Ex_Generalmatrix::save(string filename, string variable){
	save_str(self_str(0),filename,variable + "_Eu");
	save_str(self_str(1),filename,variable + "_Ed");
	save_str(dynamic_str(0),filename,variable + "_Puu_dyn");
	save_str(dynamic_str(1),filename,variable + "_Pdd_dyn");
	save_str(dynamic_str(2),filename,variable + "_Pud_dyn");
	save_str(dynamic_str(3),filename,variable + "_Xud_dyn");
	save_str(dynamic_str(4),filename,variable + "_Duu_dyn");
	save_str(dynamic_str(5),filename,variable + "_Ddd_dyn");
	save_str(dynamic_str(6),filename,variable + "_Dud_dyn");
	save_str(static_str(0),filename,variable + "_Puu_stat");
	save_str(static_str(1),filename,variable + "_Pdd_stat");
	save_str(static_str(2),filename,variable + "_Pud_stat");
	save_str(static_str(3),filename,variable + "_Xud_stat");
	save_str(static_str(4),filename,variable + "_Duu_stat");
	save_str(static_str(5),filename,variable + "_Ddd_stat");
	save_str(static_str(6),filename,variable + "_Dud_stat");
}

void Ex_Generalmatrix::load(string filename, string variable){
	self_str.resize(2);
	dynamic_str.resize(7);
	static_str.resize(7);
	load_str(self_str(0),filename,variable + "_Eu");
	load_str(self_str(1),filename,variable + "_Ed");
	load_str(dynamic_str(0),filename,variable + "_Puu_dyn");
	load_str(dynamic_str(1),filename,variable + "_Pdd_dyn");
	load_str(dynamic_str(2),filename,variable + "_Pud_dyn");
	load_str(dynamic_str(3),filename,variable + "_Xud_dyn");
	load_str(dynamic_str(4),filename,variable + "_Duu_dyn");
	load_str(dynamic_str(5),filename,variable + "_Ddd_dyn");
	load_str(dynamic_str(6),filename,variable + "_Dud_dyn");
	load_str(static_str(0),filename,variable + "_Puu_stat");
	load_str(static_str(1),filename,variable + "_Pdd_stat");
	load_str(static_str(2),filename,variable + "_Pud_stat");
	load_str(static_str(3),filename,variable + "_Xud_stat");
	load_str(static_str(4),filename,variable + "_Duu_stat");
	load_str(static_str(5),filename,variable + "_Ddd_stat");
	load_str(static_str(6),filename,variable + "_Dud_stat");
}

Ex_Generalmatrix Ex_Generalmatrix::operator+ (const Ex_Generalmatrix & gm2){
	Ex_Generalmatrix tmp;
	tmp.self_str = self_str + gm2.self_str;
	tmp.dynamic_str = dynamic_str + gm2.dynamic_str;
	tmp.static_str = static_str + gm2.static_str;
	return tmp;
}

double Ex_Generalmatrix::errnorm(double atol,double rtol, Ex_Generalmatrix &y1, Ex_Generalmatrix &y2){
	double err, err_self, err_dyn, err_stat;
	err_self=self_str.errnorm(atol,rtol,y1.self_str,y2.self_str);
	err_dyn=dynamic_str.errnorm(atol,rtol,y1.dynamic_str,y2.dynamic_str);
	err_stat=static_str.errnorm(atol,rtol,y1.static_str,y2.static_str);
	err=sqrt((err_self*err_self+err_dyn*err_dyn+err_stat*err_stat)/3.0);
	return err;
}

Ex_Generalmatrix operator*(const double &a, Ex_Generalmatrix &gm){
	Ex_Generalmatrix tmp;
	tmp.self_str = a*gm.self_str;
	tmp.dynamic_str = a*gm.dynamic_str;
	tmp.static_str = a*gm.static_str;
	return tmp;
}


#endif
