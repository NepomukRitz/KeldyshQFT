#ifndef GENERALMATRIX_17032017
#define GENERALMATRIX_17032017

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "Blockmatrix.h"
#include "Numerics.h"

using namespace std;


class Generalmatrix{
	public:
		Numerics num;
		matrix < matrix < syma < complex < double > > > > short_str;
		matrix < matrix < matrix < double > > > long_str;
		Generalmatrix(){};
		Generalmatrix(Numerics &num_in);
		void initialize(double init);
		void initialize_random();
		void resize(Numerics &num_in);
		void resize(Generalmatrix &gm);
		void save(char *filename, const char *variable_in);
		void load(char *filename, const char *variable_in);
		double errnorm(double atol,double rtol, Generalmatrix &y1, Generalmatrix &y2); 
		Generalmatrix operator+ (const Generalmatrix & gm2);
};

Generalmatrix::Generalmatrix(Numerics &num_in): num(num_in), short_str(9), long_str(7){
	short_str(0).resize(num.NfbP);
	short_str(1).resize(num.NfbP);
	short_str(2).resize(num.NfbP);
	for(int i=0; i<num.NfbP; ++i){
		short_str(0)(i).resize(num.Nges);
		short_str(1)(i).resize(num.Nges);
		short_str(2)(i).resize(num.Nges);
	}
	short_str(3).resize(num.NfbX);
	short_str(4).resize(num.NfbX);
	short_str(5).resize(num.NfbX);
	short_str(6).resize(num.NfbX);
	for(int i=0; i<num.NfbX; ++i){
		short_str(3)(i).resize(num.Nges);
		short_str(4)(i).resize(num.Nges);
		short_str(5)(i).resize(num.Nges);
		short_str(6)(i).resize(num.Nges);
	}
	short_str(7).resize(num.Nff);
	short_str(8).resize(num.Nff);
	for(int i=0; i<num.Nff; ++i){
		short_str(7)(i).resize(num.Nges);
		short_str(8)(i).resize(num.Nges);
	}
	for(int i=0;i<7;++i){
		Blockmatrix<double> block(num.L, num.N, long_str(i));
		block.resize(num.L, num.N);
	}
}

void Generalmatrix::initialize(double init){
	for(int i=0; i<num.NfbP; ++i){
		short_str(0)(i) = (complex<double>) init;
		short_str(1)(i) = (complex<double>) init;
		short_str(2)(i) = (complex<double>) init;
	}
	for(int i=0; i<num.NfbX; ++i){
		short_str(3)(i) = (complex<double>) init;
		short_str(4)(i) = (complex<double>) init;
		short_str(5)(i) = (complex<double>) init;
		short_str(6)(i) = (complex<double>) init;
	}
	for(int i=0; i<num.Nff; ++i){
		short_str(7)(i) = (complex<double>) init;
		short_str(8)(i) = (complex<double>) init;
	}
	for(int i=0; i<7; ++i){
		Blockmatrix<double> block(num.L, num.N, long_str(i));
		block.initialize(num.L, num.N, init);
	}
}

void Generalmatrix::initialize_random(){
	for(int i=0; i<num.NfbP; ++i){
		short_str(0)(i) = (complex<double>) (rand()%100)/100.;
		short_str(1)(i) = (complex<double>) (rand()%100)/100.;
		short_str(2)(i) = (complex<double>) (rand()%100)/100.;
	}
	for(int i=0; i<num.NfbX; ++i){
		short_str(3)(i) = (complex<double>) (rand()%100)/100.;
		short_str(4)(i) = (complex<double>) (rand()%100)/100.;
		short_str(5)(i) = (complex<double>) (rand()%100)/100.;
		short_str(6)(i) = (complex<double>) (rand()%100)/100.;
	}
	for(int i=0; i<num.Nff; ++i){
		short_str(7)(i) = (complex<double>) (rand()%100)/100.;
		short_str(8)(i) = (complex<double>) (rand()%100)/100.;
	}
	for(int i=0; i<7; ++i){
		Blockmatrix<double> block(num.L, num.N, long_str(i));
		block.initialize(num.L, num.N, (rand()%100)/100. );
	}
}

void Generalmatrix::resize(Numerics &num_in){
	num=num_in;
	short_str.resize(9); 
	long_str.resize(7); 
	short_str(0).resize(num.NfbP);
	short_str(1).resize(num.NfbP);
	short_str(2).resize(num.NfbP);
	for(int i=0; i<num.NfbP; ++i){
		short_str(0)(i).resize(num.Nges);
		short_str(1)(i).resize(num.Nges);
		short_str(2)(i).resize(num.Nges);
	}
	short_str(3).resize(num.NfbX);
	short_str(4).resize(num.NfbX);
	short_str(5).resize(num.NfbX);
	short_str(6).resize(num.NfbX);
	for(int i=0; i<num.NfbX; ++i){
		short_str(3)(i).resize(num.Nges);
		short_str(4)(i).resize(num.Nges);
		short_str(5)(i).resize(num.Nges);
		short_str(6)(i).resize(num.Nges);
	}
	
	short_str(7).resize(num.Nff);
	short_str(8).resize(num.Nff);
	for(int i=0; i<num.Nff; ++i){
		short_str(7)(i).resize(num.Nges);
		short_str(8)(i).resize(num.Nges);
	}
	
	for(int i=0;i<7;++i){
		Blockmatrix<double> block(num.L, num.N, long_str(i));
		block.resize(num.L, num.N);
	}
}

void Generalmatrix::resize(Generalmatrix &gm){
	resize(gm.num);
}

void Generalmatrix::save(char *filename, const char *variable_in){ 
	char variable[255]="";
	char variable2[255]="";
	strcat(variable,variable_in);
	strcpy(variable2, variable);
	strcat(variable2,"_short_str");
	short_str.save(filename,variable2);
	//strcat(variable2,"_long_str");
	//long_str.save(filename,variable2); //Caveat: The save for mmmd does not support different sizes of submatrices yet.
	strcpy(variable2, variable);
	strcat(variable2,"_Puu");
	Blockmatrix<double> Puu(num.L, num.N, long_str(0));
	Puu.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Pdd");
	Blockmatrix<double> Pdd(num.L, num.N, long_str(1));
	Pdd.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Pud");
	Blockmatrix<double> Pud(num.L, num.N, long_str(2));
	Pud.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Xud");
	Blockmatrix<double> Xud(num.L, num.N, long_str(3));
	Xud.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Duu");
	Blockmatrix<double> Duu(num.L, num.N, long_str(4));
	Duu.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Ddd");
	Blockmatrix<double> Ddd(num.L, num.N, long_str(5));
	Ddd.save(filename, variable2);
	strcpy(variable2, variable);
	strcat(variable2,"_Dud");
	Blockmatrix<double> Dud(num.L, num.N, long_str(6));
	Dud.save(filename, variable2);
}

//Caveat: This load function does not check if data are consistent with num parameters!
void Generalmatrix::load(char *filename, const char *variable_in){
	char variable[255]="";
	char variable2[255]="";
	strcat(variable,variable_in);
	strcpy(variable2, variable);
	strcat(variable2,"_short_str");
	short_str.load(filename,variable2);
	
	matrix<double> tmp;
	matrix<matrix<double> > tmp2;
	Blockmatrix<double> tmp2_block(num.L,num.N,tmp2);
	
	strcpy(variable2, variable);
	strcat(variable2,"_Puu");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(0) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Pdd");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(1) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Pud");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(2) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Xud");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(3) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Duu");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(4) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Ddd");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(5) = tmp2;
	
	strcpy(variable2, variable);
	strcat(variable2,"_Dud");
	tmp.load(filename,variable2);
	tmp2_block.convert_to_blockmatrix(num.L,num.N,tmp);
	long_str(6) = tmp2;
}


Generalmatrix Generalmatrix::operator+ (const Generalmatrix & gm2){
	Generalmatrix temp(num);
	temp.short_str=short_str+gm2.short_str;
	temp.long_str=long_str+gm2.long_str;
	return temp;
}


double Generalmatrix::errnorm(double atol,double rtol, Generalmatrix &y1, Generalmatrix &y2) {
	double err, err_short, err_long;
	err_short=short_str.errnorm(atol,rtol,y1.short_str,y2.short_str);
	err_long=long_str.errnorm(atol,rtol,y1.long_str,y2.long_str);
	err=sqrt((err_short*err_short+err_long*err_long)/2.0);
	return err;
}


Generalmatrix operator*(const double &a, Generalmatrix &gm){
	Generalmatrix temp(gm.num);
	temp.short_str=a*gm.short_str;
	temp.long_str=a*gm.long_str;
	return temp;
}








#endif
 
