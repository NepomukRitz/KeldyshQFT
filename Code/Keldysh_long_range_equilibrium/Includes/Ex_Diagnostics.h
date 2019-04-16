#ifndef EX_DIAGNOSTICS_01022019
#define EX_DIAGNOSTICS_01022019

#include <list>
#include <string>

#include "matrix.h" 
#include "Ex_functions.h" 

using namespace std;

class Ex_Diagnostics{
	public:
		list<double> Lambda_values;
		list<double> x_values;
		list<double> Precomputation_time;
		list<double> Preintegration_time;
		list<double> ePuu_bubble_time; 
		list<double> ePdd_bubble_time; 
		list<double> ePud_bubble_time; 
		list<double> ePdu_bubble_time; 
		list<double> eXuu_bubble_time; 
		list<double> eXdd_bubble_time; 
		list<double> eXud_bubble_time; 
		list<double> eXdu_bubble_time; 
		list<double> P_rhs_time;
		list<double> XD_rhs_time;
		list<double> Self_rhs_time;
		Ex_Diagnostics(){};
		void feed_P_bubble(bool spin1, bool spin2, double time);
		void feed_X_bubble(bool spin1, bool spin2, double time);
		void save(string filename);
		void load(string filename);
};

void Ex_Diagnostics::feed_P_bubble(bool spin1, bool spin2, double time){
	if(spin1==1 && spin2==1){
		ePuu_bubble_time.push_back(time);
	}
	else if(spin1==1 && spin2==0){
		ePud_bubble_time.push_back(time);
	}
	else if(spin1==0 && spin2==1){
		ePdu_bubble_time.push_back(time);
	}
	else{
		ePdd_bubble_time.push_back(time);
	}
}

void Ex_Diagnostics::feed_X_bubble(bool spin1, bool spin2, double time){
	if(spin1==1 && spin2==1){
		eXuu_bubble_time.push_back(time);
	}
	else if(spin1==1 && spin2==0){
		eXud_bubble_time.push_back(time);
	}
	else if(spin1==0 && spin2==1){
		eXdu_bubble_time.push_back(time);
	}
	else{
		eXdd_bubble_time.push_back(time);
	}
}

void Ex_Diagnostics::save(string filename){
	matrix<double> Lambda_values_matrix  = list_to_matrix(Lambda_values);
	matrix<double> x_values_matrix  = list_to_matrix(x_values);
	matrix<double> Precomputation_time_matrix  = list_to_matrix(Precomputation_time);
	matrix<double> Preintegration_time_matrix  = list_to_matrix(Preintegration_time);
	matrix<double> ePuu_bubble_time_matrix  = list_to_matrix(ePuu_bubble_time);
	matrix<double> ePdd_bubble_time_matrix  = list_to_matrix(ePdd_bubble_time);
	matrix<double> ePud_bubble_time_matrix  = list_to_matrix(ePud_bubble_time);
	matrix<double> ePdu_bubble_time_matrix  = list_to_matrix(ePdu_bubble_time);
	matrix<double> eXuu_bubble_time_matrix  = list_to_matrix(eXuu_bubble_time);
	matrix<double> eXdd_bubble_time_matrix  = list_to_matrix(eXdd_bubble_time);
	matrix<double> eXud_bubble_time_matrix  = list_to_matrix(eXud_bubble_time);
	matrix<double> eXdu_bubble_time_matrix  = list_to_matrix(eXdu_bubble_time);
	matrix<double> P_rhs_time_matrix  = list_to_matrix(P_rhs_time);
	matrix<double> XD_rhs_time_matrix  = list_to_matrix(XD_rhs_time);
	matrix<double> Self_rhs_time_matrix  = list_to_matrix(Self_rhs_time);

	string variable = "Diagnostics";
	Lambda_values_matrix.save(filename.c_str(),(variable + "_Lambda_values").c_str());
	x_values_matrix.save(filename.c_str(),(variable + "_x_values").c_str());
	Precomputation_time_matrix.save(filename.c_str(),(variable + "_Precomputation_time").c_str());
	Preintegration_time_matrix.save(filename.c_str(),(variable + "_Preintegration_time").c_str());
	ePuu_bubble_time_matrix.save(filename.c_str(),(variable + "_ePuu_bubble_time").c_str());
	ePdd_bubble_time_matrix.save(filename.c_str(),(variable + "_ePdd_bubble_time").c_str());
	ePud_bubble_time_matrix.save(filename.c_str(),(variable + "_ePud_bubble_time").c_str());
	ePdu_bubble_time_matrix.save(filename.c_str(),(variable + "_ePdu_bubble_time").c_str());
	eXuu_bubble_time_matrix.save(filename.c_str(),(variable + "_eXuu_bubble_time").c_str());
	eXdd_bubble_time_matrix.save(filename.c_str(),(variable + "_eXdd_bubble_time").c_str());
	eXud_bubble_time_matrix.save(filename.c_str(),(variable + "_eXud_bubble_time").c_str());
	eXdu_bubble_time_matrix.save(filename.c_str(),(variable + "_eXdu_bubble_time").c_str());
	P_rhs_time_matrix.save(filename.c_str(),(variable + "_P_rhs_time").c_str());
	XD_rhs_time_matrix.save(filename.c_str(),(variable + "_XD_rhs_time").c_str());
	Self_rhs_time_matrix.save(filename.c_str(),(variable + "_Self_rhs_time").c_str());
}

void Ex_Diagnostics::load(string filename){
	string variable = "Diagnostics";

	matrix<double> Lambda_values_matrix;
	Lambda_values_matrix.load(filename.c_str(),(variable + "_Lambda_values").c_str());
	Lambda_values = matrix_to_list(Lambda_values_matrix);

	matrix<double> x_values_matrix;
	x_values_matrix.load(filename.c_str(),(variable + "_x_values").c_str());
	x_values = matrix_to_list(x_values_matrix);
	
	matrix<double> Precomputation_time_matrix;
	Precomputation_time_matrix.load(filename.c_str(),(variable + "_Precomputation_time").c_str());
	Precomputation_time = matrix_to_list(Precomputation_time_matrix);
	
	matrix<double> Preintegration_time_matrix;
	Preintegration_time_matrix.load(filename.c_str(),(variable + "_Preintegration_time").c_str());
	Preintegration_time = matrix_to_list(Preintegration_time_matrix);

	matrix<double> ePuu_bubble_time_matrix;
	ePuu_bubble_time_matrix.load(filename.c_str(),(variable + "_ePuu_bubble_time").c_str());
	ePuu_bubble_time = matrix_to_list(ePuu_bubble_time_matrix);

	matrix<double> ePdd_bubble_time_matrix;
	ePdd_bubble_time_matrix.load(filename.c_str(),(variable + "_ePdd_bubble_time").c_str());
	ePdd_bubble_time = matrix_to_list(ePdd_bubble_time_matrix);

	matrix<double> ePud_bubble_time_matrix;
	ePud_bubble_time_matrix.load(filename.c_str(),(variable + "_ePud_bubble_time").c_str());
	ePud_bubble_time = matrix_to_list(ePud_bubble_time_matrix);

	matrix<double> ePdu_bubble_time_matrix;
	ePdu_bubble_time_matrix.load(filename.c_str(),(variable + "_ePdu_bubble_time").c_str());
	ePdu_bubble_time = matrix_to_list(ePdu_bubble_time_matrix);

	matrix<double> eXuu_bubble_time_matrix;
	eXuu_bubble_time_matrix.load(filename.c_str(),(variable + "_eXuu_bubble_time").c_str());
	eXuu_bubble_time = matrix_to_list(eXuu_bubble_time_matrix);

	matrix<double> eXdd_bubble_time_matrix;
	eXdd_bubble_time_matrix.load(filename.c_str(),(variable + "_eXdd_bubble_time").c_str());
	eXdd_bubble_time = matrix_to_list(eXdd_bubble_time_matrix);

	matrix<double> eXud_bubble_time_matrix;
	eXud_bubble_time_matrix.load(filename.c_str(),(variable + "_eXud_bubble_time").c_str());
	eXud_bubble_time = matrix_to_list(eXud_bubble_time_matrix);

	matrix<double> eXdu_bubble_time_matrix;
	eXdu_bubble_time_matrix.load(filename.c_str(),(variable + "_eXdu_bubble_time").c_str());
	eXdu_bubble_time = matrix_to_list(eXdu_bubble_time_matrix);
	
	matrix<double> P_rhs_time_matrix;
	P_rhs_time_matrix.load(filename.c_str(),(variable + "_P_rhs_time").c_str());
	P_rhs_time = matrix_to_list(P_rhs_time_matrix);
	
	matrix<double> XD_rhs_time_matrix;
	XD_rhs_time_matrix.load(filename.c_str(),(variable + "_XD_rhs_time").c_str());
	XD_rhs_time = matrix_to_list(XD_rhs_time_matrix);
	
	matrix<double> Self_rhs_time_matrix;
	Self_rhs_time_matrix.load(filename.c_str(),(variable + "_Self_rhs_time").c_str());
	Self_rhs_time = matrix_to_list(Self_rhs_time_matrix);
}

#endif

