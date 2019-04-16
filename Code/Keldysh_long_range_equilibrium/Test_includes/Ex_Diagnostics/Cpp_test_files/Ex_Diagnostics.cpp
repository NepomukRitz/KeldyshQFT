#define ACCURACY_P_BUB 1e-6
#define ACCURACY_X_BUB 1e-6
#define TEMPORAER_1 1
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_Diagnostics.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	Ex_Diagnostics diagnostics;
	double Lambda = 1e-2;
	double x = -15.0;
	double time_pre=5.0;
	double time_ePuu= 10.0;
	double time_ePdd= 19.0;
	double time_ePud= 15.0;
	double time_ePdu= 13.0;
	double time_eXuu= 10.7;
	double time_eXdd= 19.3;
	double time_eXud= 15.2;
	double time_eXdu= 13.1;
	double time_self = 40.1;
	double time_p = 30.1;
	double time_x = 50.1;
	diagnostics.Lambda_values.push_back(Lambda);
	diagnostics.x_values.push_back(x);
	diagnostics.Precomputation_time.push_back(time_pre);
	diagnostics.feed_P_bubble(1,1,time_ePuu);
	diagnostics.feed_P_bubble(0,0,time_ePdd);
	diagnostics.feed_P_bubble(1,0,time_ePud);
	diagnostics.feed_P_bubble(0,1,time_ePdu);
	diagnostics.feed_X_bubble(1,1,time_eXuu);
	diagnostics.feed_X_bubble(0,0,time_eXdd);
	diagnostics.feed_X_bubble(1,0,time_eXud);
	diagnostics.feed_X_bubble(0,1,time_eXdu);
	diagnostics.Self_rhs_time.push_back(time_self);
	diagnostics.P_rhs_time.push_back(time_p);
	diagnostics.XD_rhs_time.push_back(time_x);
	string filename ="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_Diagnostics"; 
	diagnostics.save(filename);

	Ex_Diagnostics diagnostics_loaded;
	diagnostics_loaded.load(filename);
	cout<<"diagnostics_loaded.Lambda_values="<<diagnostics_loaded.Lambda_values<<endl;
	cout<<"diagnostics_loaded.x_values="<<diagnostics_loaded.x_values<<endl;
	cout<<"diagnostics_loaded.Precomputation_time="<<diagnostics_loaded.Precomputation_time<<endl;
	cout<<"diagnostics_loaded.ePuu_bubble_time="<<diagnostics_loaded.ePuu_bubble_time<<endl;
	cout<<"diagnostics_loaded.ePdd_bubble_time="<<diagnostics_loaded.ePdd_bubble_time<<endl;
	cout<<"diagnostics_loaded.ePud_bubble_time="<<diagnostics_loaded.ePud_bubble_time<<endl;
	cout<<"diagnostics_loaded.ePdu_bubble_time="<<diagnostics_loaded.ePdu_bubble_time<<endl;
	cout<<"diagnostics_loaded.eXuu_bubble_time="<<diagnostics_loaded.eXuu_bubble_time<<endl;
	cout<<"diagnostics_loaded.eXdd_bubble_time="<<diagnostics_loaded.eXdd_bubble_time<<endl;
	cout<<"diagnostics_loaded.eXud_bubble_time="<<diagnostics_loaded.eXud_bubble_time<<endl;
	cout<<"diagnostics_loaded.eXdu_bubble_time="<<diagnostics_loaded.eXdu_bubble_time<<endl;
	cout<<"diagnostics_loaded.P_rhs_time="<<diagnostics_loaded.P_rhs_time<<endl;
	cout<<"diagnostics_loaded.XD_rhs_time="<<diagnostics_loaded.XD_rhs_time<<endl;
	cout<<"diagnostics_loaded.Self_rhs_time="<<diagnostics_loaded.Self_rhs_time<<endl;
	
	return 0;
}

