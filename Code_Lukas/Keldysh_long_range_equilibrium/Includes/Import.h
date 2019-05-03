#ifndef IMPORT_13042016
#define IMPORT_13042016

#include <string.h>

using namespace std;

double get_option(int inputN,const char *inputV[],char *was){
	int n;
	char option[30];
	sprintf(option,"-%s",was);
	for (n=1;n<(inputN-1);n++){
		if (strcmp(inputV[n],option)==0)
		return (double) atof(inputV[n+1]);
	}
	return 0;
}

#endif
