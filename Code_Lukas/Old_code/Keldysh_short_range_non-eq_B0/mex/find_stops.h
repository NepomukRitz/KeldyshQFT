#ifndef FIND_STOPS_DLKJHALD
#define FIND_STOPS_DLKJHALD

#include <matrix.h>
#include <basic.h>
#include <complex>
#include <omp.h>

using namespace std;

template<class interpol>
matrix<double> find_stops(interpol G, matrix<double> &positions, double taul, int frequency_grid_size=10000, int max_num=1000)
{
	matrix<double> ret(max_num);
	matrix<double> freq(frequency_grid_size);
	matrix<double> val(frequency_grid_size);
	int num=0;
	for (int pos=0; pos<positions.dim_c; pos++)
	{
		int pos_val=(int)(positions(pos)+.5);
		omp_set_num_threads(16);
		#pragma omp parallel for
		for (int i=0; i<frequency_grid_size; i++)
		{
			freq(i)=-2.*taul+4.*taul*(double)i/(double)frequency_grid_size;
			val(i) = abs(imag(G(freq(i))(pos_val, pos_val)));
		}
		for (int i=1; i<frequency_grid_size-1; i++)
		{
			if(val(i)>val(i-1) && val(i)>val(i+1) && num<max_num)
			{
				ret(num)=freq(i);
				num++;
			}
		}
	}
	if (num==0)
		num = 1;
	matrix<double> stops(num);
	if (num==0)
		stops(0) = .0;
	else
	{
		for (int i=0; i<num; i++)
			stops(i)=ret(i);
		stops.sort();
	}
	return stops;
};

#endif
