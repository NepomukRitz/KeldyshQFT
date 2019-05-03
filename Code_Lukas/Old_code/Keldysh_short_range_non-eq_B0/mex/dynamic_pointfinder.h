#ifndef DYNAMIC_POINTFINDER_67THY_905_GH
#define DYNAMIC_POINTFINDER_67THY_905_GH

#include <matrix.h>
#include <omp.h>
#include <iostream>

#define DYNAMIC_POINTFINDER_EVALUATION_INFO_TH78_RT_56 0

using namespace std;

//TODO: OPTIMIZE

template <class fun, class val, class Double>
void dynamic_pointfinder(fun &f, matrix<Double> &xi, matrix<val> &yi, double err=1e-06, int max_points=1e6)
{
	int number_of_evaluations = 0;
	int N = xi.dim_c;
	matrix<int> good_interval(N-1);   //0: interval is guaranteed to be good; 1: inverval might be good; 2: interval is bad
	matrix<int> good_intervaln(N-1);  //0: interval is guaranteed to be good; 1: inverval might be good; 2: interval is bad
	matrix<Double> yn_x(N);
	matrix<val> yn_y(N);
	matrix<int> yn_prior_new_points(N);

	int number_of_new_points;

	//Compute the function at the points xn;
	#pragma omp parallel for
	for (int i=0; i<N; i++)
	{
		yn_x(i) = xi(i);
		yn_y(i) = f(yn_x(i));
		yn_prior_new_points(i) = 0;
	}

	yi = yn_y;

	//All intervals might be good
	for (int i=0; i<N-1; i++)
		good_interval(i) = 1;

	number_of_new_points = 1;
	while (number_of_new_points != 0 && N<max_points)
	{
		//Is the interval actually good?
		//TODO: If the following line (omp stuff) is omitted, the result is always the same. otherwise there are some fluctuations in the number of points computed...
		//#pragma omp parallel for
		for (int i=0; i<N-2; i++)
		{
			if (good_interval(i)==1)
			{
				//Double x = .5*(yn_x(i+2)+yn_x(i));
				//val yintermediate = (yn_y(i)*(x-yn_x(i+2))-yn_y(i+2)*(x-yn_x(i)))/(yn_x(i)-yn_x(i+2));
				double errloc = (1./(yn_x(i)-yn_x(i+2))*(f.select(yn_y(i))*(yn_x(i+1)-yn_x(i+2))-f.select(yn_y(i+2))*(yn_x(i+1)-yn_x(i)))-f.select(yn_y(i+1))).mabs().mmax();
				if (errloc > err)
				{
					//#pragma omp critical
					good_interval(i)   = 2;
					//#pragma omp critical
					good_interval(i+1) = 2;
				}
			}
		}
		//How many new points have to be added? How much does each single_point need to be shifted?
		number_of_new_points = 0;

		yn_prior_new_points(0) = 0;
		for (int i=0; i<N-1; i++)
		{
			if (good_interval(i)==2)
			{
				number_of_new_points++;
			}
			else
			{
				good_interval(i)=0;
			}
			yn_prior_new_points(i+1) = number_of_new_points;
		}
		
		if (number_of_new_points != 0)
		{
			good_intervaln.resize(N+number_of_new_points-1);

			//Assume all new intervals will be good, pick out the ones for recomputation later on
			for (int i=0; i<N+number_of_new_points-1; i++)
			{
				good_intervaln(i) = 0;
			}

			//Insert points
			xi.resize(N+number_of_new_points);
			yi.resize(N+number_of_new_points);
			xi(0) = yn_x(0);
			yi(0) = yn_y(0);
			#pragma omp parallel for
			for (int i=0; i<N-1; i++)
			{
				if(good_interval(i)==0)
				{
					xi(i+yn_prior_new_points(i+1)+1) = yn_x(i+1);
					yi(i+yn_prior_new_points(i+1)+1) = yn_y(i+1);
				}
				else
				{
					xi(i+yn_prior_new_points(i+1))   = .5*(yn_x(i+1)+yn_x(i));
					yi(i+yn_prior_new_points(i+1))   = f(xi(i+yn_prior_new_points(i+1)));
					xi(i+yn_prior_new_points(i+1)+1) = yn_x(i+1);
					yi(i+yn_prior_new_points(i+1)+1) = yn_y(i+1);
				}
			}
			

			for (int i=0; i<N-1; i++)
			{
				if(good_interval(i)!=0)
				{
				#if DYNAMIC_POINTFINDER_EVALUATION_INFO_TH78_RT_56
					number_of_evaluations++;
				#endif
					good_intervaln(i+yn_prior_new_points(i))   = 1;
					good_intervaln(i+yn_prior_new_points(i)+1) = 1;
					if (i!=0)
					{
						good_intervaln(i+yn_prior_new_points(i)-1) = 1;
					}
					if (i!=N-2)
					{
						good_intervaln(i+yn_prior_new_points(i)+2) = 1;
					}
				}
			}
			N += number_of_new_points;
			yn_x = xi;
			yn_y = yi;
			good_interval = good_intervaln;
			yn_prior_new_points.resize(N);
		}
	}
#if DYNAMIC_POINTFINDER_EVALUATION_INFO_TH78_RT_56
	cout << number_of_evaluations << endl;
#endif
};

#endif
