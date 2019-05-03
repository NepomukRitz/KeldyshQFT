General Notation:

If not specified differently, we think of the static structure
"matrix<matrix<T> > &A"
as Blockmatrix in position space. The outer matrix structure is given by the short indices, the inner matrix structure by the long indices: A(l,k)(j,i), with -L<=l,k<=L and max(0,-l)<=j<min(Nges,Nges-l), max(0,-k)<=i<min(Nges,Nges-k).

If not specified differently, we think of the dynamic structure
"matrix<matrix<matrix<T> > > &A"
as a vector (in freq_indices) of Blockmatrices in position space. 


template<typename T> int determine_mult_count(matrix<matrix<matrix<T> > > &A)
Computes the total amount of submatrices A(i,l,k) that comprise the dynamic structure A. The name "determine_mult_count" is chosen because this number is the number of queue elements in our parallelization-scheduler for dynamic structure multiplication.

int determine_mult_count(matrix<matrix<double> > &A)
Computes the total amount of submatrices A(l,k) that comprise the static structure A. The name "determine_mult_count" is chosen because this number is the number of queue elements in our parallelization-scheduler for static structure multiplication.

template<typename T> matrix<matrix<int> > determine_mult_job_list(matrix<matrix<matrix<T> > > &A)
Generates a job_list for dynamic structure multiplication. The output of this function is fed to the parallelization-scheduler.

template<typename T> matrix<matrix<int> > determine_mult_job_list(matrix<matrix<T> > &A)
Generates a job_list for static structure multiplication. The output of this function is fed to the parallelization-scheduler.

template <typename T1, typename T2> class Ex_mult_dyn{
	public:
		T1 &A;
		T2 &B;
		Ex_mult_dyn(T1 &A_in,
		            T2 &B_in);
		matrix<complex<double> >  operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);
		int dim_c(matrix<int> &job);
		int volume(matrix<int> &job);
};
Functor class that provides the multiplication function of dynamical structures for the parallelization-scheduler.
The operator() gets an element (i,l,k) of job_list and computes the corresponding matrix element
C(i)(l,k) = A(i)(l,q)B(i)(q,k), with short indices l,q,k that run between -L(i) and L(i), where L(i) is the frequency dependent feedback length. 

template<typename T> class Ex_mult_stat{
	public:
		matrix<matrix<T> > &A;
		matrix<matrix<T> > &B;
		Ex_mult_stat(matrix<matrix<T> > &A_in,
		             matrix<matrix<T> > &B_in);
		matrix<T>  operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);
		int dim_c(matrix<int> &job);
		int volume(matrix<int> &job);
};
Functor class that provides the multiplication function of static structures for the parallelization-scheduler.
The operator() gets an element (l,k) of job_list and computes the corresponding matrix element
C(l,k) = A(l,q)B(q,k), with short indices l,k that run between -L and L, where L is the static feedback length. 

template<typename T> class Ex_exterior_mult : public Ex_mult_stat<T>{
	public:
		int L_inner;
		Ex_exterior_mult(matrix<matrix<T> > &A_in,
		                 matrix<matrix<T> > &B_in,
		                 int L_inner_in);
		matrix<T> operator()(matrix<int> &job);
};
Derived functor class that provides the multiplication function of static structures for the parallelization-scheduler,
with summation index only running over exterior short indices: -L<=q<-L(i) and L(i)<q<=L.

template<typename T1, typename T2> matrix<matrix<matrix<complex<double> > > > mult_sort(matrix<matrix<complex<double> > > &res_mpi, Ex_mult_dyn<T1,T2> &comp_obj)
Takes the resulting list of matrices C(z(i,l,k)) computed with the parallelization scheduler and transforms it to the dynamic structure C(i)(l,k).

template<typename T> matrix<matrix<T> > mult_sort(matrix<matrix<T> > &res_mpi, Ex_mult_stat<T> &comp_obj)
Takes the resulting list of matrices C(z(l,k)) computed with the parallelization scheduler and transforms it to the static structure C(l,k).

template<typename T1, typename T2> matrix<matrix<matrix<complex<double> > > > mpi_mult(T1 &A, T2 &B)
Multiplies two dynamical structures using the parallelization scheduler (which uses both mpi & openmp). Concrete implementations for T1 and T2 follow below. 
