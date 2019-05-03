EX_mpi.h: //Should be ok; maybe remove the volume property at job. 

int ex_mpi_determine_counts(matrix<int>  &job_list_volumes, int nprocs, matrix<int> &receive_count, matrix<int> &displs)

matrix<int> ex_job_displacement(int upper_boundary, int r, int nprocs, matrix<int> &job_list_volumes)

template <typename Tfund, typename Tj, typename Tc, typename Comp_obj> matrix<Tc> ex_mpi_computation(matrix<Tj> &job_list, Comp_obj &comp_obj)



Ex_Compute_bubble.h: //Organize the job lists new; Computation itself should be ok

int ex_determine_number_of_matrix_components(matrix<int> &L_structure)

matrix<matrix<double> > ex_job_list_for_dyn_bubble(matrix<int> &L_structure, matrix<double> &wb)

matrix<matrix<double> > ex_job_list_for_static_bubble(int L)

template <typename T> void ex_compute_dyn_bubble(Ex_freq_str &bubble, T &Bubble_dyn)

template <typename T> void ex_compute_static_bubble(Ex_freq_str &bubble, T &Bubble_stat)



Ex_freq_str.h: //Independent functions should be ok; Ex_freq_str could be transformed in Wrapper class

void ex_freq_str_resize(int L,
                        int N,
                        matrix<int> &L_structure,
                        matrix<matrix<matrix<complex<double> > > > &dynamic_str,
                        matrix<matrix<double> > &static_str)

void ex_freq_str_initialize(int L,
                            int N,
                            matrix<int> &L_structure,
                            matrix<matrix<matrix<complex<double> > > > &dynamic_str,
                            matrix<matrix<double> > &static_str, 
                            double init)

void ex_freq_str_save(int L,
                      int N,
                      matrix<int> &L_structure,
                      matrix<double> &wb,
                      matrix<matrix<matrix<complex<double> > > > &dynamic_str,
                      matrix<matrix<double> > &static_str,
                      string filename,
                      string variable)

void ex_freq_str_load(matrix<matrix<matrix<complex<double> > > > &dynamic_str,
                      matrix<matrix<double> > &static_str,
                      string filename,
                      string variable)

class Ex_freq_str{
	public:
		int L;
		int N;
		matrix<int> &L_structure;
		matrix<double> &wb;
		int Nges;
		int N_freq;
		matrix<matrix<matrix<complex<double> > > > dynamic_str;
		matrix<matrix<double> > static_str;
		Ex_freq_str(int L_in,
		            int N_in,
		            matrix<int> &L_structure_in,
		            matrix<double> &wb_in);
		void initialize(double init);
		void save(string filename, string variable);
		void load(string filename, string variable);
};



Ex_functions.h: //Contains the very generic functions; Should be extended if possible

template <typename T> double abs(matrix<T> A){

template <typename T> double abs(syma<T> A){



Ex_Generalmatrix.h: //Fundamental container class that is suitable for the ode-solver; Should be ok, however the save and load syntax in matrix and syma could be transformed to take strings. 

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



Ex_Vertex.h: //Wrapper class without any deep properties; Check if Blockmatrix-structure is actually used at any point. 

template <int mode> class Ex_Vertex{
	public:
		Numerics & num;
		Substitution<mode> sub;
		matrix<double> wf_subst;
		Ex_Generalmatrix &data;
		matrix<syma<complex<double> > > &ERetu;
		matrix<syma<complex<double> > > &ERetd;
		matrix<matrix<matrix<complex<double> > > > &aPuu_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aPdd_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aPud_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aXud_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aDuu_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aDdd_dynamic;
		matrix<matrix<matrix<complex<double> > > > &aDud_dynamic;
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol_subst; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol_subst; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aPuu_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aPdd_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aPud_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aXud_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aDuu_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aDdd_dynamic_ipol; 
		linear_ipol_bin<matrix<matrix<complex<double> > > > aDud_dynamic_ipol; 
		matrix<matrix< double > > &aPuu_feedback_data;
		matrix<matrix< double > > &aPdd_feedback_data;
		matrix<matrix< double > > &aPud_feedback_data;
		matrix<matrix< double > > &aXud_feedback_data;
		matrix<matrix< double > > &aDuu_feedback_data;
		matrix<matrix< double > > &aDdd_feedback_data;
		matrix<matrix< double > > &aDud_feedback_data;
		Blockmatrix<double> aPuu_feedback;
		Blockmatrix<double> aPdd_feedback;
		Blockmatrix<double> aPud_feedback;
		Blockmatrix<double> aXud_feedback;
		Blockmatrix<double> aDuu_feedback;
		Blockmatrix<double> aDdd_feedback;
		Blockmatrix<double> aDud_feedback;
		Ex_Vertex(Numerics &num_in, Substitution<mode> &sub_in, Ex_Generalmatrix &data_in);
};



Ex_multiplication.h: //Computes various products between frequency structures (mpi+omp or omp only).

int determine_mult_count(matrix<matrix<matrix<complex<double> > > > &A)

int determine_mult_count(matrix<matrix<double> > &A)

matrix<matrix<int> > determine_mult_job_list(matrix<matrix<matrix<complex<double> > > > &A)

matrix<matrix<int> > determine_mult_job_list(matrix<matrix<double> > &A){

template <typename T> class Ex_mult{
	public:
		T &A;
		T &B;
		Ex_dyn_mult(T &A_in,
		            T &B_in);
		auto operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);
		int dim_c(matrix<int> &job);
		int volume(matrix<int> &job);
};

class Ex_exterior_mult : Ex_mult<matrix<matrix<double> > {
	public:
		int L_inner;
		matrix<double> operator()(matrix<int> &job);
};

matrix<matrix<matrix<complex<double> > > > mult_sort(matrix<matrix<complex<double> > > &res_mpi, Ex_mult &comp_obj)

matrix<matrix<double> > mult_sort(matrix<matrix<double> > &res_mpi, Ex_mult &comp_obj)

template<typename T> T mpi_mult(T &A, T &B)

template<typename T> T mpi_ext_mult(T &A, T &B, int L_inner)

matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<complex<double> > > > &A, matrix<matrix<matrix<complex<double> > > > &B)

matrix<matrix<double> > omp_mult(matrix<matrix<double> > &A, matrix<matrix<double> > &B)

matrix<matrix<double> > omp_ext_mult(matrix<matrix<double> > &A, matrix<matrix<double> > &B, int L_inner)

matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<complex<double> > > > &A, matrix<matrix<double> > &B)

matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<double> > &B, matrix<matrix<matrix<complex<double> > > > &A)
