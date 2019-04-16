#ifndef EX_VERTEX_30102018
#define EX_VERTEX_30102018


#include "Ex_freq_str.h"
#include "Ex_Generalmatrix.h"
#include "Substitution.h"

template <int mode> class Ex_Vertex{
	public:
		Numerics &num;
		Substitution<mode> &sub;
		matrix<double> wf_subst;
		Ex_Generalmatrix &data;
		matrix<syma<complex<double> > > &ERetu;
		matrix<syma<complex<double> > > &ERetd;
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol_subst; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol_subst; 
		Ex_freq_str aPuu;
		Ex_freq_str aPdd;
		Ex_freq_str aPud;
		Ex_freq_str aXud;
		Ex_freq_str aDuu;
		Ex_freq_str aDdd;
		Ex_freq_str aDud;
		Ex_Vertex(Numerics &num_in, Substitution<mode> &sub_in, Ex_Generalmatrix &data_in);
		void save(string filename, string variable);
		void load(string filename, string variable);
};

template<int mode> Ex_Vertex<mode>::Ex_Vertex(Numerics &num_in,
                                              Substitution<mode> &sub_in,
                                              Ex_Generalmatrix &data_in):
                                              num(num_in),
                                              sub(sub_in),
                                              data(data_in),
                                              wf_subst(num.Nff),
                                              ERetu(data.self_str(0)),
                                              ERetd(data.self_str(1)),
                                              ERetu_ipol(num.wf, data.self_str(0)),
                                              ERetd_ipol(num.wf, data.self_str(1)),
                                              ERetu_ipol_subst(wf_subst, data.self_str(0)),
                                              ERetd_ipol_subst(wf_subst, data.self_str(1)),
                                              aPuu(num.L, num.N, num.Lp_structure, num.wbP, data.dynamic_str(0), data.static_str(0)),
                                              aPdd(num.L, num.N, num.Lp_structure, num.wbP, data.dynamic_str(1), data.static_str(1)),
                                              aPud(num.L, num.N, num.Lp_structure, num.wbP, data.dynamic_str(2), data.static_str(2)),
                                              aXud(num.L, num.N, num.Lx_structure, num.wbX, data.dynamic_str(3), data.static_str(3)),
                                              aDuu(num.L, num.N, num.Lx_structure, num.wbX, data.dynamic_str(4), data.static_str(4)),
                                              aDdd(num.L, num.N, num.Lx_structure, num.wbX, data.dynamic_str(5), data.static_str(5)),
                                              aDud(num.L, num.N, num.Lx_structure, num.wbX, data.dynamic_str(6), data.static_str(6)){
	for(int i=0; i<num.Nff; ++i){
		wf_subst(i) = sub.subst_concatenated(num.wf(i));
	}
} 

template<int mode> void Ex_Vertex<mode>::save(string filename, string variable){
	ERetu.save(filename.c_str(),(variable + "_ERetu").c_str());
	ERetd.save(filename.c_str(),(variable + "_ERetd").c_str());
	aPuu.save(filename,variable + "_aPuu");
	aPdd.save(filename,variable + "_aPdd");
	aPud.save(filename,variable + "_aPud");
	aXud.save(filename,variable + "_aXud");
	aDuu.save(filename,variable + "_aDuu");
	aDdd.save(filename,variable + "_aDdd");
	aDud.save(filename,variable + "_aDud");
}

template<int mode> void Ex_Vertex<mode>::load(string filename, string variable){
	ERetu.load(filename.c_str(),(variable + "_ERetu").c_str());
	ERetd.load(filename.c_str(),(variable + "_ERetd").c_str());
	aPuu.load(filename,variable + "_aPuu");
	aPdd.load(filename,variable + "_aPdd");
	aPud.load(filename,variable + "_aPud");
	aXud.load(filename,variable + "_aXud");
	aDuu.load(filename,variable + "_aDuu");
	aDdd.load(filename,variable + "_aDdd");
	aDud.load(filename,variable + "_aDud");
}

#endif
