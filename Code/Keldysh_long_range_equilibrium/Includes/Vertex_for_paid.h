#ifndef VERTEX_FOR_PAID_17052017
#define VERTEX_FOR_PAID_17052017

#include "Generalmatrix.h"

template <int mode> class Vertex{
	public:
		Numerics & num;
		Substitution<mode> &sub;
		matrix<double> wf_subst;
		Generalmatrix & data;
		matrix<syma<complex<double> > > &aPuu_central;
		matrix<syma<complex<double> > > &aPdd_central;
		matrix<syma<complex<double> > > &aPud_central;
		matrix<syma<complex<double> > > &aXud_central;
		matrix<syma<complex<double> > > &aDuu_central;
		matrix<syma<complex<double> > > &aDdd_central;
		matrix<syma<complex<double> > > &aDud_central;
		matrix<syma<complex<double> > > &ERetu;
		matrix<syma<complex<double> > > &ERetd;
		linear_ipol_bin<syma<complex<double> > > aPuu_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aPdd_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aPud_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aXud_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aDuu_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aDdd_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > aDud_central_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol; 
		linear_ipol_bin<syma<complex<double> > > ERetu_ipol_subst; 
		linear_ipol_bin<syma<complex<double> > > ERetd_ipol_subst; 
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
		complex<double> interpolate_componentwise(double freq, int j1, int j2, int k, matrix<double> &w, matrix<syma<complex<double> > > &a);
		complex<double> interpolate_componentwise_aPuu_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aPdd_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aPud_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aXud_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aDuu_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aDdd_central(double freq, int j1, int j2);
		complex<double> interpolate_componentwise_aDud_central(double freq, int j1, int j2);
		Vertex(Numerics &num_in, Substitution<mode> &sub_in, Generalmatrix &data_in);
};

template<int mode> Vertex<mode>::Vertex(Numerics &num_in,
                                        Substitution<mode> &sub_in,
                                        Generalmatrix &data_in):
                                        num(num_in),
                                        sub(sub_in),
                                        wf_subst(num.Nff),
                                        data(data_in),
                                        aPuu_central(data.short_str(0)),
                                        aPdd_central(data.short_str(1)),
                                        aPud_central(data.short_str(2)),
                                        aXud_central(data.short_str(3)),
                                        aDuu_central(data.short_str(4)),
                                        aDdd_central(data.short_str(5)),
                                        aDud_central(data.short_str(6)),
                                        ERetu(data.short_str(7)),
                                        ERetd(data.short_str(8)),
                                        aPuu_central_ipol(num.wbP, aPuu_central),
                                        aPdd_central_ipol(num.wbP, aPdd_central),
                                        aPud_central_ipol(num.wbP, aPud_central),
                                        aXud_central_ipol(num.wbX, aXud_central),
                                        aDuu_central_ipol(num.wbX, aDuu_central),
                                        aDdd_central_ipol(num.wbX, aDdd_central),
                                        aDud_central_ipol(num.wbX, aDud_central),
                                        ERetu_ipol(num.wf, data.short_str(7)),
                                        ERetd_ipol(num.wf, data.short_str(8)),
                                        ERetu_ipol_subst(wf_subst, data.short_str(7)),
                                        ERetd_ipol_subst(wf_subst, data.short_str(8)),
                                        aPuu_feedback_data(data.long_str(0)),
                                        aPdd_feedback_data(data.long_str(1)),
                                        aPud_feedback_data(data.long_str(2)),
                                        aXud_feedback_data(data.long_str(3)),
                                        aDuu_feedback_data(data.long_str(4)),
                                        aDdd_feedback_data(data.long_str(5)),
                                        aDud_feedback_data(data.long_str(6)),
                                        aPuu_feedback(num.L, num.N, aPuu_feedback_data),
                                        aPdd_feedback(num.L, num.N, aPdd_feedback_data),
                                        aPud_feedback(num.L, num.N, aPud_feedback_data),
                                        aXud_feedback(num.L, num.N, aXud_feedback_data),
                                        aDuu_feedback(num.L, num.N, aDuu_feedback_data),
                                        aDdd_feedback(num.L, num.N, aDdd_feedback_data),
                                        aDud_feedback(num.L, num.N, aDud_feedback_data){
	for(int i=0; i<num.Nff; ++i){
		wf_subst(i) = sub.subst_concatenated(num.wf(i));
	}
} 
		

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise(double freq, int j1, int j2, int k, matrix<double> &w, matrix<syma<complex<double> > > &a){
	if (freq==w(k)){
		return a(k)(j1,j2);
	}
	else if (freq==w(k-1)){
		return a(k-1)(j1,j2);
	}
	return ((freq-w(k-1))/(w(k)-w(k-1)))*a(k)(j1,j2)+((w(k)-freq)/(w(k)-w(k-1)))*a(k-1)(j1,j2);
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aPuu_central(double freq, int j1, int j2){
	int k= num.interpol_P(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbP, aPuu_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aPdd_central(double freq, int j1, int j2){
	int k= num.interpol_P(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbP, aPdd_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aPud_central(double freq, int j1, int j2){
	int k= num.interpol_P(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbP, aPud_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aXud_central(double freq, int j1, int j2){
	int k= num.interpol_X(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbX, aXud_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aDuu_central(double freq, int j1, int j2){
	int k= num.interpol_X(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbX, aDuu_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aDdd_central(double freq, int j1, int j2){
	int k= num.interpol_X(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbX, aDdd_central); 
}

template<int mode> complex<double> Vertex<mode>::interpolate_componentwise_aDud_central(double freq, int j1, int j2){
	int k= num.interpol_X(freq);	
	return interpolate_componentwise(freq, j1, j2, k, num.wbX, aDud_central); 
}


#endif
