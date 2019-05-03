#ifndef RPA_VERTEX_04052017
#define RPA_VERTEX_04052017


template <int mode> class RPA_vertex{
 public:
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> sub;
 matrix<syma<complex<double> > > &Bubble_central; 
 matrix<double> wbXs;
 syma<complex<double> > U_x;
 linear_ipol_bin<syma<complex<double> > > iBubble_central;
 matrix<syma<complex<double> > > aXud_central;
 linear_ipol_bin<syma<complex<double> > > iaXud_central; 
 syma<complex<double> >  aXud_central_stat;
 matrix<syma<complex<double> > > aXud_central_dyn;
 linear_ipol_bin<syma<complex<double> > > iaXud_central_dyn; 

 RPA_vertex(Numerics &num, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, matrix<syma<complex<double> > > Bubble_central_in, Barevertex &bare);


};

template <int mode> RPA_vertex<mode>::RPA_vertex(Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, matrix<syma<complex<double> > > Bubble_central_in, Barevertex &bare): num(num_in), pre(pre_in), sub(sub_in), Bubble_central(Bubble_central_in), wbXs(num.NfbX), U_x(num.Nges), iBubble_central(wbXs, Bubble_central), aXud_central(num.num_freq_pre), iaXud_central(pre.freq_pre, aXud_central), aXud_central_dyn(num.num_freq_pre), iaXud_central_dyn(pre.freq_pre, aXud_central_dyn){

 for(int i=0; i<num.NfbX; ++i){
  wbXs(i) =sub.subst_concatenated(num.wbX(i)); 
 }

 for(int i=0; i<num.Nges; ++i){
  for(int j=0; j<=i; ++j){
   U_x(i,j) = 0.5*bare(i,1,j,0,j,1,i,0);
  }
 }
 
 omp_set_num_threads(16);
 #pragma omp parallel for
 for(int i=0; i<num.num_freq_pre; ++i){
  //cout<<"i="<<i<<endl;
  double freq = pre.freq_pre(i);
  matrix<complex<double> > Tmp;
  Tmp = -0.5*U_x*iBubble_central(freq);
  for(int j=0; j<num.Nges; ++j){
	Tmp(j,j) += (complex<double>)1.0;
  }
  Tmp.inv();
  aXud_central(i)= Tmp*U_x;
 }

 aXud_central_stat = aXud_central(0);
 
 omp_set_num_threads(16);
 #pragma omp parallel for
 for(int i=0; i<num.num_freq_pre; ++i){
  aXud_central_dyn(i) = aXud_central(i) - aXud_central_stat;
 }
}

#endif
