#include <mex.h>
#include <iostream>
#include <complex>
#include <matrix.h>

void mu_cdv_ml2c(vector<complex<double> > &c,const mxArray *ml,int &nf)
{
	int k;
	double *mlr, *mli;
	nf=mxGetNumberOfElements(ml);
	c.resize(nf);
	mlr=mxGetPr(ml);
	if (mxIsComplex(ml))
	{
		mli=mxGetPi(ml);
		for (k=0;k<nf;k++)
		{
			real(c[k]) = mlr[k];
			imag(c[k]) = mli[k];
		}
	}
	else
		for (k=0;k<nf;k++)
		{
			real(c[k]) = mlr[k];
			imag(c[k]) = 0.0;
		}
}

void mu_dv_ml2c(vector<double> &c,const mxArray *ml,int &nf)
{
	int k;
	double *mlr;
	nf=mxGetNumberOfElements(ml);
	c.resize(nf);
	mlr=mxGetPr(ml);
    for (k=0;k<nf;k++)
	{
		c[k] = mlr[k];
	}
}

mxArray * operator= (mxArray &*ml,matrix<double> &c) {
 if (!mxIsEmpty(ml)) {
  mxDestroyArray(ml);
 }
 double *mlr;
 mwSize ndim=2, dims[2];
 dims[0]=c.dim_r;
 dims[1]=c.dim_c;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
 mlr=mxGetPr(ml);
 for (int i=0;i<c.dim_r;i++)
  for (int j=0;j<c.dim_c;j++)
   mlr[i+j*c.dim_r] = c[i][j].real;
 return ml;
}

void mu_c2ml(matrix<double> &c,mxArray *ml[])
{
	double *mlr;
	ml[0]=mxCreateDoubleMatrix((mwSize) c.dim_r,(mwSize) c.dim_c,mxREAL);
	mlr=mxGetPr(ml[0]);
	for (int i=0;i<c.dim_r;i++)
     for (int j=0;j<c.dim_c;j++)
	  mlr[i+j*c.dim_r] = c[i][j].real;
}

void mu_c2ml(matrix<complex<double> > &c,mxArray *ml[])
{
	double *mlr, *mli;
	ml[0]=mxCreateDoubleMatrix((mwSize) c.dim_r,(mwSize) c.dim_c,mxCOMPLEX);
	mlr=mxGetPr(ml[0]);
	mli=mxGetPi(ml[0]);
	for (int i=0;i<c.dim_r;i++)
     for (int j=0;j<c.dim_c;j++)
	 {
	  mlr[i+j*c.dim_r] = c[i][j].real;
	  mli[i+j*c.dim_r] = c[i][j].imag;
	 }
}

void mu_dv_c2ml(vector<double> &c,mxArray *ml[],int nf)
{
	double *mlr;
	ml[0]=mxCreateDoubleMatrix((mwSize) nf,(mwSize) 1,mxREAL);
	mlr=mxGetPr(ml[0]);
	for (int i=0;i<nf;i++)
	{
		mlr[i] = c[i];
	}
}

void mu_cdr2t_c2ml(mymatrix<complex<double> > &c,mxArray *ml[],int col)
{
	int k,i,j;
	double *mlr, *mli;
	mwSize nd=2,d[2];
	d[1]=col;
	d[0]=col;
	ml[0]=mxCreateNumericArray(nd,d,mxDOUBLE_CLASS,mxCOMPLEX);
	mlr=mxGetPr(ml[0]);
	mli=mxGetPi(ml[0]);
	
	for (i=0;i<col;i++)
	{
		for (j=0;j<col;j++)
		{
			mlr[i+j* col] = real(c(i,j));
			mli[i+j* col] = imag(c(i,j));
		}
	}
}

void mu_cdr3t_c2ml(vector< mymatrix<complex<double> > > &c,mxArray *ml[],int nf,int ChLe)
{
	int k,i,j;
	double *mlr, *mli;
	mwSize nd=3,d[3];
	d[2]=nf;
	d[1]=ChLe;
	d[0]=ChLe;
	ml[0]=mxCreateNumericArray(nd,d,mxDOUBLE_CLASS,mxCOMPLEX);
	mlr=mxGetPr(ml[0]);
	mli=mxGetPi(ml[0]);
	for (k=0;k<nf;k++)
	{
		for (i=0;i<ChLe;i++)
		{
			for (j=0;j<ChLe;j++)
			{
				mlr[i+j* ChLe+k*ChLe*ChLe] = real(c[k](j,i));
				mli[i+j* ChLe+k*ChLe*ChLe] = imag(c[k](j,i));
			}
		}
	}
}

//-----------------------------------
//-------------- MAIN ---------------
//-----------------------------------

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[]) 
{

	//-------------- default + declaration ---------------
	int j,k,l,m,N=100;
	double T = 0.00001, U0 = 1, mu = 0, dist = 0.001,hopping = 1.;
	vector<double> e,V,R_alpha,omega_alpha;

	int omega_n;
	if (gam==0) Gamma = 1;
	
	//-------------- errors and handling ---------------
	if (nrhs==0) mexErrMsgTxt("Proper handling: pto2_chain(T,energy discretization,U,mu,V<----vector,R_alpha,omega_alpha)");
	if (nrhs>7) mexErrMsgTxt("you need 7 arguments!");
	if (nrhs<7) mexErrMsgTxt("you need 7 arguments!");
    
	//-------------- read in from matlab ---------------
	T = (double) *mxGetPr(prhs[0]); // read T
    int e_num = (int) *mxGetPr(prhs[1]); // construct energy vector
	e_num = 61;
	e.resize(e_num);
	if (T != 0)
	{
	if (gam==0)
	{
		if (30*T>=2)
		{
			dist = 3.9999/(e_num-1.);
			for (j=0;j<e_num;j++) e[j] = -1.99995+j*dist; 
		}
		else
		{
			dist = 60*T/(e_num-1.);
			for (j=0;j<e_num;j++) e[j] = -30*T+j*dist;
		}
	}
	if (gam==1)
	{
		dist = 60*T/(e_num-1.);
		for (j=0;j<e_num;j++) e[j] = -30*T+j*dist; 
	}
	}
	if (T==0) 
	{
	e.resize(1);
	e[0] = 0;
	e_num = 1;
	}
    U0 = (double) *mxGetPr(prhs[2]); // read U
	mu = (double) *mxGetPr(prhs[3]); // read mu
	mu_dv_ml2c(V,prhs[4],N); // read V and N
	mu_dv_ml2c(R_alpha,prhs[5],omega_n); // read R_alpha
	mu_dv_ml2c(omega_alpha,prhs[6],omega_n); // read omega_alpha
	
	//-------------- calculate non-interacting Greens function G0 at energy-vector e --------------
	vector<double> b(N-1);
	vector<mymatrix<complex<double> > > G0(e_num);
	vector<double> U(N);
	if (V[N/2]!=0)
		for (j=0;j<N;j++) U[j] = U0*V[j]/V[N/2];
	else 
		for (j=0;j<N;j++) U[j] = U0;
	for (j=0;j<N-1;j++)
	{
		b[j] = hopping;
	}
	full_green_td_RealAxis(V,b,0.,G0,N,e,e_num);
	
	//-------------- calculate Self-Energy Sigma and Vertex part P in 2nd order perturbation theory --------------
	vector<mymatrix<complex<double> > > sigma_f(e_num);
	vector<mymatrix<complex<double> > > P(e_num);
	for (j=0;j<e_num;j++) 
	{
		P[j].resize(N,N);
		sigma_f[j].resize(N,N);
	}
	
	vector<double> sigma_h(N,.5);
	hartree_sigma(T,V,b,U,0.,N,sigma_h,R_alpha,omega_alpha,omega_n);
	fock_sigma(G0,T,U,mu,e,e_num,dist,N,sigma_f);
	calc_P(G0,T,U,mu,e,e_num,dist,N,P);
	

	//-------------- calculate interacting Hamiltonian and Greens function --------------
	vector<mymatrix<complex<double> > > G(e_num);
	mymatrix<complex<double> > H(N,N);
    for (j=0;j<e_num;j++) 
	{
		G[j].resize(N,N);
		for (l=0;l<N;l++){
			for (m=0;m<N;m++){
				H(l,m) = 0.;
				H(l,m) = -sigma_f[j](l,m);
			}
			H(l,l)-= sigma_h[l];
		}
		for (k=0;k<N;k++) H(k,k) = H(k,k)-V[k]+e[j];
		for (k=0;k<N-1;k++) H(k,k+1) = H(k,k+1)-b[k];
		for (k=1;k<N;k++) H(k,k-1) = H(k,k-1)-b[k-1];
		if (gam == 1)
		{
			H(0,0) = H(0,0) + I*Gamma;
			H(N-1,N-1) = H(N-1,N-1) +I*Gamma;
		}
		else
		{
			H(0,0) = H(0,0) - greenOFlead_RealAxis(e[j],mu);
			H(N-1,N-1) = H(N-1,N-1) - greenOFlead_RealAxis(e[j],mu);
		}
		m_inv_c(H,G[j],N);
//	 LWORK=-1;
//	    zsytrf_("U", &N,H[0],&N,IPIV,&WOQ,&LWORK,&INFO);
//		LWORK=(int) WOQ;
//	    WORK.rezise(LWORK);
//	    zsytrf_("U", &N,H[0],&N,IPIV,WORK,&LWORK,&INFO);
//		zsytri_("U", &N,H[0],&N,IPIV,WORK,&INFO )

	 }

	//-------------- calculate landauer transmission --------------
	vector<double> Tr_a(e_num);
	vector<double> Tr_0(e_num);
	double A;
	for (j=0;j<e_num;j++) 
	{
		if (gam==1)
		{
		Tr_a[j] = 4*Gamma*Gamma*(real(G[j](0,N-1))*real(G[j](0,N-1))+imag(G[j](0,N-1))*imag(G[j](0,N-1)));
		Tr_0[j] = 4*Gamma*Gamma*(real(G0[j](0,N-1))*real(G0[j](0,N-1))+imag(G0[j](0,N-1))*imag(G0[j](0,N-1)));
		}
		else 
		{
			if (e[j]>2 || e[j]<-2) A=0;
			else A = 4-pow(e[j],2.);
			Tr_a[j] = A*(real(G[j](0,N-1))*real(G[j](0,N-1))+imag(G[j](0,N-1))*imag(G[j](0,N-1)));
			Tr_0[j] = A*(real(G0[j](0,N-1))*real(G0[j](0,N-1))+imag(G0[j](0,N-1))*imag(G0[j](0,N-1)));
		}
	}
	
	//-------------- calculate vertex transmission --------------
	vector<complex<double> > Tr_b(e_num);
	for (j=0;j<e_num;j++) Tr_b[j] = 0;
	for (j=0;j<e_num;j++)
	{
		for (k=0;k<N;k++)
		{
			for (l=0;l<N;l++)
			{
				if (gam == 1)
				{
					Tr_b[j] += 2*Gamma*(real(G[j](0,l))-I*imag(G[j](0,l)))*P[j](l,k)*G[j](k,0);
				}
				else
				{
					Tr_b[j] -= 2*imag(greenOFlead_RealAxis(e[j],mu))*(real(G[j](0,l))-I*imag(G[j](0,l)))*P[j](l,k)*G[j](k,0);
				}
			}
		}
	}
	
	//--------------  derivative of fermi function and cond. integrand --------------
	vector<double> df(e_num);
	vector<double> f(e_num);
	vector<double> integrand_0(e_num);
	vector<double> integrand_a(e_num);
	vector<complex<double> > integrand_b(e_num);
	for (j=0;j<e_num;j++) 
	{
		df[j] = dfermi(T,mu,e[j]);
		f[j] = fermi(T,mu,e[j]);
		integrand_0[j] = df[j]*Tr_0[j];
		integrand_a[j] = df[j]*Tr_a[j];
		integrand_b[j] = df[j]*Tr_b[j];
		
	}
	
	//-------------- calculate conductance --------------
	double C_0;
	double C_a;
	complex<double> C_b;
	C_0 = extSimpson(integrand_0,e,e_num,dist);
	C_a = extSimpson(integrand_a,e,e_num,dist);
	C_b = cextSimpson(integrand_b,e,e_num,dist);
	C_0 = -C_0;
	C_a = -C_a;
	C_b = -C_b;
	vector<double> Conductance_0(1);
	vector<double> Conductance_a(1);
	vector<complex<double> > Conductance_b(1);
	Conductance_0[0] = C_0;
	Conductance_a[0] = C_a;
	Conductance_b[0] = C_b;
	
	//-------------- convert to matlab --------------
	//mu_dv_c2ml(sigma_h,plhs,N);
	//mu_cdv_c2ml(Tr_a,plhs,e_num);
	mu_dv_c2ml(Conductance_a,plhs,1);
	mu_cdv_c2ml(Conductance_b,plhs+1,1);
	//mu_cdv_c2ml(Tr_b,plhs+2,e_num);
	//mu_dv_c2ml(e,plhs+3,e_num);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

complex<double> greenOFlead_RealAxis(double e,double mu)
{
	complex<double> g;
	if (4-pow(e+mu,2)>=0) g = 0.5*(e+mu-I*sqrt(4-pow(e+mu,2)));
	//else g = 0.5*(e+mu+sqrt(pow(e+mu,2)-4));
	else g=0.;
	return g;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

complex<double> greenOFlead_ImagAxis(complex<double> z,double mu)
{
	double r,i;
	if (imag(z)<0.) return 1./2.*(z+mu+I*sqrt((complex<double>) (4.+imag(z)*imag(z))));
	else return 1./2.*(z+mu-I*sqrt((complex<double>) (4.+imag(z)*imag(z))));
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void full_green_td_RealAxis(vector<double> a,vector<double> b, double muh,
				   vector<mymatrix<complex<double> > > &G,int ChLe,vector<double> e, int e_num)
{
	int i;
    for (i=0;i<e_num;i++) G[i].resize(ChLe,ChLe);
	complex<double> g;
	vector<complex<double> > D1, D2, U1, U2;
	int n,rn,m;
	D1.resize(ChLe);
	D2.resize(ChLe);
	U1.resize(ChLe-1);
	U2.resize(ChLe-1);
	for (i=0;i<e_num;i++)
	{
		if (gam ==1)
		{
			g= -I*Gamma;
		}
		else
		{
			g = greenOFlead_RealAxis(e[i],muh);
		}
		D1[0]=e[i] + muh - a[0] - g;
		D2[ChLe-1]=e[i] + muh - a[ChLe-1] - g;
		for (n=0, rn=ChLe-2;n < (ChLe-2);n++, rn--)
		{
			U1[n]=-b[n]/D1[n];
			D1[n+1]=(e[i] + muh - a[n+1])+b[n]*U1[n];
			U2[rn]=-b[rn]/D2[rn+1];
			D2[rn]=(e[i] + muh - a[rn])+b[rn]*U2[rn];
		}
		U1[n]=-b[n]/D1[n];
		D1[n+1]=(e[i] + muh - a[n+1] - g)+b[n]*U1[n];
		U2[rn]=-b[rn]/D2[rn+1];
		D2[rn]=(e[i] + muh - a[rn] - g)+b[rn]*U2[rn];
		G[i](0,0)=1./D2[0];
		for (n=0;n < ChLe-1;n++)
		{
			for (m=n+1;m<ChLe;m++)
				G[i](m,n)=-U2[m-1]*G[i](m-1,n);
			G[i](n+1,n+1)=G[i](n,n)*D1[n]/D2[n+1];
		}
		for (n=0;n <ChLe-1;n++)
		{
			for (m=n+1;m<ChLe;m++)
			{
				G[i](n,m) = G[i](m,n);
			}
		}
	}
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void full_green_td_ImagAxis(vector<double> a,vector<double> b, double muh,
							vector<mymatrix<complex<double> > > &G,int ChLe,vector<complex<double> > e, int e_num)
{
	int i;
    for (i=0;i<e_num;i++) G[i].resize(ChLe,ChLe);
	complex<double> g;
	vector<complex<double> > D1, D2, U1, U2;
	int n,rn,m;
	D1.resize(ChLe);
	D2.resize(ChLe);
	U1.resize(ChLe-1);
	U2.resize(ChLe-1);
	for (i=0;i<e_num;i++)
	{
		if (gam ==1)
		{
			g= -I*Gamma;
		}
		else
		{
			g = Gamma*greenOFlead_ImagAxis(e[i],muh);
		}
		D1[0]=e[i] + muh - a[0] - g;
		D2[ChLe-1]=e[i] + muh - a[ChLe-1] - g;
		for (n=0, rn=ChLe-2;n < (ChLe-2);n++, rn--)
		{
			U1[n]=-b[n]/D1[n];
			D1[n+1]=(e[i] + muh - a[n+1])+b[n]*U1[n];
			U2[rn]=-b[rn]/D2[rn+1];
			D2[rn]=(e[i] + muh - a[rn])+b[rn]*U2[rn];
		}
		U1[n]=-b[n]/D1[n];
		D1[n+1]=(e[i] + muh - a[n+1] - g)+b[n]*U1[n];
		U2[rn]=-b[rn]/D2[rn+1];
		D2[rn]=(e[i] + muh - a[rn] - g)+b[rn]*U2[rn];
		G[i](0,0)=1./D2[0];
		for (n=0;n < ChLe-1;n++)
		{
			for (m=n+1;m<ChLe;m++)
				G[i](m,n)=-U2[m-1]*G[i](m-1,n);
			G[i](n+1,n+1)=G[i](n,n)*D1[n]/D2[n+1];
		}
		for (n=0;n <ChLe-1;n++)
		{
			for (m=n+1;m<ChLe;m++)
			{
				G[i](n,m) = G[i](m,n);
			}
		}
	}
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void fock_sigma(vector<mymatrix<complex<double> > > G,double T,vector<double> U0,double mu,vector<double> e, int e_num,double dist, int N, vector<mymatrix<complex<double> > > &Sigma)
{
	int i,j,k,l,m;
	vector<double> innerIntegrand(e_num);
    vector<double> outerIntegrand(e_num);
	vector<double> f(e_num);
	vector<mymatrix<double> > F(e_num);
	for (i=0;i<e_num;i++) f[i] = fermi(T,mu,e[i]);

	int minusT = (e_num-1)/3;
	int plusT = 2*(e_num-1)/3;
	
    //-------------- imag(Sigma) ---------------
	
	for (m=0;m<e_num;m++) 
	{
		for (i=0;i<N;i++)
		{
			for (j=0;j<=i;j++)
			{
				for (k=0;k<e_num;k++) 
				{
					for (l=0;l<e_num;l++) 
					{
						if (l-k+m<0 || l-k+m>e_num-1) innerIntegrand[l] = 0.;
						else innerIntegrand[l] = imag(G[k](i,j))*imag(G[l](j,i))*imag(G[l-k+m](i,j))*(f[l]*(1.-f[k])*(1.-f[l-k+m])+f[k]*f[l-k+m]*(1.-f[l]));		
					}
					outerIntegrand[k] = extSimpson(innerIntegrand,e,e_num,dist);
				}
				Sigma[m](i,j) = I*U0[i]*U0[j]/pi/pi*extSimpson(outerIntegrand,e,e_num,dist); // assign imaginary part
			}
		}
	}
	
	//-------------- real(Sigma) --------------

	for (m=0;m<e_num;m++) 
	{
		for (i=0;i<N;i++)
		{
			for (j=0;j<=i;j++)
			{
				for (k=0;k<e_num;k++) 
				{
					for (l=0;l<e_num;l++)
					{
						if (l-k+m<0 || l-k+m>e_num-1) innerIntegrand[l] = 0.;
						else innerIntegrand[l] = real(G[k](i,j))*imag(G[l](j,i))*imag(G[l-k+m](i,j))*f[l]*(1.-f[l-k+m])-imag(G[k](i,j))*imag(G[l](j,i))*real(G[l-k+m](i,j))*f[l]*f[k]-imag(G[k](i,j))*real(G[l](j,i))*imag(G[l-k+m](i,j))*f[k]*f[l-k+m];
					}
					outerIntegrand[k] = extSimpson(innerIntegrand,e,e_num,dist);
				}
				Sigma[m](i,j) += U0[i]*U0[j]/pi/pi*extSimpson(outerIntegrand,e,e_num,dist); // add real part
				Sigma[m](j,i) = Sigma[m](i,j); 
			}
		}
	}
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------


void calc_P(vector<mymatrix<complex<double> > > G,double T,vector<double> U0,double mu,vector<double> e, int e_num,double dist, int N, vector<mymatrix<complex<double> > > &P)
{
	int i,j,k,l,m;
	vector<complex<double> > innerIntegrand(e_num);
    vector<complex<double> > outerIntegrand(e_num);
	vector<double> f(e_num);
	for (i=0;i<e_num;i++) f[i] = fermi(T,mu,e[i]);
	for (m=(e_num-1)/3;m<2*(e_num-1)/3+1;m++)
	{
		for (i=0;i<N;i++)
		{
			for (j=0;j<N;j++)
			{
				for (k=0;k<e_num;k++)
				{
					for (l=0;l<e_num;l++)
					{
						if (gam == 1)
						{
						if (l-k+m<0 || l-k+m>e_num-1) innerIntegrand[l] = 0.;
						else innerIntegrand[l] = 2*Gamma*(real(G[k](i,N-1))-I*imag(G[k](i,N-1)))*G[k](N-1,j)*imag(G[l](j,i))*imag(G[l-k+m](i,j))*(f[l]*(1.-f[k])*(1.-f[l-k+m])+f[k]*f[l-k+m]*(1.-f[l]));
						}
						else
						{
							if (l-k+m<0 || l-k+m>e_num-1) innerIntegrand[l] = 0.;
							else innerIntegrand[l] = Gamma*sqrt(4-pow(e[m],2.))*(real(G[k](i,N-1))-I*imag(G[k](i,N-1)))*G[k](N-1,j)*imag(G[l](j,i))*imag(G[l-k+m](i,j))*(f[l]*(1.-f[k])*(1.-f[l-k+m])+f[k]*f[l-k+m]*(1.-f[l]));

						}
					}
					outerIntegrand[k] = cextSimpson(innerIntegrand,e,e_num,dist);
				}
				P[m](i,j) = U0[i]*U0[j]/pi/pi*cextSimpson(outerIntegrand,e,e_num,dist); // assign vertex 
			}
		}
	}
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

double extSimpson(vector<double> f, vector<double> xn, int e_num,double h)
{
	int i;
	double F;
	if (e_num>2)
	{
		if (e_num%2 == 1)
		{
			F = f[0]+f[e_num-1];
			for (i=1;i<e_num-1;i++)
			{
				if (i%2==1) F = F + 4*f[i];
				else F = F + 2*f[i];
			}
			F = F*h/3.;
			return F;
		}
		if (e_num%2 == 0)
		{
			double F_end1;
			double F_end;
			F = f[0]+f[e_num-2];
			for (i=1;i<e_num-2;i++)
			{
				if (i%2==1) F = F + 4*f[i];
				else F = F + 2*f[i];
			}
			F = F*h/3.;
			F = F + 0.5*(f[e_num-1]+f[e_num-2])/h;
			return F;
		}
	}
	if (e_num==2) F = 0.5*(f[0]+f[1])/h;
	if (e_num==1) F = f[0]/h;
}

complex<double> cextSimpson(vector<complex<double> > f, vector<double> xn, int e_num,double h)
{
	complex<double> F;
	if ((int) e_num%2 == 1)
	{
		int i;
		F = f[0]+f[e_num-1];
		for (i=1;i<e_num-1;i++)
		{
			if (i%2==1) F = F + 4.*f[i];
			else F = F + 2.*f[i];
		}
		F = F*h/3.;
		return F;
	}
	if ((int) e_num%2 == 0)
	{
		int i;
		F = f[0]+f[e_num-2];
		for (i=1;i<e_num-3;i++)
		{
			if (i%2==1) F = F + 4.*f[i];
			else F = F + 2.*f[i];
		}
		F = F*h/3.;
		F = F + 0.5*(f[e_num-1]+f[e_num-2])/h;
		return F;
	}
	
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

double fermi(double T, double mu, double e)
{
	return 1./(exp((e-mu)/T)+1);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

double fermi0(double e)
{
    if (e<0) return 1.;
	if (e==0) return 0.5;
	if (e>0) return 0;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

double dfermi(double T, double mu, double e)
{
	return -1/T*(exp((e-mu)/T))/pow((exp((e-mu)/T)+1),2.);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void CgreenDext(vector<double> Sigma,
				vector<complex<double> > &S,vector<complex<double> > &D1, 
				vector<complex<double> > &D2,vector<complex<double> > &U1,
				vector<complex<double> > &U2,vector<complex<double> > &G,int N)
{
	vector<complex<double> > Q1(N), Q2(N);
	int n,rn;
	Q1[0]=G[0]*Sigma[0]*G[0];
	Q2[N-1]=G[N-1]*Sigma[N-1]*G[N-1];
	for (n=0,rn=N-2;n<(N-1);n++,rn--)
	{
		Q1[n+1]=U2[n]*Q1[n]*U2[n]+G[n+1]*Sigma[n+1]*G[n+1];
		Q2[rn]=U1[rn]*Q2[rn+1]*U1[rn]+G[rn]*Sigma[rn]*G[rn];
	}
	for (n=0;n<(N);n++) S[n]=(Q1[n]-G[n]*Sigma[n]*G[n]+Q2[n]);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void greenDiag_ImagAxis(vector<double> a,vector<double> b,complex<double> z,double mu,int N,
						vector<complex<double> > &G,vector<complex<double> > &D1, 
						vector<complex<double> > &D2,vector<complex<double> > &U1,
						vector<complex<double> > &U2)
{
	complex<double> g;
	g=greenOFlead_ImagAxis(z,mu);
	
	D1[0]=z + mu - a[0] - g;
	D2[N-1]=z + mu - a[N-1] - g;
	int n, rn;
	for (n=0, rn=N-2;n < (N-2);n++, rn--)
	{
		U1[n]=-b[n]/D1[n];
		D1[n+1]=(z + mu - a[n+1])+b[n]*U1[n];
		U2[rn]=-b[rn]/D2[rn+1];
		D2[rn]=(z + mu - a[rn])+b[rn]*U2[rn];
	}
	U1[n]=-b[n]/D1[n];
	D1[n+1]=(z + mu - a[n+1] - g)+b[n]*U1[n];
	U2[rn]=-b[rn]/D2[rn+1];
	D2[rn]=(z + mu - a[rn] - g)+b[rn]*U2[rn];
	G[0]=1./D2[0];
	for (n=0;n < N-1;n++)
	{
		G[n+1]=G[n]*D1[n]/D2[n+1];
	}
}
	
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------

void hartree_sigma(double T,vector<double> a,vector<double> b,vector<double> U,double mu, int N, vector<double> &dens,
			 vector<double> R_alpha, vector<double> omega_alpha,int omega_n)
{
	int i,j;
	vector<complex<double> > G(N);
	vector<double> E2(N);
	vector<double> E3(N);
	vector<complex<double> > S(N);
	vector<complex<double> > D1(N);
	vector<complex<double> > D2(N);
	vector<complex<double> > U1(N-1);
	vector<complex<double> > U2(N-1);
	complex<double> z;
	for (i=0;i<omega_n;i++)
	{
		greenDiag_ImagAxis(a,b,I*omega_alpha[i]*T,mu,N,G,D1,D2,U1,U2);
		for (j=0;j<N;j++)
			dens[j] += 2.*T*R_alpha[i]*real(G[j]); 
	}
	for (j=0;j<N;j++)
		dens[j] *= U[j];
	for (i=0;i<omega_n;i++)
	{
		z = I*pi*(2.*i+1.)*T;
		greenDiag_ImagAxis(a,b,I*omega_alpha[i]*T,mu,N,G,D1,D2,U1,U2);
		CgreenDext(dens,S,D1,D2,U1,U2,G,N);
		for (j=0;j<N;j++)
		{
			E2[j] += 2.*T*R_alpha[i]*real(S[j]); 
		}
	}
	for (j=0;j<N;j++)
		dens[j] += U[j]*E2[j];
}
