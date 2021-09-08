#ifndef FPP_MFRG_INTERPOLATOR1D_H
#define FPP_MFRG_INTERPOLATOR1D_H

/// source:    https://kluge.in-chemnitz.de/opensource/SplineK1/SplineK1.h
/// alternative: https://github.com/igmhub/likely/blob/master/likely/TriCubicInterpolator.h
/*
 * SplineK1.h
 *
 * simple cubic SplineK1 interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../vertex_data.h"
#ifdef HAVE_SSTREAM
#include <sstream>
#include <string>
#endif // HAVE_SSTREAM

// not ideal but disable unused-function warnings
// (we get them because we have implementations in the header file,
// and this is because we want to be able to quickly separate them
// into a cpp file if necessary)

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
//namespace
//{


    // SplineK1 interpolation
    template <class DataContainer, typename Q>
    class SplineK1 : public DataContainer
    {
    public:
    // SplineK1 types
    enum spline_type {
        linear = 10,            // linear interpolation
        cspline = 30,           // cubic splines (classical C^2)
        cspline_hermite = 31    // cubic hermite splines (local, only C^1)
    };

    // boundary condition type for the SplineK1 end-points
    enum bd_type {
        first_deriv = 1,    /// known first derivative
        second_deriv = 2,    /// known second derivative
        third_deriv = 3    /// known third derivative
    };

    protected:
    std::vector<double> m_x = DataContainer::frequencies_K1.b.ts;
    //vec<Q> K1;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i, or else it won't go through grid points
    vec<Q> m_b,m_c,m_d;        // SplineK1 coefficients
    Q m_c0;                            // for left extrapolation
    spline_type m_type = cspline_hermite;         /// set const?
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    bool m_made_monotonic = false;
    void set_coeffs_from_b();               // calculate c_i, d_i from b_i
    /// given in FrequencyGrid already:
    size_t find_closest(double x) const;    // closest idx so that m_x[idx]<=x

    public:
    // default constructor: set boundary condition to be zero curvature
    // at both ends, i.e. natural splines
    /// Do I need this?
    //SplineK1(): m_type(cspline),
    //            m_left(second_deriv), m_right(second_deriv),
    //            m_left_value(0.0), m_right_value(0.0), m_made_monotonic(false)
    //{
    //    ;
    //}
    explicit SplineK1(double Lambda)
        :   DataContainer(Lambda)
    {
        this->initializeK1();
    }


    // modify boundary conditions: if called it must be before initializeK1()
    //void set_boundary(bd_type left, Q left_value,
    //                  bd_type right, Q right_value);

    void initializeK1();
    // set all data points (cubic_spline=false means linear interpolation)
    //void initializeK1(const std::vector<double>& x,
    //                const vec<Q>& y,
    //                spline_type type=cspline_hermite);

    // adjust coefficients so that the SplineK1 becomes piecewise monotonic
    // where possible
    //   this is done by adjusting slopes at grid points by a non-negative
    //   factor and this will break C^2
    //   this can also break boundary conditions if adjustments need to
    //   be made at the boundary points
    // returns false if no adjustments have been made, true otherwise
    bool make_monotonic();

    // evaluates the SplineK1 at point x
    Q interpolK1 (double x) const;
    Q deriv(int order, double x) const;

    // returns the input data points
    //std::vector<double> get_x() const { return m_x; }
    //vec<Q> get_y() const { return DataContainer::K1; }
    //Q get_x_min() const { assert(!m_x.empty()); return m_x.front(); }
    //Q get_x_max() const { assert(!m_x.empty()); return m_x.back(); }

    #ifdef HAVE_SSTREAM
    // SplineK1 info string, i.e. SplineK1 type, boundary conditions etc.
    std::string info() const;
    #endif // HAVE_SSTREAM

    };


    /// Do I need this? --> solves matrix, Too expensive?
    namespace internal
    {

    // band matrix solver
    template <typename Q>
    class band_matrix
    {
    private:
        vec< vec<Q> > m_upper;  // upper band
        vec< vec<Q> > m_lower;  // lower band
    public:
        band_matrix() {};                             // constructor
        band_matrix(int dim, int n_u, int n_l);       // constructor
        ~band_matrix() {};                            // destructor
        void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
        int dim() const;                             // matrix dimension
        int num_upper() const
        {
            return (int)m_upper.size()-1;
        }
        int num_lower() const
        {
            return (int)m_lower.size()-1;
        }
        // access operator
        Q & operator () (int i, int j);            // write
        Q   operator () (int i, int j) const;      // read
        // we can store an additional diagonal (in m_lower)
        Q& saved_diag(int i);
        Q  saved_diag(int i) const;
        void lu_decompose();
        vec<Q> r_solve(const vec<Q>& b) const;
        vec<Q> l_solve(const vec<Q>& b) const;
        vec<Q> lu_solve(const vec<Q>& b, bool is_lu_decomposed=false);

    };

    } // namespace internal




    // ---------------------------------------------------------------------
    // implementation part, which could be separated into a cpp file
    // ---------------------------------------------------------------------

    // SplineK1 implementation
    // -----------------------

    //template <class DataContainer, typename Q>
    //void SplineK1<DataContainer,Q>::set_boundary(SplineK1::bd_type left, Q left_value,
    //                               SplineK1::bd_type right, Q right_value)
    //{
    //assert(m_x.size()==0);          // initializeK1() must not have happened yet
    //m_left=left;
    //m_right=right;
    //m_left_value=left_value;
    //m_right_value=right_value;
    //}

    template <class DataContainer, typename Q>
    void SplineK1<DataContainer,Q>::set_coeffs_from_b()
    {
    assert(m_x.size() == (DataContainer::K1).size());
    assert(m_x.size()==m_b.size());
    assert(m_x.size()>2);
    size_t n=m_b.size();
    if(m_c.size()!=n)
        m_c.resize(n);
    if(m_d.size()!=n)
        m_d.resize(n);

    for(size_t i=0; i<n-1; i++) {       /// i=n-1 not treated (only used for extrapolation to the right)
        const double h  = m_x[i+1]-m_x[i];      /// spacing
        // from continuity and differentiability condition
        m_c[i] = (3.0 * (DataContainer::K1[i + 1] - DataContainer::K1[i]) / h - (2.0 * m_b[i] + m_b[i + 1]) ) / h;   /// checked
        // from differentiability condition
        m_d[i] = ( (m_b[i+1]-m_b[i])/(3.0*h) - 2.0/3.0*m_c[i] ) / h;
    }

    // for left extrapolation coefficients
    m_c0 = (m_left==first_deriv) ? 0.0 : m_c[0];
    }

    template <class DataContainer, typename Q>
    void SplineK1<DataContainer,Q>::initializeK1()
    {
    //assert(x.size()==y.size());
    //assert(x.size()>2);
    //m_type=type;
    //m_made_monotonic=false;
    //m_x=x;
    //    DataContainer::K1=y;
    int n = (int) m_x.size();         /// replace with info on dims?
    // check strict monotonicity of input vector x
    //for(int i=0; i<n-1; i++) {
    //    assert(m_x[i]<m_x[i+1]);
    //}


    if(m_type==linear) {
        // linear interpolation
        m_d.resize(n);
        m_c.resize(n);
        m_b.resize(n);
        for(int i=0; i<n-1; i++) {
            m_d[i]=0.0;
            m_c[i]=0.0;
            m_b[i]= (DataContainer::K1[i + 1] - DataContainer::K1[i]) / (m_x[i + 1] - m_x[i]);
        }
        // ignore boundary conditions, set slope equal to the last segment
        m_b[n-1]=m_b[n-2];
        m_c[n-1]=0.0;
        m_d[n-1]=0.0;
    }
    /*else if(m_type==cspline) {
        // classical cubic splines which are C^2 (twice cont differentiable)
        // this requires solving an equation system

        // setting up the matrix and right hand side of the equation system
        // for the parameters b[]
        internal::band_matrix<Q> A(n,1,1);
        std::vector<double>  rhs(n);
        for(int i=1; i<n-1; i++) {
            A(i,i-1)=1.0/3.0*(m_x[i]-m_x[i-1]);
            A(i,i)=2.0/3.0*(m_x[i+1]-m_x[i-1]);
            A(i,i+1)=1.0/3.0*(m_x[i+1]-m_x[i]);
            rhs[i]=(DataContainer::K1[i+1]-DataContainer::K1[i])/(m_x[i+1]-m_x[i]) - (DataContainer::K1[i]-DataContainer::K1[i-1])/(m_x[i]-m_x[i-1]);
        }
        // boundary conditions
        if(m_left == SplineK1::second_deriv) {
            // 2*c[0] = f''
            A(0,0)=2.0;
            A(0,1)=0.0;
            rhs[0]=m_left_value;
        } else if(m_left == SplineK1::first_deriv) {
            // b[0] = f', needs to be re-expressed in terms of c:
            // (2c[0]+c[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
            A(0,0)=2.0*(m_x[1]-m_x[0]);
            A(0,1)=1.0*(m_x[1]-m_x[0]);
            rhs[0]=3.0*((DataContainer::K1[1]-DataContainer::K1[0])/(m_x[1]-m_x[0])-m_left_value);
        } else {
            assert(false);
        }
        if(m_right == SplineK1::second_deriv) {
            // 2*c[n-1] = f''
            A(n-1,n-1)=2.0;
            A(n-1,n-2)=0.0;
            rhs[n-1]=m_right_value;
        } else if(m_right == SplineK1::first_deriv) {
            // b[n-1] = f', needs to be re-expressed in terms of c:
            // (c[n-2]+2c[n-1])(x[n-1]-x[n-2])
            // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
            A(n-1,n-1)=2.0*(m_x[n-1]-m_x[n-2]);
            A(n-1,n-2)=1.0*(m_x[n-1]-m_x[n-2]);
            rhs[n-1]=3.0*(m_right_value-(DataContainer::K1[n-1]-DataContainer::K1[n-2])/(m_x[n-1]-m_x[n-2]));
        } else {
            assert(false);
        }

        // solve the equation system to obtain the parameters c[]
        m_c=A.lu_solve(rhs);

        // calculate parameters b[] and d[] based on c[]
        m_d.resize(n);
        m_b.resize(n);
        for(int i=0; i<n-1; i++) {
            m_d[i]=1.0/3.0*(m_c[i+1]-m_c[i])/(m_x[i+1]-m_x[i]);
            m_b[i]=(DataContainer::K1[i+1]-DataContainer::K1[i])/(m_x[i+1]-m_x[i])
                   - 1.0/3.0*(2.0*m_c[i]+m_c[i+1])*(m_x[i+1]-m_x[i]);
        }
        // for the right extrapolation coefficients (zero cubic term)
        // f_{n-1}(x) = y_{n-1} + b*(x-x_{n-1}) + c*(x-x_{n-1})^2
        double h=m_x[n-1]-m_x[n-2];
        // m_c[n-1] is determined by the boundary condition
        m_d[n-1]=0.0;
        m_b[n-1]=3.0*m_d[n-2]*h*h+2.0*m_c[n-2]*h+m_b[n-2];   // = f'_{n-2}(x_{n-1})
        if(m_right==first_deriv)
            m_c[n-1]=0.0;   // force linear extrapolation

    }
    */
    else if(m_type==cspline_hermite) {
        // hermite cubic splines which are C^1 (cont. differentiable)
        // and derivatives are specified on each grid point
        // (here we use 3-point finite differences)
        m_b.resize(n);
        m_c.resize(n);
        m_d.resize(n);
        // set b to match 1st order derivative finite difference
        for(int i=1; i<n-1; i++) {
            const double h  = m_x[i+1]-m_x[i];
            const double hl = m_x[i]-m_x[i-1];
            m_b[i] = -h / (hl*(hl+h)) * DataContainer::K1[i - 1] + (h - hl) / (hl * h) * DataContainer::K1[i]
                     + hl / (h*(hl+h)) * DataContainer::K1[i + 1];
        }
        // boundary conditions determine b[0] and b[n-1]
        if(m_left==first_deriv) {
            m_b[0]=m_left_value;
        } else if(m_left==second_deriv) {
            const double h = m_x[1]-m_x[0];
            m_b[0]=0.5*(-m_b[1]-0.5*m_left_value*h+ 3.0 * (DataContainer::K1[1] - DataContainer::K1[0]) / h);  /// checked
        } else if (m_left==third_deriv) {
            const double h = m_x[1]-m_x[0];
            m_b[0]=-m_b[1]+m_left_value/6.*h*h+ 2.0 * (DataContainer::K1[1] - DataContainer::K1[0]) / h;  /// added by me
        } else {
            assert(false);
        }
        if(m_right==first_deriv) {
            m_b[n-1]=m_right_value;
            m_c[n-1]=0.0;
        } else if(m_right==second_deriv) {
            const double h = m_x[n-1]-m_x[n-2];
            m_b[n-1]=0.5*(-m_b[n-2]+0.5*m_right_value*h+ 3.0 * (DataContainer::K1[n - 1] - DataContainer::K1[n - 2]) / h); /// checked
            m_c[n-1]=0.5*m_right_value; /// m_d[n-1] is set to 0. Is this correct/necessary?
        } else if (m_right==third_deriv) {
            const double h = m_x[n-1]-m_x[n-2];
            m_b[n-1]=-m_b[n-2]-m_right_value/6.*h*h+ 2.0 * (DataContainer::K1[n - 1] - DataContainer::K1[n - 2]) / h;  /// added by me
            m_d[n-1]= DataContainer::K1[n - 1] - DataContainer::K1[n - 2] - m_b[n - 2]; /// ???
        } else {
            assert(false);
        }
        m_d[n-1]=0.0;

        // parameters c and d are determined by continuity and differentiability
        set_coeffs_from_b();

    } else {
        assert(false);
    }

    // for left extrapolation coefficients
    m_c0 = (m_left==first_deriv) ? 0.0 : m_c[0];
    }

    template <class DataContainer, typename Q>
    bool SplineK1<DataContainer,Q>::make_monotonic()
    {
    assert(m_x.size() == DataContainer::K1.size());
    assert(m_x.size()==m_b.size());
    assert(m_x.size()>2);
    bool modified = false;
    const int n=(int)m_x.size();
    // make sure: input data monotonic increasing --> b_i>=0
    //            input data monotonic decreasing --> b_i<=0
    for(int i=0; i<n; i++) {
        int im1 = std::max(i-1, 0);
        int ip1 = std::min(i+1, n-1);
        if(((DataContainer::K1[im1] <= DataContainer::K1[i]) && (DataContainer::K1[i] <= DataContainer::K1[ip1]) && m_b[i] < 0.0) ||
           ((DataContainer::K1[im1] >= DataContainer::K1[i]) && (DataContainer::K1[i] >= DataContainer::K1[ip1]) && m_b[i] > 0.0) ) {
            modified=true;
            m_b[i]=0.0;
        }
    }
    // if input data is monotonic (b[i], b[i+1], avg have all the same sign)
    // ensure a sufficient criteria for monotonicity is satisfied:
    //     sqrt(b[i]^2+b[i+1]^2) <= 3 |avg|, with avg=(y[i+1]-y[i])/h,
    for(int i=0; i<n-1; i++) {
        double h = m_x[i+1]-m_x[i];
        double avg = (DataContainer::K1[i + 1] - DataContainer::K1[i]) / h;
        if( avg==0.0 && (m_b[i]!=0.0 || m_b[i+1]!=0.0) ) {
            modified=true;
            m_b[i]=0.0;
            m_b[i+1]=0.0;
        } else if( (m_b[i]>=0.0 && m_b[i+1]>=0.0 && avg>0.0) ||
                   (m_b[i]<=0.0 && m_b[i+1]<=0.0 && avg<0.0) ) {
            // input data is monotonic
            double r = sqrt(m_b[i]*m_b[i]+m_b[i+1]*m_b[i+1])/std::fabs(avg);
            if(r>3.0) {
                // sufficient criteria for monotonicity: r<=3
                // adjust b[i] and b[i+1]
                modified=true;
                m_b[i]   *= (3.0/r);
                m_b[i+1] *= (3.0/r);
            }
        }
    }

    if(modified==true) {
        set_coeffs_from_b();
        m_made_monotonic=true;
    }

    return modified;
    }

    // return the closest idx so that m_x[idx] <= x (return 0 if x<m_x[0])
    template <class DataContainer, typename Q>
    size_t SplineK1<DataContainer,Q>::find_closest(double x) const
    {
    //std::vector<double>::const_iterator it;
    //it=std::upper_bound(m_x.begin(),m_x.end(),x);       // *it > x
    //size_t idx = std::max( int(it-m_x.begin())-1, 0);   // m_x[idx] <= x
    //return idx;
       return DataContainer::frequencies_K1.b.fconv(x);
    }

    template <class DataContainer, typename Q>
    Q SplineK1<DataContainer,Q>::interpolK1 (double x) const
    {
    // polynomial evaluation using Horner's scheme
    // TODO: consider more numerically accurate algorithms, e.g.:
    //   - Clenshaw
    //   - Even-Odd method by A.C.R. Newbery
    //   - Compensated Horner Scheme
    size_t n=m_x.size();
    size_t idx=find_closest(x);

    double h=DataContainer::frequencies_K1.b.grid_transf(x) - m_x[idx];
    Q interpol;
    /// don't support extrapolation
    //if(x<m_x[0]) {
    //    // extrapolation to the left
    //    interpol=(m_c0*h + m_b[0])*h + DataContainer::K1[0];
    //} else if(x>m_x[n-1]) {
    //    // extrapolation to the right
    //    interpol=(m_c[n-1]*h + m_b[n-1])*h + DataContainer::K1[n - 1];
    //} else {
        // interpolation
        interpol=((m_d[idx]*h + m_c[idx])*h + m_b[idx])*h + DataContainer::K1[idx];
    //}
    assert(isfinite(interpol));
    return interpol;
    }

    template <class DataContainer, typename Q>
    Q SplineK1<DataContainer,Q>::deriv(int order, double x) const
    {
    assert(order>0);
    size_t n=m_x.size();
    size_t idx = find_closest(x);

    double h=x-m_x[idx];
    Q interpol;
    /// We don't need extrapolation
    if(x<m_x[0]) {
        // extrapolation to the left
        switch(order) {
            case 1:
                interpol=2.0*m_c0*h + m_b[0];
                break;
            case 2:
                interpol=2.0*m_c0;
                break;
            default:
                interpol=0.0;
                break;
        }
    } else if(x>m_x[n-1]) {
        // extrapolation to the right
        switch(order) {
            case 1:
                interpol=2.0*m_c[n-1]*h + m_b[n-1];
                break;
            case 2:
                interpol=2.0*m_c[n-1];
                break;
            default:
                interpol=0.0;
                break;
        }
    } else {
        // interpolation
        switch(order) {
            case 1:
                interpol=(3.0*m_d[idx]*h + 2.0*m_c[idx])*h + m_b[idx];
                break;
            case 2:
                interpol=6.0*m_d[idx]*h + 2.0*m_c[idx];
                break;
            case 3:
                interpol=6.0*m_d[idx];
                break;
            default:
                interpol=0.0;
                break;
        }
    }
    return interpol;
    }

    #ifdef HAVE_SSTREAM
    std::string SplineK1::info() const
    {
    std::stringstream ss;
    ss << "type " << m_type << ", left boundary deriv " << m_left << " = ";
    ss << m_left_value << ", right boundary deriv " << m_right << " = ";
    ss << m_right_value << std::endl;
    if(m_made_monotonic) {
    ss << "(SplineK1 has been adjusted for piece-wise monotonicity)";
    }
    return ss.str();
    }
    #endif // HAVE_SSTREAM


    namespace internal
    {

        // band_matrix implementation
        // -------------------------
        template <typename Q>
        band_matrix<Q>::band_matrix(int dim, int n_u, int n_l)
        {
            resize(dim, n_u, n_l);
        }
            // -------------------------
        template <typename Q>
        void band_matrix<Q>::resize(int dim, int n_u, int n_l)
        {
            assert(dim>0);
            assert(n_u>=0);
            assert(n_l>=0);
            m_upper.resize(n_u+1);
            m_lower.resize(n_l+1);
            for(size_t i=0; i<m_upper.size(); i++) {
                m_upper[i].resize(dim);
            }
            for(size_t i=0; i<m_lower.size(); i++) {
                m_lower[i].resize(dim);
            }
        }
        template <typename Q>
        int band_matrix<Q>::dim() const
        {
            if(m_upper.size()>0) {
                return m_upper[0].size();
            } else {
                return 0;
            }
        }


        // defines the new operator (), so that we can access the elements
        // by A(i,j), index going from i=0,...,dim()-1
        template <typename Q>
        Q & band_matrix<Q>::operator () (int i, int j)
        {
            int k=j-i;       // what band is the entry
            assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
            assert( (-num_lower()<=k) && (k<=num_upper()) );
            // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
            if(k>=0)    return m_upper[k][i];
            else        return m_lower[-k][i];
        }
        template <typename Q>
        Q band_matrix<Q>::operator () (int i, int j) const
        {
            int k=j-i;       // what band is the entry
            assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
            assert( (-num_lower()<=k) && (k<=num_upper()) );
            // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
            if(k>=0)    return m_upper[k][i];
            else        return m_lower[-k][i];
        }
        // second diag (used in LU decomposition), saved in m_lower
        template <typename Q>
        Q band_matrix<Q>::saved_diag(int i) const
        {
            assert( (i>=0) && (i<dim()) );
            return m_lower[0][i];
        }
        template <typename Q>
        Q & band_matrix<Q>::saved_diag(int i)
        {
            assert( (i>=0) && (i<dim()) );
            return m_lower[0][i];
        }

        // LR-Decomposition of a band matrix
        template <typename Q>
        void band_matrix<Q>::lu_decompose()
        {
            int  i_max,j_max;
            int  j_min;
            Q x;

            // preconditioning
            // normalize column i so that a_ii=1
            for(int i=0; i<this->dim(); i++) {
                assert(this->operator()(i,i)!=0.0);
                this->saved_diag(i)=1.0/this->operator()(i,i);
                j_min=std::max(0,i-this->num_lower());
                j_max=std::min(this->dim()-1,i+this->num_upper());
                for(int j=j_min; j<=j_max; j++) {
                    this->operator()(i,j) *= this->saved_diag(i);
                }
                this->operator()(i,i)=1.0;          // prevents rounding errors
            }

            // Gauss LR-Decomposition
            for(int k=0; k<this->dim(); k++) {
                i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
                for(int i=k+1; i<=i_max; i++) {
                    assert(this->operator()(k,k)!=0.0);
                    x=-this->operator()(i,k)/this->operator()(k,k);
                    this->operator()(i,k)=-x;                         // assembly part of L
                    j_max=std::min(this->dim()-1,k+this->num_upper());
                    for(int j=k+1; j<=j_max; j++) {
                        // assembly part of R
                        this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
                    }
                }
            }
        }
        // solves Ly=b
        template <typename Q>
        vec<Q> band_matrix<Q>::l_solve(const vec<Q>& b) const
        {
            assert( this->dim()==(int)b.size() );
            vec<Q> x(this->dim());
            int j_start;
            Q sum;
            for(int i=0; i<this->dim(); i++) {
                sum=0;
                j_start=std::max(0,i-this->num_lower());
                for(int j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
                x[i]=(b[i]*this->saved_diag(i)) - sum;
            }
            return x;
        }
        // solves Rx=y
        template <typename Q>
        vec<Q> band_matrix<Q>::r_solve(const vec<Q>& b) const
        {
            assert( this->dim()==(int)b.size() );
            vec<Q> x(this->dim());
            int j_stop;
            Q sum;
            for(int i=this->dim()-1; i>=0; i--) {
                sum=0;
                j_stop=std::min(this->dim()-1,i+this->num_upper());
                for(int j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
                x[i]=( b[i] - sum ) / this->operator()(i,i);
            }
            return x;
        }
        template <typename Q>
        vec<Q> band_matrix<Q>::lu_solve(const vec<Q>& b,
                                                  bool is_lu_decomposed)
        {
            assert( this->dim()==(int)b.size() );
            vec<Q>  x,y;
            if(is_lu_decomposed==false) {
                this->lu_decompose();
            }
            y=this->l_solve(b);
            x=this->r_solve(y);
            return x;
        }

    } // namespace internal



//} // namespace





#endif //FPP_MFRG_INTERPOLATOR1D_H
