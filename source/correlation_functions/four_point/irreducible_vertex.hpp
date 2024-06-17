#ifndef KELDYSH_MFRG_IRREDUCIBLE_VERTEX_HPP
#define KELDYSH_MFRG_IRREDUCIBLE_VERTEX_HPP

/**
 * The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
 * @tparam Q Type of the data.
 */
template <class Q>
class irreducible{
    friend State<state_datatype,false> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);

    using buffer_type_bare = multidimensional::multiarray<Q,2>;
    buffer_type_bare empty_bare() {
        if (KELDYSH) return buffer_type_bare ({16, n_in});
        else return buffer_type_bare ({1, n_in});
    }

    using freqGrid_type_K3 = bufferFrequencyGrid<k3>;
    using buffer_type_NRG = dataBuffer<Q, k3, K3_config.rank, K3_config.num_freqs, K3_config.position_first_freq_index, freqGrid_type_K3, INTERPOLATION>;

    mutable buffer_type_bare bare;

    bool has_NRG_input = false;
public:
    mutable buffer_type_NRG NRG_rest;
    // shall hold R + \sum_r K_{3, r} - \Gamma_0 -> Do we need two of those, one for each spin component?

    /**
     * Standard constructor for just the bare vertex.
     */
    irreducible() {
        bare = empty_bare();
    };

    /**
     * Used if NRG input for the rest term and K3 is to be used.
     * Needs some parameters to set the frequency grid.
     * @param lambda    regulator
     * @param config    parameters
     */
    void initialize_NRG_input(const double lambda, const fRG_config config){
        bare = empty_bare();
        NRG_rest = buffer_type_NRG(lambda, K3_config.dims, config);
        has_NRG_input = true;
    };


    /**
     * Returns the value of the irreducible vertex. Just constants, if it is the bare vertex.
     * @tparam result_type Type of the result. Typically =Q.
     * @param iK Keldysh index.
     * @param i_in Internal index. Always = 0 for the SIAM.
     * @param spin Spin index.
     * @return Value of the irreducible vertex.
     */
    template<typename result_type=Q> auto val(my_index_t iK, my_index_t i_in, my_index_t spin) const -> result_type;

    /**
     * Read out the irreducible vertex if it contains more than just the bare vertex,
     * e.g. if it holds non-trivial data from NRG. Smoothly interpolates along frequencies.
     * @tparam result_type Type of the result. Typically = Q.
     * @param input VertexInput specifying where to read out the vertex.
     * @return
     */
    template<typename result_type=Q> auto valsmooth(const VertexInput& input) const -> result_type;

    auto acc(int i) const -> Q;
    void direct_set(int i,Q value);

    /**
     * Set the value of the bare interaction to Q.
     * @param iK Keldysh index.
     * @param i_in Internal index
     */
    void setvert(int iK, int i_in, Q);

    /**
     *
     * @param input
     * @param val
     */

    /**
     * Set the value of the dynamical NRG rest term (+K3)
     * @param iK        Keldysh index
     * @param i_spin    Spin index
     * @param iw_t      Bosonic transfer frequency index in the t-channel
     * @param iv_t      Fermionic transfer frequency index in the t-channel
     * @param ivp_t     Fermionic transfer frequency index in the t-channel
     * @param val       Value of the vertex.
     */
    void set_NRG_rest(const int iK, const int i_spin,
                      const int iw_t, const int iv_t, const int ivp_t, const Q val);

    /**
     * Initialize the irreducible vertex.
     * @param val Value of the bare interaction.
     */
    void initialize(Q val);

    buffer_type_bare get_vec() const {return bare;}
    void set_vec(const buffer_type_bare& bare_in) { bare = bare_in;}

    // Various operators for the irreducible vertex
    auto operator+= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare +=vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator+(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator-= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare -=vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator-(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs -= rhs; return lhs;
    }
    auto operator+= (const double& alpha) -> irreducible<Q> {
        this->bare +=alpha;
        return *this;
    }
    friend irreducible<Q> operator+(irreducible<Q> lhs, const double& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator*= (const double& alpha) -> irreducible<Q> {
        this->bare *=alpha;
        return *this;
    }
    friend irreducible<Q> operator*(irreducible<Q> lhs, const double& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator*= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare *= vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator*(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator/= (const irreducible<Q>& vertex) -> irreducible<Q> {
        //his->bare /= vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator/(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        //lhs /= rhs;
        return lhs;
    }
};

/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> template<typename result_type> auto irreducible<Q>::val(const my_index_t iK, const my_index_t i_in, const my_index_t spin) const -> result_type {
    if constexpr(std::is_same_v<result_type,Q>) {
        switch (spin) {
            case 0:
                return bare.at(iK, i_in);
                break;
            case 1:
                return -bare.at(iK, i_in);
                break;
            case 2:
                return 0.;
                break;
            default:
                utils::print("Problems in irred.val. Abort.");
                assert(false);
                return 0.;
        }
    }
    else {
        result_type result;
        constexpr int rows = result_type::RowsAtCompileTime;
        switch (spin) {
            case 0:
                result = bare.template at_vectorized<0,0,rows>(iK,i_in);
                break;
            case 1:
                result = -bare.template at_vectorized<0,0,rows>(iK,i_in);
                break;
            case 2:
                result = myzero<result_type>();
                break;
            default:
                utils::print("Problems in irred.val. Abort.");
                assert(false);
        }
        return result;
    }
}

template <typename Q> template<typename result_type> auto irreducible<Q>::valsmooth(const VertexInput& input) const -> result_type {
    // read out bare vertex
    result_type bare_part = val(input.iK, input.i_in, input.spin);

    if (not has_NRG_input) return bare_part;

    // read out NRG rest term
    result_type dynamical_part = NRG_rest.interpolate(input) ; // TODO: How to specify the channel parametrization?

    return bare_part + dynamical_part;
}


template <typename Q> auto irreducible<Q>::acc(int i) const -> Q {
    assert(i>=0 && i<bare.size());
    return bare.flat_at(i);
}

template <typename Q> void irreducible<Q>::direct_set(int i, Q value) {
    assert(i>=0 && i<bare.size());
    bare.flat_at(i)=value;
}

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value) {
    bare.at(iK, i_in) = value;
}

template <typename Q> void irreducible<Q>::set_NRG_rest(const int iK, const int i_spin,
        const int iw_t, const int iv_t, const int ivp_t, const Q val) {
    assert(has_NRG_input);
    NRG_rest.setvert(val, i_spin, iw_t, iv_t, ivp_t, iK, 0);
}

template <typename Q> void irreducible<Q>::initialize(Q val) {
    if (KELDYSH){
        if (CONTOUR_BASIS != 1) {
            // Keldysh basis:
            for (auto i:odd_Keldysh) {
                for (int i_in=0; i_in<n_in; ++i_in) {
                    this->setvert(i, i_in, val);
                }
            }
        }
        else {
            // Contour basis:
            for (int i_in=0; i_in<n_in; ++i_in) {
                this->setvert( 0, i_in, val); // for forward contour
                this->setvert(15, i_in,-val); // for backward contour
            }
        }
    }
    else{
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setvert(0, i_in, val);
        }
    }
}




#endif //KELDYSH_MFRG_IRREDUCIBLE_VERTEX_HPP
