#ifndef KELDYSH_MFRG_ANDERSON_ACCELERATION_HPP
#define KELDYSH_MFRG_ANDERSON_ACCELERATION_HPP


#include <deque>
#include <Eigen/Dense>

#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/state.hpp"

namespace anderson_impl
{

    template<typename Q>
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> colwise_difference(Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> m)
    {
        const auto cols = m.cols();
        return m.rightCols(cols - 1) - m.leftCols(cols - 1);
    }


    /**
     * @brief Solves an ordinary least-squares problem.
     *
     * Return an $\vec x$ that minimizes the square norm of the
     * error $\vec y - A\vec x$.
     * assume that A consists out of two parts A = (a, a^*)^T and same for y = (b, b^*)^T.
     * Hence we have A^dagger.A = a^\dagger.a + (a^\dagger.a)^*
     *           and A^dagger.y = a^\dagger.b + (a^\dagger.b)^*
     * in the Keldysh formalism we don't store \Sigma^A => Take real part (see above).
     *
     * @param y Vector
     * @param A Matrix
     * @return auto Vector
     */
     template<typename Q>
    Eigen::MatrixXd ordinary_least_squares(const Eigen::Matrix< Q, Eigen::Dynamic, 1> y, const Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> A)
    {

        Eigen::MatrixXd ATA = (A.transpose().conjugate() * A).real(); // A^dagger A
        return Eigen::Inverse(ATA) * (A.transpose().conjugate()* y).real();
    }

} // namespace anderson_impl

/**
 * @brief Perform an anderson mixing update
 *
 * @tparam num_components
 * @tparam lattice_type
 * @param new_state Newly calculated state.
 * @param iter_states Previous Anderson steps.
 * @param selfenergy_evals Previous selfenergy evaluations.
 * @return state<num_components, lattice_type>
 */
template <typename Q, bool diff>
State<Q, diff> anderson_update(
        const std::deque<State<Q, diff>> &rhs_evals,
        const std::deque<State<Q, diff>> &iteration_steps,
        double zeta = 1.0)
{
    using namespace anderson_impl;
    using myMatrixXQ = Eigen::Matrix< Q, Eigen::Dynamic, Eigen::Dynamic>;
    using myVectorXQ = Eigen::Matrix< Q, Eigen::Dynamic, 1>;

    int m = rhs_evals.size();
    Eigen::VectorXd alpha(m);

    if (m == 0)
    {
        throw std::invalid_argument(
        "Cannot do Anderson update without initial value.");
    }
    else if (m == 1)
    {
        alpha( 0) = 1.0;
    }
    else
    {
        size_t sigma_size = getFlatSize(rhs_evals[0].selfenergy.Sigma.get_dims());

    myMatrixXQ error_matrix(sigma_size, m);
    for (int i = 0; i < m; ++i)
    {
        error_matrix.col(i) =
                rhs_evals[i].selfenergy.get_selfenergy_vector_incl_hartree() -
                iteration_steps[i].selfenergy.get_selfenergy_vector_incl_hartree();
    }


        myMatrixXQ error_differences = colwise_difference(error_matrix);

    // Parametrization of alpha, to achieve $\sum_i \alpha_i = 1$.
    myVectorXQ y = -error_matrix.col(0);
    Eigen::VectorXd gamma = ordinary_least_squares(
            y,
            error_differences);

    alpha(0) = 1. - gamma(0);
    for (int i = 1; i < m - 1; ++i)
    {
    alpha(i) = gamma(i - 1) - gamma(i);
    }
    alpha(m - 1) = gamma(m - 2);
    }

    // const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // logging::log_message_single("Anderson iteration mixing vector: ", alpha.transpose().format(CleanFmt));
    if (mpi_world_rank() == 0){
        utils::print("Anderson acceleration gives the coefficients: \t [");
        for (int j = 0; j < alpha.size(); j++) {std::cout << alpha[j] << ", ";}
        std::cout << "]\n";

    }

    State<Q, diff> result =
            alpha.coeff(0) *
            (zeta * rhs_evals[0] + (1 - zeta) * iteration_steps[0]);
    for (int i = 1; i < m; ++i)
    {
        result += alpha.coeff(i) *
        (zeta * rhs_evals[i] + (1 - zeta) * iteration_steps[i]);
    }
    return result;
}


#endif //KELDYSH_MFRG_ANDERSON_ACCELERATION_HPP
