#ifndef KELDYSH_MFRG_ANDERSON_ACCELERATION_HPP
#define KELDYSH_MFRG_ANDERSON_ACCELERATION_HPP


#include <deque>
#include <Eigen/Dense>

#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/state.hpp"

namespace anderson_impl
{
    Eigen::MatrixXd colwise_difference(Eigen::MatrixXd m);

    /**
     * @brief Solves an ordinary least-squares problem.
     *
     * Return an $\vec x$ that minimizes the square norm of the
     * error $\vec y - A\vec x$.
     *
     * @param y Vector
     * @param A Matrix
     * @return auto Vector
     */
    Eigen::VectorXd ordinary_least_squares(Eigen::VectorXd y, Eigen::MatrixXd A);

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

    Eigen::MatrixXd error_matrix(sigma_size, m);
    for (int i = 0; i < m; ++i)
    {
        error_matrix.col(i) =
                rhs_evals[i].selfenergy.get_selfenergy_vector_incl_hartree() -
                iteration_steps[i].selfenergy.get_selfenergy_vector_incl_hartree();
    }

    Eigen::MatrixXd error_differences = colwise_difference(error_matrix);

    // Parametrization of alpha, to achieve $\sum_i \alpha_i = 1$.
    Eigen::VectorXd gamma = ordinary_least_squares(
            -error_matrix.col(0),
            error_differences);

    alpha(0) = 1 - gamma(0);
    for (int i = 1; i < m - 1; ++i)
    {
    alpha(i) = gamma(i - 1) - gamma(i);
    }
    alpha(m - 1) = gamma(m - 2);
    }

    // const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // logging::log_message_single("Anderson iteration mixing vector: ", alpha.transpose().format(CleanFmt));
    if (mpi_world_rank() == 0){
        utils::print("Anderson acceleration gives: \t [");
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
