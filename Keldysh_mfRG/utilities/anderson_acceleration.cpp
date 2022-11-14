#include "./anderson_acceleration.hpp"



namespace anderson_impl
{
    Eigen::MatrixXd colwise_difference(Eigen::MatrixXd m)
    {
        const auto cols = m.cols();
        return m.rightCols(cols - 1) - m.leftCols(cols - 1);
    }


    Eigen::VectorXd ordinary_least_squares(Eigen::VectorXd y, Eigen::MatrixXd A)
    {
        Eigen::MatrixXd ATA = A.transpose() * A;
        return Eigen::Inverse(ATA) * A.transpose() * y;
    }
} // namespace anderson_impl
