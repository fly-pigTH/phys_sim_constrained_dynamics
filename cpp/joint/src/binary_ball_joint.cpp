#include "joint/include/binary_ball_joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

const VectorXr BinaryBallJoint::phi(const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "BinaryBallJoint::phi");

    return q[0]->ToWorldPoint(first_attach_point_)
        - q[1]->ToWorldPoint(second_attach_point_);
}

const std::vector<MatrixX6r> BinaryBallJoint::Jphi(
    const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "BinaryBallJoint::Jphi");

    std::vector<MatrixX6r> J(this->n(), MatrixX6r::Zero(this->m(), 6));

    J[0] = q[0]->ComputePointJacobian(first_attach_point_);
    J[1] = -q[1]->ComputePointJacobian(second_attach_point_);

    return J;
}

void BinaryBallJoint::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "joint::BinaryBallJoint::InitializeDerived";

    for (const std::string& key :
        { "first_attach_point", "second_attach_point" }) {
        CheckCondition(opt.GetVectorOptionSize(key) == 3, error_location,
            "Incompatible " + key + " size.");
    }
    first_attach_point_ = opt.vector_option().at("first_attach_point");
    second_attach_point_ = opt.vector_option().at("second_attach_point");
}

}
}