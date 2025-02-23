#include "joint/include/unary_ball_joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

const VectorXr UnaryBallJoint::phi(const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryBallJoint::phi");

    return q[0]->ToWorldPoint(attach_point_) - anchor_point_;
}

const std::vector<MatrixX6r> UnaryBallJoint::Jphi(
    const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryBallJoint::Jphi");

    std::vector<MatrixX6r> J(this->n(), MatrixX6r::Zero(this->m(), 6));
    J[0] = q[0]->ComputePointJacobian(attach_point_);

    return J;
}

void UnaryBallJoint::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "joint::UnaryBallJoint::InitializeDerived";

    CheckCondition(opt.GetVectorOptionSize("attach_point") == 3,
        error_location, "Incompatible attach_point size.");
    CheckCondition(opt.GetVectorOptionSize("anchor_point") == 3,
        error_location, "Incompatible anchor_point size.");
    attach_point_ = opt.vector_option().at("attach_point");
    anchor_point_ = opt.vector_option().at("anchor_point");
}

}
}