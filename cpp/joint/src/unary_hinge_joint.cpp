#include "joint/include/unary_hinge_joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {


// Q: phi的表达形式随意？只要等于零就保持了约束？
const VectorXr UnaryHingeJoint::phi(const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryHingeJoint::phi");

    Vector6r point_diff; point_diff.setZero();
    point_diff.head(3) = q[0]->ToWorldPoint(attach_point_) - anchor_point_;
    point_diff.tail(3) = q[0]->ToWorldPoint(attach_point_ + attach_direction_)
        - (anchor_point_ + anchor_direction_);
    return point_diff;
}

const std::vector<MatrixX6r> UnaryHingeJoint::Jphi(
    const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryHingeJoint::Jphi");

    std::vector<MatrixX6r> J(this->n(), MatrixX6r::Zero(this->m(), 6));

    J[0].topRows(3) = q[0]->ComputePointJacobian(attach_point_);
    J[0].bottomRows(3) = q[0]->ComputePointJacobian(
        attach_point_ + attach_direction_);

    return J;
}

void UnaryHingeJoint::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "joint::UnaryHingeJoint::InitializeDerived";

    for (const std::string& key : { "attach_point", "anchor_point",
        "attach_direction", "anchor_direction" }) {
        CheckCondition(opt.GetVectorOptionSize(key) == 3, error_location,
            "Incompatible " + key + " size.");
    }
    attach_point_ = opt.vector_option().at("attach_point");
    anchor_point_ = opt.vector_option().at("anchor_point");

    attach_direction_ =
        opt.vector_option().at("attach_direction").stableNormalized();
    CheckCondition(IsClose(attach_direction_.stableNorm(), 1,
        ToReal(1e-6), ToReal(1e-3)), error_location,
        "Too small attach_direction.");
    anchor_direction_ =
        opt.vector_option().at("anchor_direction").stableNormalized();
    CheckCondition(IsClose(anchor_direction_.stableNorm(), 1,
        ToReal(1e-6), ToReal(1e-3)), error_location,
        "Too small anchor_direction.");
}

}
}