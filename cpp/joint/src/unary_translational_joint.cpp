#include "joint/include/unary_translational_joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {


// Core Function: transfer the attach frame (in relative frame, constant, I guess) to the world-attach frame,
// and calulate the error between this new frame and the anchor frame
const VectorXr UnaryTranslationalJoint::phi(const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryTranslationalJoint::phi");

    Matrix3r world_frame; world_frame.setZero();
    for (integer i = 0; i < 3; ++i)
        world_frame.col(i) = q[0]->ToWorldVector(attach_frame_.col(i));
    return (world_frame - anchor_frame_).reshaped(); // to X x 1 vector (9*1)
}

const std::vector<MatrixX6r> UnaryTranslationalJoint::Jphi(
    const link::LinkGroup& q) const {
    this->CheckLinkGroupSize(q, "UnaryTranslationalJoint::Jphi");

    std::vector<MatrixX6r> J(this->n(), MatrixX6r::Zero(this->m(), 6));    // m = 9
    // 直接计算了partial phi/partial q, 不过对phi进行了拼接, reshaped 方法采用了列优先重拍，和这里拼接相符合
    for (integer i = 0; i < 3; ++i)
        J[0].middleRows<3>(3 * i) = // This is a special usage, i.e. 3*i 往后数三个
            q[0]->ComputeVectorJacobian(attach_frame_.col(i));
    return J;
}

void UnaryTranslationalJoint::InitializeDerived(const Options& opt) {
    const std::string error_location =
        "joint::UnaryTranslationalJoint::InitializeDerived";

    for (const std::string& key : { "attach_frame", "anchor_frame" }) {
        CheckCondition(opt.GetMatrixOptionRows(key) == 3 &&
            opt.GetMatrixOptionCols(key) == 3, error_location,
            "Incompatible matrix size for key " + key + ".");
    }
    attach_frame_ = opt.matrix_option().at("attach_frame");
    anchor_frame_ = opt.matrix_option().at("anchor_frame");

    CheckCondition(IsRotationMatrix(attach_frame_) &&
        IsRotationMatrix(anchor_frame_), error_location,
        "The rotation matrix is not orthogonal.");
}

}
}