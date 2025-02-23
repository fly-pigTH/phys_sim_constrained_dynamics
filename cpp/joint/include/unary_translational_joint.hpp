#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_TRANSLATIONAL_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_TRANSLATIONAL_JOINT

#include "joint/include/joint.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

class UnaryTranslationalJoint : public Joint {
public:
    UnaryTranslationalJoint() : Joint(JointType::kUnaryTranslational, 9, 1),
        attach_frame_(Matrix3r::Zero()), anchor_frame_(Matrix3r::Zero()) {}
    ~UnaryTranslationalJoint() {}

    const VectorXr phi(const link::LinkGroup& q) const override;

    const std::vector<MatrixX6r> Jphi(const link::LinkGroup& q) const override;

    const Matrix3r attach_frame() const { return attach_frame_; }
    const Matrix3r anchor_frame() const { return anchor_frame_; }

protected:
    void InitializeDerived(const Options& opt) override;

private:
    Matrix3r attach_frame_;
    Matrix3r anchor_frame_;
};

}
}

#endif