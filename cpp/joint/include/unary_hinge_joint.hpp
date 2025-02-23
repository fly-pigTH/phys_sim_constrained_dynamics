#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_HINGE_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_HINGE_JOINT

#include "joint/include/joint.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

class UnaryHingeJoint : public Joint {
public:
    UnaryHingeJoint() : Joint(JointType::kUnaryHinge, 6, 1),
        attach_point_(Vector3r::Zero()), anchor_point_(Vector3r::Zero()),
        attach_direction_(Vector3r::Zero()),
        anchor_direction_(Vector3r::Zero()) {}
    ~UnaryHingeJoint() {}

    const VectorXr phi(const link::LinkGroup& q) const override;

    const std::vector<MatrixX6r> Jphi(const link::LinkGroup& q) const override;

    const Vector3r& attach_point() const { return attach_point_; }
    const Vector3r& anchor_point() const { return anchor_point_; }
    const Vector3r& attach_direction() const { return attach_direction_; }
    const Vector3r& anchor_direction() const { return anchor_direction_; }

protected:
    void InitializeDerived(const Options& opt) override;

private:
    Vector3r attach_point_;
    Vector3r anchor_point_;

    // The directions below are normalized during initialization.
    Vector3r attach_direction_;
    Vector3r anchor_direction_;
};

}
}

#endif