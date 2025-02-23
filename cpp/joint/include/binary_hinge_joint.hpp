#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_BINARY_HINGE_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_BINARY_HINGE_JOINT

#include "joint/include/joint.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

class BinaryHingeJoint : public Joint {
public:
    BinaryHingeJoint() : Joint(JointType::kBinaryHinge, 6, 2),
        first_attach_point_(Vector3r::Zero()),
        second_attach_point_(Vector3r::Zero()),
        first_attach_direction_(Vector3r::Zero()),
        second_attach_direction_(Vector3r::Zero()) {}
    ~BinaryHingeJoint() {}

    const VectorXr phi(const link::LinkGroup& q) const override;

    const std::vector<MatrixX6r> Jphi(const link::LinkGroup& q) const override;

    const Vector3r& first_attach_point() const {
        return first_attach_point_;
    }
    const Vector3r& second_attach_point() const {
        return second_attach_point_;
    }
    const Vector3r& first_attach_direction() const {
        return first_attach_direction_;
    }
    const Vector3r& second_attach_direction() const {
        return second_attach_direction_;
    }

protected:
    void InitializeDerived(const Options& opt) override;

private:
    Vector3r first_attach_point_;
    Vector3r second_attach_point_;
    // The two directions below are normalized during initialization.
    Vector3r first_attach_direction_;
    Vector3r second_attach_direction_;
};

}
}

#endif