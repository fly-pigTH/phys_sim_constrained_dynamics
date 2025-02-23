#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_BINARY_BALL_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_BINARY_BALL_JOINT

#include "joint/include/joint.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

class BinaryBallJoint : public Joint {
public:
    BinaryBallJoint() : Joint(JointType::kBinaryBall, 3, 2),
        first_attach_point_(Vector3r::Zero()),
        second_attach_point_(Vector3r::Zero()) {}
    ~BinaryBallJoint() {}

    const VectorXr phi(const link::LinkGroup& q) const override;

    const std::vector<MatrixX6r> Jphi(const link::LinkGroup& q) const override;

    const Vector3r& first_attach_point() const { return first_attach_point_; }
    const Vector3r& second_attach_point() const { return second_attach_point_; }

protected:
    void InitializeDerived(const Options& opt) override;

private:
    Vector3r first_attach_point_;
    Vector3r second_attach_point_;
};

}
}

#endif