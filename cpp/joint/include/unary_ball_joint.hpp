#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_BALL_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_UNARY_BALL_JOINT

#include "joint/include/joint.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

class UnaryBallJoint : public Joint {
public:
    UnaryBallJoint() : Joint(JointType::kUnaryBall, 3, 1),
        attach_point_(Vector3r::Zero()), anchor_point_(Vector3r::Zero()) {}
    ~UnaryBallJoint() {}

    const VectorXr phi(const link::LinkGroup& q) const override;

    const std::vector<MatrixX6r> Jphi(const link::LinkGroup& q) const override;

    const Vector3r& attach_point() const { return attach_point_; }
    const Vector3r& anchor_point() const { return anchor_point_; }

protected:
    void InitializeDerived(const Options& opt) override;

private:
    Vector3r attach_point_;
    Vector3r anchor_point_;
};

}
}

#endif