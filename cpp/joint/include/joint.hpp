#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_JOINT
#define PHYS_SIM_CONSTRAINED_DYNAMICS_JOINT_JOINT

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"
#include "link/include/link.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

enum class JointType {
    kUnaryBall = 0,
    kBinaryBall,
    kUnaryHinge,
    kBinaryHinge,
    kUnaryTranslational,
    kTotalNum
};

class Joint {
public:
    // Each constraint describes a constraint function phi on link states:
    //
    // phi(q0, q1, ...) = 0.
    Joint(const JointType type, const integer constraint_num,
        const integer link_num);
    virtual ~Joint() {}

    // Initialize customized data members in the derived class.
    // The argument compliance = inverse stiffness. A small compliance (e.g.,
    // -> 0) means a hard constraint.
    void Initialize(const real compliance, const Options& opt);

    // Compute the constraint function \phi(q0, q1, ...) = 0, where q0, q1, ...,
    // are the DoFs (translations and rotations) of each link. The returned
    // vector's size is constraint_num_.
    const VectorXr ComputeConstraint(const link::LinkGroup& link_group) const {
        return phi(link_group);
    }
    virtual const VectorXr phi(const link::LinkGroup& q) const = 0;

    // Compute Jphi(q0, q1, ...). The returned vector's size is link_num_.
    // Each element in the vector is a matrix of size constraint_num_ x (
    // translation dofs + rotation dofs). Here, translation_dofs = 3 and
    // rotation dofs = 3 in 3D.
    const std::vector<MatrixX6r> ComputeJacobian(
        const link::LinkGroup& link_group) const {
        return Jphi(link_group);
    }
    // NOTE 抽象类的纯虚函数，要求子类给出覆盖实现
    virtual const std::vector<MatrixX6r> Jphi(
        const link::LinkGroup& q) const = 0;

    const JointType type() const { return type_; }
    const integer constraint_num() const { return constraint_num_; }
    const integer m() const { return constraint_num_; }
    const integer link_num() const { return link_num_; }
    const integer n() const { return link_num_; }

    const real compliance() const { return compliance_; }
    const real c() const { return compliance_; }

protected:
    virtual void InitializeDerived(const Options& opt) = 0;

    void CheckConstraintIndex(const integer constraint_index,
        const std::string& error_location) const;
    void CheckLinkGroupSize(const link::LinkGroup& group,
        const std::string& error_location) const;
    void CheckConstraintForceSize(const VectorXr& constraint_force,
        const std::string& error_location) const;

    JointType type_;
    integer constraint_num_;
    integer link_num_;

    // Compliance.
    real compliance_;
};

const std::string ToString(const JointType type);

}
}

#endif