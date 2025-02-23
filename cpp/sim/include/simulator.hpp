#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_SIM_SIMULATOR
#define PHYS_SIM_CONSTRAINED_DYNAMICS_SIM_SIMULATOR

#include "link/include/link.hpp"
#include "joint/include/joint.hpp"
#include "contact/include/contact_proxy.hpp"
#include "contact/include/collision_detector.hpp"

namespace phys_sim_constrained_dynamics {
namespace sim {

struct JointInfo {
public:
    std::shared_ptr<joint::Joint> joint;

    integer offset;
    integer link_num;
    std::vector<integer> link_indices;
};

using JointInfoGroup = std::vector<JointInfo>;

struct ContactProxyPairInfo {
public:
    std::array<integer, 2> contact_proxy_indices;
    real compliance;
    real frictional_coefficient;
};

using ContactProxyPairInfoGroup = std::vector<ContactProxyPairInfo>;

using ContactJacobian = Eigen::Matrix<real, 3, 6>;

struct ContactInfo {
public:
    contact::CollisionEvent event;
    std::array<std::shared_ptr<link::Link>, 2> links;
    std::array<ContactJacobian, 2> jacobians;
    real compliance;
    real frictional_coefficient;
};

using ContactInfoGroup = std::vector<ContactInfo>;

struct StepInfo {
public:
    // For debugging purposes.
    ContactInfoGroup contacts;
};

class Simulator {
public:
    Simulator();
    ~Simulator() {}

    const integer AddLink(const link::LinkType type,
        const Vector3r& translation, const Matrix3r& rotation,
        const Vector3r& linear_velocity, const Vector3r& angular_velocity,
        const Options& opt);

    const integer AddJoint(const joint::JointType type,
        const std::vector<integer>& link_indices,
        const real compliance, const Options& opt);

    const integer AddLinkContactProxy(const integer link_index,
        const contact::ContactProxyType type, const Options& opt);
    const integer AddEnvironmentalContactProxy(
        const contact::ContactProxyType type, const Options& opt);

    const integer AddContactProxyPair(const integer first_contact_proxy_index,
        const integer second_contact_proxy_index, const real compliance,
        const real frictional_coefficient);

    void ClearForce(const integer link_index);
    void ClearForces();
    void ClearTorque(const integer link_index);
    void ClearTorques();

    void ApplyForce(const integer link_index, const Vector3r& attach_point,
        const Vector3r& force);
    void ApplyTorque(const integer link_index, const Vector3r& torque);

    // Reset a link's position and velocity.
    void ResetLinkPositionsAndVelocities(const integer link_index,
        const Vector3r& translation,
        const Matrix3r& rotation,
        const Vector3r& linear_velocity,
        const Vector3r& angular_velocity);

    const StepInfo Step(const real time_step, const Options& opt);

    // Get data member states.
    const integer link_num() const { return link_num_; }
    const link::LinkGroup& links() const { return links_; }
    const std::shared_ptr<link::Link> link(const integer link_index) const;
    const integer joint_num() const { return joint_num_; }
    const JointInfoGroup& joints() const { return joints_; }
    const JointInfo& joint(const integer joint_index) const;

    const integer contact_proxy_num() const { return contact_proxy_num_; }
    const contact::ContactProxyGroup contact_proxies() const {
        return contact_proxies_;
    }
    const std::shared_ptr<contact::ContactProxy> contact_proxy(
        const integer contact_proxy_index) const;

    const integer contact_proxy_pair_num() const {
        return contact_proxy_pair_num_;
    }
    const ContactProxyPairInfoGroup& contact_proxy_pairs() const {
        return contact_proxy_pairs_;
    }
    const ContactProxyPairInfo& contact_proxy_pair(
        const integer contact_proxy_pair_index) const;

private:
    void CheckLinkIndex(const integer link_index,
        const std::string& error_location) const;
    void CheckJointIndex(const integer joint_index,
        const std::string& error_location) const;
    void CheckContactProxyIndex(const integer contact_proxy_index,
        const std::string& error_location) const;
    void CheckContactProxyPairIndex(const integer contact_proxy_pair_index,
        const std::string& error_location) const;

    const std::shared_ptr<link::Link> CreateLink(const link::LinkType type,
        const integer link_index, const Options& opt);
    const std::shared_ptr<joint::Joint> CreateJoint(const joint::JointType type,
        const real compliance, const Options& opt) const;
    const std::shared_ptr<contact::ContactProxy> CreateContactProxy(
        const contact::ContactProxyType type, const Options& opt) const;

    // Dynamics-related functions.
    const SparseMatrixXr ComputeMomentumAndJointConstraintMatrix(
        const real time_step, const std::vector<Matrix3r>& inertia,
        const Options& opt) const;

    const SparseMatrixXr ComputeMomentumMatrix(const real time_step,
        const std::vector<Matrix3r>& inertia, const Options& opt) const;
    const VectorXr ComputeMomentumVector(const real time_step,
        const std::vector<Matrix3r>& inertia, const Options& opt) const;

    const SparseMatrixXr ComputeJointConstraintJacobianMatrix(
        const Options& opt) const;
    const VectorXr ComputeJointConstraintComplianceVector(
        const real time_step, const Options& opt) const;
    const VectorXr ComputeJointConstraintVector(const real time_step,
        const Options& opt) const;

    const ContactInfoGroup DetectCollisions() const;

    const SparseMatrixXr ComputeContactJacobianMatrix(
        const ContactInfoGroup& contacts) const;
    const VectorXr ComputeContactComplianceVector(const real time_step,
        const ContactInfoGroup& contacts) const;
    const VectorXr ComputeContactConstraintVector(const real time_step,
        const ContactInfoGroup& contacts) const;

    void SolveDynamicsWithoutComplementarity(
        const SparseMatrixXr& momentum_and_joint_constraint_matrix,
        const VectorXr& momentum_and_joint_constraint_vector,
        const Options& opt,
        Matrix6Xr& new_velocities, VectorXr& new_constraint_impulses) const;
    // The core LCP solver.
    void SolveDynamicsWithComplementarity(
        const SparseMatrixXr& momentum_and_joint_constraint_matrix,
        const VectorXr& momentum_and_joint_constraint_vector,
        const ContactInfoGroup& contacts,
        const SparseMatrixXr& contact_jacobian_matrix,
        const VectorXr& contact_compliance_vector,
        const VectorXr& contact_constraint_vector,
        const Options& opt,
        Matrix6Xr& new_velocities, VectorXr& new_constraint_impulses) const;

    // See the comments in SolveDynamicsWithContact.
    const Vector3r SolveLocalFrictionalBlcp(
        const Matrix3r& M,
        const Eigen::HouseholderQR<Matrix3r>& M_inverse,
        const Vector3r& b, const real neg_phi_inv_h,
        const real mu) const;

    void StepLinksAndJoints(const real time_step,
        const Matrix6Xr& new_velocities, const Options& opt);

    // Links are rigid bodies.
    integer link_num_;
    link::LinkGroup links_;

    // Joints are equality constraints.
    integer joint_num_;
    JointInfoGroup joints_;

    // Contact proxies.
    integer contact_proxy_num_;
    contact::ContactProxyGroup contact_proxies_;
    // The list of all candidate contact pairs.
    integer contact_proxy_pair_num_;
    ContactProxyPairInfoGroup contact_proxy_pairs_;
};

}
}

#endif