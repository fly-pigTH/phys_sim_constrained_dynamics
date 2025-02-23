#include "sim/include/simulator.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "link/include/box_link.hpp"
#include "link/include/sphere_link.hpp"
#include "joint/include/unary_ball_joint.hpp"
#include "joint/include/binary_ball_joint.hpp"
#include "joint/include/unary_hinge_joint.hpp"
#include "joint/include/binary_hinge_joint.hpp"
#include "joint/include/unary_translational_joint.hpp"
#include "contact/include/point_cloud_contact_proxy.hpp"
#include "contact/include/sphere_contact_proxy.hpp"
#include "contact/include/box_contact_proxy.hpp"
#include "contact/include/plane_contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace sim {

Simulator::Simulator()
    : link_num_(0), links_(0), joint_num_(0), joints_(0),
    contact_proxy_num_(0), contact_proxies_(0),
    contact_proxy_pair_num_(0), contact_proxy_pairs_(0) {}

const integer Simulator::AddLink(const link::LinkType type,
    const Vector3r& translation, const Matrix3r& rotation,
    const Vector3r& linear_velocity, const Vector3r& angular_velocity,
    const Options& opt) {

    const std::string error_location = "sim::Simulator::AddLink";

    const integer link_idx = link_num_;
    ++link_num_;
    std::shared_ptr<link::Link> link = CreateLink(type, link_idx, opt);

    link->set_t(translation);
    link->set_R(rotation);
    link->set_v(linear_velocity);
    link->set_omega(angular_velocity);

    links_.push_back(link);
    // Allocate places for contact proxies.
    ++contact_proxy_num_;
    contact_proxies_.push_back(nullptr);
    CheckCondition(link_num_ == contact_proxy_num_, error_location,
        "Do not call AddEnvironmentalContactProxy before finishing "
        "all Addlink calls.");

    return link_idx;
}

const integer Simulator::AddJoint(const joint::JointType type,
    const std::vector<integer>& link_indices,
    const real compliance, const Options& opt) {

    const std::string error_location = "sim::Simulator::AddJoint";
    CheckCondition(compliance > 0, error_location, "Non-positive compliance. "
        "In particular, please refrain from using zero compliance and use a "
        "small epsilon instead to improve numerical robustness.");

    const std::shared_ptr<joint::Joint> joint = CreateJoint(type,
        compliance, opt);
    JointInfo joint_info;
    joint_info.joint = joint;
    // Set up offset.
    if (joint_num_ == 0) {
        // The very first joint.
        joint_info.offset = 0;
    } else {
        joint_info.offset = joints_.back().offset + joints_.back().joint->m();
    }
    // Sanity check link indices.
    CheckCondition(static_cast<integer>(link_indices.size()) == joint->n(),
        error_location, "The size of link_indices mismatches the joint.");
    for (const integer i : link_indices)
        CheckLinkIndex(i, "AddJoint");
    joint_info.link_num = static_cast<integer>(link_indices.size());
    joint_info.link_indices = link_indices;

    const integer joint_idx = joint_num_;
    joints_.push_back(joint_info);
    ++joint_num_;

    return joint_idx;
}

const integer Simulator::AddLinkContactProxy(const integer link_index,
    const contact::ContactProxyType type, const Options& opt) {
    CheckLinkIndex(link_index, "AddLinkContactProxy");

    contact_proxies_[link_index] = CreateContactProxy(type, opt);
    return link_index;
}

const integer Simulator::AddEnvironmentalContactProxy(
    const contact::ContactProxyType type, const Options& opt) {

    const integer idx = contact_proxy_num_;
    contact_proxies_.push_back(CreateContactProxy(type, opt));
    ++contact_proxy_num_;

    return idx;
}

const integer Simulator::AddContactProxyPair(
    const integer first_contact_proxy_index,
    const integer second_contact_proxy_index,
    const real compliance, const real frictional_coefficient) {

    const std::string error_location =
        "sim::Simulator::AddContactProxyPair";

    ContactProxyPairInfo info;
    CheckContactProxyIndex(first_contact_proxy_index, "AddContactProxyPair");
    CheckContactProxyIndex(second_contact_proxy_index, "AddContactProxyPair");
    info.contact_proxy_indices = {
        first_contact_proxy_index, second_contact_proxy_index };

    CheckCondition(compliance > 0, error_location, "Non-positive compliance. "
        "In particular, please refrain from using zero compliance and use a "
        "small epsilon instead to improve numerical robustness.");
    info.compliance = compliance;

    CheckCondition(0 <= frictional_coefficient && frictional_coefficient <= 1,
        error_location, "Frictional coefficient must be in [0, 1].");
    info.frictional_coefficient = frictional_coefficient;

    const integer idx = contact_proxy_pair_num_;
    contact_proxy_pairs_.push_back(info);
    ++contact_proxy_pair_num_;

    return idx;
}

void Simulator::ResetLinkPositionsAndVelocities(const integer link_index,
    const Vector3r& translation, const Matrix3r& rotation,
    const Vector3r& linear_velocity, const Vector3r& angular_velocity) {
    CheckLinkIndex(link_index, "ResetLinkPositionsAndVelocities");

    links_[link_index]->set_t(translation);
    links_[link_index]->set_R(rotation);
    links_[link_index]->set_v(linear_velocity);
    links_[link_index]->set_omega(angular_velocity);
}

// Get data member states.
const std::shared_ptr<link::Link> Simulator::link(
    const integer link_index) const {
    CheckLinkIndex(link_index, "link");
    return links_[link_index];
}

const JointInfo& Simulator::joint(const integer joint_index) const {
    CheckJointIndex(joint_index, "joint");
    return joints_[joint_index];
}

const std::shared_ptr<contact::ContactProxy> Simulator::contact_proxy(
    const integer contact_proxy_index) const {
    CheckContactProxyIndex(contact_proxy_index, "contact_proxy");
    return contact_proxies_[contact_proxy_index];
}

const ContactProxyPairInfo& Simulator::contact_proxy_pair(
    const integer contact_proxy_pair_index) const {
    CheckContactProxyPairIndex(contact_proxy_pair_index, "contact_proxy_pair");
    return contact_proxy_pairs_[contact_proxy_pair_index];
}

void Simulator::CheckLinkIndex(const integer link_index,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Out-of-bound link_index (link_num_): " << link_index << " (" <<
        link_num_ << ").";
    CheckCondition(0 <= link_index && link_index < link_num_,
        "sim::Simulator::" + error_location, ss.str());
}

void Simulator::CheckJointIndex(const integer joint_index,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Out-of-bound joint_index (joint_num_): "
        << joint_index << " (" << joint_num_ << ").";
    CheckCondition(0 <= joint_index && joint_index < joint_num_,
        "sim::Simulator::" + error_location, ss.str());
}

void Simulator::CheckContactProxyIndex(const integer contact_proxy_index,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Out-of-bound contact_proxy_index (contact_proxy_num): "
        << contact_proxy_index << " (" << contact_proxy_num_ << ").";
    CheckCondition(0 <= contact_proxy_index
        && contact_proxy_index < contact_proxy_num_,
        "sim::Simulator::" + error_location, ss.str());   
}

void Simulator::CheckContactProxyPairIndex(
    const integer contact_proxy_pair_index,
    const std::string& error_location) const {

    std::stringstream ss;
    ss << "Out-of-bound contact_proxy_pair_index (contact_proxy_pair_num): "
        << contact_proxy_pair_index << " (" << contact_proxy_pair_num_ << ").";
    CheckCondition(0 <= contact_proxy_pair_index
        && contact_proxy_pair_index < contact_proxy_pair_num_,
        "sim::Simulator::" + error_location, ss.str());
}

const std::shared_ptr<link::Link> Simulator::CreateLink(
    const link::LinkType type, const integer link_index, const Options& opt) {
    const std::string error_location = "sim::Simulator::CreateLink";

    CheckLinkIndex(link_index, "CreateLink");

    std::shared_ptr<link::Link> link(nullptr);
    switch (type) {
        case link::LinkType::kBox: {
            link = std::make_shared<link::BoxLink>();
            break;
        }
        case link::LinkType::kSphere: {
            link = std::make_shared<link::SphereLink>();
            break;
        }
        default: {
            CheckCondition(false, error_location, "Unsupported link type: " +
                ToString(type) + ".");
            break;
        }
    }
    link->Initialize(link_index, opt);

    return link;
}

const std::shared_ptr<joint::Joint> Simulator::CreateJoint(
    const joint::JointType type, const real compliance,
    const Options& opt) const {

    const std::string error_location = "sim::Simulator::CreateJoint";

    std::shared_ptr<joint::Joint> joint(nullptr);
    switch (type) {
        case joint::JointType::kUnaryBall: {
            joint = std::make_shared<joint::UnaryBallJoint>();
            break;
        }
        case joint::JointType::kBinaryBall: {
            joint = std::make_shared<joint::BinaryBallJoint>();
            break;
        }
        case joint::JointType::kUnaryHinge: {
            joint = std::make_shared<joint::UnaryHingeJoint>();
            break;
        }
        case joint::JointType::kBinaryHinge: {
            joint = std::make_shared<joint::BinaryHingeJoint>();
            break;
        }
        case joint::JointType::kUnaryTranslational: {
            joint = std::make_shared<joint::UnaryTranslationalJoint>();
            break;
        }
        default: {
            CheckCondition(false, error_location, "Unsupported joint type: " +
                ToString(type) + ".");
            break;
        }
    }

    joint->Initialize(compliance, opt);

    return joint;
}

const std::shared_ptr<contact::ContactProxy> Simulator::CreateContactProxy(
    const contact::ContactProxyType type, const Options& opt) const {

    const std::string error_location = "sim::Simulator::CreateContactProxy";

    std::shared_ptr<contact::ContactProxy> contact_proxy(nullptr);
    switch (type) {
        case contact::ContactProxyType::kPointCloud: {
            contact_proxy = std::make_shared<contact::PointCloudContactProxy>();
            break;
        }
        case contact::ContactProxyType::kSphere: {
            contact_proxy = std::make_shared<contact::SphereContactProxy>();
            break;
        }
        case contact::ContactProxyType::kBox: {
            contact_proxy = std::make_shared<contact::BoxContactProxy>();
            break;
        }
        case contact::ContactProxyType::kPlane: {
            contact_proxy = std::make_shared<contact::PlaneContactProxy>();
            break;
        }
        default: {
            CheckCondition(false, error_location,
                "Unsupported contact proxy type: " + ToString(type) + ".");
            break;
        }
    }

    contact_proxy->Initialize(opt);

    return contact_proxy;
}

}
}