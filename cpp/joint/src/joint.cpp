#include "joint/include/joint.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace joint {

Joint::Joint(const JointType type,
    const integer constraint_num, const integer link_num)
    : type_(type), constraint_num_(constraint_num),
    link_num_(link_num), compliance_(0) {
    CheckCondition(constraint_num > 0 && link_num > 0,
        "joint::Joint::Joint", "Empty constraint or link.");
}

// Initialize customized data members in the derived class.
// The argument compliance = inverse stiffness. A small compliance (e.g.,
// -> 0) means a hard constraint.
void Joint::Initialize(const real compliance, const Options& opt) {
// maybe this is a check
    CheckCondition(compliance >= 0, "joint::Joint::Initialize",
        "Compliance must be >= 0.");
    compliance_ = compliance;

    InitializeDerived(opt);
}

void Joint::CheckConstraintIndex(const integer constraint_index,
    const std::string& error_location) const {
    std::stringstream ss;
    ss << "Incorrect constraint_index (): "
        << constraint_index << " (" << constraint_num_ << ").";
    CheckCondition(0 <= constraint_index && constraint_index < constraint_num_,
        "joint::Joint::" + error_location, ss.str());
}

void Joint::CheckLinkGroupSize(const link::LinkGroup& group,
    const std::string& error_location) const {
    std::stringstream ss;
    ss << "Incorrect group size (): " << group.size()
        << " (" << link_num_ << ").";
    CheckCondition(static_cast<integer>(group.size()) == link_num_,
        "link::" + error_location, ss.str());
}

void Joint::CheckConstraintForceSize(const VectorXr& constraint_force,
    const std::string& error_location) const {
    std::stringstream ss;
    ss << "Incorrect constraint force size (): " << constraint_force.size()
        << " (" << constraint_num_ << ").";
    CheckCondition(static_cast<integer>(constraint_force.size())
        == constraint_num_, "link::" + error_location, ss.str());
}

const std::string ToString(const JointType type) {
    switch (type) {
        case JointType::kUnaryBall:
            return "unary_ball";
        case JointType::kBinaryBall:
            return "binary_ball";
        case JointType::kUnaryHinge:
            return "unary_hinge";
        case JointType::kBinaryHinge:
            return "binary_hinge";
        case JointType::kUnaryTranslational:
            return "unary_translational";
        default:
            return "unsupported type";
    }
}

}
}