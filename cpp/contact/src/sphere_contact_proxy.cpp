#include "contact/include/sphere_contact_proxy.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

SphereContactProxy::SphereContactProxy()
    : ContactProxy(ContactProxyType::kSphere),
    center_(Vector3r::Zero()), radius_(0) {}

void SphereContactProxy::Initialize(const Options& opt) {
    const std::string error_location =
        "contact::SphereContactProxy::Initialize";

    CheckCondition(opt.GetVectorOptionSize("center") == 3, error_location,
        "Incompatible center size.");
    center_ = opt.vector_option().at("center");
    radius_ = opt.real_option().at("radius");
    CheckCondition(radius_ > 0, error_location, "The radius must be positive.");
}

}
}