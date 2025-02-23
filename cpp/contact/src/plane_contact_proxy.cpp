#include "contact/include/plane_contact_proxy.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

PlaneContactProxy::PlaneContactProxy()
    : ContactProxy(ContactProxyType::kPlane),
    local_frame_(Matrix3r::Zero()), offset_(0) {}

void PlaneContactProxy::Initialize(const Options& opt) {
    const std::string error_location =
        "contact::PlaneContactProxy::Initialize";

    CheckCondition(opt.GetVectorOptionSize("normal") == 3, error_location,
        "Incompatible normal size.");
    const Vector3r normal =
        opt.vector_option().at("normal").stableNormalized();
    CheckCondition(IsClose(normal.stableNorm(), 1, ToReal(1e-6), ToReal(1e-3)),
        error_location, "Too small normal.");
    local_frame_ = BuildFrameFromUnitNormal(normal);

    offset_ = opt.real_option().at("offset");
}

}
}