#include "contact/include/contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

const std::string ToString(const ContactProxyType type) {
    switch (type) {
        case ContactProxyType::kPointCloud:
            return "point_cloud";
        case ContactProxyType::kSphere:
            return "sphere";
        case ContactProxyType::kBox:
            return "box";
        case ContactProxyType::kPlane:
            return "plane";
        default:
            return "unsupported type";
    }
}

const integer ToInt(const ContactProxyType type) {
    return static_cast<integer>(type);
}

}
}