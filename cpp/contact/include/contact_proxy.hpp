#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_CONTACT_PROXY
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_CONTACT_PROXY

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

enum class ContactProxyType {
    kPointCloud = 0,
    kSphere,
    kBox,
    kPlane,
    kTotalNum
};

class ContactProxy {
public:
    ContactProxy(const ContactProxyType type) : type_(type) {}
    virtual ~ContactProxy() {}

    virtual void Initialize(const Options& opt) = 0;

    const ContactProxyType type() const { return type_; }

private:
    ContactProxyType type_;
};

using ContactProxyGroup = std::vector<std::shared_ptr<ContactProxy>>;

const std::string ToString(const ContactProxyType type);
const integer ToInt(const ContactProxyType type);

}
}

#endif