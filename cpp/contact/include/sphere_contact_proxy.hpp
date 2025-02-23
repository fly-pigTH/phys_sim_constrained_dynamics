#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_SPHERE_CONTACT_PROXY
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_SPHERE_CONTACT_PROXY

#include "contact/include/contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

class SphereContactProxy : public ContactProxy {
public:
    SphereContactProxy();

    void Initialize(const Options& opt) override;

    const Vector3r& center() const { return center_; }
    const Vector3r& c() const { return center_; }
    const real radius() const { return radius_; }
    const real r() const { return radius_; }

private:
    Vector3r center_;
    real radius_;
};

}
}

#endif