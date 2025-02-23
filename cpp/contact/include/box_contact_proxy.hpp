#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_BOX_CONTACT_PROXY
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_BOX_CONTACT_PROXY

#include "contact/include/contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

class BoxContactProxy : public ContactProxy {
public:
    BoxContactProxy();

    void Initialize(const Options& opt) override;

    const Vector3r& size() const { return size_; }
    const Vector3r& s() const { return size_; }

    const Vector3r& translation() const { return translation_; }
    const Vector3r& t() const { return translation_; }

    const Matrix3r& rotation() const { return rotation_; }
    const Matrix3r& R() const { return rotation_; }

private:
    Vector3r size_;

    // Position and orientation.
    Vector3r translation_;
    Matrix3r rotation_;
};

}
}

#endif