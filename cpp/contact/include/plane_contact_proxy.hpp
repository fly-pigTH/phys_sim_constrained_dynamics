#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_PLANE_CONTACT_PROXY
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_PLANE_CONTACT_PROXY

#include "contact/include/contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

class PlaneContactProxy : public ContactProxy {
public:
    PlaneContactProxy();

    void Initialize(const Options& opt) override;

    const Matrix3r& local_frame() const { return local_frame_; }
    const Matrix3r& R() const { return local_frame_; }
    const Vector3r normal() const { return local_frame_.col(2); }
    const Vector3r n() const { return local_frame_.col(2); }
    const real offset() const { return offset_; }
    const real d() const { return offset_; }

private:
    Matrix3r local_frame_;
    real offset_;
};

}
}

#endif