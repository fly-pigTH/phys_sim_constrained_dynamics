#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_POINT_CLOUD_CONTACT_PROXY
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_POINT_CLOUD_CONTACT_PROXY

#include "contact/include/contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

class PointCloudContactProxy : public ContactProxy {
public:
    PointCloudContactProxy();

    void Initialize(const Options& opt) override;

    const integer vertex_num() const { return vertex_num_; }
    const Matrix3Xr& vertices() const { return vertices_; }

private:
    integer vertex_num_;
    Matrix3Xr vertices_;
};

}
}

#endif