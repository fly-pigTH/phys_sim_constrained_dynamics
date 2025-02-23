#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_SPHERE_LINK
#define PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_SPHERE_LINK

#include "link/include/link.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

// A sphere centered at the origin in its local frame.
class SphereLink : public Link {
public:
    SphereLink() : Link(LinkType::kSphere), radius_(0) {}

    void Initialize(const integer index, const Options& opt) override;

    const real radius() const { return radius_; }

private:
    real radius_;
};

const real ComputeSphereVolume(const real radius);

const Vector3r ComputeSphereInertia(const real mass, const real radius);

// Create a sphere with radius = 1.
const std::pair<Matrix3Xr, Matrix3Xi> CreateSphere();

}
}

#endif 