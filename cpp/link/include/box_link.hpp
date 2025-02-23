#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_BOX_LINK
#define PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_BOX_LINK

#include "link/include/link.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

// An axis-aligned box centered at the origin in the local frame.
class BoxLink : public Link {
public:
    BoxLink() : Link(LinkType::kBox), size_(Vector3r::Zero()) {}

    void Initialize(const integer index, const Options& opt) override;

    const Vector3r size() const { return size_; }

private:
    Vector3r size_;
};

const Vector3r ComputeBoxInertia(const real mass, const Vector3r& size);

// Create a [-1, 1]^3 box.
const std::pair<Matrix3Xr, Matrix3Xi> CreateBox();

}
}

#endif