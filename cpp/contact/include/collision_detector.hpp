#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_COLLISION_DETECTOR
#define PHYS_SIM_CONSTRAINED_DYNAMICS_CONTACT_COLLISION_DETECTOR

#include "contact/include/contact_proxy.hpp"
#include "contact/include/point_cloud_contact_proxy.hpp"
#include "contact/include/sphere_contact_proxy.hpp"
#include "contact/include/box_contact_proxy.hpp"
#include "contact/include/plane_contact_proxy.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

struct CollisionEvent {
public:
    CollisionEvent() : local_frame(Matrix3r::Zero()),
        collision_point(Vector3r::Zero()), distance(0) {}

    // The first (dim - 1) column(s) is local tangents.
    // The last column is normal.
    // The normal is always set to point from the second body to the first body.
    // In other words, the normal is the direction of the contact force offered
    // by the second object applied to the first object.
    Matrix3r local_frame;
    Vector3r collision_point;
    // We use <= 0 distance to indicate penetration. When a collision is
    // detected, this number should always be <= 0.
    real distance;
};

using Isometry = std::pair<Matrix3r, Vector3r>;
using CollisionEventGroup = std::vector<CollisionEvent>;

class CollisionDetector {
public:
    CollisionDetector() { InitializeCollisionDetectionTable(); }
    ~CollisionDetector() {}

    const CollisionEventGroup DetectCollisions(
        const std::shared_ptr<ContactProxy>& first_contact_proxy,
        const Isometry& first_transformation,
        const std::shared_ptr<ContactProxy>& second_contact_proxy,
        const Isometry& second_transformation) const;

private:
    typedef const CollisionEventGroup(*DetectCollisionPointer)(
        const std::shared_ptr<ContactProxy>&, const Isometry&,
        const std::shared_ptr<ContactProxy>&, const Isometry&);

    void InitializeCollisionDetectionTable();

    std::array<std::array<DetectCollisionPointer,
        static_cast<integer>(ContactProxyType::kTotalNum)>,
        static_cast<integer>(ContactProxyType::kTotalNum)>
        collision_detection_table_;
};

// Four types of contact proxies: points, sphere, box, and plane.
const CollisionEventGroup DetectPointCloudPlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation);

const CollisionEventGroup DetectSphereSphereCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation);

const CollisionEventGroup DetectSphereBoxCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation);

const CollisionEventGroup DetectSpherePlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation);

const CollisionEventGroup DetectBoxPlaneCollisions(
    const std::shared_ptr<ContactProxy>& first,
    const Isometry& first_transformation,
    const std::shared_ptr<ContactProxy>& second,
    const Isometry& second_transformation);

}
}

#endif