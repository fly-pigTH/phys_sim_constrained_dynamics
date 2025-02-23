#include "link/include/sphere_link.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

void SphereLink::Initialize(const integer index, const Options& opt) {
    const std::string error_location = "link::SphereLink::Initialize";

    this->set_index(index);
    radius_ = opt.real_option().at("radius");
    CheckCondition(radius_ > 0, error_location, "Radius must be positive.");
    const real density = opt.real_option().at("density");
    CheckCondition(density > 0, error_location, "Density should be positive.");

    this->set_mass_and_inv_mass(ComputeSphereVolume(radius_) * density);

    // Inertia of a sphere:
    this->set_inertia_and_inv_inertia(ComputeSphereInertia(this->mass(),
        radius_));

    const auto sphere = CreateSphere();
    this->set_vertices(sphere.first * radius_);
    this->set_elements(sphere.second);
}

const real ComputeSphereVolume(const real radius) {
    return ToReal(4) * Pi() * radius * radius * radius / ToReal(3);
}

const Vector3r ComputeSphereInertia(const real mass, const real radius) {
    const Vector3r inertia = Vector3r::Constant(
        ToReal(2) * mass * radius * radius / ToReal(5));
    return inertia;
}

const std::pair<Matrix3Xr, Matrix3Xi> CreateSphere() {
    // Load obj mesh file.
    return LoadObjMesh("./asset/mesh/sphere.obj");
}

}
}