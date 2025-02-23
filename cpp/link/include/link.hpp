#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_LINK
#define PHYS_SIM_CONSTRAINED_DYNAMICS_LINK_LINK

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

enum class LinkType {
    kBox = 0,
    kSphere,
    kTotalNum
};

class Link {
public:
    Link(const LinkType type);
    virtual ~Link() {}

    // A few possibilities in the derived classes:
    // - Take density and geometric parameters (e.g., box size), and compute
    //   mass, inertia, vertices, and elements;
    // - Take density and a triangle mesh (loaded by a file name), and compute
    //   mass and inertia;
    // - Take mass, inertia, and a triangle mesh.
    virtual void Initialize(const integer index, const Options& opt) = 0;

    // Force and torque.
    void ClearForce() { force_.setZero(); }
    void ClearTorque() { torque_.setZero(); }
    void ApplyForce(const Vector3r& attach_point, const Vector3r& force);
    void ApplyTorque(const Vector3r& torque) { torque_ += torque; }

    void Step(const real time_step, const Options& opt);

    // Getters.
    const LinkType type() const { return type_; }
    const integer index() const { return index_; }
    const real mass() const { return mass_; }
    const real inv_mass() const { return inv_mass_; }
    const Vector3r& inertia() const { return inertia_; }
    const Vector3r& inv_inertia() const { return inv_inertia_; }
    const Matrix3r ToWorldInertia() const;
    const Matrix3r ToInvWorldInertia() const;
    const integer vertex_num() const { return vertex_num_; }
    const Matrix3Xr& vertices() const { return vertices_; }
    const integer element_num() const { return element_num_; }
    const Matrix3Xi& elements() const { return elements_; }

    const Vector3r& translation() const { return translation_; }
    const Vector3r& t() const { return translation_; }
    const Matrix3r& rotation() const { return rotation_; }
    const Matrix3r& R() const { return rotation_; }
    const Vector3r& linear_velocity() const { return linear_velocity_; }
    const Vector3r& v() const { return linear_velocity_; }
    const Vector3r& angular_velocity() const { return angular_velocity_; }
    const Vector3r& omega() const { return angular_velocity_; }
    const Vector3r& force() const { return force_; }
    const Vector3r& f() const { return force_; }
    const Vector3r& torque() const { return torque_; }
    const Vector3r& tau() const { return torque_; }

    // Setters.
    void set_translation(const Vector3r& translation) {
        translation_ = translation;
    }
    void set_t(const Vector3r& t) { translation_ = t; }
    void set_rotation(const Matrix3r& rotation) { set_R(rotation); }
    void set_R(const Matrix3r& R);
    void set_linear_velocity(const Vector3r& linear_velocity) {
        linear_velocity_ = linear_velocity;
    }
    void set_v(const Vector3r& v) { linear_velocity_ = v; }
    void set_angular_velocity(const Vector3r& angular_velocity) {
        angular_velocity_ = angular_velocity;
    }
    void set_omega(const Vector3r& omega) { angular_velocity_ = omega; }
    void set_force(const Vector3r& force) { force_ = force; }
    void set_f(const Vector3r& f) { force_ = f; }
    void set_torque(const Vector3r& torque) { torque_ = torque; }
    void set_tau(const Vector3r& tau) { torque_ = tau; }

    const Vector3r ToWorldPoint(const Vector3r& point) const {
        return rotation_ * point + translation_;
    }
    const Matrix3Xr ToWorldPoints(const Matrix3Xr& points) const {
        return (rotation_ * points).colwise() + translation_;
    }
    const Vector3r ToWorldVector(const Vector3r& vector) const {
        return rotation_ * vector;
    }
    const Vector3r ToBodyPoint(const Vector3r& point) const {
        return rotation_.transpose() * (point - translation_);
    }
    const Vector3r ToBodyVector(const Vector3r& vector) const {
        return rotation_.transpose() * vector;
    }

    // Useful Jacobian functions.
    const Eigen::Matrix<real, 3, 6> ComputePointJacobian(
        const Vector3r& point) const;
    const Eigen::Matrix<real, 3, 6> ComputeVectorJacobian(
        const Vector3r& vector) const;

protected:
    void set_index(const integer index);
    void set_mass_and_inv_mass(const real mass);
    void set_inertia_and_inv_inertia(const Vector3r& inertia);
    void set_vertices(const Matrix3Xr& vertices);
    void set_elements(const Matrix3Xi& elements);

private:
    LinkType type_;
    integer index_;

    // Material parameters.
    real mass_;
    real inv_mass_;
    Vector3r inertia_;
    Vector3r inv_inertia_;

    // Geometry -- We use triangle meshes for both 2D and 3D.
    integer vertex_num_;
    Matrix3Xr vertices_;
    integer element_num_;
    Matrix3Xi elements_;

    // Dynamics.
    Vector3r translation_;
    Matrix3r rotation_;
    Vector3r linear_velocity_;
    Vector3r angular_velocity_;

    Vector3r force_;
    Vector3r torque_;
};

using LinkGroup = std::vector<std::shared_ptr<Link>>;

const std::string ToString(const LinkType type);

}
}

#endif