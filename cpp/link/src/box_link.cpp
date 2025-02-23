#include "link/include/box_link.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace link {

void BoxLink::Initialize(const integer index, const Options& opt) {
    const std::string error_location = "link::BoxLink::Initialize";

    this->set_index(index);

    CheckCondition(opt.GetVectorOptionSize("size") == 3, error_location,
        "Incompatible size.");
    size_ = opt.vector_option().at("size");
    CheckCondition(size_.minCoeff() > 0, error_location, "Non-positive size.");

    const real density = opt.real_option().at("density");
    CheckCondition(density > 0, error_location, "Density should be positive.");
    this->set_mass_and_inv_mass(size_.prod() * density);

    // Inertia of a box:
    this->set_inertia_and_inv_inertia(ComputeBoxInertia(this->mass(), size_));

    const auto box = CreateBox();
    this->set_vertices(box.first.cwiseProduct(
        size_ / 2 * RowVectorXr::Ones(box.first.cols())));
    this->set_elements(box.second);
}

const Vector3r ComputeBoxInertia(const real mass, const Vector3r& size) {
    Vector3r inertia;

    const real x = size.x(), y = size.y(), z = size.z();

    inertia(0) = mass * (y * y + z * z) / 12;
    inertia(1) = mass * (x * x + z * z) / 12;
    inertia(2) = mass * (x * x + y * y) / 12;

    return inertia;
}

const std::pair<Matrix3Xr, Matrix3Xi> CreateBox() {

    Matrix3Xr vertices(3, 8);
    vertices << -1,  -1,  -1,  -1,  1,  1,  1,  1,
                -1,  -1,  1,  1,  -1,  -1,  1,  1,
                -1,  1,  -1,  1,  -1,  1,  -1,  1;
    Matrix3Xi elements(3, 12);
    elements << 0,  0,  0,  0,  4,  4,  2,  2,  1,  1,  0,  0,
                1,  3,  4,  5,  6,  7,  7,  3,  5,  7,  6,  2,
                3,  2,  5,  1,  7,  5,  6,  7,  7,  3,  4,  6;

    return { vertices, elements };
}

}
}