#include "contact/include/box_contact_proxy.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

BoxContactProxy::BoxContactProxy()
    : ContactProxy(ContactProxyType::kBox),
    size_(Vector3r::Zero()), translation_(Vector3r::Zero()),
    rotation_(Matrix3r::Zero()) {}

void BoxContactProxy::Initialize(const Options& opt) {
    const std::string error_location = "contact::BoxContactProxy::Initialize";

    for (const std::string& key : { "size", "translation" }) {
        CheckCondition(opt.GetVectorOptionSize(key) == 3, error_location,
            "Incompatible " + key + " size.");
    }
    size_ = opt.vector_option().at("size");
    CheckCondition(size_.minCoeff() > 0, error_location, "The box size must be "
        "strictly positive.");
    translation_ = opt.vector_option().at("translation");
    const MatrixXr R = opt.matrix_option().at("rotation");
    CheckCondition(static_cast<integer>(R.rows()) == 3 &&
        R.rows() == R.cols(), error_location, "The rotation is not square.");

    CheckCondition(IsRotationMatrix(R), error_location,
        "The rotation matrix is not orthogonal.");
    rotation_ = R;
}

}
}