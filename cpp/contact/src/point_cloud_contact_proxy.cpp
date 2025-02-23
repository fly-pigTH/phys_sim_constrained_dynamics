#include "contact/include/point_cloud_contact_proxy.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_constrained_dynamics {
namespace contact {

PointCloudContactProxy::PointCloudContactProxy()
    : ContactProxy(ContactProxyType::kPointCloud), vertex_num_(0),
    vertices_(3, 0) {}

void PointCloudContactProxy::Initialize(const Options& opt) {
    const std::string error_location =
        "contact::PointCloudContactProxy::Initialize";
    const MatrixXr vertices = opt.matrix_option().at("vertices");
    CheckCondition(static_cast<integer>(vertices.rows()) == 3 &&
        static_cast<integer>(vertices.cols()) > 0, error_location,
        "Invalid vertices size.");

    vertex_num_ = static_cast<integer>(vertices.cols());
    vertices_ = vertices;
}

}
}