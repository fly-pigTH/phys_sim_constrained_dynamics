#include "sim/include/simulator.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "contact/include/collision_detector.hpp"

namespace phys_sim_constrained_dynamics {
namespace sim {

const ContactInfoGroup Simulator::DetectCollisions() const {
    const std::string error_location =
        "sim::Simulator::DetectCollision";
    // Loop over all contact proxy pairs and return a list of detected contact
    // events, including their rigid bodies, local frames, and so on.
    ContactInfoGroup info;
    contact::CollisionDetector detector;
    const auto get_contact_proxy_isometry = [&](const integer idx)
        -> const contact::Isometry {
        if (0 <= idx && idx < link_num_) {
            // This is a link contact proxy. Use its current translation and
            // rotation to construct the isometry.
            return { links_[idx]->R(), links_[idx]->t() };
        } else if (link_num_ <= idx && idx < contact_proxy_num_) {
            // This is an environment contact proxy.
            return { Matrix3r::Identity(), Vector3r::Zero() };
        } else {
            // This should never happen, but just to be safe...
            CheckCondition(false, error_location, "Out-of-bound idx.");
            return { Matrix3r::Identity(), Vector3r::Zero() };
        }
    };
    const auto is_link = [&](const integer idx) -> const bool {
        return 0 <= idx && idx < link_num_;
    };
    const ContactJacobian zero = ContactJacobian::Zero();
    for (const ContactProxyPairInfo& pair : contact_proxy_pairs_) {
        const integer first_idx = pair.contact_proxy_indices[0];
        const integer second_idx = pair.contact_proxy_indices[1];
        const contact::CollisionEventGroup events = detector.DetectCollisions(
            contact_proxies_[first_idx],
            get_contact_proxy_isometry(first_idx),
            contact_proxies_[second_idx],
            get_contact_proxy_isometry(second_idx));
        for (const contact::CollisionEvent& e : events) {
            ContactInfo contact_info;
            contact_info.event = e;
            contact_info.links = { nullptr, nullptr };
            contact_info.jacobians = { zero, zero };
            if (is_link(first_idx)) {
                contact_info.links[0] = links_[first_idx];
                contact_info.jacobians[0] = e.local_frame.transpose() *
                    links_[first_idx]->ComputePointJacobian(
                        links_[first_idx]->ToBodyPoint(e.collision_point)
                    );
            }
            if (is_link(second_idx)) {
                contact_info.links[1] = links_[second_idx];
                contact_info.jacobians[1] = e.local_frame.transpose() *
                    links_[second_idx]->ComputePointJacobian(
                        links_[second_idx]->ToBodyPoint(e.collision_point)
                    );
            }
            contact_info.compliance = pair.compliance;
            contact_info.frictional_coefficient = pair.frictional_coefficient;
            info.push_back(contact_info);
        }
    }

    return info;
}

const SparseMatrixXr Simulator::ComputeContactJacobianMatrix(
    const ContactInfoGroup& contacts) const {
    const std::string error_location =
        "sim::Simulator::ComputeContactJacobianMatrix";
    const integer contact_num = static_cast<integer>(contacts.size());

    std::vector<Eigen::Triplet<real>> nonzeros;
    const auto add_block = [&](const integer begin_row_idx,
        const integer begin_col_idx, const ContactJacobian& jacobian) -> void {
        for (integer i = 0; i < 3; ++i)
            for (integer j = 0; j < 6; ++j)
                nonzeros.push_back(Eigen::Triplet<real>(
                    begin_row_idx + i, begin_col_idx + j, jacobian(i, j)));
    };

    for (integer i = 0; i < contact_num; ++i) {
        const ContactInfo& c = contacts[i];
        // Sanity check.
        CheckCondition(!(c.links[0] == nullptr
            && c.links[1] == nullptr), error_location,
            "Contact between two static objects.");
        for (integer j = 0; j < 2; ++j) {
            if (c.links[j] == nullptr) continue;
            // Assemble the Jacobian.
            // Note that the second object needs a minus sign to compute the
            // relative velocity in the local frame.
            // The local frame's normal points from the second object to the
            // first one. The velocity computed from J * v is the first object'
            // velocity w.r.t. the second object in this local frame.
            add_block(i * 3, c.links[j]->index() * 6,
                (j == 0 ? 1 : -1) * c.jacobians[j]);
        }
    }

    return FromTriplet(contact_num * 3, link_num_ * 6, nonzeros);
}

const VectorXr Simulator::ComputeContactComplianceVector(
    const real time_step, const ContactInfoGroup& contacts) const {
    const std::string error_location =
        "sim::Simulator::ComputeContactComplianceVector";

    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real inv_h = ToReal(1) / time_step;

    const integer contact_num = static_cast<integer>(contacts.size());
    VectorXr neg_c_inv_h_sqr = VectorXr::Zero(contact_num * 3);
    for (integer i = 0; i < contact_num; ++i) {
        neg_c_inv_h_sqr(3 * i + 2) = -contacts[i].compliance * inv_h * inv_h;
    }
    return neg_c_inv_h_sqr;
}

const VectorXr Simulator::ComputeContactConstraintVector(
    const real time_step, const ContactInfoGroup& contacts) const {
    const std::string error_location =
        "sim::Simulator::ComputeContactConstraintVector";

    CheckCondition(time_step > 0, error_location, "The time_step argument "
        "must be positive.");
    const real inv_h = ToReal(1) / time_step;

    const integer contact_num = static_cast<integer>(contacts.size());
    VectorXr neg_phi_inv_h = VectorXr::Zero(contact_num * 3);
    for (integer i = 0; i < contact_num; ++i) {
        // Distance is the violation along the normal direction.
        neg_phi_inv_h(3 * i + 2) = -contacts[i].event.distance * inv_h;
    }
    return neg_phi_inv_h;
}

const Vector3r Simulator::SolveLocalFrictionalBlcp(
    const Matrix3r& M, const Eigen::HouseholderQR<Matrix3r>& M_inverse,
    const Vector3r& b, const real neg_phi_inv_h, const real mu) const {

    // v[i] = M[ii] * neg_mu[i] + b.
    // The BLCP condition:
    // - Contact:
    //   * v.z = -phi / h, neg_mu.z > 0.
    //   * v.x > 0 => neg_mu.x = -mu * neg_mu.z.
    //     v.x < 0 => neg_mu.x = mu * neg_mu.z.
    //     v.x = 0 => -mu * neg_mu.z <= neg_mu.x <= mu * neg_mu.z.
    // - No contact:
    //   * neg_mu == 0 && v.z >= -phi / h.
    if (b.z() >= neg_phi_inv_h) {
        // Case 1: no contact.
        return Vector3r::Zero();
    } else {
        // Case 2: with contact. We have nine cases to consider.
        Matrix3r A;
        const Vector3r c = -b + Vector3r(0, 0, neg_phi_inv_h);
        for (const integer dx : { -1, 0, 1 })
            for (const integer dy : { -1, 0, 1 }) {
                A.setZero();
                // Assemble the first column.
                if (dx == 0) A.col(0) = M.col(0);
                else A.col(0) = -Vector3r::UnitX();
                // Assemble the second column.
                if (dy == 0) A.col(1) = M.col(1);
                else A.col(1) = -Vector3r::UnitY();
                // Assemble the last column.
                A.col(2) = M.col(0) * dx * mu + M.col(1) * dy * mu + M.col(2);

                // Solve the linear system of equations.
                const Vector3r unknown = A.householderQr().solve(c);
                const real neg_mu_z = unknown(2);
                // Check the contact impulse.
                if (neg_mu_z < 0) continue;

                // Fetch the impulse.
                const real max_f = mu * neg_mu_z;
                Vector3r neg_mu(0, 0, neg_mu_z);
                neg_mu.x() = (dx == 0) ? unknown.x() : dx * max_f;
                neg_mu.y() = (dy == 0) ? unknown.y() : dy * max_f;
                const Vector3r vel = M * neg_mu + b;
                // Check the LCP condition on x.
                const bool static_x = (dx == 0) && -max_f <= neg_mu.x() &&
                    neg_mu.x() <= max_f;
                const bool dynamic_x = (dx != 0) && dx * vel.x() <= 0;
                if (!(static_x || dynamic_x)) continue;
                // Check the LCP condition on y.
                const bool static_y = (dy == 0) && -max_f <= neg_mu.y() &&
                    neg_mu.y() <= max_f;
                const bool dynamic_y = (dy != 0) && dy * vel.y() <= 0;
                if (!(static_y || dynamic_y)) continue;
                // This is a valid impulse.
                return neg_mu;
            }
        CheckCondition(false,
            "sim::Simulator::SolveLocalFrictionalBlcp",
            "Local BLCP failed. Please consider reducing the time step.");
        return Vector3r::Zero();
    }
}

}
}