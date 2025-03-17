#include "sim/include/simulator.hpp"
#include "basic/include/log.hpp"
#include "basic/include/math.hpp"
#include "joint/include/binary_hinge_joint.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

namespace phys_sim_constrained_dynamics {

// A simple PID controller.
class PidController {
public:
    PidController() : p_coeff_(0), i_coeff_(0), d_coeff_(0),
        error_integral_(0) {}

    void Initialize(const real p_coeff, const real i_coeff,
        const real d_coeff) {
        p_coeff_ = p_coeff;
        i_coeff_ = i_coeff;
        d_coeff_ = d_coeff;
        error_integral_ = 0;
    }

    void Reset() { error_integral_ = 0; }

    const real ComputeControlSignal(const real error,
        const real error_derivative) {
        error_integral_ += error;
        return p_coeff_ * error + d_coeff_ * error_derivative
            + i_coeff_ * error_integral_;
    }

private:
    real p_coeff_;
    real i_coeff_;
    real d_coeff_;

    real error_integral_;
};

void TestGripper() {
    PrintInfo("phys_sim_constrained_dynamics::sim::TestGripper",
        "TestGripper begins...");

    const real link_length = 0.2;
    const real link_width = 0.05;
    const real link_density = 1e3;
    const integer finger_num = 6;
    const integer finger_link_num = 3;
    const real initial_gripper_height = 0.6;

    const real sphere_radius = 0.2;
    const real sphere_density = 1e2;

    const Vector3r gravitational_acceleration = Vector3r(0, 0, -9.81);
    const real joint_compliance = 1e-6;
    const real contact_compliance = 1e-6;
    const real frictional_coefficient = 0.1;
    // An impulse tolerance controls mass * delta v.
    const real blcp_ip_tol = link_density * link_length * link_width
        * link_width * ToReal(0.01);

    integer substep_num = 20;
    integer max_blcp_iter = 50;
    float time_step = 0.017f;

    // A simple PID position controller.
    PidController x_ctrl; x_ctrl.Initialize(500.0, 0.0, 100.0);
    PidController y_ctrl; y_ctrl.Initialize(500.0, 0.0, 100.0);
    PidController z_ctrl; z_ctrl.Initialize(500.0, 1.0, 100.0);
    const Vector3r init_target_pos(0, 0,
        initial_gripper_height + link_length / 2);
    Vector3r target_pos = init_target_pos;
    // A simple hinge angle controller.
    std::vector<PidController> hinge_angle_ctrl(finger_link_num * finger_num);
    for (auto& c : hinge_angle_ctrl) c.Initialize(10.0, 1e-3, 1.0);
    real target_angle = 0;
    std::vector<real> last_angle(finger_link_num * finger_num, 0);

    // Initialize the simulator.
    sim::Simulator sim;
    // Add base.
    Options link_opt;
    link_opt.vector_option()["size"] =
        Vector3r(link_length, link_width, link_width);
    link_opt.real_option()["density"] = link_density;
    sim.AddLink(link::LinkType::kBox,
        Vector3r(0, 0, initial_gripper_height + link_length / 2),
        Eigen::AngleAxis<real>(Pi() / 2, Vector3r::UnitY()).toRotationMatrix(),
        Vector3r::Zero(), Vector3r::Zero(), link_opt);
    // Add fingers.
    for (integer d = 0; d < finger_num; ++d) {
        const Matrix3r R = Eigen::AngleAxis<real>(d * 2 * Pi() / finger_num,
            Vector3r::UnitZ()).toRotationMatrix();
        for (integer i = 0; i < finger_link_num; ++i) {
            const integer link_idx = sim.AddLink(link::LinkType::kBox, R *
                Vector3r((i + 0.5) * link_length, 0, initial_gripper_height),
                R, Vector3r::Zero(), Vector3r::Zero(), link_opt);
        }
    }
    // Add sphere.
    link_opt.real_option()["radius"] = sphere_radius;
    link_opt.real_option()["density"] = sphere_density;
    const integer sphere_idx = sim.AddLink(link::LinkType::kSphere,
        Vector3r(0, 0, sphere_radius), Matrix3r::Identity(), Vector3r::Zero(),
        Vector3r::Zero(), link_opt);
    // Add all hinge joints.
    Options joint_opt;
    for (integer d = 0; d < finger_num; ++d) {
        joint_opt.vector_option()["first_attach_point"] =
            Vector3r(link_length / 2, 0, 0);
        joint_opt.vector_option()["second_attach_point"] =
            Vector3r(-link_length / 2, 0, 0);
        joint_opt.vector_option()["first_attach_direction"] =
            Eigen::AngleAxis<real>(d * 2 * Pi() / finger_num,
                -Vector3r::UnitX()) * Vector3r::UnitY();
        joint_opt.vector_option()["second_attach_direction"] =
            Vector3r::UnitY();
        sim.AddJoint(joint::JointType::kBinaryHinge,
            { 0, 1 + d * finger_link_num }, joint_compliance, joint_opt);
        for (integer i = 0; i < finger_link_num - 1; ++i) {
            joint_opt.vector_option()["first_attach_direction"] =
                Vector3r::UnitY();
            sim.AddJoint(joint::JointType::kBinaryHinge, {
                1 + d * finger_link_num + i, 1 + d * finger_link_num + i + 1 },
                joint_compliance, joint_opt);
        }
    }
    // Add link contact proxies.
    Options contact_proxy_opt;
    for (integer d = 0; d < finger_num; ++d) {
        for (integer i = 0; i < finger_link_num; ++i) {
            const integer link_idx = 1 + d * finger_link_num + i;
            contact_proxy_opt.vector_option()["size"] =
                Vector3r(link_length, link_width, link_width);
            contact_proxy_opt.vector_option()["translation"] = Vector3r::Zero();
            contact_proxy_opt.matrix_option()["rotation"] =
                Matrix3r::Identity();
            sim.AddLinkContactProxy(link_idx, contact::ContactProxyType::kBox,
                contact_proxy_opt);
        }
    }
    contact_proxy_opt.real_option()["radius"] = sphere_radius;
    contact_proxy_opt.vector_option()["center"] = Vector3r::Zero();
    sim.AddLinkContactProxy(sphere_idx, contact::ContactProxyType::kSphere,
        contact_proxy_opt);
    // Add environmental contact proxy.
    contact_proxy_opt.real_option()["offset"] = 0;
    contact_proxy_opt.vector_option()["normal"] = Vector3r::UnitZ();
    const integer ground_idx = sim.AddEnvironmentalContactProxy(
        contact::ContactProxyType::kPlane, contact_proxy_opt);
    // Add contact proxy pairs.
    for (integer d = 0; d < finger_num; ++d) {
        for (integer i = 0; i < finger_link_num; ++i) {
            sim.AddContactProxyPair(sphere_idx, 1 + d * finger_link_num + i,
                contact_compliance, frictional_coefficient);
        }
    }
    sim.AddContactProxyPair(sphere_idx, ground_idx, contact_compliance,
        frictional_coefficient);
    // Compute gripper mass.
    real gripper_mass = 0;
    for (integer i = 0; i < 1 + finger_num * finger_link_num; ++i)
        gripper_mass += sim.link(i)->mass();

    // A helper function for computing hinge joint angles.
    const auto get_hinge_angle = [&](const integer joint_index) -> const real {
        const sim::JointInfo& joint_info = sim.joint(joint_index);
        const std::shared_ptr<joint::BinaryHingeJoint> joint =
            std::dynamic_pointer_cast<joint::BinaryHingeJoint>(
                joint_info.joint);
        const integer first_index = joint_info.link_indices[0];
        const integer second_index = joint_info.link_indices[1];
        const Vector3r c1_p = -sim.link(first_index)->ToWorldVector(
            joint->first_attach_point());
        const Vector3r c2_p = -sim.link(second_index)->ToWorldVector(
            joint->second_attach_point());
        const Vector3r d = sim.link(first_index)->ToWorldVector(
            joint->first_attach_direction());
        const Vector3r n1 = d.cross(c1_p).stableNormalized();
        const Vector3r n2 = d.cross(c2_p).stableNormalized();
        // n1 -> x.
        // d -> z.
        const Vector3r& x = n1;
        const Vector3r& z = d;
        const Vector3r y = z.cross(x);
        const real angle = std::atan2(n2.dot(y), n2.dot(x));
        return (angle > 0) ? (angle - Pi()) : (angle + Pi());
    };

    // Collect initial states --- this is useful for resetting the sim.
    const integer link_num = sim.link_num();
    std::vector<Vector3r> init_t(link_num);
    std::vector<Matrix3r> init_R(link_num);
    std::vector<Vector3r> init_v(link_num);
    std::vector<Vector3r> init_omega(link_num);
    for (integer i = 0; i < link_num; ++i) {
        init_t[i] = sim.link(i)->t();
        init_R[i] = sim.link(i)->R();
        init_v[i] = sim.link(i)->v();
        init_omega[i] = sim.link(i)->omega();
    }

    // Register visualizer.
    polyscope::init();
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::YFront);
    polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Free);
    polyscope::view::lookAt({ 0.0, -3.0, 0.4 }, { 0.0, 0.0, 0.1 });
    // For taking screenshots.
    polyscope::options::screenshotExtension = ".jpg";

    // Scene extent.
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = 1.;
    polyscope::state::boundingBox =
        std::tuple<glm::vec3, glm::vec3>{ { -2.0, -2.0, 0.0 },
            { 2.0, 2.0, 1.0 } };

    polyscope::removeAllStructures();
    // Register surface mesh structures.
    for (integer i = 0; i < link_num; ++i) {
        polyscope::registerSurfaceMesh("link_" + std::to_string(i),
            sim.link(i)->ToWorldPoints(sim.link(i)->vertices()).transpose(),
            sim.link(i)->elements().transpose());
        polyscope::getSurfaceMesh("link_" + std::to_string(i))
            ->setSurfaceColor({ 0x4C / 255.0, 0xC9 / 255.0, 0xFE / 255.0 });
        polyscope::getSurfaceMesh("link_" + std::to_string(i))
            ->setEdgeWidth(1.0);
    }

    // The following UI design is borrowed from the SIGGRAPH course:
    // https://github.com/siggraphcontact/rigidBodyTutorial.
    bool paused = true;
    bool step_once = false;
    Options sim_opt;
    const auto callback = [&]() -> void {
        // Key presed.
        if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_LeftArrow)) {
            target_pos.x() -= 0.01;
        } else if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_RightArrow)) {
            target_pos.x() += 0.01;
        }
        if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_DownArrow)) {
            target_pos.z() -= 0.01;
        } else if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_UpArrow)) {
            target_pos.z() += 0.01;
        }
        if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Q)) {
            target_angle -= 1.0 / 180 * Pi();
        } else if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_W)) {
            target_angle += 1.0 / 180 * Pi();
        }

        ImGui::PushItemWidth(100);

        ImGui::SliderFloat("Time step", &time_step, 0.001f, 0.1f, "%.3f");
        ImGui::SliderInt("Number of substeps", &substep_num, 1, 20, "%u");
        ImGui::SliderInt("BLCP iterations", &max_blcp_iter, 1, 100, "%u");
        if (ImGui::Button("Reset")) {
            sim.ClearForces();
            sim.ClearTorques();
            for (integer i = 0; i < link_num; ++i) {
                sim.ResetLinkPositionsAndVelocities(i, init_t[i],
                    init_R[i], init_v[i], init_omega[i]);
            }
            paused = true;
            step_once = false;
            x_ctrl.Reset();
            y_ctrl.Reset();
            z_ctrl.Reset();
            target_pos = init_target_pos;
            for (auto& c : hinge_angle_ctrl) c.Reset();
            target_angle = 0;
            for (real& a : last_angle) a = 0;
        }
        ImGui::SameLine();
        if (ImGui::Button("Step once")) step_once = true;

        ImGui::Checkbox("Pause", &paused);

        Vector3r base_pos = sim.link(0)->t();
        Vector3r base_vel = sim.link(0)->v();
        sim::StepInfo step_info;
        if (!paused || step_once) {
            Tic();
            // Keep simulating.
            sim_opt.integer_option()["max_blcp_iter"] = max_blcp_iter;
            sim_opt.real_option()["blcp_impulse_abs_tol"] = blcp_ip_tol;
            for (integer i = 0; i < substep_num; ++i) {
                sim.ClearForces();
                sim.ClearTorques();
                // Apply control force.
                base_pos = sim.link(0)->t();
                base_vel = sim.link(0)->v();
                Vector3r ctrl_force(0, 0, 0);
                ctrl_force.x() = x_ctrl.ComputeControlSignal(
                    target_pos.x() - base_pos.x(), -base_vel.x()
                );
                ctrl_force.y() = y_ctrl.ComputeControlSignal(
                    target_pos.y() - base_pos.y(), -base_vel.y()
                );
                ctrl_force.z() = z_ctrl.ComputeControlSignal(
                    target_pos.z() - base_pos.z(), -base_vel.z()
                );
                sim.ApplyForce(0, Vector3r::Zero(), ctrl_force);
                // Apply torque.
                for (integer j = 0; j < finger_num * finger_link_num; ++j) {
                    const sim::JointInfo& joint_info = sim.joint(j);
                    const std::shared_ptr<joint::BinaryHingeJoint> joint =
                        std::dynamic_pointer_cast<joint::BinaryHingeJoint>(
                            joint_info.joint);
                    const integer first_index = joint_info.link_indices[0];
                    const integer second_index = joint_info.link_indices[1];
                    const real angle = get_hinge_angle(j);
                    const Vector3r tq_dir =
                        sim.link(first_index)->ToWorldVector(
                            joint->first_attach_direction());
                    const real target_angle_j = (j % finger_link_num == 0) ?
                        -target_angle : target_angle;
                    const real tq = hinge_angle_ctrl[j].ComputeControlSignal(
                        target_angle_j - angle,
                        -(angle - last_angle[j]) / time_step);
                    sim.ApplyTorque(first_index, -tq * tq_dir);
                    sim.ApplyTorque(second_index, tq * tq_dir);
                }
                // Apply gravity.
                for (integer j = 0; j < link_num; ++j) {
                    sim.ApplyForce(j, Vector3r::Zero(),
                        sim.link(j)->mass() * gravitational_acceleration);
                }
                step_info = sim.Step(ToReal(time_step) / substep_num, sim_opt);
                // Record last angles.
                for (integer i = 0; i < finger_num * finger_link_num; ++i) {
                    last_angle[i] = get_hinge_angle(i);
                }
            }

            step_once = false;
        }

        ImGui::Text("Target pos: (%.2f, %.2f, %.2f).",
            target_pos.x(), target_pos.y(), target_pos.z());
        ImGui::Text("Base pos: (%.2f, %.2f, %.2f).",
            base_pos.x(), base_pos.y(), base_pos.z());
        ImGui::Text("Target hinge angle: %.0f degrees.",
            180.0 / Pi() * target_angle);
        for (integer i = 0; i < finger_num * finger_link_num; ++i) {
            ImGui::Text("Hinge joint %d angle: %.0f degrees.", i,
                180.0 / Pi() * last_angle[i]);
        }

        // Visualization.
        // The robot itself:
        for (integer i = 0; i < link_num; ++i) {
            polyscope::getSurfaceMesh("link_" + std::to_string(i)
                )->updateVertexPositions(sim.link(i)->ToWorldPoints(
                    sim.link(i)->vertices()).transpose());
        }
        // The contact information.
        if (!step_info.contacts.empty()) {
            std::vector<Vector3r> contact_points;
            std::vector<Vector3r> contact_normals;
            std::vector<Vector3r> contact_fronts;
            for (const auto& c : step_info.contacts) {
                contact_points.push_back(c.event.collision_point);
                contact_normals.push_back(c.event.local_frame.col(2));
                contact_fronts.push_back(c.event.local_frame.col(0));
            }
            polyscope::removePointCloud("contacts");
            auto pc = polyscope::registerPointCloud("contacts", contact_points);
            pc->setPointColor({ 1.0, 0.0, 0.0 });
            pc->setPointRadius(0.01);
            pc->addVectorQuantity("normal", contact_normals)->setVectorColor(
                { 1.0, 1.0, 0.0 })->setVectorLengthScale(
                    0.1)->setEnabled(true);
            pc->addVectorQuantity("front", contact_fronts)->setVectorColor(
                { 0.0, 1.0, 0.0 })->setVectorLengthScale(
                    0.1)->setEnabled(true);
        }

        ImGui::PopItemWidth();

        // Call the line below to take a screenshot.
        polyscope::screenshot();        // open here
    };

    polyscope::state::userCallback = callback;
    polyscope::show();

    PrintInfo("phys_sim_constrained_dynamics::sim::TestGripper",
        "TestGripper ends...");
}

}

int main() {
    phys_sim_constrained_dynamics::TestGripper();
    return 0;
}