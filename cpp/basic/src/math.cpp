#include "basic/include/math.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_constrained_dynamics {

const real ToReal(const float value) {
    return static_cast<real>(value);
}

const real ToReal(const double value) {
    return static_cast<real>(value);
}

const real ToReal(const integer value) {
    return static_cast<real>(value);
}

const double ToDouble(const real value) {
    return static_cast<double>(value);
}

const float ToFloat(const real value) {
    return static_cast<float>(value);
}

const real Pi() {
    return ToReal(3.141592653589793238462643383);
}

const bool IsClose(const real a, const real b, const real abs_tol,
    const real rel_tol) {
    return std::abs(a - b) <= std::abs(b) * rel_tol + abs_tol;
}

const bool IsClose(const VectorXr& a, const VectorXr& b, const real abs_tol,
    const real rel_tol) {
    CheckCondition(a.size() == b.size(), "basic::math::IsClose",
        "Sizes of two vectors are different.");
    return ((a - b).cwiseAbs().array() <= b.cwiseAbs().array() * rel_tol
        + abs_tol).all();
}

const Matrix3r BuildFrameFromUnitNormal(const Vector3r& n) {
    real max_len = -1;
    Matrix3r R; R.setZero();
    R.col(2) = n;
    for (integer i = 0; i < 3; ++i) {
        const Vector3r x = n.cross(Vector3r::Unit(i));
        const real x_len = x.norm();
        if (x_len > max_len) {
            max_len = x_len;
            R.col(0) = x / x_len;
        }
    }
    R.col(1) = n.cross(R.col(0));
    return R;
}

const Matrix3r ToCrossProductMatrix(const Vector3r& a) {
    Matrix3r A;
    A << 0, -a(2), a(1),
        a(2), 0, -a(0),
        -a(1), a(0), 0;
    return A;
}

const Matrix3r ToRotationMatrix(const Vector3r& angle) {
    // Define
    //     a_0 = cos(theta).
    //     a_1 = sin(theta) / theta.
    //     a_2 = (1 - cos(theta)) / theta^2.
    // Then, R = I + a_1(theta) A + a_2(theta) A^2. Here A = [axis_].
    // To see this, consider any vector v:
    // Rv = v + sin(theta) * unit_axis x v +
    //      (1 - cos(theta)) unit_axis x (unit_axis x v).
    const real theta = angle.stableNorm();
    const real theta_sqr = theta * theta;
    const real a0 = std::cos(theta);
    // 1 rad = 57 degrees.
    // 2e-3 rad = 0.1 degrees, which should be sufficiently small for motions
    // we care about.
    const real tol = ToReal(2e-3);
    real a1 = 0;
    real a2 = 0;
    if (theta < tol) {
        // sinx = x - x3 / 3! + x5 / 5! - ...
        // sinx / x = 1 - x2 / 3! + x4 / 5! - ...
        a1 = 1 - theta_sqr / 6;
        // sinx / x = a1 + O(theta^4).
        // cosx = 1 - x2 / 2! + x4 / 4! - x6 / 6! + ...
        // (1 - cosx) / x2 = 1 / 2! - x2 / 4! + x4 / 6! - ...
        a2 = ToReal(0.5) - theta_sqr / 24;
        // (1 - cosx) / x2 = a2 + O(theta^4).
    } else {
        a1 = std::sin(theta) / theta;
        a2 = (1 - std::cos(theta)) / theta_sqr;
    }
    const Matrix3r A = ToCrossProductMatrix(angle);
    return Matrix3r::Identity() + a1 * A + a2 * A * A;
}

// Sparse matrices.
const SparseMatrixXr FromTriplet(const integer row_num, const integer col_num,
    const std::vector<Eigen::Triplet<real>>& nonzeros) {
    SparseMatrixXr mat(row_num, col_num);
    mat.setFromTriplets(nonzeros.begin(), nonzeros.end());
    mat.makeCompressed();
    return mat;
}

const std::vector<Eigen::Triplet<real>> ToTriplet(const SparseMatrixXr& mat) {
    SparseMatrixXr mat_compressed = mat;
    mat_compressed.makeCompressed();
    std::vector<Eigen::Triplet<real>> nonzeros;
    const integer outer_size = static_cast<integer>(mat_compressed.outerSize());
    for (integer k = 0; k < outer_size; ++k)
        for (SparseMatrixXr::InnerIterator it(mat_compressed, k); it; ++it) {
            nonzeros.push_back(Eigen::Triplet<real>(it.row(), it.col(),
                it.value()));
        }
    return nonzeros;
}

const SparseMatrixXr FromDiagonal(const VectorXr& diagonal) {
    CheckCondition(diagonal.size() > 0, "basic::FromDiagonal",
        "Empty diagonal.");
    const integer matrix_size = static_cast<integer>(diagonal.size());
    std::vector<Eigen::Triplet<real>> nonzeros;
    for (integer i = 0; i < matrix_size; ++i) {
        nonzeros.push_back(Eigen::Triplet<real>(i, i, diagonal(i)));
    }
    return FromTriplet(matrix_size, matrix_size, nonzeros);
}

const bool IsRotationMatrix(const Matrix3r& R) {
    const Matrix3r I = Matrix3r::Identity();
    return IsClose(R.determinant(), 1, ToReal(1e-3), 0) &&
        IsClose((R * R.transpose() - I).squaredNorm() / (3 * 3), 0,
            ToReal(1e-6), 0);
}

}