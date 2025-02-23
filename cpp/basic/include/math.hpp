#ifndef PHYS_SIM_CONSTRAINED_DYNAMICS_BASIC_MATH
#define PHYS_SIM_CONSTRAINED_DYNAMICS_BASIC_MATH

#include "basic/include/config.hpp"
#include "basic/include/options.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_constrained_dynamics {

const real ToReal(const float value);
const real ToReal(const double value);
const real ToReal(const integer value);
const double ToDouble(const real value);
const float ToFloat(const real value);

const real Pi();

// |a - b| <= abs_tol + rel_tol * |b|.
const bool IsClose(const real a, const real b, const real abs_tol,
    const real rel_tol);

// For all i, |ai - bi| <= abs_tol + rel_tol * |bi|.
const bool IsClose(const VectorXr& a, const VectorXr& b, const real abs_tol,
    const real rel_tol);

// Build a local frame.
const Matrix3r BuildFrameFromUnitNormal(const Vector3r& n);

const Matrix3r ToCrossProductMatrix(const Vector3r& a);

const Matrix3r ToRotationMatrix(const Vector3r& angle);

// Sparse matrices.
const SparseMatrixXr FromTriplet(const integer row_num, const integer col_num,
    const std::vector<Eigen::Triplet<real>>& nonzeros);
const std::vector<Eigen::Triplet<real>> ToTriplet(const SparseMatrixXr& mat);
const SparseMatrixXr FromDiagonal(const VectorXr& diagonal);

// Rotation matrix.
const bool IsRotationMatrix(const Matrix3r& R);

}

#endif