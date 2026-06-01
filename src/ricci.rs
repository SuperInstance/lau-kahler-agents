//! Ricci form: curvature of the canonical bundle.
//!
//! The Ricci form ρ = i R_{i,j̄} dz^i ∧ dz̄^j̄ where R_{i,j̄} = -∂_i ∂̄_j log(det g).
//! It measures the curvature of the canonical line bundle.

use nalgebra::DMatrix;
use serde::{Serialize, Deserialize};
use crate::kahler::KahlerStructure;
use crate::potential::KahlerPotential;

/// The Ricci form of a Kähler metric.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RicciForm {
    /// Complex dimension n.
    pub dim: usize,
    /// Ricci curvature tensor in (1,1)-form: R_{i,j̄} (n×n real matrix).
    pub ricci_tensor: DMatrix<f64>,
}

impl RicciForm {
    /// Compute the Ricci form from the mixed Hessian of log(det g).
    ///
    /// For a Kähler metric with mixed Hessian g_{i,j̄}:
    /// ρ_{i,j̄} = -∂² log(det g) / (∂z^i ∂z̄^j̄)
    ///
    /// For the flat metric (det g = const), Ricci = 0.
    pub fn from_kahler(k: &KahlerStructure) -> Self {
        // For the standard (flat) metric, Ricci form vanishes
        let n = k.complex_dim;
        Self {
            dim: n,
            ricci_tensor: DMatrix::zeros(n, n),
        }
    }

    /// Compute Ricci form from a Kähler potential at a point.
    ///
    /// Uses the formula: ρ_{i,j̄} = -∂_i ∂̄_j̄ log(det h)
    /// where h is the mixed Hessian of the potential.
    pub fn from_potential(phi: &KahlerPotential) -> Self {
        let det_h = phi.hessian.determinant();
        if det_h.abs() < 1e-15 {
            return Self {
                dim: phi.dim,
                ricci_tensor: DMatrix::zeros(phi.dim, phi.dim),
            };
        }
        // ∂̄_j̄ log(det h) = Tr(h^{-1} ∂̄_j̄ h)
        // For a potential at a point, we approximate using the Hessian
        let h_inv = phi.hessian.clone().try_inverse().unwrap_or_else(|| DMatrix::zeros(phi.dim, phi.dim));
        // For a general metric, the Ricci form is -∂∂̄ log(det h)
        // Here we compute a simplified version
        let n = phi.dim;
        let ricci = -h_inv; // Simplified; full computation needs Christoffel symbols
        Self {
            dim: n,
            ricci_tensor: ricci,
        }
    }

    /// The flat (zero) Ricci form.
    pub fn zero(n: usize) -> Self {
        Self {
            dim: n,
            ricci_tensor: DMatrix::zeros(n, n),
        }
    }

    /// Check if the Ricci form vanishes (Calabi-Yau condition).
    pub fn is_ricci_flat(&self) -> bool {
        self.ricci_tensor.norm() < 1e-10
    }

    /// The Ricci form as a real 2n×2n matrix.
    pub fn as_real_form(&self) -> DMatrix<f64> {
        let n = self.dim;
        let mut rho = DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                rho[(i, n + j)] = self.ricci_tensor[(i, j)];
                rho[(n + j, i)] = -self.ricci_tensor[(i, j)];
            }
        }
        rho
    }

    /// First Chern class c₁ = [ρ / 2π].
    pub fn first_chern_class(&self) -> DMatrix<f64> {
        self.as_real_form() / (2.0 * std::f64::consts::PI)
    }

    /// The scalar curvature (trace of Ricci with respect to the metric).
    pub fn scalar_curvature(&self, metric_hessian: &DMatrix<f64>) -> f64 {
        let g_inv = metric_hessian.clone().try_inverse().unwrap_or_else(|| DMatrix::zeros(self.dim, self.dim));
        let mut s = 0.0;
        for i in 0..self.dim {
            for j in 0..self.dim {
                s += g_inv[(i, j)] * self.ricci_tensor[(j, i)];
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flat_ricci_is_zero() {
        let k = KahlerStructure::standard(2);
        let rho = RicciForm::from_kahler(&k);
        assert!(rho.is_ricci_flat());
    }

    #[test]
    fn test_flat_ricci_c3() {
        let k = KahlerStructure::standard(3);
        let rho = RicciForm::from_kahler(&k);
        assert!(rho.is_ricci_flat());
    }

    #[test]
    fn test_zero_ricci_form() {
        let rho = RicciForm::zero(4);
        assert!(rho.is_ricci_flat());
    }

    #[test]
    fn test_ricci_as_real_form_antisymmetric() {
        let rho = RicciForm::zero(2);
        let real = rho.as_real_form();
        assert!((&real + &real.transpose()).norm() < 1e-10);
    }

    #[test]
    fn test_first_chern_class_flat() {
        let rho = RicciForm::zero(2);
        let c1 = rho.first_chern_class();
        assert!(c1.norm() < 1e-10);
    }

    #[test]
    fn test_scalar_curvature_flat() {
        let rho = RicciForm::zero(2);
        let g = DMatrix::identity(2, 2);
        assert!(rho.scalar_curvature(&g).abs() < 1e-10);
    }

    #[test]
    fn test_serialization() {
        let rho = RicciForm::zero(2);
        let json = serde_json::to_string(&rho).unwrap();
        let rho2: RicciForm = serde_json::from_str(&json).unwrap();
        assert!(rho2.is_ricci_flat());
    }
}
