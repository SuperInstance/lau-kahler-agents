//! Kähler potential: ω = i∂∂̄Φ.
//!
//! The Kähler form can locally be written as ω = i ∂∂̄ Φ
//! for a real-valued function Φ called the Kähler potential.

use nalgebra::{DMatrix, DVector, Complex};
use serde::{Serialize, Deserialize};

/// A Kähler potential Φ on C^n.
///
/// The Kähler form is recovered as ω = i ∂∂̄ Φ.
/// In coordinates: ω_{i,j̄} = ∂²Φ / ∂z^i ∂z̄^j̄.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KahlerPotential {
    /// Complex dimension n.
    pub dim: usize,
    /// The mixed Hessian ∂²Φ/∂z^i∂z̄^j (n×n real matrix, since Φ is real).
    pub hessian: DMatrix<f64>,
}

impl KahlerPotential {
    /// Create from a positive-definite mixed Hessian matrix.
    pub fn new(hessian: DMatrix<f64>) -> Self {
        let dim = hessian.nrows();
        Self { dim, hessian }
    }

    /// The flat potential Φ = |z|² = Σ|z_i|², giving the standard Kähler form.
    pub fn flat(n: usize) -> Self {
        Self {
            dim: n,
            hessian: DMatrix::identity(n, n),
        }
    }

    /// Fubini-Study potential on CP^n (restricted to affine chart):
    /// Φ = log(1 + |z|²).
    /// Mixed Hessian: h_{i,j̄} = (δ_{ij}(1+|z|²) - z_i z̄_j) / (1+|z|²)²
    pub fn fubini_study(n: usize, z: &DVector<Complex<f64>>) -> Self {
        let zsq: f64 = z.iter().map(|zi| zi.norm_sqr()).sum();
        let denom = (1.0 + zsq).powi(2);
        let mut h = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let kron = if i == j { 1.0 } else { 0.0 };
                h[(i, j)] = (kron * (1.0 + zsq) - (z[i].re * z[j].re + z[i].im * z[j].im)) / denom;
            }
        }
        Self { dim: n, hessian: h }
    }

    /// Evaluate the potential at a point (for flat: Φ = |z|²).
    pub fn evaluate(&self, z: &DVector<Complex<f64>>) -> f64 {
        z.iter().map(|zi| zi.norm_sqr()).sum()
    }

    /// Recover the Kähler form ω from the potential.
    /// ω = i Σ h_{i,j̄} dz^i ∧ dz̄^j̄
    /// As a real 2n×2n matrix:
    pub fn kahler_form(&self) -> DMatrix<f64> {
        let n = self.dim;
        let mut omega = DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                // The block structure of ω from the mixed Hessian
                omega[(i, n + j)] = self.hessian[(i, j)];
                omega[(n + j, i)] = -self.hessian[(i, j)];
            }
        }
        omega
    }

    /// The associated metric g(u,v) = ω(u, Jv).
    /// For ω from potential, g_{i,j̄} = h_{i,j̄}.
    pub fn metric(&self) -> DMatrix<f64> {
        let n = self.dim;
        let mut g = DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                g[(i, j)] = self.hessian[(i, j)];
                g[(n + i, n + j)] = self.hessian[(i, j)];
            }
        }
        g
    }

    /// Check positive definiteness of the mixed Hessian.
    pub fn is_positive_definite(&self) -> bool {
        let eigenvals = self.hessian.symmetric_eigenvalues();
        eigenvals.iter().all(|&e| e > -1e-10)
    }

    /// The volume element: det(h_{i,j̄}).
    pub fn volume_element(&self) -> f64 {
        self.hessian.determinant().max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flat_potential() {
        let phi = KahlerPotential::flat(2);
        assert!(phi.is_positive_definite());
        let omega = phi.kahler_form();
        // Standard symplectic form ω = [[0,I],[-I,0]] = J^T
        let j = crate::complex::ComplexStructure::standard(2);
        assert!((omega - j.matrix.transpose()).norm() < 1e-10);
    }

    #[test]
    fn test_flat_potential_metric() {
        let phi = KahlerPotential::flat(3);
        let g = phi.metric();
        assert!((g - DMatrix::identity(6, 6)).norm() < 1e-10);
    }

    #[test]
    fn test_fubini_study_potential() {
        let z = DVector::from_vec(vec![Complex::new(0.0, 0.0); 2]);
        let phi = KahlerPotential::fubini_study(2, &z);
        // At origin, FS potential gives identity Hessian
        assert!((phi.hessian - DMatrix::identity(2, 2)).norm() < 1e-10);
    }

    #[test]
    fn test_fubini_study_positive_definite() {
        let z = DVector::from_vec(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 1.0)]);
        let phi = KahlerPotential::fubini_study(2, &z);
        assert!(phi.is_positive_definite());
    }

    #[test]
    fn test_evaluate_potential() {
        let phi = KahlerPotential::flat(2);
        let z = DVector::from_vec(vec![Complex::new(3.0, 4.0), Complex::new(0.0, 0.0)]);
        assert!((phi.evaluate(&z) - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_kahler_form_antisymmetric() {
        let phi = KahlerPotential::flat(3);
        let omega = phi.kahler_form();
        assert!((&omega + &omega.transpose()).norm() < 1e-10);
    }

    #[test]
    fn test_volume_element_flat() {
        let phi = KahlerPotential::flat(2);
        assert!((phi.volume_element() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_serialization() {
        let phi = KahlerPotential::flat(2);
        let json = serde_json::to_string(&phi).unwrap();
        let phi2: KahlerPotential = serde_json::from_str(&json).unwrap();
        assert!(phi2.is_positive_definite());
    }
}
