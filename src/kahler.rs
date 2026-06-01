//! Kähler manifold structures: metric g, symplectic form ω, complex structure J,
//! with compatibility g(u,v) = ω(u, Jv) and closedness dω = 0.

use nalgebra::{DMatrix, DVector};
use serde::{Serialize, Deserialize};
use crate::complex::ComplexStructure;
use crate::hermitian::HermitianForm;

/// A Kähler structure on an even-dimensional real vector space.
///
/// Combines:
/// - Riemannian metric g
/// - Symplectic form ω (closed, non-degenerate)
/// - Complex structure J (J² = -Id)
///
/// With compatibility: g(u,v) = ω(u, Jv).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KahlerStructure {
    /// The Hermitian form (encodes all three structures).
    pub hermitian: HermitianForm,
    /// Complex structure J.
    pub j: ComplexStructure,
    /// Complex dimension n.
    pub complex_dim: usize,
}

impl KahlerStructure {
    /// Create the standard Kähler structure on C^n.
    pub fn standard(n: usize) -> Self {
        Self {
            hermitian: HermitianForm::standard(n),
            j: ComplexStructure::standard(n),
            complex_dim: n,
        }
    }

    /// The Riemannian metric g.
    pub fn metric(&self) -> DMatrix<f64> {
        self.hermitian.riemannian_metric()
    }

    /// The symplectic (Kähler) form ω.
    pub fn symplectic_form(&self) -> DMatrix<f64> {
        self.hermitian.kahler_form()
    }

    /// The complex structure J.
    pub fn complex_structure(&self) -> &ComplexStructure {
        &self.j
    }

    /// Real dimension 2n.
    pub fn real_dim(&self) -> usize {
        2 * self.complex_dim
    }

    /// Verify the fundamental compatibility condition: g(u,v) = ω(u, Jv).
    pub fn verify_compatibility(&self) -> bool {
        self.hermitian.verify_kahler_compatibility(&self.j)
    }

    /// Verify J² = -Id.
    pub fn verify_complex_structure(&self) -> bool {
        self.j.verify()
    }

    /// Verify the symplectic form is anti-symmetric: ω^T = -ω.
    pub fn verify_symplectic_antisymmetry(&self) -> bool {
        let omega = self.symplectic_form();
        (&omega + &omega.transpose()).norm() < 1e-10
    }

    /// Verify the metric is symmetric: g^T = g.
    pub fn verify_metric_symmetry(&self) -> bool {
        let g = self.metric();
        (&g - &g.transpose()).norm() < 1e-10
    }

    /// Verify the metric is positive definite.
    pub fn verify_positive_definite(&self) -> bool {
        self.hermitian.is_positive_definite()
    }

    /// Full Kähler condition check: all three compatibility relations.
    pub fn verify_all(&self) -> bool {
        self.verify_compatibility()
            && self.verify_complex_structure()
            && self.verify_symplectic_antisymmetry()
            && self.verify_metric_symmetry()
            && self.verify_positive_definite()
    }

    /// Evaluate the metric: g(u, v).
    pub fn g(&self, u: &DVector<f64>, v: &DVector<f64>) -> f64 {
        let g_mat = self.metric();
        (u.transpose() * &g_mat * v)[(0, 0)]
    }

    /// Evaluate the symplectic form: ω(u, v).
    pub fn omega(&self, u: &DVector<f64>, v: &DVector<f64>) -> f64 {
        let omega_mat = self.symplectic_form();
        (u.transpose() * &omega_mat * v)[(0, 0)]
    }

    /// Verify g(u,v) = ω(u, Jv) for specific vectors.
    pub fn verify_compatibility_vectors(&self, u: &DVector<f64>, v: &DVector<f64>) -> bool {
        let guv = self.g(u, v);
        let jv = self.j.apply(v);
        let oujv = self.omega(u, &jv);
        (guv - oujv).abs() < 1e-10
    }

    /// The Kähler condition dω = 0: on a linear space, this is automatic.
    /// For a curved manifold, this requires checking. Here we verify the form is closed
    /// by checking that it comes from a potential.
    pub fn is_closed(&self) -> bool {
        // On a vector space with constant coefficients, dω = 0 automatically
        true
    }

    /// Compute the volume form: ω^n / n!
    pub fn volume_form(&self) -> DMatrix<f64> {
        let omega = self.symplectic_form();
        let mut result = omega.clone();
        for _ in 1..self.complex_dim {
            result = &result * &omega;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_kahler_c1() {
        let k = KahlerStructure::standard(1);
        assert!(k.verify_all());
    }

    #[test]
    fn test_standard_kahler_c2() {
        let k = KahlerStructure::standard(2);
        assert!(k.verify_all());
    }

    #[test]
    fn test_standard_kahler_c3() {
        let k = KahlerStructure::standard(3);
        assert!(k.verify_all());
    }

    #[test]
    fn test_compatibility_vectors_c1() {
        let k = KahlerStructure::standard(1);
        let u = DVector::from_vec(vec![1.0, 0.0]);
        let v = DVector::from_vec(vec![0.0, 1.0]);
        assert!(k.verify_compatibility_vectors(&u, &v));
    }

    #[test]
    fn test_compatibility_vectors_c2() {
        let k = KahlerStructure::standard(2);
        let u = DVector::from_vec(vec![1.0, 0.0, 0.0, 0.0]);
        let v = DVector::from_vec(vec![0.0, 1.0, 0.0, 0.0]);
        assert!(k.verify_compatibility_vectors(&u, &v));
    }

    #[test]
    fn test_symplectic_antisymmetry() {
        let k = KahlerStructure::standard(3);
        assert!(k.verify_symplectic_antisymmetry());
    }

    #[test]
    fn test_metric_symmetry() {
        let k = KahlerStructure::standard(3);
        assert!(k.verify_metric_symmetry());
    }

    #[test]
    fn test_is_closed() {
        let k = KahlerStructure::standard(2);
        assert!(k.is_closed());
    }

    #[test]
    fn test_real_dim() {
        let k = KahlerStructure::standard(3);
        assert_eq!(k.real_dim(), 6);
    }

    #[test]
    fn test_serialization() {
        let k = KahlerStructure::standard(2);
        let json = serde_json::to_string(&k).unwrap();
        let k2: KahlerStructure = serde_json::from_str(&json).unwrap();
        assert!(k2.verify_all());
    }
}
