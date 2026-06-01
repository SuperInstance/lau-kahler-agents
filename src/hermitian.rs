//! Hermitian form: a complex bilinear (sesquilinear) form on a complex vector space.

use nalgebra::{DMatrix, DVector, Complex};
use serde::{Serialize, Deserialize};
use crate::complex::ComplexStructure;

/// A Hermitian form h on C^n, represented as an n×n Hermitian matrix.
///
/// h(u, v) = u^* H v, with H = H^* (conjugate transpose equals self).
/// The real part gives the Riemannian metric, the imaginary part gives the symplectic form.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HermitianForm {
    /// Hermitian matrix H (complex n×n).
    pub matrix: DMatrix<Complex<f64>>,
    /// Complex dimension.
    pub dim: usize,
}

impl HermitianForm {
    /// Create from a complex Hermitian matrix.
    pub fn new(matrix: DMatrix<Complex<f64>>) -> Self {
        let dim = matrix.nrows();
        Self { matrix, dim }
    }

    /// Create the standard Hermitian form (identity matrix).
    pub fn standard(n: usize) -> Self {
        let mut mat = DMatrix::from_element(n, n, Complex::new(0.0, 0.0));
        for i in 0..n {
            mat[(i, i)] = Complex::new(1.0, 0.0);
        }
        Self {
            matrix: mat,
            dim: n,
        }
    }

    /// Verify the matrix is actually Hermitian: H = H^*.
    pub fn is_hermitian(&self) -> bool {
        let conj_t = self.matrix.adjoint();
        (self.matrix.clone() - conj_t).norm() < 1e-10
    }

    /// Evaluate h(u, v) = u^* H v.
    pub fn evaluate(&self, u: &DVector<Complex<f64>>, v: &DVector<Complex<f64>>) -> Complex<f64> {
        let hu = self.matrix.transpose() * u;
        let mut result = Complex::new(0.0, 0.0);
        for i in 0..self.dim {
            result += hu[i].conj() * v[i];
        }
        result
    }

    /// Extract the Riemannian metric g(u,v) = Re(h(u,v)).
    /// Returns a 2n×2n real matrix.
    pub fn riemannian_metric(&self) -> DMatrix<f64> {
        let n = self.dim;
        let mut g = DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                let h_ij = self.matrix[(i, j)];
                // g = [[Re(H), -Im(H)], [Im(H), Re(H)]]
                g[(i, j)] = h_ij.re;
                g[(n + i, n + j)] = h_ij.re;
                g[(i, n + j)] = -h_ij.im;
                g[(n + i, j)] = h_ij.im;
            }
        }
        g
    }

    /// Extract the symplectic (Kähler) form ω(u,v) = Im(h(u,v)).
    /// Returns a 2n×2n real matrix.
    pub fn kahler_form(&self) -> DMatrix<f64> {
        let n = self.dim;
        let mut omega = DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                let h_ij = self.matrix[(i, j)];
                // ω = [[Im(H), Re(H)], [-Re(H), Im(H)]]
                omega[(i, j)] = h_ij.im;
                omega[(n + i, n + j)] = h_ij.im;
                omega[(i, n + j)] = h_ij.re;
                omega[(n + i, j)] = -h_ij.re;
            }
        }
        omega
    }

    /// Verify positive definiteness: h(v,v) > 0 for all nonzero v.
    /// Checks via eigenvalues of the real metric.
    pub fn is_positive_definite(&self) -> bool {
        let g = self.riemannian_metric();
        let eigenvals = g.symmetric_eigenvalues();
        eigenvals.iter().all(|&e| e > -1e-10)
    }

    /// Compatibility: g(u,v) = ω(u, Jv).
    pub fn verify_kahler_compatibility(&self, j: &ComplexStructure) -> bool {
        let g = self.riemannian_metric();
        let omega = self.kahler_form();
        let gjv = &omega * &j.matrix;
        (g - gjv).norm() < 1e-10
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_hermitian_is_hermitian() {
        let h = HermitianForm::standard(3);
        assert!(h.is_hermitian());
    }

    #[test]
    fn test_standard_hermitian_positive_definite() {
        let h = HermitianForm::standard(3);
        assert!(h.is_positive_definite());
    }

    #[test]
    fn test_hermitian_evaluate() {
        let h = HermitianForm::standard(2);
        let u = DVector::from_vec(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);
        let v = DVector::from_vec(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);
        let val = h.evaluate(&u, &v);
        assert!((val.re - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_riemannian_metric_identity() {
        let h = HermitianForm::standard(2);
        let g = h.riemannian_metric();
        let expected = DMatrix::identity(4, 4);
        assert!((g - expected).norm() < 1e-10);
    }

    #[test]
    fn test_kahler_form_standard() {
        let h = HermitianForm::standard(2);
        let omega = h.kahler_form();
        // Standard symplectic form ω = [[0,I],[-I,0]] = -J = J^T
        let j = ComplexStructure::standard(2);
        assert!((omega - j.matrix.transpose()).norm() < 1e-10);
    }

    #[test]
    fn test_kahler_compatibility() {
        let h = HermitianForm::standard(2);
        let j = ComplexStructure::standard(2);
        assert!(h.verify_kahler_compatibility(&j));
    }

    #[test]
    fn test_hermitian_serialization() {
        let h = HermitianForm::standard(2);
        let json = serde_json::to_string(&h).unwrap();
        let h2: HermitianForm = serde_json::from_str(&json).unwrap();
        assert!(h2.is_hermitian());
    }

    #[test]
    fn test_custom_hermitian() {
        // H = [[2, 1+i], [1-i, 3]]
        let h = HermitianForm::new(DMatrix::from_row_slice(2, 2, &[
            Complex::new(2.0, 0.0), Complex::new(1.0, 1.0),
            Complex::new(1.0, -1.0), Complex::new(3.0, 0.0),
        ]));
        assert!(h.is_hermitian());
        assert!(h.is_positive_definite());
    }
}
