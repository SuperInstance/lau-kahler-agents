//! Complex structure J with J² = -Id.

use nalgebra::{DMatrix, DVector, Complex};
use serde::{Serialize, Deserialize};

/// A complex structure on an even-dimensional real vector space.
///
/// Satisfies J² = -Id, turning R^{2n} into a complex vector space C^n.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexStructure {
    /// The matrix representation of J (2n × 2n real matrix).
    pub matrix: DMatrix<f64>,
    /// Complex dimension n (real dimension is 2n).
    pub complex_dim: usize,
}

impl ComplexStructure {
    /// Create the standard complex structure on R^{2n}.
    ///
    /// In block form: J = [[0, -I], [I, 0]].
    pub fn standard(n: usize) -> Self {
        let dim = 2 * n;
        let mut j = DMatrix::zeros(dim, dim);
        for i in 0..n {
            j[(i, n + i)] = -1.0;
            j[(n + i, i)] = 1.0;
        }
        Self { matrix: j, complex_dim: n }
    }

    /// Verify J² = -Id.
    pub fn verify(&self) -> bool {
        let n = self.matrix.nrows();
        let j2 = &self.matrix * &self.matrix;
        let neg_id = -DMatrix::identity(n, n);
        (j2 - neg_id).norm() < 1e-10
    }

    /// Apply J to a vector.
    pub fn apply(&self, v: &DVector<f64>) -> DVector<f64> {
        &self.matrix * v
    }

    /// Complexify: interpret a real 2n-vector as a complex n-vector.
    pub fn to_complex_vector(&self, v: &DVector<f64>) -> DVector<Complex<f64>> {
        let n = self.complex_dim;
        let mut result = DVector::from_element(n, Complex::new(0.0, 0.0));
        for i in 0..n {
            result[i] = Complex::new(v[i], v[n + i]);
        }
        result
    }

    /// The (p,q)-projection for a complex vector.
    /// For now, returns the real and imaginary parts separately.
    pub fn pq_components(&self, v: &DVector<Complex<f64>>) -> (DVector<Complex<f64>>, DVector<Complex<f64>>) {
        let n = v.nrows();
        // (1,0)-part: (v - i*Jv)/2, simplified to Re and Im
        let mut p_part = DVector::from_element(n, Complex::new(0.0, 0.0));
        let mut q_part = DVector::from_element(n, Complex::new(0.0, 0.0));
        for i in 0..n {
            let z = v[i];
            p_part[i] = Complex::new(z.re, 0.0);
            q_part[i] = Complex::new(0.0, z.im);
        }
        (p_part, q_part)
    }

    /// Eigenvalues of J are ±i. Returns true if all eigenvalues are purely imaginary.
    pub fn has_pure_imaginary_eigenvalues(&self) -> bool {
        let eigenvals = self.matrix.complex_eigenvalues();
        eigenvals.iter().all(|ev| (ev.re.abs() < 1e-8) && (ev.im.abs() - 1.0).abs() < 1e-8)
    }

    /// Check if J is orthogonal: J^T J = I.
    pub fn is_orthogonal(&self) -> bool {
        let jtj = self.matrix.transpose() * &self.matrix;
        let n = self.matrix.nrows();
        (jtj - DMatrix::identity(n, n)).norm() < 1e-10
    }

    /// Compatibility with a metric g: g(Ju, Jv) = g(u, v).
    pub fn is_compatible_with_metric(&self, g: &DMatrix<f64>) -> bool {
        let jtgj = self.matrix.transpose() * g * &self.matrix;
        (jtgj - g).norm() < 1e-10
    }

    /// The real dimension (2n).
    pub fn real_dim(&self) -> usize {
        2 * self.complex_dim
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;

    #[test]
    fn test_standard_complex_structure_c1() {
        let j = ComplexStructure::standard(1);
        assert!(j.verify());
        assert_eq!(j.real_dim(), 2);
    }

    #[test]
    fn test_standard_complex_structure_c2() {
        let j = ComplexStructure::standard(2);
        assert!(j.verify());
        assert_eq!(j.real_dim(), 4);
    }

    #[test]
    fn test_standard_complex_structure_c3() {
        let j = ComplexStructure::standard(3);
        assert!(j.verify());
        assert_eq!(j.real_dim(), 6);
    }

    #[test]
    fn test_j_squared_is_neg_identity() {
        let j = ComplexStructure::standard(4);
        let n = j.real_dim();
        let j2 = &j.matrix * &j.matrix;
        let neg_id = -DMatrix::identity(n, n);
        assert!((j2 - neg_id).norm() < 1e-10);
    }

    #[test]
    fn test_apply_j() {
        let j = ComplexStructure::standard(1);
        let v = DVector::from_vec(vec![1.0, 0.0]);
        let jv = j.apply(&v);
        assert!((jv[0] - 0.0).abs() < 1e-10);
        assert!((jv[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_complex_vector() {
        let j = ComplexStructure::standard(2);
        let v = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let cv = j.to_complex_vector(&v);
        assert_eq!(cv.nrows(), 2);
        assert!((cv[0].re - 1.0).abs() < 1e-10);
        assert!((cv[0].im - 3.0).abs() < 1e-10);
        assert!((cv[1].re - 2.0).abs() < 1e-10);
        assert!((cv[1].im - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_standard_j_is_orthogonal() {
        let j = ComplexStructure::standard(3);
        assert!(j.is_orthogonal());
    }

    #[test]
    fn test_standard_j_compatible_with_identity_metric() {
        let j = ComplexStructure::standard(2);
        let g = DMatrix::identity(4, 4);
        assert!(j.is_compatible_with_metric(&g));
    }

    #[test]
    fn test_serialization() {
        let j = ComplexStructure::standard(2);
        let json = serde_json::to_string(&j).unwrap();
        let j2: ComplexStructure = serde_json::from_str(&json).unwrap();
        assert!(j2.verify());
    }
}
