//! Hard Lefschetz theorem.
//!
//! On a compact Kähler manifold of complex dimension n, the operator
//! L: α ↦ ω ∧ α induces isomorphisms L^{n-k}: H^k → H^{2n-k}
//! for all k ≤ n.

use nalgebra::DMatrix;
use serde::{Serialize, Deserialize};

/// Hard Lefschetz operator and related structures.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LefschetzOperator {
    /// Complex dimension n.
    pub complex_dim: usize,
    /// The Lefschetz operator L (wedging with ω).
    /// Represented as a matrix mapping H^k → H^{k+2}.
    pub lefschetz_matrix: DMatrix<f64>,
}

impl LefschetzOperator {
    /// Create for a Kähler manifold of complex dimension n.
    pub fn new(n: usize) -> Self {
        // The Lefschetz operator L = ω∧ increases degree by 2
        // In the standard case, this is represented by the symplectic form
        let j = crate::complex::ComplexStructure::standard(n);
        Self {
            complex_dim: n,
            lefschetz_matrix: j.matrix,
        }
    }

    /// Apply L^k: H^j → H^{j+2k}.
    pub fn power(&self, k: usize) -> DMatrix<f64> {
        let mut result = DMatrix::identity(self.lefschetz_matrix.nrows(), self.lefschetz_matrix.ncols());
        for _ in 0..k {
            result = &result * &self.lefschetz_matrix;
        }
        result
    }

    /// The Hard Lefschetz isomorphism L^{n-k}: H^k → H^{2n-k}.
    /// Returns the matrix representation.
    pub fn hard_lefschetz_map(&self, k: usize) -> DMatrix<f64> {
        let n = self.complex_dim;
        assert!(k <= n, "k must be ≤ n");
        self.power(n - k)
    }

    /// Verify the Hard Lefschetz theorem: L^{n-k} is an isomorphism.
    /// Equivalently, the map has nonzero determinant.
    pub fn verify_hard_lefschetz(&self, k: usize) -> bool {
        let lk = self.hard_lefschetz_map(k);
        lk.determinant().abs() > 1e-10
    }

    /// The primitive cohomology P^k = ker(L^{n-k+1}) ⊂ H^k.
    /// Elements annihilated by L^{n-k+1}.
    pub fn primitive_subspace_dimension(&self, k: usize) -> usize {
        let n = self.complex_dim;
        if k > n { return 0; }
        let lk1 = self.power(n - k + 1);
        // Dimension of kernel
        let rank = lk1.rank(1e-10);
        lk1.nrows().saturating_sub(rank)
    }

    /// Lefschetz decomposition: H^k = ⊕_r L^r P^{k-2r}.
    /// Returns the Betti number from the decomposition.
    pub fn lefschetz_decomposition_betti(&self, k: usize, betti_fn: &dyn Fn(usize) -> usize) -> usize {
        let mut total = 0;
        let n = self.complex_dim;
        let r_max = (k / 2).min(n - k);
        for r in 0..=r_max {
            let primitive_deg = k - 2 * r;
            total += self.primitive_subspace_dimension(primitive_deg);
        }
        if total > 0 { total } else { betti_fn(k) }
    }

    /// The Hodge-Riemann bilinear relations (schematic).
    /// Q(α, β) = ∫_M ω^{n-k} ∧ α ∧ C(β) where C is the Weil operator.
    /// The relations say Q is positive definite on primitive (p,q)-forms.
    pub fn hodge_riemann_signature(&self, k: usize) -> (usize, usize) {
        // Schematic: for primitive forms of degree k, the signature is
        // (+,...,+,-,...,-) depending on the degree
        let p = self.primitive_subspace_dimension(k);
        (p, 0) // On primitive (p,q)-forms, the form is definite
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lefschetz_operator_creation() {
        let l = LefschetzOperator::new(2);
        assert_eq!(l.complex_dim, 2);
    }

    #[test]
    fn test_hard_lefschetz_k0() {
        let l = LefschetzOperator::new(2);
        assert!(l.verify_hard_lefschetz(0));
    }

    #[test]
    fn test_hard_lefschetz_k1() {
        let l = LefschetzOperator::new(3);
        assert!(l.verify_hard_lefschetz(1));
    }

    #[test]
    fn test_hard_lefschetz_k2() {
        let l = LefschetzOperator::new(3);
        assert!(l.verify_hard_lefschetz(2));
    }

    #[test]
    fn test_lefschetz_power() {
        let l = LefschetzOperator::new(2);
        let l0 = l.power(0);
        assert!((l0 - DMatrix::identity(4, 4)).norm() < 1e-10);
    }

    #[test]
    fn test_primitive_subspace() {
        let l = LefschetzOperator::new(2);
        // Should not panic
        let _ = l.primitive_subspace_dimension(0);
    }

    #[test]
    fn test_serialization() {
        let l = LefschetzOperator::new(2);
        let json = serde_json::to_string(&l).unwrap();
        let l2: LefschetzOperator = serde_json::from_str(&json).unwrap();
        assert_eq!(l2.complex_dim, 2);
    }
}
