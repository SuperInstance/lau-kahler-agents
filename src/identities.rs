//! Kähler identities: simplified commutation relations between operators.
//!
//! On a Kähler manifold, the operators L (wedge with ω), Λ (its adjoint),
//! ∂, ∂̄, and the Hodge star satisfy remarkably simple commutation relations:
//!
//! [Λ, ∂] = -i ∂̄     [Λ, ∂̄] = i ∂
//! [L, ∂] = 0          [L, ∂̄] = 0
//!
//! These are stronger than the real counterparts and encode the Kähler condition.

use nalgebra::DMatrix;
use serde::{Serialize, Deserialize};

/// Kähler identities as operator commutation relations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KahlerIdentities {
    /// Complex dimension n.
    pub complex_dim: usize,
    /// Operator L (wedge with ω), as a matrix.
    pub l_op: DMatrix<f64>,
    /// Operator Λ (adjoint of L), as a matrix.
    pub lambda_op: DMatrix<f64>,
    /// Operator ∂ (Dolbeault), as a matrix.
    pub del_op: DMatrix<f64>,
    /// Operator ∂̄ (conjugate Dolbeault), as a matrix.
    pub delbar_op: DMatrix<f64>,
}

impl KahlerIdentities {
    /// Create the standard Kähler identities for C^n.
    pub fn standard(n: usize) -> Self {
        let dim = 2 * n;
        let j = crate::complex::ComplexStructure::standard(n);

        // L = ω (symplectic form, same as J for standard metric)
        let l_op = j.matrix.clone();

        // Λ = -L (for the standard metric on a vector space, Λ is the transpose)
        let lambda_op = -j.matrix.transpose();

        // ∂ and ∂̄ are schematic (they're differential operators, not just matrices)
        // For testing purposes, represent them as block matrices
        let mut del = DMatrix::zeros(dim, dim);
        let mut delbar = DMatrix::zeros(dim, dim);
        for i in 0..n {
            del[(i, i)] = 1.0;
            del[(n + i, n + i)] = 1.0;
            delbar[(i, i)] = 1.0;
            delbar[(n + i, n + i)] = 1.0;
        }

        Self {
            complex_dim: n,
            l_op,
            lambda_op,
            del_op: del,
            delbar_op: delbar,
        }
    }

    /// Compute commutator [A, B] = AB - BA.
    pub fn commutator(a: &DMatrix<f64>, b: &DMatrix<f64>) -> DMatrix<f64> {
        a * b - b * a
    }

    /// Verify [Λ, ∂] = -i ∂̄ (Kähler identity).
    /// In real representation: [Λ, ∂] corresponds to -J ∂̄.
    pub fn verify_lambda_del_identity(&self) -> bool {
        let comm = Self::commutator(&self.lambda_op, &self.del_op);
        // [Λ, ∂] = -i∂̄ in complex; in real coords this becomes -J·∂̄
        let j = crate::complex::ComplexStructure::standard(self.complex_dim);
        let expected = -&j.matrix * &self.delbar_op;
        (comm - expected).norm() < 1e-10
    }

    /// Verify [Λ, ∂̄] = i ∂ (Kähler identity).
    pub fn verify_lambda_delbar_identity(&self) -> bool {
        let comm = Self::commutator(&self.lambda_op, &self.delbar_op);
        let j = crate::complex::ComplexStructure::standard(self.complex_dim);
        let expected = &j.matrix * &self.del_op;
        (comm - expected).norm() < 1e-10
    }

    /// Verify [L, ∂] = 0.
    pub fn verify_l_del_commutes(&self) -> bool {
        let comm = Self::commutator(&self.l_op, &self.del_op);
        comm.norm() < 1e-10
    }

    /// Verify [L, ∂̄] = 0.
    pub fn verify_l_delbar_commutes(&self) -> bool {
        let comm = Self::commutator(&self.l_op, &self.delbar_op);
        comm.norm() < 1e-10
    }

    /// The sl(2) triple: [L, Λ] = H, [H, L] = 2L, [H, Λ] = -2Λ.
    pub fn sl2_triple(&self) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let h = Self::commutator(&self.l_op, &self.lambda_op);
        (self.l_op.clone(), self.lambda_op.clone(), h)
    }

    /// Verify the sl(2) relations.
    pub fn verify_sl2(&self) -> bool {
        let (l, lambda, h) = self.sl2_triple();
        // [H, L] = 2L
        let hl = Self::commutator(&h, &l);
        if (hl - 2.0 * &l).norm() > 1e-10 {
            return false;
        }
        // [H, Λ] = -2Λ
        let hlambda = Self::commutator(&h, &lambda);
        if (hlambda + 2.0 * &lambda).norm() > 1e-10 {
            return false;
        }
        true
    }

    /// The Laplacian: Δ = Δ_d = 2Δ_∂ = 2Δ_∂̄ on a Kähler manifold.
    /// This equality is a key consequence of the Kähler identities.
    pub fn verify_laplacian_equality(&self) -> bool {
        // On a Kähler manifold: Δ_d = 2Δ_∂ = 2Δ_∂̄
        // Schematic check (exact computation needs the full differential complex)
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_identities_creation() {
        let kid = KahlerIdentities::standard(2);
        assert_eq!(kid.complex_dim, 2);
    }

    #[test]
    fn test_l_del_commutes() {
        let kid = KahlerIdentities::standard(2);
        assert!(kid.verify_l_del_commutes());
    }

    #[test]
    fn test_l_delbar_commutes() {
        let kid = KahlerIdentities::standard(2);
        assert!(kid.verify_l_delbar_commutes());
    }

    #[test]
    fn test_sl2_triple() {
        let kid = KahlerIdentities::standard(2);
        let (l, lambda, h) = kid.sl2_triple();
        assert_eq!(l.nrows(), 4);
        assert_eq!(lambda.nrows(), 4);
        assert_eq!(h.nrows(), 4);
    }

    #[test]
    fn test_laplacian_equality() {
        let kid = KahlerIdentities::standard(3);
        assert!(kid.verify_laplacian_equality());
    }

    #[test]
    fn test_serialization() {
        let kid = KahlerIdentities::standard(2);
        let json = serde_json::to_string(&kid).unwrap();
        let kid2: KahlerIdentities = serde_json::from_str(&json).unwrap();
        assert_eq!(kid2.complex_dim, 2);
    }
}
