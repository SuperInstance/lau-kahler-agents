//! Hodge decomposition on Kähler manifolds.
//!
//! On a Kähler manifold, the Hodge decomposition is finer than the real one:
//! H^k(M) = ⊕_{p+q=k} H^{p,q}(M)
//! with H^{p,q} = H̄^{q,p} (complex conjugation symmetry).

use serde::{Serialize, Deserialize};

/// Hodge diamond for a Kähler manifold of complex dimension n.
///
/// The (p,q) entry counts the dimension of H^{p,q}.
/// The diamond is (2n+1) rows high, with the middle row being the widest.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HodgeDiamond {
    /// Complex dimension n.
    pub complex_dim: usize,
    /// Hodge numbers h^{p,q} stored as diamond[p][q].
    pub diamond: Vec<Vec<usize>>,
}

impl HodgeDiamond {
    /// Create the Hodge diamond for CP^n.
    ///
    /// CP^n has h^{p,p} = 1 for 0 ≤ p ≤ n, and all other h^{p,q} = 0.
    pub fn cp_n(n: usize) -> Self {
        let mut diamond = vec![vec![0; n + 1]; n + 1];
        for (p, row) in diamond.iter_mut().enumerate() {
            if p <= n {
                row[p] = 1;
            }
        }
        Self { complex_dim: n, diamond }
    }

    /// Create the Hodge diamond for a Calabi-Yau 3-fold.
    ///
    /// Generic form with h^{1,1} and h^{2,1} as parameters.
    pub fn cy3(h11: usize, h21: usize) -> Self {
        // h^{0,0}=1, h^{1,0}=0, h^{2,0}=0, h^{3,0}=1
        // h^{1,1}=h11, h^{2,1}=h21
        // Serre duality: h^{p,q} = h^{3-p,3-q}
        // h^{2,2}=h^{1,1}=h11, h^{3,1}=h^{0,2}=0, h^{1,3}=h^{2,0}=0
        // h^{3,2}=h^{0,1}=0, h^{2,3}=h^{1,0}=0
        let mut diamond = vec![vec![0; 4]; 4];
        diamond[0][0] = 1;
        diamond[1][1] = h11;
        diamond[2][2] = h11;
        diamond[3][3] = 1;
        diamond[2][1] = h21;
        diamond[1][2] = h21;
        diamond[3][0] = 1;
        diamond[0][3] = 1;
        Self { complex_dim: 3, diamond }
    }

    /// Get h^{p,q}.
    pub fn get(&self, p: usize, q: usize) -> usize {
        if p < self.diamond.len() && q < self.diamond[p].len() {
            self.diamond[p][q]
        } else {
            0
        }
    }

    /// The k-th Betti number: b_k = Σ_{p+q=k} h^{p,q}.
    pub fn betti_number(&self, k: usize) -> usize {
        let mut sum = 0;
        for p in 0..=k.min(self.complex_dim) {
            let q = k - p;
            if q <= self.complex_dim {
                sum += self.get(p, q);
            }
        }
        sum
    }

    /// Verify the symmetry h^{p,q} = h^{q,p}.
    pub fn verify_symmetry(&self) -> bool {
        for p in 0..=self.complex_dim {
            for q in 0..=self.complex_dim {
                if self.get(p, q) != self.get(q, p) {
                    return false;
                }
            }
        }
        true
    }

    /// Verify the Serre duality: h^{p,q} = h^{n-p,n-q}.
    pub fn verify_serre_duality(&self) -> bool {
        let n = self.complex_dim;
        for p in 0..=n {
            for q in 0..=n {
                if self.get(p, q) != self.get(n - p, n - q) {
                    return false;
                }
            }
        }
        true
    }

    /// The Hodge decomposition: H^k = ⊕_{p+q=k} H^{p,q}.
    /// Returns the (p,q) decomposition of the k-th cohomology.
    pub fn decomposition(&self, k: usize) -> Vec<(usize, usize, usize)> {
        let mut result = vec![];
        for p in 0..=k.min(self.complex_dim) {
            let q = k - p;
            if q <= self.complex_dim {
                let h = self.get(p, q);
                if h > 0 {
                    result.push((p, q, h));
                }
            }
        }
        result
    }

    /// The Euler characteristic: χ = Σ_k (-1)^k b_k.
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi = 0i64;
        for k in 0..=2 * self.complex_dim {
            chi += if k % 2 == 0 { 1 } else { -1 } * self.betti_number(k) as i64;
        }
        chi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cp1_hodge_diamond() {
        let hd = HodgeDiamond::cp_n(1);
        assert_eq!(hd.get(0, 0), 1);
        assert_eq!(hd.get(1, 1), 1);
        assert_eq!(hd.get(0, 1), 0);
    }

    #[test]
    fn test_cp2_hodge_diamond() {
        let hd = HodgeDiamond::cp_n(2);
        assert_eq!(hd.get(0, 0), 1);
        assert_eq!(hd.get(1, 1), 1);
        assert_eq!(hd.get(2, 2), 1);
    }

    #[test]
    fn test_cp_n_symmetry() {
        let hd = HodgeDiamond::cp_n(3);
        assert!(hd.verify_symmetry());
    }

    #[test]
    fn test_cp_n_serre_duality() {
        let hd = HodgeDiamond::cp_n(3);
        assert!(hd.verify_serre_duality());
    }

    #[test]
    fn test_cp1_betti_numbers() {
        let hd = HodgeDiamond::cp_n(1);
        // For CP^1: h^{0,0}=1, h^{1,1}=1
        // b_0 = h^{0,0} = 1, b_1 = h^{1,0}+h^{0,1} = 0, b_2 = h^{1,1} = 1
        assert_eq!(hd.betti_number(0), 1);
        assert_eq!(hd.betti_number(1), 0);
        assert_eq!(hd.betti_number(2), 1);
    }

    #[test]
    fn test_cy3_hodge_diamond() {
        let hd = HodgeDiamond::cy3(5, 10);
        assert_eq!(hd.get(0, 0), 1);
        assert_eq!(hd.get(1, 1), 5);
        assert_eq!(hd.get(2, 1), 10);
        assert_eq!(hd.get(3, 3), 1);
    }

    #[test]
    fn test_cy3_symmetry() {
        let hd = HodgeDiamond::cy3(5, 10);
        assert!(hd.verify_symmetry());
    }

    #[test]
    fn test_cy3_euler() {
        let hd = HodgeDiamond::cy3(5, 10);
        // χ = 2(h^{1,1} - h^{2,1}) = 2(5-10) = -10
        assert_eq!(hd.euler_characteristic(), -10);
    }

    #[test]
    fn test_decomposition_cp2() {
        let hd = HodgeDiamond::cp_n(2);
        let d2 = hd.decomposition(2);
        // For CP^2: h^{2,0}=0, h^{1,1}=1, h^{0,2}=0
        assert!(d2.iter().any(|(p, q, h)| *p == 1 && *q == 1 && *h == 1));
    }

    #[test]
    fn test_serialization() {
        let hd = HodgeDiamond::cp_n(2);
        let json = serde_json::to_string(&hd).unwrap();
        let hd2: HodgeDiamond = serde_json::from_str(&json).unwrap();
        assert_eq!(hd2.get(1, 1), 1);
    }
}
