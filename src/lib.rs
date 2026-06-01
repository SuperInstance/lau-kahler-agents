//! # lau-kahler-agents
//!
//! Kähler geometry for agents — where Riemannian, symplectic, and complex geometry unify.
//!
//! A Kähler manifold carries three simultaneous structures: a Riemannian metric `g`,
//! a symplectic form `ω`, and a complex structure `J`, satisfying the compatibility
//! condition `g(u, v) = ω(u, Jv)`.

pub mod complex;
pub mod hermitian;
pub mod kahler;
pub mod potential;
pub mod ricci;
pub mod calabi_yau;
pub mod hodge;
pub mod lefschetz;
pub mod identities;

pub mod prelude {
    pub use crate::complex::*;
    pub use crate::hermitian::*;
    pub use crate::kahler::*;
    pub use crate::potential::*;
    pub use crate::ricci::*;
    pub use crate::calabi_yau::*;
    pub use crate::hodge::*;
    pub use crate::lefschetz::*;
    pub use crate::identities::*;
}
