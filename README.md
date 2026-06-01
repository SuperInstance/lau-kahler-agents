# lau-kahler-agents

**Kähler geometry for agents — where Riemannian, symplectic, and complex geometry unify.**

A Rust library implementing the core structures of Kähler geometry: complex structures (`J`), Hermitian forms, symplectic forms (`ω`), Kähler potentials, Ricci curvature, Calabi–Yau manifolds, Hodge diamonds, Hard Lefschetz, and the Kähler identities. Every structure is verified computationally and serializable via serde.

[![74 tests passing](https://img.shields.io/badge/tests-74%20passing-brightgreen)]()

---

## Table of Contents

- [What This Does](#what-this-does)
- [Key Idea](#key-idea)
- [Install](#install)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
  - [ComplexStructure](#complexstructure)
  - [HermitianForm](#hermitianform)
  - [KahlerStructure](#kahlerstructure)
  - [KahlerPotential](#kahlerpotential)
  - [RicciForm](#ricciform)
  - [CalabiYau](#calabiyau)
  - [HodgeDiamond](#hodgediamond)
  - [LefschetzOperator](#lefschetzoperator)
  - [KahlerIdentities](#kahleridentities)
- [How It Works](#how-it-works)
- [The Math](#the-math)
  - [Complex Structures](#complex-structures)
  - [Hermitian Forms and the Triple Structure](#hermitian-forms-and-the-triple-structure)
  - [Kähler Potentials](#kähler-potentials)
  - [Ricci Curvature and Calabi–Yau](#ricci-curvature-and-calabi-yau)
  - [Hodge Theory on Kähler Manifolds](#hodge-theory-on-kähler-manifolds)
  - [Hard Lefschetz and sl(2) Representations](#hard-lefschetz-and-sl2-representations)
  - [Kähler Identities](#kähler-identities-1)
- [License](#license)

---

## What This Does

This library provides computational tools for Kähler geometry — the remarkable class of manifolds where three fundamental geometric structures coexist compatibly:

1. **Complex structure** `J` — turns `ℝ²ⁿ` into `ℂⁿ`, satisfying `J² = -Id`
2. **Riemannian metric** `g` — measures lengths and angles, symmetric and positive-definite
3. **Symplectic form** `ω` — closed, non-degenerate 2-form for measuring areas

The Kähler compatibility condition `g(u, v) = ω(u, Jv)` makes these three structures determine each other, leading to extraordinary mathematical consequences: Hodge decomposition, Hard Lefschetz, the Kähler identities, and the existence of Ricci-flat Calabi–Yau metrics.

Modules:

- **complex** — Complex structure `J` with `J² = -Id`, orthogonality checks, metric compatibility
- **hermitian** — Hermitian forms `h`, extraction of the Riemannian metric and symplectic form
- **kahler** — Full Kähler structure combining all three, with verification of compatibility
- **potential** — Kähler potentials `Φ` with `ω = i∂∂̄Φ`, flat and Fubini–Study potentials
- **ricci** — Ricci curvature form, first Chern class, scalar curvature
- **calabi_yau** — Calabi–Yau manifolds (Ricci-flat Kähler), holomorphic volume forms, Hodge numbers
- **hodge** — Hodge diamonds, Betti numbers from Hodge decomposition, Serre duality verification
- **lefschetz** — Hard Lefschetz operator, primitive cohomology, Lefschetz decomposition
- **identities** — Kähler identities as operator commutation relations, sl(2) representation

---

## Key Idea

A **Kähler manifold** is a manifold that is simultaneously complex, Riemannian, and symplectic — and these three structures are compatible. This compatibility is so rigid that it forces deep consequences:

- The three operators `∂`, `∂̄`, and the Lefschetz operator `L = ω∧` satisfy the **Kähler identities**: simple commutation relations that are stronger than on a generic Riemannian or complex manifold.
- Cohomology splits: `Hᵏ(M) = ⊕_{p+q=k} H^{p,q}(M)` (**Hodge decomposition**), with symmetry `h^{p,q} = h^{q,p}`.
- The **Hard Lefschetz theorem** gives isomorphisms `Lⁿ⁻ᵏ: Hᵏ → H²ⁿ⁻ᵏ`.
- The **Kähler class** (cohomology class of `ω`) controls the geometry through the potential `Φ`.
- **Calabi–Yau manifolds** (Ricci-flat Kähler) exist whenever the first Chern class vanishes (Yau's theorem), and their Hodge numbers determine topological invariants like the Euler characteristic.

This library makes all of these structures computable, verifiable, and serializable.

---

## Install

```toml
[dependencies]
lau-kahler-agents = "0.1.0"
```

Or:

```bash
cargo add lau-kahler-agents
```

Requires Rust 2021+. Dependencies: `nalgebra` (with serde), `serde`, `serde_json`.

---

## Quick Start

### Build and verify a Kähler structure

```rust
use lau_kahler_agents::prelude::*;

// Standard Kähler structure on C²
let k = KahlerStructure::standard(2);

// Verify all compatibility conditions
assert!(k.verify_all()); // J²=-Id, g symmetric, ω antisymmetric, g=ω∘J, positive definite

// Evaluate the metric and symplectic form on vectors
let u = nalgebra::DVector::from_vec(vec![1.0, 0.0, 0.0, 0.0]);
let v = nalgebra::DVector::from_vec(vec![0.0, 1.0, 0.0, 0.0]);
println!("g(u,v) = {}", k.g(&u, &v));
println!("ω(u,v) = {}", k.omega(&u, &v));
```

### Kähler potential and Fubini–Study metric

```rust
use lau_kahler_agents::prelude::*;
use nalgebra::{DVector, Complex};

// Flat potential: Φ = |z|²
let phi_flat = KahlerPotential::flat(2);
assert!(phi_flat.is_positive_definite());

// Fubini–Study potential on CP² at a point
let z = DVector::from_vec(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 1.0)]);
let phi_fs = KahlerPotential::fubini_study(2, &z);
assert!(phi_fs.is_positive_definite());
```

### Hodge diamond for a Calabi–Yau 3-fold

```rust
use lau_kahler_agents::prelude::*;

// Quintic CY3: h^{1,1} = 1, h^{2,1} = 101
let hd = HodgeDiamond::cy3(1, 101);
println!("h^{1,1} = {}", hd.get(1, 1)); // 1
println!("h^{2,1} = {}", hd.get(2, 1)); // 101
println!("χ = {}", hd.euler_characteristic()); // -200
assert!(hd.verify_symmetry());     // h^{p,q} = h^{q,p}
assert!(hd.verify_serre_duality()); // h^{p,q} = h^{n-p,n-q}
```

### Calabi–Yau and Ricci-flatness

```rust
use lau_kahler_agents::prelude::*;

let cy = CalabiYau::flat(3);
assert!(cy.verify_ricci_flat());
println!("SU(n) holonomy dim = {}", cy.su_n_holonomy_dimension()); // 8 for SU(3)
println!("χ(CY3) = {}", CalabiYau::euler_characteristic_cy3(5, 10)); // -10
```

### Kähler identities and sl(2) representation

```rust
use lau_kahler_agents::prelude::*;

let kid = KahlerIdentities::standard(2);
assert!(kid.verify_l_del_commutes());     // [L, ∂] = 0
assert!(kid.verify_l_delbar_commutes());   // [L, ∂̄] = 0
assert!(kid.verify_laplacian_equality());  // Δ_d = 2Δ_∂ = 2Δ_∂̄
```

---

## API Reference

### ComplexStructure

```rust
pub struct ComplexStructure {
    pub matrix: DMatrix<f64>,
    pub complex_dim: usize,
}
```

A complex structure `J` on `ℝ²ⁿ` with `J² = -Id`.

| Method | Description |
|--------|-------------|
| `standard(n)` | Standard `J = [[0,-I],[I,0]]` on `ℂⁿ` |
| `verify()` | Check `J² = -Id` |
| `apply(v)` | Apply `J` to a vector |
| `to_complex_vector(v)` | Interpret `ℝ²ⁿ` vector as `ℂⁿ` |
| `pq_components(v)` | Decompose into (1,0) and (0,1) parts |
| `has_pure_imaginary_eigenvalues()` | Eigenvalues of `J` are `±i` |
| `is_orthogonal()` | `JᵀJ = I` |
| `is_compatible_with_metric(g)` | `g(Ju, Jv) = g(u,v)` |
| `real_dim()` | `2n` |

### HermitianForm

```rust
pub struct HermitianForm {
    pub matrix: DMatrix<Complex<f64>>,
    pub dim: usize,
}
```

An `n×n` Hermitian matrix `H = H*` representing the sesquilinear form `h(u,v) = u*Hv`.

| Method | Description |
|--------|-------------|
| `new(matrix)` | Create from Hermitian matrix |
| `standard(n)` | Identity matrix (standard inner product on `ℂⁿ`) |
| `is_hermitian()` | Verify `H = H*` |
| `evaluate(u, v)` | `u*Hv` |
| `riemannian_metric()` | `g = Re(h)` as `2n×2n` real matrix |
| `kahler_form()` | `ω = Im(h)` as `2n×2n` real matrix |
| `is_positive_definite()` | All eigenvalues of `g` non-negative |
| `verify_kahler_compatibility(j)` | `g(u,v) = ω(u, Jv)` |

### KahlerStructure

```rust
pub struct KahlerStructure {
    pub hermitian: HermitianForm,
    pub j: ComplexStructure,
    pub complex_dim: usize,
}
```

The full Kähler structure: metric `g`, symplectic form `ω`, complex structure `J`.

| Method | Description |
|--------|-------------|
| `standard(n)` | Standard Kähler on `ℂⁿ` |
| `metric()` | Riemannian metric `g` |
| `symplectic_form()` | Kähler form `ω` |
| `complex_structure()` | Reference to `J` |
| `real_dim()` | `2n` |
| `verify_compatibility()` | `g = ω ∘ J` |
| `verify_complex_structure()` | `J² = -Id` |
| `verify_symplectic_antisymmetry()` | `ωᵀ = -ω` |
| `verify_metric_symmetry()` | `gᵀ = g` |
| `verify_positive_definite()` | `g` is positive definite |
| `verify_all()` | All five conditions |
| `g(u, v)` | Evaluate metric on vectors |
| `omega(u, v)` | Evaluate symplectic form |
| `verify_compatibility_vectors(u, v)` | `g(u,v) = ω(u,Jv)` for specific vectors |
| `is_closed()` | `dω = 0` (automatic on vector spaces) |
| `volume_form()` | `ωⁿ/n!` |

### KahlerPotential

```rust
pub struct KahlerPotential {
    pub dim: usize,
    pub hessian: DMatrix<f64>,
}
```

A Kähler potential `Φ` with `ω = i∂∂̄Φ`, stored as the mixed Hessian `∂²Φ/∂zⁱ∂z̄ʲ`.

| Method | Description |
|--------|-------------|
| `new(hessian)` | From mixed Hessian |
| `flat(n)` | `Φ = |z|²` (standard metric) |
| `fubini_study(n, z)` | `Φ = log(1 + |z|²)` on `ℂℂPⁿ` |
| `evaluate(z)` | `Φ(z) = |z|²` for flat |
| `kahler_form()` | Recover `ω` from Hessian |
| `metric()` | Associated metric `g` |
| `is_positive_definite()` | Hessian positive definite |
| `volume_element()` | `det(h_{i,j̄})` |

### RicciForm

```rust
pub struct RicciForm {
    pub dim: usize,
    pub ricci_tensor: DMatrix<f64>,
}
```

The Ricci form `ρ = iR_{i,j̄} dzⁱ ∧ dz̄ʲ`, measuring curvature of the canonical bundle.

| Method | Description |
|--------|-------------|
| `from_kahler(k)` | Compute from Kähler structure (flat → zero) |
| `from_potential(phi)` | Compute from potential's Hessian |
| `zero(n)` | Zero Ricci form |
| `is_ricci_flat()` | `ρ = 0` |
| `as_real_form()` | `2n×2n` real matrix |
| `first_chern_class()` | `c₁ = ρ/(2π)` |
| `scalar_curvature(g)` | Trace of Ricci w.r.t. metric |

### CalabiYau

```rust
pub struct CalabiYau {
    pub kahler: KahlerStructure,
    pub ricci: RicciForm,
}
```

A Calabi–Yau manifold: Ricci-flat Kähler.

| Method | Description |
|--------|-------------|
| `flat(n)` | Flat Calabi–Yau on `ℂⁿ` |
| `verify_ricci_flat()` | Check `ρ = 0` |
| `holomorphic_volume_form()` | `Ω = dz¹ ∧ ... ∧ dzⁿ` |
| `su_n_holonomy_dimension()` | `n² - 1` (dimension of SU(n)) |
| `hodge_numbers_cy3(h11, h21)` | Hodge numbers for a CY3 |
| `euler_characteristic_cy3(h11, h21)` | `χ = 2(h^{1,1} - h^{2,1})` |
| `verify_stability()` | Bogomolov stability |
| `volume()` | Volume from potential |

### HodgeDiamond

```rust
pub struct HodgeDiamond {
    pub complex_dim: usize,
    pub diamond: Vec<Vec<usize>>,
}
```

The Hodge diamond `h^{p,q}` for a Kähler manifold.

| Method | Description |
|--------|-------------|
| `cp_n(n)` | Hodge diamond for `ℂℂPⁿ` |
| `cy3(h11, h21)` | Hodge diamond for a CY3-fold |
| `get(p, q)` | Get `h^{p,q}` |
| `betti_number(k)` | `βₖ = Σ_{p+q=k} h^{p,q}` |
| `verify_symmetry()` | `h^{p,q} = h^{q,p}` |
| `verify_serre_duality()` | `h^{p,q} = h^{n-p,n-q}` |
| `decomposition(k)` | `(p,q,h)` triples summing to `Hᵏ` |
| `euler_characteristic()` | `χ = Σ(-1)ᵏβₖ` |

### LefschetzOperator

```rust
pub struct LefschetzOperator {
    pub complex_dim: usize,
    pub lefschetz_matrix: DMatrix<f64>,
}
```

The Hard Lefschetz operator `L = ω∧`.

| Method | Description |
|--------|-------------|
| `new(n)` | Create for dimension `n` |
| `power(k)` | `Lᵏ` as matrix |
| `hard_lefschetz_map(k)` | `Lⁿ⁻ᵏ: Hᵏ → H²ⁿ⁻ᵏ` |
| `verify_hard_lefschetz(k)` | Check `Lⁿ⁻ᵏ` is an isomorphism |
| `primitive_subspace_dimension(k)` | `dim Pᵏ = ker(Lⁿ⁻ᵏ⁺¹)` |
| `lefschetz_decomposition_betti(k, betti_fn)` | Betti from Lefschetz decomposition |
| `hodge_riemann_signature(k)` | Hodge–Riemann bilinear form signature |

### KahlerIdentities

```rust
pub struct KahlerIdentities {
    pub complex_dim: usize,
    pub l_op: DMatrix<f64>,
    pub lambda_op: DMatrix<f64>,
    pub del_op: DMatrix<f64>,
    pub delbar_op: DMatrix<f64>,
}
```

The Kähler identity operators: `L`, `Λ`, `∂`, `∂̄`.

| Method | Description |
|--------|-------------|
| `standard(n)` | Standard identities for `ℂⁿ` |
| `commutator(a, b)` | `[A, B] = AB - BA` |
| `verify_lambda_del_identity()` | `[Λ, ∂] = -i∂̄` |
| `verify_lambda_delbar_identity()` | `[Λ, ∂̄] = i∂` |
| `verify_l_del_commutes()` | `[L, ∂] = 0` |
| `verify_l_delbar_commutes()` | `[L, ∂̄] = 0` |
| `sl2_triple()` | `(L, Λ, H)` where `H = [L, Λ]` |
| `verify_sl2()` | `[H,L] = 2L`, `[H,Λ] = -2Λ` |
| `verify_laplacian_equality()` | `Δ_d = 2Δ_∂ = 2Δ_∂̄` |

---

## How It Works

The library builds up Kähler geometry from the ground up:

1. **Start with a complex structure** `J` — a matrix satisfying `J² = -Id`. This turns real `2n`-dimensional space into complex `n`-dimensional space.

2. **Add a Hermitian form** `h` — a sesquilinear inner product on `ℂⁿ`. Its real part gives the Riemannian metric `g`, its imaginary part gives the symplectic form `ω`.

3. **Verify Kähler compatibility** — check that `g(u,v) = ω(u, Jv)` and `g(Ju, Jv) = g(u,v)`. These conditions make the three structures mutually determined.

4. **Express through a potential** — the Kähler form comes from a potential function: `ω = i∂∂̄Φ`. The mixed Hessian `∂²Φ/∂zⁱ∂z̄ʲ` is positive-definite and encodes the metric.

5. **Compute curvature** — the Ricci form is `ρ = -i∂∂̄ log(det g)`. If `ρ = 0`, the manifold is **Calabi–Yau**.

6. **Build the Hodge diamond** — on a Kähler manifold, cohomology splits: `Hᵏ = ⊕_{p+q=k} H^{p,q}` with symmetry `h^{p,q} = h^{q,p}` and Serre duality `h^{p,q} = h^{n-p,n-q}`.

7. **Apply Hard Lefschetz** — the operator `L = ω∧` gives isomorphisms `Lⁿ⁻ᵏ: Hᵏ ≅ H²ⁿ⁻ᵏ`. The primitive cohomology `Pᵏ = ker(Lⁿ⁻ᵏ⁺¹)` generates everything via the Lefschetz decomposition.

8. **Verify the Kähler identities** — the operators `L, Λ, ∂, ∂̄` satisfy `sl(2)` commutation relations: `[Λ,∂] = -i∂̄`, `[Λ,∂̄] = i∂`, `[L,∂] = 0`, `[L,∂̄] = 0`. These imply the Laplacian equality `Δ_d = 2Δ_∂ = 2Δ_∂̄`.

---

## The Math

### Complex Structures

A **complex structure** on a `2n`-dimensional real vector space is a linear map `J: V → V` with `J² = -Id`. This makes `V` into an `n`-dimensional complex vector space by defining `i·v = Jv`.

The eigenvalues of `J` are `±i`, and `J` is orthogonal when `JᵀJ = I`. A metric `g` is **compatible** with `J` if `g(Ju, Jv) = g(u,v)`.

### Hermitian Forms and the Triple Structure

A **Hermitian form** on `ℂⁿ` is a sesquilinear map `h(u,v) = u*Hv` with `H = H*` (Hermitian matrix). It splits into:

- **Riemannian metric**: `g(u,v) = Re h(u,v)` — symmetric, positive-definite
- **Symplectic form**: `ω(u,v) = Im h(u,v)` — antisymmetric, non-degenerate

The **Kähler compatibility** condition is:

```
g(u, v) = ω(u, Jv)
```

This means knowing any two of `(g, ω, J)` determines the third. The whole geometry is encoded in the single Hermitian form `h`.

### Kähler Potentials

On a Kähler manifold, the Kähler form is locally exact:

```
ω = i ∂∂̄ Φ
```

where `Φ` is the **Kähler potential** (a real-valued function). In coordinates:

```
ω_{i,j̄} = ∂²Φ / ∂zⁱ ∂z̄ʲ
```

**Flat metric**: `Φ = |z|² = Σ|zᵢ|²` gives `ω_{i,j̄} = δ_{ij}`.

**Fubini–Study metric** on `ℂℂPⁿ`: `Φ = log(1 + |z|²)`. At the origin, this reduces to the flat metric. The Fubini–Study metric has positive Ricci curvature and is the canonical metric on projective space.

### Ricci Curvature and Calabi–Yau

The **Ricci form** measures the curvature of the canonical line bundle:

```
ρ_{i,j̄} = -∂² log(det g) / ∂zⁱ ∂z̄ʲ
```

The **first Chern class** is `c₁ = [ρ/(2π)]` ∈ `H²(M, ℝ)`.

**Calabi–Yau manifolds** are Kähler manifolds with `c₁ = 0`. By **Yau's theorem** (proving the Calabi conjecture), every Kähler manifold with vanishing first Chern class admits a Ricci-flat Kähler metric. This is equivalent to the holonomy being contained in `SU(n)`.

For a Calabi–Yau 3-fold, the topological invariants are the **Hodge numbers** `h^{1,1}` and `h^{2,1}`, and the Euler characteristic is:

```
χ = 2(h^{1,1} - h^{2,1})
```

### Hodge Theory on Kähler Manifolds

On a Kähler manifold, cohomology has a **Hodge decomposition**:

```
Hᵏ(M, ℂ) = ⊕_{p+q=k} H^{p,q}(M)
```

with the symmetry `H^{p,q} = H̄^{q,p}` (complex conjugation).

The **Hodge diamond** arranges `h^{p,q} = dim H^{p,q}` in a diamond pattern. It satisfies:

- **Symmetry**: `h^{p,q} = h^{q,p}` (complex conjugation)
- **Serre duality**: `h^{p,q} = h^{n-p,n-q}`
- **Kähler package**: the entire diamond is determined by the primitive cohomology via the Lefschetz decomposition

Betti numbers are recovered: `βₖ = Σ_{p+q=k} h^{p,q}`.

Examples:
- **ℂℂP¹** (Riemann sphere): Diamond `[[1],[0,0],[1]]`, betti `[1,0,1]`, χ = 2.
- **ℂℂPⁿ**: `h^{p,p} = 1` for `0 ≤ p ≤ n`, all others zero.
- **Calabi–Yau 3-fold**: `h^{0,0}=1, h^{1,1}, h^{2,1}, h^{3,0}=1` and their reflections.

### Hard Lefschetz and sl(2) Representations

The **Lefschetz operator** `L = ω∧` maps `Hᵏ → Hᵏ⁺²` by wedging with the Kähler form.

**Hard Lefschetz Theorem**: For `k ≤ n`, the map `Lⁿ⁻ᵏ: Hᵏ → H²ⁿ⁻ᵏ` is an isomorphism.

This gives the **Lefschetz decomposition**:

```
Hᵏ = ⊕ᵣ Lʳ Pᵏ⁻²ʳ
```

where `Pʲ = ker(Lⁿ⁻ʲ⁺¹) ⊂ Hʲ` is the **primitive cohomology**.

The operators `{L, Λ, H}` form an **sl(2) representation**:

```
[L, Λ] = H,  [H, L] = 2L,  [H, Λ] = -2Λ
```

where `Λ` is the adjoint of `L` and `H` acts by `(n-k)` on `Hᵏ`.

### Kähler Identities

The Kähler condition implies commutation relations between the Lefschetz operators and the Dolbeault operators:

```
[Λ, ∂] = -i∂̄       [Λ, ∂̄] = i∂
[L, ∂]  = 0         [L, ∂̄]  = 0
```

These are **stronger** than the identities on a generic Riemannian or complex manifold. Their key consequence is the **Laplacian equality**:

```
Δ_d = 2Δ_∂ = 2Δ_∂̄
```

This means the de Rham Laplacian, the Dolbeault Laplacian, and the conjugate Dolbeault Laplacian all agree (up to a factor of 2). This is why harmonic forms on a Kähler manifold decompose by type `(p,q)` — the foundation of Hodge theory.

---

## License

MIT
