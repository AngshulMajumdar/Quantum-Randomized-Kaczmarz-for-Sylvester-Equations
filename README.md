# QRK-Sylv — Rust Solvers and Baseline Comparisons

**Structure-Preserving Quantum Randomized Kaczmarz for General Sylvester Equations: Convergence and End-to-End Resource Analysis**

This repository contains Rust implementations of QRK-Sylv, its explicit-Kronecker comparator, four named matrix-equation baseline **mathematical kernels**, a classical direct reference, complex finite-shot measurement emulation, the exact small-system Julia clock-dilation state vector, and the paper's instance generators. All numerical code is Rust, with no external runtime crates.

## Solvers

| CLI method | Rust entry point | Implementation | Domain |
|---|---|---|---|
| QRK-Sylv | `solvers::qrk` | Structured norm-proportional Kaczmarz with frozen sampling sequence, active-branch normalization and success probability | Any uniquely solvable complex Sylvester equation; caller supplies a valid solution-norm bound |
| Generic Kronecker RK | `solvers::generic_kaczmarz` | Explicit $K=I\otimes A+B^T\otimes I$ Kaczmarz on the **same** selected rows | General, small dense matrices |
| Somma et al. | `solvers::somma` | Numerical spectral solution-map emulator: complex diagonal / symmetric functional calculus or positive-real Laplace quadrature | Implemented diagonal, real-symmetric and positive-real regimes |
| Wang–Liu | `solvers::wang_liu` | Classical matrix-sign embedding and Newton sign iterations | Certified same-side strict half-plane inputs |
| Clayton et al. | `solvers::clayton` | Deterministic Lyapunov evolution-integral emulator | $B=A^\dagger$, Hurwitz $A$, positive-semidefinite $-C$ |
| Benedetti et al. | `solvers::benedetti` | Expected mixed-state semigroup channel emulator | Same Lyapunov domain plus normal $A$ |
| Classical direct | `Sylvester::direct_solve` | Pivoted complex Gaussian elimination on small Kronecker systems | Small square uniquely solvable systems |

**Scope of the baselines:** The four named baselines implement numerical evaluations of the mathematical solution maps underlying the cited quantum algorithms. They are **not gate-by-gate implementations** of Somma's block encoding, Wang–Liu's quantum shifted-resolvent multiplexing, Clayton's quantum policy gradient, or Benedetti's measured/restarted quantum process. Consequently the CPU timings describe the Rust emulators, not quantum circuit runtime, quantum speedups, or the exact emulator timings reported in the manuscript. An `N/A` means that the configured emulator declines an input; its applicability checks may be more conservative than the corresponding published algorithm's full hypotheses.

## Run

Install a recent Rust stable toolchain, then from this repository:

```bash
cargo run --release -- --family all --size 4 --instances 30 --sweeps 50 --shots 8192 --seed 123
```

For the manuscript's four application families, run:

```bash
cargo run --release -- --family screened --instances 30
cargo run --release -- --family convection --instances 30
cargo run --release -- --family wave --instances 30
cargo run --release -- --family control --instances 30
```

Two additional diagnostic families are `complex` and `normal-lyapunov`. The latter lets the Benedetti emulator run on inputs satisfying its normality assumption; it is not one of the manuscript's four main families.

```bash
cargo run --release -- --family normal-lyapunov --size 4 --instances 10
cargo run --release -- --family complex --size 2 --instances 10
cargo test
```

The CLI prints a Markdown table including relative Frobenius error, phase-sensitive finite-shot error, and emulator CPU milliseconds. QRK samples are drawn once and frozen for all measurement emulations; the generic Kronecker RK gets the exact same row sequence. The finite-shot module uses $8192$ shots by default, visibility $0.975$, and symmetric readout error $0.015$ with a simple explicit mitigation model. The CLI's classical reference is only used to evaluate errors and supply the valid $R$ bound in this diagnostic harness; a deployed quantum algorithm must obtain $R$ independently. The table of CPU times is environment dependent.

`clock::simulate` separately propagates the explicit complex Julia dilation over the active and distinct failure-clock sectors and checks state normalization and active-branch agreement. It is intentionally small-system-only and does not replace circuit synthesis or quantum-hardware execution.

## Source layout

- `src/complex.rs` — complex doubles and Hermitian inner products.
- `src/matrix.rs` — Sylvester geometry, implicit actions, exact row normals, direct small-system reference.
- `src/rng.rs` — deterministic pseudorandom sampling and Gaussian values.
- `src/sampling.rs` — alias/rejection norm-proportional sampling and frozen sequences.
- `src/linalg.rs` — complex matrix inverse, exponential, Hermitian tests, and symmetric eigensolver.
- `src/solvers.rs` — QRK, all four named numerical baselines, generic RK, finite-shot estimates.
- `src/clock.rs` — full small-system failure-clock state vectors.
- `src/fixtures.rs` — paper's four matrices and two extra diagnostic families.
- `src/main.rs` — command-line table output.
- `tests/algorithms.rs` — mathematical invariant and domain tests.

## Tables from the manuscript

These are the **reported manuscript tables**, not numerical output claimed for this Rust implementation. `N/A` entries and their reasons have been retained. Values and stated uncertainties have not been changed.

### Table 1. Structural applicability

| Method | Independent A, B | Nonnormal or defective | Additional condition | Full class |
|---|:---:|:---:|---|:---:|
| Somma et al. | ✓ | R / × | Diagonalizability + listed efficient regimes | × |
| Wang–Liu | ✓ | ✓ | Half-plane separation + FoV/resolvent | × |
| Clayton et al. | × | ✓ | Lyapunov form, Hurwitz stability, positive forcing | × |
| Benedetti et al. | × | × | Lyapunov form, normality, stability, PSD forcing | × |
| **QRK-Sylv** | ✓ | ✓ | None | ✓ |

### Table 2. Screened diffusion

| $n$ | Method | Ideal error | Finite-shot error | Emulator ms |
|---|---|---|---|---|
| 2 | **QRK-Sylv** | $(5.83\pm17.90)\times10^{-9}$ | $0.0310\pm0.0081$ | $6.06\pm0.20$ |
| 2 | Somma et al. | $(2.80\pm1.97)\times10^{-4}$ | $0.0329\pm0.0073$ | $16.77\pm9.52$ |
| 2 | Wang–Liu | $(1.07\pm0.31)\times10^{-3}$ | $0.0318\pm0.0071$ | $20.13\pm8.65$ |
| 2 | Clayton et al. | N/A | N/A | -- |
| 2 | Benedetti et al. | N/A | N/A | -- |
| 4 | **QRK-Sylv** | $(3.70\pm3.65)\times10^{-4}$ | $0.0517\pm0.0075$ | $27.06\pm4.75$ |
| 4 | Somma et al. | $(1.67\pm1.07)\times10^{-3}$ | $0.0445\pm0.0074$ | $28.88\pm17.86$ |
| 4 | Wang–Liu | $(2.01\pm0.45)\times10^{-3}$ | $0.0448\pm0.0069$ | $19.74\pm9.34$ |
| 4 | Clayton et al. | N/A | N/A | -- |
| 4 | Benedetti et al. | N/A | N/A | -- |
| 8 | **QRK-Sylv** | $(3.45\pm2.24)\times10^{-3}$ | $0.0913\pm0.0098$ | $99.30\pm5.01$ |
| 8 | Somma et al. | $(2.82\pm1.14)\times10^{-3}$ | $0.0608\pm0.0060$ | $107.20\pm14.10$ |
| 8 | Wang–Liu | $(3.66\pm0.42)\times10^{-3}$ | $0.0617\pm0.0057$ | $21.93\pm9.20$ |
| 8 | Clayton et al. | N/A | N/A | -- |
| 8 | Benedetti et al. | N/A | N/A | -- |

### Table 3. Convection–diffusion

| $n$ | Method | Ideal error | Finite-shot error | Emulator ms |
|---|---|---|---|---|
| 2 | **QRK-Sylv** | $(7.03\pm9.18)\times10^{-9}$ | $0.0325\pm0.0095$ | $5.89\pm0.04$ |
| 2 | Somma et al. | $(1.04\pm0.28)\times10^{-12}$ | $0.0281\pm0.0060$ | $16.09\pm11.00$ |
| 2 | Wang–Liu | $(1.08\pm0.37)\times10^{-3}$ | $0.0318\pm0.0076$ | $18.98\pm8.32$ |
| 2 | Clayton et al. | N/A | N/A | -- |
| 2 | Benedetti et al. | N/A | N/A | -- |
| 4 | **QRK-Sylv** | $(8.61\pm7.65)\times10^{-4}$ | $0.0510\pm0.0074$ | $24.42\pm3.05$ |
| 4 | Somma et al. | $(1.14\pm0.19)\times10^{-12}$ | $0.0413\pm0.0078$ | $26.01\pm9.81$ |
| 4 | Wang–Liu | $(1.90\pm0.38)\times10^{-3}$ | $0.0438\pm0.0072$ | $23.75\pm11.45$ |
| 4 | Clayton et al. | N/A | N/A | -- |
| 4 | Benedetti et al. | N/A | N/A | -- |
| 8 | **QRK-Sylv** | $(8.88\pm9.45)\times10^{-3}$ | $0.0909\pm0.0072$ | $106.64\pm10.85$ |
| 8 | Somma et al. | $(1.12\pm0.16)\times10^{-12}$ | $0.0602\pm0.0066$ | $100.36\pm12.79$ |
| 8 | Wang–Liu | $(3.43\pm0.59)\times10^{-3}$ | $0.0625\pm0.0071$ | $22.80\pm10.62$ |
| 8 | Clayton et al. | N/A | N/A | -- |
| 8 | Benedetti et al. | N/A | N/A | -- |

### Table 4. Indefinite coupled-wave family

Wang–Liu is N/A because its required half-plane spectral separation is absent.

| $n$ | Method | Ideal error | Finite-shot error | Emulator ms |
|---|---|---|---|---|
| 2 | **QRK-Sylv** | $(4.01\pm8.30)\times10^{-4}$ | $0.0318\pm0.0084$ | $6.07\pm0.21$ |
| 2 | Somma et al. | $(1.89\pm0.79)\times10^{-12}$ | $0.0302\pm0.0068$ | $16.58\pm9.58$ |
| 2 | Wang–Liu | N/A | N/A | -- |
| 2 | Clayton et al. | N/A | N/A | -- |
| 2 | Benedetti et al. | N/A | N/A | -- |
| 4 | **QRK-Sylv** | $(8.53\pm13.39)\times10^{-4}$ | $0.0480\pm0.0078$ | $24.94\pm1.96$ |
| 4 | Somma et al. | $(2.40\pm0.54)\times10^{-12}$ | $0.0423\pm0.0072$ | $26.39\pm12.56$ |
| 4 | Wang–Liu | N/A | N/A | -- |
| 4 | Clayton et al. | N/A | N/A | -- |
| 4 | Benedetti et al. | N/A | N/A | -- |
| 8 | **QRK-Sylv** | $(9.28\pm9.66)\times10^{-4}$ | $0.0928\pm0.0061$ | $104.79\pm9.32$ |
| 8 | Somma et al. | $(2.45\pm0.38)\times10^{-12}$ | $0.0602\pm0.0070$ | $103.16\pm20.62$ |
| 8 | Wang–Liu | N/A | N/A | -- |
| 8 | Clayton et al. | N/A | N/A | -- |
| 8 | Benedetti et al. | N/A | N/A | -- |

### Table 5. Nonnormal stable-control Lyapunov family

Somma is N/A because the coefficient is defective; Benedetti is N/A because it is nonnormal.

| $n$ | Method | Ideal error | Finite-shot error | Emulator ms |
|---|---|---|---|---|
| 2 | **QRK-Sylv** | $(4.04\pm13.99)\times10^{-11}$ | $0.0280\pm0.0068$ | $6.17\pm0.37$ |
| 2 | Somma et al. | N/A | N/A | -- |
| 2 | Wang–Liu | $(6.01\pm2.08)\times10^{-4}$ | $0.0301\pm0.0070$ | $16.81\pm5.48$ |
| 2 | Clayton et al. | $(1.53\pm0.62)\times10^{-6}$ | $0.0304\pm0.0069$ | $2.17\pm0.66$ |
| 2 | Benedetti et al. | N/A | N/A | -- |
| 4 | **QRK-Sylv** | $(3.31\pm6.52)\times10^{-5}$ | $0.0495\pm0.0082$ | $24.43\pm1.17$ |
| 4 | Somma et al. | N/A | N/A | -- |
| 4 | Wang–Liu | $(8.48\pm2.26)\times10^{-4}$ | $0.0447\pm0.0082$ | $17.40\pm3.57$ |
| 4 | Clayton et al. | $(4.23\pm3.62)\times10^{-5}$ | $0.0434\pm0.0080$ | $2.70\pm0.90$ |
| 4 | Benedetti et al. | N/A | N/A | -- |
| 8 | **QRK-Sylv** | $(7.12\pm16.75)\times10^{-4}$ | $0.0940\pm0.0085$ | $98.14\pm5.55$ |
| 8 | Somma et al. | N/A | N/A | -- |
| 8 | Wang–Liu | $(1.06\pm0.21)\times10^{-3}$ | $0.0708\pm0.0075$ | $21.91\pm8.96$ |
| 8 | Clayton et al. | $(4.03\pm4.36)\times10^{-4}$ | $0.0725\pm0.0093$ | $5.93\pm6.86$ |
| 8 | Benedetti et al. | N/A | N/A | -- |

### Table 6. Complex clock-dilation validation

Mean ± sample SD over 30 independent complex instances; the active-branch error compares the simulated and analytical branch.

| Quantity | $n=2,\ T=200$ | $n=4,\ T=800$ |
|---|---|---|
| Active-branch error | $(3.37\pm2.48)\times10^{-15}$ | $(6.71\pm4.34)\times10^{-15}$ |
| Active probability | $0.7904\pm0.0223$ | $0.9403\pm0.0061$ |
| Final RK relative error | $(1.37\pm0.75)\times10^{-16}$ | $(2.64\pm4.12)\times10^{-11}$ |
| Data + clock qubits | 11 | 15 |
| Theorem lower bound on active probability | $e^{-1}$ | $e^{-1}$ |

### Table 7. Qubit scaling

| $n$ | $N=n^2$ | $T$ | Data/reference qubits | Clock qubits |
|---|---|---|---|---|
| 2 | 4 | 200 | 3 | 8 |
| 4 | 16 | 800 | 5 | 10 |
| 8 | 64 | 3200 | 7 | 12 |
| 16 | 256 | 12800 | 9 | 14 |
| 32 | 1024 | 51200 | 11 | 16 |
| 64 | 4096 | 204800 | 13 | 18 |
| 128 | 16384 | 819200 | 15 | 20 |

## Algorithm notes

The complex row normal has entries $\overline{A_{ik}}$ and $\overline{B_{\ell j}}$, and squared norm

$$d_{ij}=\|A_{i,:}\|_2^2+\|B_{:,j}\|_2^2+2\operatorname{Re}(A_{ii}\overline{B_{jj}}).$$

For general complex inputs, the conjugation in the cross term follows directly from the definition of the Hermitian row normal $\Gamma_{ij}$ and is necessary for the row to have the stated squared norm. The source manuscript displays the cross term without the conjugation; the Rust implementation uses the algebraically consistent norm. This distinction disappears for real-valued coefficients.

The sampler uses $p_{ij}=d_{ij}/\sum_{r,s}d_{rs}$ without materializing the Kronecker operator. The exact clock module uses $\sigma=\sqrt{T}R$, $\alpha_t=\sqrt{1+|\beta_t/\sigma|^2}$ and a distinct orthogonal failure sector at every step. The expected probability of the active branch is at least $e^{-1}$ when $R$ is a valid bound.

The paper's tables report **historical experiment values**, while the Rust executable computes new numbers using the explicit numerical kernels and runtime of the local machine. Comparing their milliseconds or finite-shot standard deviations as identical procedures would be misleading. The paper's structural comparisons and every `N/A` statement are retained above.

## References for the algorithms

- Somma, Low, Berry, and Babbush, [Quantum algorithm for linear matrix equations](https://arxiv.org/abs/2508.02822), 2025/2026.
- Wang and Liu, [Sign Embedding Quantum Algorithms for Matrix Equations and Matrix Functions](https://arxiv.org/abs/2604.25333), 2026.
- Clayton et al., [Differentiable Quantum Computing for Large-scale Linear Control](https://arxiv.org/abs/2411.01391), 2024.
- Benedetti, Rosmanis, and Rosenkranz, [Probabilistic quantum algorithm for Lyapunov equations and matrix inversion](https://arxiv.org/abs/2508.04689), 2025/2026.
