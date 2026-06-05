# Quantum Chemistry Method Decision Matrix

## Which Method Should I Use?

| Problem | Recommended Method | Notes |
|---|---|---|
| Ground-state geometry | DFT | B3LYP, PBE0, or ωB97X-D; cheap and reliable for most organic/inorganic systems |
| Single-point benchmark energy | CCSD(T) | "Gold standard"; use with large basis (cc-pVTZ or better); CBS extrapolation recommended |
| Excited states (single-reference) | EOM-CCSD | Good for low-lying valence and Rydberg states; EOM-CC3 for higher accuracy |
| Excited states (charge-transfer / large systems) | TD-DFT | Use range-separated functionals (CAM-B3LYP, ωB97X); cheaper than EOM-CC |
| Strong correlation / near-degeneracy | CASSCF / CASPT2 | Choose active space carefully; NEVPT2 as alternative to CASPT2 |
| Strong correlation + large systems | DMRG, FCIQMC, pair CCD | pair CCD (pCCD) scales better; DMRG for 1D-like systems |
| Heavy atoms (4th row and beyond) | X2C or DKH | X2C (exact two-component) preferred; pair with relativistic ECPs for very heavy atoms |
| Magnetic coupling (exchange constants) | EOM-CCSD, NEVPT2, broken-symmetry DFT | BS-DFT is fastest; EOM-CCSD and NEVPT2 for higher accuracy |
| Thermochemistry (high accuracy) | W4, HEAT, or focal-point analysis | Composite methods; include ZPVE, core–valence, SO, and scalar relativistic corrections |
| Dispersion-dominated non-covalent interactions | DFT-D3/D4, DLPNO-CCSD(T) | Always add dispersion correction to DFT; SAPT for decomposition |
| Large systems (100s of atoms) | DFT, semi-empirical (GFN2-xTB), ONIOM | xTB for pre-screening; ONIOM for active-site accuracy with MM environment |
| Open-shell / radical systems | UHF-based CCSD(T), ROHF-CCSD(T) | Check ⟨S²⟩ for spin contamination; DLPNO for large open-shell systems |
| Periodic solids / band structure | DFT (plane-wave, HSE06/PBE0) | Use hybrid functionals for band gaps; GW for quasiparticle corrections |
| NMR chemical shifts | GIAO-DFT (B3LYP, PBE0) | Gauge-including atomic orbitals; pcS-n or cc-pVTZ basis sets |
| Reaction dynamics / PES sampling | AIMD (Born–Oppenheimer or CP), RPMD | AIMD at DFT level; enhanced sampling (metadynamics, umbrella) for rare events |

---

## Quick Accuracy vs. Cost Guide

```
Accuracy  ↑  CCSD(T)/CBS  ──  FCI
           │  CCSD(T)
           │  CCSD
           │  MP2
           │  DFT (hybrid)
           │  DFT (GGA)
           │  HF
           │  semi-empirical (xTB, PM7)
Cost      ↑  (scales with system size)
```

| Method | Formal Scaling | Practical Limit (atoms) |
|---|---|---|
| semi-empirical (xTB) | O(N) – O(N²) | thousands |
| DFT | O(N³) | ~500 (conventional), ~10,000 (linear-scaling) |
| MP2 | O(N⁵) | ~100–200 |
| CCSD | O(N⁶) | ~50–80 |
| CCSD(T) | O(N⁷) | ~30–50 |
| CASSCF (n,m) | exponential in active space | ≤ 18 electrons / 18 orbitals (conventional) |

---

## Basis Set Cheat Sheet

| Purpose | Suggested Basis |
|---|---|
| Quick geometry / screening | def2-SVP, 6-31G* |
| Production geometry | def2-TZVP, cc-pVTZ |
| Single-point / benchmark | def2-TZVPP, cc-pVTZ → cc-pVQZ (CBS) |
| Anions / diffuse functions needed | aug-cc-pVTZ, def2-TZVPD |
| Heavy atoms (no ECP) | x2c-TZVPall (+ X2C Hamiltonian) |
| Heavy atoms (with ECP) | def2-TZVP + ECP (e.g., Stuttgart ECPs) |
| NMR | pcS-2, cc-pVTZ-J |

---

## Common Software by Task

| Task | Software Options |
|---|---|
| DFT / general QC | ORCA, Gaussian, Q-Chem, PySCF |
| Coupled cluster | ORCA, CFOUR, MRCC, Psi4 |
| Multireference | MOLPRO, OpenMolcas, NEVPT2 in ORCA |
| Periodic DFT | VASP, Quantum ESPRESSO, CP2K |
| AIMD | CP2K, VASP, i-PI + ORCA |
| Semi-empirical | xTB (Grimme group), MOPAC |
| DMRG | Block2, CheMPS2, ITensor |
