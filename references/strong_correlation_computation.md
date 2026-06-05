# Workflows for Strong Correlation: CASSCF/CASPT2 and pCCD-Based Methods

> Covers static and dynamic electron correlation for strongly correlated systems.
> Methods: CASSCF, CASPT2, NEVPT2, pCCD, pCCD-LCCSD, pCCD+PT2.

---

## Table of Contents

1. [Scientific Goal](#1-scientific-goal)
2. [Recommended Methods](#2-recommended-methods)
3. [Accuracy / Cost Tradeoffs](#3-accuracy--cost-tradeoffs)
4. [Example Python Code](#4-example-python-code)
5. [Visualization Steps](#5-visualization-steps)
6. [Common Pitfalls](#6-common-pitfalls)
7. [Validation Strategy](#7-validation-strategy)

---

## 1. Scientific Goal

Strong (static) correlation arises when a single Slater determinant is an inadequate zero-order description. Canonical examples include:

- **Bond breaking / potential energy surfaces (PES)** — homolytic dissociation opens near-degenerate configurations.
- **Transition-metal complexes** — partially filled d/f shells with many low-lying states.
- **Biradicals and diradicaloids** — singlet–triplet splittings, conical intersections.
- **Extended π-systems** — polyacenes, graphene nanoflakes with strong radical character.
- **Excited states** — charge-transfer, double-excitation, dark states inaccessible to TDDFT.

**The two-layer challenge:**

| Layer | Physical origin | Must capture? |
|-------|----------------|--------------|
| Static correlation | Near-degeneracy, multi-reference character | Yes — wrong qualitative picture otherwise |
| Dynamic correlation | Short-range electron–electron cusp | Yes — 10–50 kcal/mol errors if missing |

A correct workflow addresses **both layers** in a balanced, size-consistent way.

---

## 2. Recommended Methods

### 2.1 CASSCF — Complete Active Space SCF

CASSCF performs a full CI within a chosen *active space* (n electrons in m orbitals, written CAS(n,m)) while simultaneously optimizing all orbitals. It captures **static correlation exactly within the active space** but contains **no dynamic correlation**.

**Active space selection rules:**
- Include all orbitals involved in bond breaking or near-degeneracy.
- Include bonding + antibonding pairs for each breaking bond.
- Include all d/f orbitals in metal complexes.
- Upper practical limit: CAS(18,18) with conventional FCI; DMRG or FCIQMC extend this.

**Active space diagnostic:** natural orbital occupation numbers (NOONs) outside [0.02, 1.98] signal strongly correlated orbitals that belong in the active space.

### 2.2 CASPT2 — Complete Active Space Perturbation Theory (2nd order)

CASPT2 adds **dynamic correlation** on top of CASSCF via multi-reference second-order perturbation theory. Uses the CASSCF wavefunction as the reference.

Key parameters:
- **IPEA shift** (default 0.25 a.u. in OpenMolcas/Molpro) — corrects the ionization potential / electron affinity imbalance. Set to 0.0 for ground-state PES; 0.25 for excitation energies.
- **Level shift** (typically 0.1–0.5 a.u.) — regularizes near-zero denominators (intruder states). Use the smallest value that eliminates intruder states.
- **SS-CASPT2 vs. MS-CASPT2 vs. XMS-CASPT2** — multi-state variants for near-degenerate electronic states; XMS-CASPT2 is invariant to active-space rotations and preferred for conical intersections.

### 2.3 NEVPT2 — N-Electron Valence State PT2

Alternative to CASPT2 using a different zero-order Hamiltonian (Dyall Hamiltonian). Key advantages:
- **Intruder-state free** by construction — no level shift needed.
- Slightly more expensive than CASPT2 but more robust.
- Strongly contracted (SC) and partially contracted (PC) variants.

Use NEVPT2 as a cross-check when CASPT2 requires large level shifts.

### 2.4 pCCD — Pair Coupled Cluster Doubles

pCCD (also called AP1roG) restricts the coupled cluster amplitudes to electron-pair excitations: each electron pair moves from occupied orbital *i* to virtual orbital *a* together. This gives a **geminal wavefunction** that:

- Scales as **O(o·v)** in amplitude optimization — dramatically cheaper than CASSCF for large active spaces.
- Captures strong (static) correlation in the pair channel with near-CASSCF quality for single-bond breaking.
- Is **size-consistent** and **size-extensive** by construction.
- Has an exact orbital gradient — gradient-driven orbital optimization is straightforward.

**Limitations:** neglects inter-pair correlation (seniority-2 sector only). Breaks down for multi-bond breaking without correction.

### 2.5 pCCD-LCCSD — pCCD with Linearized CC Correction

Adds linearized single and double excitations on top of pCCD to recover inter-pair dynamic correlation. Captures:
- All seniority sectors (0, 2) at low cost.
- Better than CCSD for strongly correlated systems where CCSD diverges.

### 2.6 pCCD+PT2 — pCCD with Perturbative Dynamic Correlation

Analogous to CASPT2 but using the pCCD reference wavefunction. Currently the most practical route to both static + dynamic correlation for **large strongly correlated systems** (>30 correlated electrons) where CASPT2 is intractable.

Two flavors:
- **pCCD+MP2**: simple, fast, works well near equilibrium.
- **pCCD+MRPT2**: proper multi-reference PT2 using pCCD density matrices.

### 2.7 Method Selection Flowchart

```
Start
  │
  ├─ Active space ≤ (18,18)?
  │     ├─ Yes → CASSCF → CASPT2 or NEVPT2         [gold standard]
  │     └─ No  ──────────────────────────────────┐
  │                                              │
  ├─ Primarily single-bond breaking?             │
  │     ├─ Yes → pCCD+PT2                        │
  │     └─ No  → pCCD-LCCSD+PT2 or DMRG-PT2     │
  │                                              │
  └─ >50 correlated electrons? ──────────────────┘
        └─ DMRG-CASPT2 / FCIQMC-PT2  [cutting edge]
```

---

## 3. Accuracy / Cost Tradeoffs

### Computational Scaling

| Method | Formal Scaling | Memory | Notes |
|--------|---------------|--------|-------|
| CASSCF(n,m) | O(m! / (n/2)!(m−n/2)!) for FCI step | O(m²) for CI | Exponential in active space |
| CASPT2 | O(N⁵) | O(N⁴) for AO integrals | N = number of basis functions |
| NEVPT2 (SC) | O(N⁵) | O(N⁴) | Slightly more expensive than CASPT2 |
| pCCD | O(o·v·N_iter) ≈ O(N³) | O(o·v) | o = occupied, v = virtual orbitals |
| pCCD-LCCSD | O(o²v²) | O(o²v²) | Similar to CCSD |
| pCCD+PT2 | O(N⁵) | O(N⁴) | PT2 step dominates |
| DMRG | O(M³·k³) | O(M²·k²) | M = bond dimension, k = sites |

### Accuracy Benchmarks (Typical Errors)

| Property | Method | MAE vs. FCI/Experiment |
|----------|--------|----------------------|
| Equilibrium geometry | CASSCF | 0.01–0.02 Å (bond breaking region worse) |
| Reaction energy | CASPT2 | 1–3 kcal/mol |
| Reaction energy | NEVPT2 | 1–3 kcal/mol |
| Excitation energy | XMS-CASPT2 | 0.1–0.2 eV |
| Singlet–triplet gap | CASPT2 | 0.1–0.3 eV |
| PES (single-bond breaking) | pCCD | ~2 kcal/mol vs. CASSCF |
| PES (single-bond breaking) | pCCD+PT2 | ~1 kcal/mol vs. CASPT2 |
| PES (multi-bond breaking) | pCCD-LCCSD | 2–5 kcal/mol |

### Practical Guidance

**Use CASSCF/CASPT2 when:**
- Active space fits in ≤ (16,16) orbitals.
- You need quantitative excitation energies or spin-state splittings.
- High-accuracy reaction energies for a relatively small system.

**Use pCCD-based methods when:**
- Active space would exceed (18,18) — e.g., polynuclear metal clusters, long polyenes.
- You need PES scans across many geometries (pCCD is much cheaper per point).
- Screening many candidate molecules quickly before CASPT2 refinement.
- The correlation is dominated by electron-pair breaking.

---

## 4. Example Python Code

### 4.1 CASSCF + CASPT2 with PySCF

```python
"""
CASSCF -> CASPT2 workflow for N2 dissociation curve.
Requires: pyscf >= 2.4, mrpt (pyscf extension) or pyscf-forge
"""

import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf, mcscf
from pyscf.mcscf import casci

# ── Helper: build molecule at a given bond distance ────────────────────────
def build_n2(r_angstrom: float) -> gto.Mole:
    mol = gto.Mole()
    mol.atom = f"N 0 0 0; N 0 0 {r_angstrom}"
    mol.basis = "cc-pvdz"
    mol.symmetry = True   # D2h for N2
    mol.spin = 0
    mol.charge = 0
    mol.verbose = 3
    mol.build()
    return mol


# ── Active space: CAS(6,6) — bonding + antibonding σ, πx, πy ──────────────
def run_casscf_caspt2(r: float) -> dict:
    mol = build_n2(r)
    
    # 1. RHF guess
    mf = scf.RHF(mol).run()
    
    # 2. CASSCF(6,6)
    mc = mcscf.CASSCF(mf, ncas=6, nelecas=6)
    mc.fcisolver.nroots = 1          # ground state only
    mc.max_cycle_macro = 200
    mc.conv_tol = 1e-10
    # Sort active orbitals by energy; refine with NOON inspection afterwards
    mc.kernel()
    
    casscf_energy = mc.e_tot
    
    # 3. Natural orbital occupations (NOONs)
    noons, natorbs = mcscf.addons.make_natural_orbitals(mc)
    
    # 4. CASPT2 via pyscf mrpt module  
    try:
        from pyscf.mrpt import NEVPT2   # fallback: NEVPT2 if CASPT2 unavailable
        pt2 = NEVPT2(mc)
        pt2.kernel()
        pt2_energy = mc.e_tot + pt2.e_corr
        method_used = "NEVPT2"
    except ImportError:
        # Try OpenMolcas via pyscf interface if available
        pt2_energy = None
        method_used = "CASSCF only"
    
    return {
        "r": r,
        "E_HF": mf.e_tot,
        "E_CASSCF": casscf_energy,
        "E_PT2": pt2_energy,
        "method_PT2": method_used,
        "NOONs": noons,
    }


# ── Scan the PES ───────────────────────────────────────────────────────────
bond_lengths = np.linspace(0.9, 3.5, 27)   # Angstrom
results = [run_casscf_caspt2(r) for r in bond_lengths]

E_HF      = np.array([d["E_HF"]      for d in results])
E_CASSCF  = np.array([d["E_CASSCF"]  for d in results])
E_PT2     = np.array([d["E_PT2"] if d["E_PT2"] else np.nan for d in results])
```

### 4.2 Active Space Diagnostics

```python
"""
Automatic active-space selection helpers.
"""

def analyze_noons(mc, threshold_lower=0.02, threshold_upper=1.98):
    """
    Return orbital indices with fractional occupation (strongly correlated).
    NOONs outside [threshold_lower, threshold_upper] are candidates
    for the active space.
    """
    noons, _ = mcscf.addons.make_natural_orbitals(mc)
    n_active = mc.ncas
    n_core   = mc.ncore
    
    strongly_correlated = []
    for i, occ in enumerate(noons):
        if threshold_lower < occ < threshold_upper:
            orbital_idx = n_core + i  # absolute index
            strongly_correlated.append((orbital_idx, occ))
    
    print(f"  Active space NOON analysis ({n_active} orbitals):")
    for idx, occ in strongly_correlated:
        flag = "  ← strongly correlated" if 0.1 < occ < 1.9 else ""
        print(f"    Orbital {idx:3d}: NOON = {occ:.4f}{flag}")
    
    return strongly_correlated


def t1_diagnostic(cc_object):
    """
    T1 diagnostic from a CCSD object (pyscf).
    > 0.02 signals potential multi-reference character.
    """
    t1 = cc_object.t1
    n_occ = t1.shape[0]
    t1_diag = np.sqrt(np.sum(t1**2) / n_occ)
    print(f"  T1 diagnostic: {t1_diag:.4f} "
          f"({'WARNING: MR character' if t1_diag > 0.02 else 'OK'})")
    return t1_diag
```

### 4.3 pCCD Workflow (using PyBESTLab / AP1roG)

```python
"""
pCCD and pCCD+PT2 workflow.
Requires: PyBEST (https://pypi.org/project/pybest/) or the
          AP1roG module from the Ayers group.

PyBEST installation:  pip install pybest
"""

from pybest import Molecule, RHF, pCCD, pCCDpt2
import numpy as np

def run_pccd_scan(distances: np.ndarray, basis: str = "cc-pvdz") -> dict:
    """
    Scan PES with pCCD and pCCD+PT2 for a diatomic (here: H2 for illustration).
    """
    energies_pccd  = []
    energies_pt2   = []

    for r in distances:
        # Build molecule
        mol = Molecule(
            atoms=[("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, r))],
            basis=basis,
            charge=0,
            multiplicity=1,
        )

        # RHF reference
        rhf = RHF(mol)
        rhf.kernel()

        # pCCD with orbital optimization
        pccd = pCCD(rhf)
        pccd.max_iter = 500
        pccd.conv_tol = 1e-9
        pccd.orbital_optimization = True   # crucial for strong correlation
        pccd.kernel()
        energies_pccd.append(pccd.energy_total)

        # pCCD+PT2 dynamic correction
        pt2 = pCCDpt2(pccd)
        pt2.kernel()
        energies_pt2.append(pt2.energy_total)

        print(f"  r = {r:.3f} Å | pCCD = {pccd.energy_total:.6f} "
              f"| pCCD+PT2 = {pt2.energy_total:.6f} Eh")

    return {
        "distances": distances,
        "E_pCCD":  np.array(energies_pccd),
        "E_pCCDpt2": np.array(energies_pt2),
    }


# ── pCCD amplitude inspection ──────────────────────────────────────────────
def inspect_pccd_amplitudes(pccd_object, threshold: float = 0.05):
    """
    Print dominant pair amplitudes t_ia.
    Large |t_ia| > 0.3 indicates strong correlation in that pair channel.
    """
    t_amp = pccd_object.t_amplitudes   # shape (n_occ, n_virt)
    print(f"\n  Dominant pCCD amplitudes (|t| > {threshold}):")
    for i in range(t_amp.shape[0]):
        for a in range(t_amp.shape[1]):
            if abs(t_amp[i, a]) > threshold:
                corr_label = " ← strong" if abs(t_amp[i, a]) > 0.3 else ""
                print(f"    t[{i},{a}] = {t_amp[i,a]:+.4f}{corr_label}")
```

### 4.4 Orbital Optimization Utility

```python
"""
Orbital selection and ordering utilities for CASSCF.
"""

def select_active_space_from_mp2(mf, n_active_electrons: int, n_active_orbs: int):
    """
    Use UNO (unrestricted natural orbital) or MP2 NO as initial active space guess.
    Works for closed-shell systems at equilibrium; use with care near dissociation.
    """
    from pyscf import mp

    # MP2 natural orbitals as initial guess
    mp2 = mp.MP2(mf).run()
    noons, natorbs = mp2.make_rdm1()   # density matrix → NOs

    # Sort by fractional occupation
    occ_sorted = np.sort(noons)[::-1]
    print("  MP2 natural orbital occupations (top 10):")
    for i, occ in enumerate(occ_sorted[:10]):
        print(f"    NO {i:3d}: {occ:.4f}")

    # Hand these NOs to CASSCF
    mc = mcscf.CASSCF(mf, ncas=n_active_orbs, nelecas=n_active_electrons)
    mc.mo_coeff = natorbs   # start from MP2 NOs
    return mc


def freeze_core_electrons(mc, n_frozen: int):
    """Freeze innermost n_frozen occupied orbitals in CASPT2 step."""
    mc.frozen = list(range(n_frozen))
    return mc
```

---

## 5. Visualization Steps

### 5.1 Potential Energy Surface Comparison

```python
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def plot_pes(bond_lengths, energies: dict, ref_key: str = None, title: str = "PES"):
    """
    Plot absolute energies or energies relative to the minimum of ref_key.
    
    energies: dict mapping label → np.ndarray of energies (Hartree)
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    styles = ["-o", "-s", "-^", "-D", "-v"]

    for ax in axes:
        for (label, E), col, sty in zip(energies.items(), colors, styles):
            E_plot = E.copy()
            if ref_key and ax == axes[1]:      # right panel: relative energies
                E_min = np.nanmin(energies[ref_key])
                E_plot = (E - E_min) * 627.509   # Eh → kcal/mol
                ylabel = "Relative energy (kcal/mol)"
            else:
                ylabel = "Total energy (Hartree)"
            ax.plot(bond_lengths, E_plot, sty, label=label, color=col,
                    linewidth=1.6, markersize=4)
        ax.set_xlabel("Bond distance (Å)", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.legend(frameon=False)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.grid(True, linestyle=":", alpha=0.4)

    axes[0].set_title(f"{title} — Absolute energies")
    axes[1].set_title(f"{title} — Relative to {ref_key} minimum")
    plt.tight_layout()
    plt.savefig("pes_comparison.png", dpi=200)
    plt.show()


# Example call:
# plot_pes(
#     bond_lengths,
#     {"RHF": E_HF, "CASSCF": E_CASSCF, "CASPT2": E_PT2,
#      "pCCD": E_pCCD, "pCCD+PT2": E_pCCDpt2},
#     ref_key="CASPT2",
#     title="N₂ Dissociation"
# )
```

### 5.2 NOON Heatmap

```python
def plot_noon_heatmap(noon_data: dict, title: str = "Natural Orbital Occupations"):
    """
    noon_data: dict mapping geometry label → array of NOONs.
    Useful to track how active-space character evolves along a PES.
    """
    labels    = list(noon_data.keys())
    n_geoms   = len(labels)
    n_orbs    = max(len(v) for v in noon_data.values())

    matrix = np.full((n_geoms, n_orbs), np.nan)
    for i, (label, noons) in enumerate(noon_data.items()):
        matrix[i, :len(noons)] = noons

    fig, ax = plt.subplots(figsize=(max(8, n_orbs * 0.6), max(4, n_geoms * 0.4)))
    im = ax.imshow(matrix, cmap="RdYlGn", vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(n_orbs))
    ax.set_xticklabels([f"NO {i}" for i in range(n_orbs)], rotation=45, ha="right")
    ax.set_yticks(range(n_geoms))
    ax.set_yticklabels(labels)
    plt.colorbar(im, ax=ax, label="NOON")

    # Annotate strongly correlated orbitals
    for i in range(n_geoms):
        for j in range(n_orbs):
            val = matrix[i, j]
            if not np.isnan(val) and 0.1 < val < 1.9:
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=7, color="black", fontweight="bold")

    ax.set_title(title)
    plt.tight_layout()
    plt.savefig("noon_heatmap.png", dpi=200)
    plt.show()
```

### 5.3 pCCD Amplitude Bar Chart

```python
def plot_pccd_amplitudes(t_amplitudes: np.ndarray, threshold: float = 0.05):
    """Visualize pCCD t_ia amplitudes sorted by magnitude."""
    pairs, amps = [], []
    for i in range(t_amplitudes.shape[0]):
        for a in range(t_amplitudes.shape[1]):
            if abs(t_amplitudes[i, a]) > threshold:
                pairs.append(f"({i}→{a})")
                amps.append(t_amplitudes[i, a])

    order  = np.argsort(np.abs(amps))[::-1]
    pairs  = [pairs[k] for k in order]
    amps   = [amps[k]  for k in order]
    colors = ["#d62728" if a < 0 else "#1f77b4" for a in amps]

    fig, ax = plt.subplots(figsize=(max(8, len(pairs) * 0.5), 4))
    ax.bar(pairs, amps, color=colors, edgecolor="k", linewidth=0.5)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.axhline( 0.3, color="gray", linewidth=0.8, linestyle="--", label="|t|=0.3")
    ax.axhline(-0.3, color="gray", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Pair excitation (occ→virt)")
    ax.set_ylabel("pCCD amplitude $t_{ia}$")
    ax.set_title("Dominant pCCD Amplitudes")
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("pccd_amplitudes.png", dpi=200)
    plt.show()


### 5.4 Correlation Energy Pie Chart

def plot_correlation_breakdown(E_HF, E_static, E_dynamic, labels=None):
    """
    Visualize the partitioning of correlation energy.
    E_static  = E_CASSCF - E_HF
    E_dynamic = E_CASPT2 - E_CASSCF  (at one geometry)
    """
    if labels is None:
        labels = ["Static (CASSCF)", "Dynamic (CASPT2)"]
    
    values = [abs(E_static), abs(E_dynamic)]
    colors = ["#ff7f0e", "#1f77b4"]
    explode = (0.05, 0.05)

    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(
        values, labels=labels, colors=colors, explode=explode,
        autopct=lambda p: f"{p:.1f}%\n({p/100*sum(values)*627.5:.1f} kcal/mol)",
        startangle=90, textprops={"fontsize": 11}
    )
    ax.set_title("Correlation Energy Breakdown")
    plt.tight_layout()
    plt.savefig("correlation_breakdown.png", dpi=200)
    plt.show()
```

---

## 6. Common Pitfalls

### 6.1 Active Space Selection Errors

**Wrong orbital choice:**
- Forgetting the bonding/antibonding counterpart of an included orbital violates "orbital pairing" and gives an unbalanced description.
- Including only σ but not π bonds for triple-bond breaking.
- Fix: always include pairs; use NOON analysis to verify.

**Active space too small:**
- CASSCF with too few orbitals underestimates static correlation; CASPT2 then tries to compensate, leading to intruder states.
- Fix: Expand CAS until NOONs of the outermost included orbitals are close to 0 or 2.

**Active space too large:**
- Weakly correlated orbitals (NOONs ≈ 0 or 2) in the active space waste computational effort and can cause MCSCF convergence problems.
- Fix: prune orbitals with NOON < 0.005 or > 1.995.

### 6.2 CASSCF Convergence Problems

- **State averaging** is required for excited states or near-degenerate states. Use `mc.state_average_([0.5, 0.5])` for equal weighting of two states.
- **Orbital rotation invariance:** CASSCF is invariant to rotations within the active space. This can cause the optimizer to stall. Fix: use SuperCI or augmented Hessian microiterations.
- **Root flipping:** during optimization the ground state and excited state energies may swap. Monitor state characters at each geometry.
- **Wrong initial orbitals:** starting from canonical HF orbitals often fails for transition metals. Better: localized orbitals, PM localization, or UHF natural orbitals.

### 6.3 CASPT2 / NEVPT2 Pitfalls

- **Intruder states** appear as discontinuities in the PES or anomalously large PT2 corrections. Diagnose with: `max(|E[2]_K / (E_ref - E_K)|) >> 1`. Fix: add a level shift (0.1–0.3 a.u.) and verify the PES is smooth.
- **IPEA shift confusion:** using IPEA = 0.25 a.u. for ground-state energetics can systematically lower energies. For consistent PES scans, use IPEA = 0.0 or document your choice.
- **Basis set:** CASPT2 energies converge slowly with basis set. Use at least cc-pVTZ; ANO-RCC-VTZP for heavy elements.
- **Frozen-core inconsistency:** always freeze the same core orbitals across geometries. Changing the frozen core mid-scan introduces discontinuities.

### 6.4 pCCD-Specific Pitfalls

- **Non-pairing orbitals:** pCCD assumes a pairing structure between occupied and virtual orbitals. Without orbital optimization, the default pairing (Aufbau) is often wrong. Always enable orbital optimization.
- **Amplitude divergence:** near degeneracy, pCCD amplitudes can diverge (analogous to CCSD divergence in strongly correlated cases). This is less severe than CCSD but monitor `max|t_ia|`.
- **Missing inter-pair correlation:** pCCD handles same-pair excitations but not excitations that involve electrons from different pairs. For systems with multiple interacting bonds (e.g., N₂ triple bond), use pCCD-LCCSD or pCCD+PT2 to recover inter-pair correlation.
- **Symmetry breaking:** pCCD may spontaneously break spatial or spin symmetry. Check that the final wavefunction respects the molecular symmetry.

### 6.5 General Multi-Reference Pitfalls

- **Size inconsistency of truncated CI:** CISD is not size-consistent; CASSCF is. Always verify that E(A+B at infinite separation) = E(A) + E(B).
- **Spin contamination:** CASSCF preserves spin symmetry; pCCD does not by default for open-shell systems. For biradicals, use a spin-pure formulation or check ⟨S²⟩.
- **Unbalanced treatment of states:** when computing excitation energies, use the same active space for all states (state-averaged CASSCF).

---

## 7. Validation Strategy

### 7.1 Multi-Level Consistency Checks

Run a hierarchy of methods and verify energy differences are consistent:

```
RHF → MP2 → CCSD → CCSD(T)       [single-reference ladder]
           ↓
     T1 diagnostic > 0.02?
           ↓ Yes
CASSCF(small) → CASSCF(large) → CASPT2 → NEVPT2
pCCD → pCCD+PT2 → pCCD-LCCSD
```

Agreement between CASPT2 and NEVPT2 within 1–2 kcal/mol is a good sign. Agreement between pCCD+PT2 and CASPT2 within 2–3 kcal/mol validates the pCCD approach.

### 7.2 Active Space Sensitivity Check

Run CASSCF with (n−2, m−2), (n, m), and (n+2, m+2). If the energy changes by less than 1 kcal/mol, the active space is converged.

```python
def active_space_sensitivity(mf, nelec_center: int, norbs_center: int,
                              delta: int = 2) -> dict:
    results = {}
    for dn, dm in [(-delta, -delta), (0, 0), (delta, delta)]:
        ne = nelec_center + dn
        nm = norbs_center + dm
        if ne < 2 or nm < ne // 2:
            continue
        mc = mcscf.CASSCF(mf, ncas=nm, nelecas=ne)
        mc.kernel()
        label = f"CAS({ne},{nm})"
        results[label] = mc.e_tot
        print(f"  {label}: E = {mc.e_tot:.8f} Eh")
    
    energies = list(results.values())
    spread = (max(energies) - min(energies)) * 627.509  # kcal/mol
    print(f"\n  Active space spread: {spread:.2f} kcal/mol "
          f"({'converged' if spread < 1.0 else 'NOT converged — expand CAS'})")
    return results
```

### 7.3 Basis Set Convergence

```python
def basis_set_convergence(geometry_dict: dict, bases: list[str]) -> dict:
    """
    Check CASPT2 energy convergence with basis set.
    geometry_dict: {"atoms": [...], "charge": 0, "spin": 0}
    bases: ["cc-pvdz", "cc-pvtz", "cc-pvqz"]
    """
    results = {}
    for basis in bases:
        mol = gto.Mole()
        mol.atom   = geometry_dict["atoms"]
        mol.basis  = basis
        mol.charge = geometry_dict["charge"]
        mol.spin   = geometry_dict["spin"]
        mol.build()
        # ... run CASSCF/CASPT2 and collect E_tot
        results[basis] = None   # replace with actual energy
    return results
```

Expected convergence pattern:

| Basis | Expected error vs. CBS |
|-------|----------------------|
| cc-pVDZ | 5–15 kcal/mol |
| cc-pVTZ | 1–3 kcal/mol |
| cc-pVQZ | 0.3–1 kcal/mol |
| CBS extrapolation (TZ+QZ) | < 0.1 kcal/mol |

### 7.4 Property Validation Checklist

Before trusting results, verify:

- [ ] **NOON analysis** — all strongly correlated orbitals are in the active space.
- [ ] **T1 diagnostic** — if T1 > 0.02, single-reference methods are unreliable; confirm multi-reference treatment is needed.
- [ ] **Active space sensitivity** — energy stable to ±(2,2) expansion.
- [ ] **PES smoothness** — no discontinuities along the scan (intruder states, root flipping, or orbital reordering artifacts).
- [ ] **Size consistency** — fragments at large separation give sum of isolated fragment energies.
- [ ] **Spin purity** — ⟨S²⟩ matches the target spin state.
- [ ] **Basis set convergence** — TZ and QZ agree within 1 kcal/mol for the property of interest.
- [ ] **Cross-method consistency** — CASPT2 and NEVPT2 agree; or pCCD+PT2 and CASPT2 agree within 2–3 kcal/mol.
- [ ] **Experimental benchmark** — compare to known experimental geometries, IPs, or thermochemistry for related systems.

### 7.5 Automated Validation Script

```python
def full_validation_suite(mf, nelec: int, norbs: int) -> None:
    """Run all validation checks and print a summary report."""
    print("=" * 60)
    print("STRONG CORRELATION VALIDATION SUITE")
    print("=" * 60)

    # T1 diagnostic
    from pyscf import cc
    ccsd = cc.CCSD(mf).run()
    t1 = t1_diagnostic(ccsd)

    # CASSCF
    mc = mcscf.CASSCF(mf, ncas=norbs, nelecas=nelec).run()
    noons, _ = mcscf.addons.make_natural_orbitals(mc)
    print(f"\n  CASSCF energy:  {mc.e_tot:.8f} Eh")

    # NOON analysis
    strongly_corr = analyze_noons(mc)
    frac_count = sum(1 for _, occ in strongly_corr if 0.1 < occ < 1.9)
    print(f"\n  Strongly correlated orbitals: {frac_count}")

    # Active space sensitivity
    print("\n  Active space sensitivity:")
    active_space_sensitivity(mf, nelec, norbs)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    if t1 > 0.05:
        print("  ⚠  T1 > 0.05: severe multi-reference character. CCSD unreliable.")
    elif t1 > 0.02:
        print("  ⚠  T1 > 0.02: moderate multi-reference character.")
    else:
        print("  ✓  T1 ≤ 0.02: single-reference methods may be adequate.")

    if frac_count > 0:
        print(f"  ✓  {frac_count} fractionally occupied NO(s) confirm strong correlation.")
    else:
        print("  ✓  No strongly fractional NOONs — static correlation may be mild.")
```

---

## References and Software

| Software | Methods | URL |
|----------|---------|-----|
| **PySCF** | CASSCF, NEVPT2, CCSD, MP2 | https://pyscf.org |
| **OpenMolcas** | CASSCF, CASPT2, XMS-CASPT2 | https://gitlab.com/Molcas/OpenMolcas |
| **ORCA** | CASSCF, NEVPT2, DLPNO-NEVPT2 | https://orcaforum.kofo.mpg.de |
| **PyBEST** | pCCD, pCCD+PT2, pCCD-LCCSD | https://pypi.org/project/pybest/ |
| **Block2 / DMRG** | DMRG-CASSCF, DMRG-NEVPT2 | https://github.com/block-hczhai/block2-preview |
| **Dice / FCIQMC** | SHCI, semistochastic | https://sanshar.github.io/Dice/ |

---

*Document version: 1.0 | Prepared for the Quantum Chemistry Brief project*
