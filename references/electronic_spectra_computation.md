# Electronic Spectra Workflow: UV-Vis Excitation Energies

## Scientific Goal

Compute UV-Vis excitation energies and oscillator strengths to simulate or interpret electronic absorption spectra. The primary targets are vertical excitation energies (Franck-Condon point), oscillator strengths, and the nature of each excited state (local vs. charge-transfer, bright vs. dark, singlet vs. triplet).

---

## Recommended Methods

### TD-DFT — medium systems (10–200 atoms)

The workhorse for UV-Vis spectra. Fast and generally reliable for valence excitations of closed-shell molecules. Accuracy is ~0.2–0.5 eV for well-behaved states using hybrid functionals.

- **B3LYP / PBE0**: good for local excitations
- **CAM-B3LYP / ωB97X-D**: range-separated; required for charge-transfer (CT) states — standard hybrids severely underestimate CT energies
- **M06-2X**: reasonable for Rydberg states

### EOM-CCSD — high accuracy, small/medium systems (< ~30 heavy atoms)

Gold standard for singly excited states. Typical accuracy 0.1–0.2 eV. Use when TD-DFT is unreliable or a benchmark is needed. Cost scales as O(N⁶).

### CASSCF / CASPT2 — multireference and strongly correlated states

Required when the ground or excited state has strong multireference character (diradicals, bond breaking, conical intersections, transition metal complexes). CASSCF provides qualitative wavefunctions; CASPT2 adds dynamic correlation for quantitative energies. Active space selection is the critical human judgment step.

---

## Accuracy / Cost Tradeoffs

| Method       | Typical error (eV) | Scaling | Best for                          |
|--------------|--------------------|---------|-----------------------------------|
| TD-DFT/B3LYP | 0.2–0.5            | O(N³)   | Valence excitations, large systems |
| TD-DFT/CAM-B3LYP | 0.2–0.4        | O(N³)   | CT states, Rydberg                |
| EOM-CCSD     | 0.1–0.2            | O(N⁶)   | Benchmark, singly excited states  |
| CASPT2       | 0.1–0.3            | Exp(N_active) | Multireference, conical intersections |

Use TD-DFT for a first pass. Upgrade to EOM-CCSD or CASPT2 only where TD-DFT is known to fail or where high accuracy is required.

---

## Example Python Code (PySCF TD-DFT)

```python
from pyscf import gto, dft, tddft
import numpy as np

# Define molecule
mol = gto.M(
    atom="N 0 0 0; N 0 0 1.1",
    basis="def2-SVP",
    verbose=4
)

# Ground-state DFT
mf = dft.RKS(mol)
mf.xc = "B3LYP"
mf.kernel()

# TD-DFT excited states
td = tddft.TDDFT(mf)
td.nstates = 10
td.kernel()

# Extract excitation energies and oscillator strengths
energies_eV = td.e * 27.2114  # Hartree → eV
oscillator_strengths = td.oscillator_strength()

print(f"{'State':>6}  {'Energy (eV)':>12}  {'Osc. Strength':>14}")
for i, (e, f) in enumerate(zip(energies_eV, oscillator_strengths)):
    print(f"{i+1:>6}  {e:>12.4f}  {f:>14.6f}")
```

For charge-transfer systems, swap the functional:

```python
mf.xc = "CAM-B3LYP"
```

For EOM-CCSD with PySCF:

```python
from pyscf import cc

mycc = cc.RCCSD(mf).run()
myeom = mycc.eomee_ccsd_singlet(nroots=5)
# myeom[0] contains excitation energies in Hartree
```

---

## Visualization Steps

1. **Broaden stick spectrum** — convolve stick spectrum (delta functions at each excitation energy, weighted by oscillator strength) with a Gaussian or Lorentzian of width σ ≈ 0.1–0.3 eV to mimic experimental broadening.

```python
import matplotlib.pyplot as plt

# Stick spectrum
wavelengths_nm = 1239.84 / energies_eV  # eV → nm

# Broadened spectrum
wl_grid = np.linspace(150, 600, 1000)
spectrum = np.zeros_like(wl_grid)
sigma = 10  # nm

for wl, f in zip(wavelengths_nm, oscillator_strengths):
    spectrum += f * np.exp(-0.5 * ((wl_grid - wl) / sigma) ** 2)

plt.figure()
plt.plot(wl_grid, spectrum, label="Broadened spectrum")
plt.stem(wavelengths_nm, oscillator_strengths, linefmt="C1-",
         markerfmt="C1o", basefmt=" ", label="Stick spectrum")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Oscillator strength / Intensity")
plt.legend()
plt.tight_layout()
plt.savefig("uv_vis_spectrum.png", dpi=150)
```

2. **Natural transition orbitals (NTOs)** — visualize the hole/particle pair for each excited state to determine the excitation character. PySCF provides `td.nto_coeff()`.

3. **Difference density plots** — ρ_excited − ρ_ground shows which regions gain/lose electron density. Useful for identifying CT direction.

4. **Overlap with experiment** — align computed spectrum by a rigid shift if needed; shifts > 0.5 eV suggest wrong functional or solvent effects are important.

---

## Interpretation

**Oscillator strengths**: proportional to absorption intensity. States with f < 0.001 are effectively dark (symmetry-forbidden or spin-forbidden). Strong bands have f > 0.1.

**Bright vs. dark states**: dark states (f ≈ 0) are spectroscopically invisible in one-photon absorption but accessible in fluorescence, two-photon absorption, or phosphorescence. Report them — they often govern photophysics.

**Charge-transfer states**: characterized by large spatial separation between hole and particle NTOs. TD-DFT with pure or standard hybrid functionals severely underestimates CT excitation energies due to the self-interaction error. Always use a range-separated functional (CAM-B3LYP, ωB97X-D) for CT-active systems.

**Spin-forbidden transitions**: singlet → triplet excitations have f = 0 in the non-relativistic limit. They become weakly allowed through spin-orbit coupling (SOC). If phosphorescence or intersystem crossing is relevant, include SOC corrections explicitly.

---

## Common Pitfalls

**Wrong functional for CT states.** Standard B3LYP artificially lowers CT excitation energies by 0.5–2 eV. Symptom: implausibly low-energy states with large spatial NTO separation. Fix: switch to CAM-B3LYP or ωB97X-D.

**Insufficient basis set.** Rydberg states require diffuse functions. def2-SVP is a reasonable starting point for valence excitations; add augmented functions (aug-cc-pVDZ, ma-def2-TZVP) when Rydberg character is suspected.

**Triplet instability.** If the RHF/RKS ground state has a triplet instability (a low-lying triplet below the singlet), TD-DFT excitation energies become unphysical. Check with `mf.stability()` and use UKS if needed.

**Ignoring solvent.** Gas-phase excitation energies can shift by 0.1–0.5 eV in polar solvents. Use PCM or SMD for comparison with solution-phase UV-Vis data.

**Too few states.** Request enough states (nstates = 15–20) to capture all relevant transitions in the wavelength window of interest — a bright state can appear at higher index than expected.

**Excited-state geometry vs. vertical excitation.** TD-DFT vertical energies correspond to absorption maxima from the ground-state geometry. If comparing to emission, optimize the excited-state geometry first.

---

## Validation Strategy

1. **Benchmark against experiment** on a small, well-characterized model compound before applying to your target system. Aim for systematic (not random) error.

2. **Functional sensitivity check**: recompute with at least two functionals (e.g., B3LYP and CAM-B3LYP). If excitation energies shift by > 0.3 eV, CT character is likely and the range-separated result should be trusted.

3. **Basis set convergence**: compare def2-SVP and def2-TZVP results. For valence excitations, energies should converge to < 0.05 eV. If they do not, Rydberg diffuse functions may be needed.

4. **EOM-CCSD spot check**: for the 2–3 lowest states, run EOM-CCSD on the same geometry. Agreement within 0.2 eV confirms TD-DFT is reliable for this system.

5. **NTO analysis**: confirm the character of each excited state. A state assigned as π→π* should have NTOs localized on the π system — unexpected delocalization signals a problem.

6. **Oscillator strength sum rule**: the integrated oscillator strength is related to the number of electrons (Thomas-Reiche-Kuhn sum rule). Large deviations indicate incomplete state sampling or numerical issues.
