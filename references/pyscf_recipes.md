# PySCF Recipes

All scripts assume `pip install pyscf` (≥ 2.5). For GPU acceleration: `pip install gpu4pyscf`.

---

## 1. Boilerplate molecule setup

```python
from pyscf import gto, scf

mol = gto.Mole()
mol.atom = """
C  0.0000  0.0000  0.0000
H  0.6276  0.6276  0.6276
H -0.6276 -0.6276  0.6276
H -0.6276  0.6276 -0.6276
H  0.6276 -0.6276 -0.6276
"""
mol.basis  = "cc-pVTZ"
mol.charge = 0
mol.spin   = 0          # 2S; 0 = singlet, 2 = triplet
mol.verbose = 4         # 0=quiet … 9=debug
mol.build()
```

---

## 2. Restricted Hartree-Fock / DFT

```python
# HF
mf = scf.RHF(mol)
mf.conv_tol = 1e-10
e_hf = mf.kernel()

# DFT — swap in any functional
from pyscf import dft
mf = dft.RKS(mol)
mf.xc = "B3LYP"
mf.grids.level = 4          # integration grid quality (0–9)
e_dft = mf.kernel()
```

**Convergence tips**:
- `mf.diis_space = 14` — bigger DIIS history
- `mf.level_shift = 0.3` — stabilises open-shell / small-gap systems
- `mf.damp = 0.5` combined with `mf.max_cycle = 200` for metals

---

## 3. CCSD and CCSD(T)

```python
from pyscf import cc

mf = scf.RHF(mol).run()

mycc = cc.CCSD(mf)
mycc.conv_tol      = 1e-8
mycc.conv_tol_normt = 1e-6
e_corr, t1, t2 = mycc.kernel()

# Perturbative triples
e_t = mycc.ccsd_t()

e_total_ccsd_t = mf.e_tot + e_corr + e_t
print(f"CCSD(T) total energy: {e_total_ccsd_t:.8f} Ha")
print(f"T1 diagnostic: {mycc.get_t1_diagnostic():.4f}")
# T1 > 0.02 → caution; T1 > 0.05 → results unreliable
```

---

## 4. Singlet-triplet gap

```python
from pyscf import gto, scf, cc

def energy(spin, charge=0):
    mol = gto.Mole(
        atom="C 0 0 0; C 0 0 1.2",   # e.g. vinylidene
        basis="aug-cc-pVTZ",
        spin=spin,
        charge=charge,
        verbose=0,
    ).build()
    mf = scf.UHF(mol).run()
    mycc = cc.UCCSD(mf)
    mycc.kernel()
    e_t_corr = mycc.ccsd_t()
    return mf.e_tot + mycc.e_corr + e_t_corr

e_s = energy(spin=0)
e_t = energy(spin=2)
gap_ev = (e_t - e_s) * 27.2114
print(f"S–T gap: {gap_ev:.3f} eV  ({'triplet GS' if gap_ev < 0 else 'singlet GS'})")
```

For **benzene** use `spin=0` (singlet) and `spin=2` (lowest triplet T1):
```python
benzene = """
C  1.2124  0.7000  0.0
C  1.2124 -0.7000  0.0
C  0.0000 -1.4000  0.0
C -1.2124 -0.7000  0.0
C -1.2124  0.7000  0.0
C  0.0000  1.4000  0.0
H  2.1560  1.2445  0.0
H  2.1560 -1.2445  0.0
H  0.0000 -2.4890  0.0
H -2.1560 -1.2445  0.0
H -2.1560  1.2445  0.0
H  0.0000  2.4890  0.0
"""
```
Use `UCCSD(T)` for both states; geometry-relax each state independently for
adiabatic gap.

---

## 5. CASSCF and NEVPT2

```python
from pyscf import mcscf, mrpt

mf = scf.RHF(mol).run()

# Active space: (n_electrons, n_orbitals)
# Rule of thumb: include all nearly-degenerate orbitals near HOMO/LUMO
mc = mcscf.CASSCF(mf, ncas=8, nelecas=8)
mc.max_cycle_macro = 200
mc.conv_tol = 1e-9
e_casscf = mc.kernel()[0]

# NEVPT2 on top of CASSCF (preferred over CASPT2 in PySCF)
e_nevpt2 = mrpt.NEVPT2(mc).kernel()
print(f"CASSCF:  {e_casscf:.8f} Ha")
print(f"NEVPT2:  {e_casscf + e_nevpt2:.8f} Ha")
```

**Active space selection guide**:
1. Run RHF, inspect MO energies near HOMO/LUMO
2. Check natural orbital occupancies from a cheap CASSCF(2,2)
3. Include all orbitals with occupancy between 0.02 and 1.98
4. For TM complexes: d-orbitals + σ-bonding counterparts ("double-shell")

---

## 6. EOM-CCSD (excited states)

```python
from pyscf import cc

mf = scf.RHF(mol).run()
mycc = cc.RCCSD(mf).run()

# EOM-EE: neutral excitations (UV-Vis)
eom = mycc.eomee_ccsd()
# Ask for nroots + 2–3 buffer states
nroots = 5
energies, vecs = eom.kernel(nroots=nroots)

ha2ev = 27.2114
for i, e in enumerate(energies):
    print(f"State {i+1}: {e*ha2ev:.3f} eV")
```

For ionisation potentials: `eom_ip_ccsd()`.
For electron affinities: `eom_ea_ccsd()`.

---

## 7. TD-DFT (UV-Vis spectra)

```python
from pyscf import dft, tddft
import numpy as np

mf = dft.RKS(mol)
mf.xc = "CAM-B3LYP"          # range-separated: better for CT states
mf.basis = "aug-cc-pVDZ"     # augmented basis crucial for Rydberg states
mf.run()

td = tddft.TDDFT(mf)
td.nstates = 10
td.kernel()

# Print spectrum
print(f"{'State':>6}  {'E (eV)':>10}  {'f (osc.str.)':>14}")
for i, (e, f) in enumerate(zip(td.e, td.oscillator_strength())):
    print(f"{i+1:>6}  {e*27.2114:>10.3f}  {f:>14.6f}")
```

**Broadened spectrum**:
```python
import numpy as np
import matplotlib.pyplot as plt

energies_ev = np.array(td.e) * 27.2114
osc = td.oscillator_strength()

# Gaussian broadening
x = np.linspace(2, 10, 500)
sigma = 0.2   # eV FWHM / 2.355
spectrum = sum(f * np.exp(-0.5*((x - e)/sigma)**2)
               for e, f in zip(energies_ev, osc))

plt.plot(x, spectrum)
plt.xlabel("Energy (eV)")
plt.ylabel("Intensity (arb.)")
plt.title("Simulated UV-Vis")
plt.savefig("uvvis.png", dpi=150)
plt.show()
```

---

## 8. Geometry scan (reaction coordinate)

```python
import numpy as np
from pyscf import gto, scf, cc

distances = np.linspace(0.8, 3.0, 30)   # Å
energies  = []

for d in distances:
    mol = gto.Mole(
        atom=f"H 0 0 0; H 0 0 {d}",
        basis="aug-cc-pVTZ",
        verbose=0,
    ).build()
    mf  = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    e_t  = mycc.ccsd_t()
    energies.append(mf.e_tot + mycc.e_corr + e_t)

energies = np.array(energies)
# Relative energies in kcal/mol
rel = (energies - energies[0]) * 627.509

import matplotlib.pyplot as plt
plt.plot(distances, rel)
plt.xlabel("Bond distance (Å)")
plt.ylabel("Relative energy (kcal/mol)")
plt.title("H₂ dissociation — CCSD(T)/aug-cc-pVTZ")
plt.savefig("scan.png", dpi=150)
plt.show()
```

**For multireference bond breaking** (e.g. N₂, Cr₂): switch to CASSCF for
the dissociation leg — RHF collapses beyond ~2× equilibrium bond length.

---

## 9. Molden file export (for orbital-viz)

```python
from pyscf.tools import molden

mf = scf.RHF(mol).run()
with open("molecule.molden", "w") as f:
    molden.header(mol, f)
    molden.orbital_coeff(mol, f, mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)
```

Then visualise — see `references/visualization.md`.

---

## 10. Cube file export

```python
from pyscf.tools import cubegen

mf = scf.RHF(mol).run()

# Export a specific MO (0-indexed from lowest)
cubegen.orbital(mol, "homo.cube", mf.mo_coeff[:, mol.nelectron//2 - 1])
cubegen.orbital(mol, "lumo.cube", mf.mo_coeff[:, mol.nelectron//2])

# Export electron density
cubegen.density(mol, "density.cube", mf.make_rdm1())
```
