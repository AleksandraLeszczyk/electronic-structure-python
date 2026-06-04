# PyBEST Recipes

PyBEST (Python-Based Electronic Structure Tools) is a coupled-cluster package
with particular strengths in seniority-based methods, pair-coupled cluster,
and FCIDUMP-based workflows. Install from the project's distribution channel
(not on PyPI by default — check the group's GitLab/documentation).

---

## 1. Basic EOM-CCSD in PyBEST

```python
from pybest import Molecule, HF, CCSD, EOMCCSD

mol = Molecule.from_xyz("molecule.xyz", basis="cc-pVTZ")

# Ground-state HF
hf = HF(mol)
hf.run()

# CCSD amplitudes
ccsd = CCSD(hf)
ccsd.run()

# EOM-CCSD for excited states
eom = EOMCCSD(ccsd)
eom.nroots = 6          # request 6 roots (ask for a few more than needed)
eom.run()

ha2ev = 27.2114
for i, state in enumerate(eom.states):
    print(f"State {i+1}: {state.excitation_energy * ha2ev:.3f} eV  "
          f"f = {state.oscillator_strength:.4f}")
```

**Note on roots**: always request `nroots` 2–3 higher than you actually need —
spurious or convergence-failed roots sometimes appear in the list.

---

## 2. IP-EOM-CCSD and EA-EOM-CCSD

```python
from pybest import IPCCEOM, EACCEOM

# Ionisation potentials
ip = IPCCEOM(ccsd)
ip.nroots = 4
ip.run()
for i, state in enumerate(ip.states):
    print(f"IP state {i+1}: {state.energy * ha2ev:.3f} eV")

# Electron affinities
ea = EACCEOM(ccsd)
ea.nroots = 4
ea.run()
```

---

## 3. Pair CCSD (pCCSD) — for strongly correlated systems

pCCSD is a seniority-zero truncation of CCSD, useful for actinides and
systems with strong geminal correlations. It scales more favourably than
CCSD for pair-dominated correlation.

```python
from pybest import HF, pCCSD

hf = HF(mol)
hf.run()

pcc = pCCSD(hf)
pcc.conv_tol = 1e-9
pcc.run()

print(f"pCCSD correlation energy: {pcc.e_corr:.8f} Ha")
print(f"pCCSD total energy:       {pcc.e_tot:.8f} Ha")
```

---

## 4. FCIDUMP interface (relativistic or CASSCF active-space Hamiltonian)

PyBEST can ingest FCIDUMP files produced by, e.g., PySCF + X2C or NEVPT2
active-space transformations. This is the recommended route for relativistic
frozen-pair CCSD on actinides.

```python
# Step 1 — produce FCIDUMP with PySCF + X2C (see references/relativistic.md)
# Step 2 — run PyBEST on the FCIDUMP

from pybest import FCIDUMPMolecule, pCCSD

mol = FCIDUMPMolecule("actinide_active.FCIDUMP")
pcc = pCCSD(mol)
pcc.conv_tol = 1e-8
pcc.run()

print(f"Relativistic pCCSD energy: {pcc.e_tot:.8f} Ha")
```

---

## 5. Exporting MO coefficients for orbital-viz

PyBEST can dump MO coefficients as numpy arrays for use with `orbital-viz`.
The recommended route is to export a Molden file from the underlying HF
object and then load it with `orbital-viz`.

```python
# Export Molden from PyBEST's HF object
hf.write_molden("molecule.molden")

# Then visualise (see references/visualization.md):
from orbital_viz import parse_molden_to_dict, BasisGTO, read_molden_c_matrix
from orbital_viz import plot_molecular_orbital

data   = parse_molden_to_dict("molecule.molden")
basis  = BasisGTO(**data)
C      = read_molden_c_matrix("molecule.molden", basis.n_basis)

fig = plot_molecular_orbital(basis, C, n=basis.n_basis // 2 - 1)  # HOMO
fig.show()
```

If Molden export is unavailable, dump numpy arrays directly:

```python
import numpy as np
# hf.mo_coeff: shape (n_basis, n_mo)
np.save("mo_coefficients.npy", hf.mo_coeff)

# Then build BasisGTO manually — see references/visualization.md §2
```

---

## 6. CCSD(T) in PyBEST

```python
from pybest import CCSD_T

ccsd_t = CCSD_T(hf)
ccsd_t.run()
print(f"CCSD(T) energy: {ccsd_t.e_tot:.8f} Ha")
```

---

## Convergence tips for PyBEST

- If CCSD amplitudes oscillate, reduce `pcc.diis_start = 3` and increase `pcc.diis_space = 10`
- For EOM-CCSD, use `eom.davidson_max_space = 100` (large Davidson subspace) for dense spectra
- If roots collapse onto each other, increase `eom.nroots` and look for degeneracies
