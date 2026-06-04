# Relativistic Effects in PySCF

Scalar relativistic effects become important for elements beyond Kr (Z > 36),
and essential for 5d/5f elements (lanthanides, actinides, heavy main-group).
PySCF supports X2C (exact two-component), DKH (Douglas-Kroll-Hess), and
PP (effective core potentials / pseudopotentials).

---

## 1. Scalar relativistic: X2C (recommended)

X2C is the modern standard — it folds relativistic effects into a one-component
Hamiltonian with minimal overhead. Use for any system containing 4d+ elements.

```python
from pyscf import gto, scf, cc
from pyscf.x2c import x2c

mol = gto.Mole()
mol.atom = """
U  0.0  0.0  0.0
O  0.0  0.0  1.8
O  0.0  0.0 -1.8
"""
# Use an all-electron relativistic basis — Dyall, Faegri, or ANO-RCC
mol.basis = {"U": "dyall.cv2z", "O": "cc-pVTZ"}
mol.charge = 0
mol.spin   = 0
mol.build()

# Wrap SCF with X2C
mf = x2c.RHF(mol)
mf.conv_tol = 1e-9
mf.run()

# CCSD(T) on X2C reference
mycc = cc.CCSD(mf).run()
e_t  = mycc.ccsd_t()
print(f"X2C-CCSD(T): {mf.e_tot + mycc.e_corr + e_t:.8f} Ha")
```

---

## 2. Scalar relativistic: DKH2

```python
from pyscf import gto, scf

mol = gto.Mole(
    atom="Pt 0 0 0; Cl 0 0 2.3; Cl 0 0 -2.3",
    basis={"Pt": "ano-rcc-vdz", "Cl": "cc-pVTZ"},
    verbose=4,
).build()

mf = scf.RHF(mol)
mf.with_x2c = False
# Enable DKH order-2
mf = scf.sfx2c1e.sfx2c1e(mf)   # DKH-equivalent scalar correction
mf.run()
```

---

## 3. Effective Core Potentials (ECPs / pseudopotentials)

ECPs are cheaper than all-electron relativistic methods and often sufficient for
valence properties. Stuttgart ECPs are built into PySCF.

```python
mol = gto.Mole()
mol.atom = "Au 0 0 0; Cl 0 0 2.2"
mol.basis = {
    "Au": "def2-tzvpp",    # Stuttgart ECP is automatic with def2 on heavy atoms
    "Cl": "def2-tzvpp",
}
mol.ecp = {"Au": "def2-tzvpp"}   # explicit ECP specification
mol.build()

mf = scf.RHF(mol).run()
```

---

## 4. Spin-orbit coupling (2-component / 4-component)

For spin-orbit effects (intersystem crossing, heavy-atom effect on phosphorescence):

```python
from pyscf.x2c import tdscf_grad     # TD-DFT with SOC
from pyscf import mcscf
from pyscf.socutils import soc_2c    # spin-orbit CASSCF

# Two-step: CASSCF → spin-orbit CI
mc = mcscf.CASSCF(mf, 8, 8).run()
# Then pass mc to soc_2c for state interaction
```

Full 4-component Dirac calculations are available via `pyscf.dirac`.

---

## 5. FCIDUMP with relativistic integrals (for PyBEST pCCSD)

Recommended workflow for actinide pCCSD:

```python
from pyscf import gto, scf, ao2mo
from pyscf.tools import fcidump
from pyscf.x2c import x2c

# 1. Build X2C-HF reference
mol = gto.Mole(atom="...", basis="dyall.cv2z", ...).build()
mf  = x2c.RHF(mol).run()

# 2. Run CASSCF to define active space
from pyscf import mcscf
mc = mcscf.CASSCF(mf, ncas=12, nelecas=10).run()

# 3. Dump active-space Hamiltonian as FCIDUMP
fcidump.from_mcscf(mc, "actinide_active.FCIDUMP")

# 4. Feed FCIDUMP to PyBEST pCCSD — see references/pybest_recipes.md §4
```

---

## 6. Relativistic basis sets — quick guide

| Element range | Recommended basis | Notes |
|---------------|-------------------|-------|
| 1st/2nd row | cc-pVXZ, def2-TZVPP | Non-relativistic fine |
| 3d TM | def2-TZVPP + ECP | or all-electron cc-pVTZ-DK with DKH |
| 4d/5d TM | cc-pVXZ-PP (Stuttgart ECP) | or all-electron with X2C |
| Lanthanides | ANO-RCC + X2C | Dyall cVXZ for core properties |
| Actinides | Dyall cVXZ or ANO-RCC + X2C | 5f occupation critical |

---

## 7. Interpretation

- Scalar relativistic effects contract/stabilise s and p orbitals → shorter bonds, higher IPs for 6s elements
- Gold anomaly, mercury being a liquid, Pb being stable +2 rather than +4 — all relativistic
- X2C total energies are NOT directly comparable with non-relativistic energies
- T1 diagnostics from relativistic CCSD remain valid: > 0.02 suggests multireference character
