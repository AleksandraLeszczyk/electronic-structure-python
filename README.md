# electronic-structure-python

A Claude skill for quantum chemistry in Python. Ask it anything about electronic structure calculations and get a working script alongside a clear explanation — ready to run, with method recommendations and convergence tips included.

---

## What it does

When you ask a quantum chemistry question, the skill:

1. **Classifies your system** using a decision tree (size, open/closed-shell, ground vs. excited state, heavy atoms, etc.)
2. **Recommends a method** with justification
3. **Writes a runnable Python script** you can execute immediately
4. **Adds convergence guidance** for the chosen method
5. **Explains what to look for** in the results

---

## Packages covered

| Package | What it's for |
|---------|---------------|
| [PySCF](https://pyscf.org) | HF, DFT, MP2, CCSD(T), CASSCF, EOM-CCSD, TD-DFT, periodic systems, GPU acceleration |
| [PyBEST](http://fizyka.umk.pl/~pybest/pybest-v2.1.0/) | pCCD, OO-pCCD, fpCCSD, DMRG-tCC, EOM-pCCD — unique strongly correlated methods |
| [orbital-viz](https://github.com/orbitalviz/orbital-viz) | Interactive orbital and density visualization from Molden / cube files |

---

## Methods covered

### Ground-state wavefunction methods
HF (RHF / UHF / ROHF), MP2, CCSD, CCSD(T), CCSDT, pCCD, OO-pCCD, fpCCSD, CASCI, CASSCF, CASPT2, NEVPT2, MRCI, DMRG-CASSCF, DMRG-tailored CC, FCI

### Excited-state methods
EOM-CCSD (EE / IP / EA), EOM-pCCD, EOM-pCCD+S, EOM-pCCD-LCCSD, TD-DFT (Casida / TDA)

### DFT
LDA, GGA (PBE, BLYP), meta-GGA (TPSS, r²SCAN), hybrid (B3LYP, PBE0, M06-2X), range-separated (CAM-B3LYP, ωB97X-D), double-hybrid (B2-PLYP), DFT-D3/D4

### Special topics
Relativistic effects (X2C, DKH, 4c-DHF, spin–orbit coupling), magnetic coupling constants (BS-DFT, NEVPT2), periodic systems (PBC, k-points, band structure, DOS), GPU acceleration via gpu4pyscf

---

## Example questions

```
How do I compute a CCSD(T) single-point energy for water?
How do I run CASSCF with a (6,6) active space in PySCF?
How do I get UV-Vis absorption spectra with TD-DFT?
How do I compute the singlet-triplet gap of a diradical?
How do I scan a bond-dissociation curve with OO-pCCD in PyBEST?
How do I visualize HOMO/LUMO orbitals from a Molden file?
How do I include scalar relativistic effects for a gold complex?
How do I compute magnetic exchange coupling constants?
```

---

## Reference files

Detailed code templates and background are organized in `references/`:

| File | Contents |
|------|----------|
| `quantum_chemistry_basics.md` | Hamiltonian, MOs, Slater determinants, basis sets, PES — core concepts |
| `electronic_structure_methods.md` | Full method table: scaling, accuracy, use cases, original references |
| `method_decision_matrix.md` | Quick-pick table: which method for which problem; basis set cheat sheet |
| `pyscf_overview.md` | Complete PySCF API reference with runnable code for every major method |
| `pyscf_recipes.md` | Task-oriented PySCF scripts (geometry scans, S-T gaps, NTO analysis, …) |
| `pybest_overview.md` | Full PyBEST reference: pCCD, fpCCSD, DMRG-tCC, EOM-pCCD, orbital entanglement |
| `pybest_recipes.md` | Step-by-step PyBEST workflows |
| `visualization.md` | orbital-viz usage; Molden and cube file generation; interactive 3D rendering |
| `basis_sets.md` | Basis set families, selection guide, custom / mixed basis sets |
| `electronic_spectra_computation.md` | UV-Vis spectra: TD-DFT, EOM-CCSD, oscillator strengths, NTOs |
| `relativistic.md` | X2C, DKH, ZORA, 4-component DHF, spin–orbit coupling in PySCF |
| `magnetic_properties.md` | J-coupling constants, broken-symmetry DFT, NEVPT2 spin-state manifolds |
| `strong_correlation_computation.md` | CASSCF active-space selection, DMRG, pCCD workflows |
| `file_formats.md` | Molden, cube, FCIDUMP, HDF5 — read/write in PySCF and PyBEST |
| `pitfalls.md` | Common mistakes and how to avoid them |

---

## Method selection at a glance

```
Small molecule, high accuracy?   → CCSD(T) / cc-pVTZ
Multireference / bond-breaking?  → CASSCF + NEVPT2, or OO-pCCD + fpCCSD
Strongly correlated, large?      → DMRG-CASSCF or pCCD-based (PyBEST)
Excited states, small molecule?  → EOM-CCSD (or EOM-pCCD for MR ground state)
Excited states, large molecule?  → TD-DFT / ωB97X-D or CAM-B3LYP
> 50 heavy atoms?                → DFT (+ D3BJ dispersion)
Heavy atoms (row 4+)?            → Add X2C scalar relativity
Magnetic coupling?               → BS-DFT (fast) or NEVPT2 (accurate)
```

### Accuracy vs. cost

```
Higher accuracy
│  FCI / CCSDTQ
│  CCSD(T)  ← gold standard
│  CCSD / CASPT2 / NEVPT2
│  MP2 / EOM-CCSD / pCCD+fpCCSD
│  Hybrid DFT (B3LYP, PBE0)
│  GGA DFT (PBE, BLYP)
│  Semiempirical (GFN2-xTB)
└─────────────────────────────── Larger systems
```

---

## Unit conversions

| From | To | Factor |
|------|----|--------|
| Hartree | eV | × 27.2114 |
| Hartree | kcal/mol | × 627.509 |
| Hartree | kJ/mol | × 2625.50 |
| Bohr | Å | × 0.52918 |

---

## Script conventions

Every generated script follows these rules so you can always find what you need:

```python
# pip install pyscf          ← install line always first

mol_geom = """               ← geometry as a named variable, never buried
O  0.000  0.000  0.117
H  0.000  0.757 -0.469
"""

basis      = "cc-pVTZ"       ← tunable parameters named at the top
conv_tol   = 1e-9
max_cycle  = 150

print(f"CCSD(T) energy: {e:.8f} Ha")   ← results printed with units
print(f"S-T gap: {gap * 27.2114:.3f} eV")
```

---

## Evals

Correctness is tracked in `evals/evals.json`. Each entry is a question–answer pair that the skill should handle correctly; run evals before modifying reference files.

---

## License

See `LICENSE`.
