# Psi4 Library Specifications

Psi4 is an open-source quantum chemistry package designed for high-performance ab initio electronic structure computations. This document covers five core capability areas.

---

## 1. Symmetry-Adapted Perturbation Theory (SAPT)

SAPT decomposes intermolecular interaction energies into physically meaningful components: electrostatics, exchange, induction, and dispersion. Psi4 supports several SAPT levels.

**Available levels:**

- `SAPT0` — lowest cost; uses HF monomer wavefunctions. Suitable for large systems and rapid screening.
- `SAPT2`, `SAPT2+`, `SAPT2+(3)`, `SAPT2+3` — progressively include higher-order MP2/MP3 corrections for induction and dispersion.
- `sSAPT0` — scaled SAPT0 with empirical dispersion corrections for improved accuracy at low cost.
- `FISAPT0` — functional-group SAPT; partitions interaction energy onto atomic or functional-group contributions.

**Basic input example:**

```python
molecule dimer {
  0 1
  O  0.000  0.000  0.117
  H  0.000  0.757 -0.468
  H  0.000 -0.757 -0.468
  --
  0 1
  O  0.000  0.000  3.500
  H  0.000  0.757  3.083
  H  0.000 -0.757  3.083
}

set {
  basis aug-cc-pVDZ
}

energy('sapt2+')
```

**Key settings:**

| Option | Description |
|--------|-------------|
| `basis` | Augmented basis sets (aug-cc-pVXZ) strongly recommended |
| `df_basis_sapt` | Auxiliary DF basis; set automatically when using standard basis |
| `freeze_core true` | Freeze core orbitals to reduce cost |
| `sapt_dft_functional` | Functional for SAPT(DFT) variant |

**Output components** (in kcal/mol by default):

- `Elst` — electrostatics
- `Exch` — exchange-repulsion
- `Ind` / `Exch-Ind` — induction and its exchange counterpart
- `Disp` / `Exch-Disp` — dispersion and its exchange counterpart
- `Total SAPT` — sum of all components

**Notes:** SAPT requires the dimer geometry to be specified with the `--` separator between monomers. Charge and multiplicity must be given for each monomer individually as well as for the dimer supermolecule.

---

## 2. Vibrational Analysis

Psi4 computes harmonic vibrational frequencies by finite-difference or analytic second derivatives of the energy with respect to nuclear coordinates.

**Triggering a frequency calculation:**

```python
molecule h2o {
  O
  H 1 0.96
  H 1 0.96 2 104.5
}

set {
  basis cc-pVTZ
}

# Optimize first, then compute frequencies
optimize('b3lyp')
frequencies('b3lyp')
```

Alternatively, `energy`, `gradient`, and `hessian` can be called directly. The `frequency()` wrapper calls `hessian()` internally.

**Output quantities:**

- Harmonic vibrational frequencies (cm⁻¹)
- Zero-point vibrational energy (ZPVE)
- Thermochemical corrections at a specified temperature and pressure (default: 298.15 K, 1 atm)
- Enthalpy `H`, entropy `S`, Gibbs free energy `G`
- IR intensities (when dipole derivatives are available)

**Thermochemistry settings:**

```python
set {
  t 298.15       # temperature in Kelvin
  p 101325       # pressure in Pa
  scale_freq 1.0 # frequency scaling factor
}
```

**Finite-difference control:**

```python
set findif {
  points 5          # 3 or 5-point stencil
  disp_size 0.005   # displacement size in Bohr
}
```

**Method support:** Analytic Hessians are available for HF, DFT, MP2 (with DF-MP2), and CCSD. Finite-difference Hessians are available for any method that provides energies or gradients.

**Symmetry:** Psi4 automatically uses molecular symmetry to reduce the number of displacements required. Frequencies are labeled by irreducible representation.

---

## 3. Symmetry Handling

Psi4 uses the Schoenflies notation for point group symmetry and applies it to reduce computational cost across all methods.

**Automatic detection:** Psi4 detects and assigns the highest-order Abelian subgroup compatible with the molecule's geometry. Non-Abelian groups (e.g., C₃ᵥ, T_d) are reduced to their Abelian subgroups (C_s, D₂, etc.) because the underlying integral and CI codes work in Abelian groups only.

**Controlling symmetry:**

```python
molecule benzene {
  symmetry d2h   # enforce a specific point group
  ...
}
```

To disable symmetry entirely:

```python
molecule {
  symmetry c1
  ...
}
```

**Supported point groups (Abelian):**

`C1`, `Ci`, `C2`, `Cs`, `C2h`, `C2v`, `D2`, `D2h`

**Symmetry-adapted basis functions:** Psi4 constructs symmetry-adapted linear combinations (SALCs) of atomic basis functions. Integrals, MO coefficients, and density matrices are stored in blocks by irreducible representation, reducing storage and diagonalization cost significantly.

**Orbital symmetry labels:** In output, molecular orbitals are labeled by irreducible representation and sequential index within that irrep (e.g., `2A1`, `1B2`).

**Symmetry breaking:** Some calculations (e.g., open-shell DFT, TDDFT on degenerate states, geometry optimizations near symmetry-equivalent structures) may require lowering symmetry manually to avoid artifacts or convergence problems.

**PSIFOUR symmetry API (Python):**

```python
mol = psi4.geometry("""...""")
mol.update_geometry()
print(mol.schoenflies_symbol())  # e.g., 'c2v'
print(mol.point_group().symbol())
```

---

## 4. DFT Benchmarking

Psi4 supports a broad array of density functionals and is well suited to systematic benchmarking studies.

**Functional categories available:**

| Class | Examples |
|-------|---------|
| LDA | `SVWN`, `SVWN5` |
| GGA | `BLYP`, `PBE`, `BP86` |
| meta-GGA | `M06-L`, `TPSS`, `SCAN` |
| Hybrid GGA | `B3LYP`, `PBE0`, `HSE06` |
| Hybrid meta-GGA | `M06`, `M06-2X`, `TPSSh` |
| Range-separated | `CAM-B3LYP`, `ωB97X-D`, `LC-ωPBE` |
| Double hybrid | `B2PLYP`, `DSD-PBEPBE-D3BJ` |
| DFT-D3 / DFT-D4 | any functional + `-D3BJ`, `-D3M`, `-D4` suffix |

**Running a benchmark series:**

```python
import psi4

molecules = [...]   # list of psi4.core.Molecule objects
methods   = ['b3lyp/aug-cc-pVTZ', 'pbe0/aug-cc-pVTZ', 'wb97x-d/aug-cc-pVTZ']

for mol in molecules:
    psi4.set_molecule(mol)
    for method in methods:
        e = psi4.energy(method)
        print(f"{mol.name()}, {method}: {e:.6f} Eh")
```

**Dispersion corrections:**

```python
set {
  dft_dispersion_parameters [1.0, 0.722, 1.217, 14.0]  # manual D3 parameters
}
energy('b3lyp-d3bj')
```

Psi4 interfaces with `DFTD3` and `DFTD4` programs automatically when they are installed and on `PATH`.

**DFT grid settings:**

```python
set {
  dft_spherical_points 590    # Lebedev angular grid points
  dft_radial_points    99     # Euler-Maclaurin radial points
  dft_pruning_scheme   robust # none, robust, treutler
}
```

Shorthand: `set dft_grid_name sg1` selects the SG-1 standard grid.

**Density fitting (DF-DFT / RI-DFT):**

```python
set {
  scf_type df
  df_basis_scf def2-universal-jkfit
}
```

DF reduces the formal O(N⁴) Coulomb/exchange build to approximately O(N³), enabling large-system DFT.

**Useful benchmark databases accessible via Psi4's `driver.database` module:**

- `S22`, `S66` — non-covalent interaction energies
- `A24` — small molecule interaction energies (high accuracy)
- `BH76` — barrier heights
- `W4-11` — atomization energies

---

## 5. Citation for Software

When publishing results obtained with Psi4, cite the primary software paper and, where applicable, the method-specific references below.

### Primary Psi4 Citation

> **Psi4 1.x:**
> D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish, M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio, A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer, R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni, J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein, B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Y. Sokolov, K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King, F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
> *J. Chem. Phys.* **2020**, *152*, 184108.
> DOI: [10.1063/5.0006002](https://doi.org/10.1063/5.0006002)

> **Psi4 legacy (1.1):**
> R. M. Parrish, L. A. Burns, D. G. A. Smith, A. C. Simmonett, A. E. DePrince III, E. G. Hohenstein, U. Bozkaya, A. Y. Sokolov, R. Di Remigio, R. M. Richard, J. F. Gonthier, A. M. James, H. R. McAlexander, A. Kumar, M. Saitow, X. Wang, B. P. Pritchard, P. Verma, H. F. Schaefer III, K. Patkowski, R. A. King, E. F. Valeev, F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
> *J. Chem. Theory Comput.* **2017**, *13*, 3185–3197.
> DOI: [10.1021/acs.jctc.7b00174](https://doi.org/10.1021/acs.jctc.7b00174)

### Method-Specific Citations

**SAPT:**
> E. G. Hohenstein and C. D. Sherrill, *WIREs Comput. Mol. Sci.* **2012**, *2*, 304–326.
> DOI: [10.1002/wcms.84](https://doi.org/10.1002/wcms.84)

**SAPT(DFT):**
> A. Hesselmann and G. Jansen, *Chem. Phys. Lett.* **2002**, *362*, 319–325.
> DOI: [10.1016/S0009-2614(02)01097-7](https://doi.org/10.1016/S0009-2614(02)01097-7)

**FISAPT:**
> E. E. Lao and J. M. Herbert, *J. Chem. Theory Comput.* **2014**, *10*, 4425–4436.
> DOI: [10.1021/ct5005593](https://doi.org/10.1021/ct5005593)

**DFT-D3:**
> S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* **2010**, *132*, 154104.
> DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)

**DF-MP2 / RI-MP2:**
> F. Weigend and M. Häser, *Theor. Chem. Acc.* **1997**, *97*, 331–340.
> DOI: [10.1007/s002140050269](https://doi.org/10.1007/s002140050269)

### BibTeX Entries

```bibtex
@article{psi4_2020,
  author  = {Smith, Daniel G. A. and Burns, Lori A. and Simmonett, Andrew C.
             and others},
  title   = {{Psi4 1.4}: An Open-Source Electronic Structure Program},
  journal = {J. Chem. Phys.},
  year    = {2020},
  volume  = {152},
  pages   = {184108},
  doi     = {10.1063/5.0006002}
}

@article{psi4_2017,
  author  = {Parrish, Robert M. and Burns, Lori A. and Smith, Daniel G. A.
             and others},
  title   = {{Psi4 1.1}: An Open-Source Electronic Structure Program
             Emphasizing Automation, Advanced Libraries, and Interoperability},
  journal = {J. Chem. Theory Comput.},
  year    = {2017},
  volume  = {13},
  pages   = {3185--3197},
  doi     = {10.1021/acs.jctc.7b00174}
}

@article{sapt_review,
  author  = {Hohenstein, Edward G. and Sherrill, C. David},
  title   = {Wavefunction Methods for Noncovalent Interactions},
  journal = {WIREs Comput. Mol. Sci.},
  year    = {2012},
  volume  = {2},
  pages   = {304--326},
  doi     = {10.1002/wcms.84}
}
```

### Acknowledging Psi4 in Text

A typical acknowledgment sentence for a methods section:

> All quantum chemical calculations were performed using Psi4 1.x \[ref\]. Interaction energies were computed at the SAPT2+/aug-cc-pVDZ level \[ref\]. Dispersion corrections were applied using the D3(BJ) scheme \[ref\].

---

## Quick Reference: Common Input Keywords

```python
# Globals
psi4.set_memory('8 GB')
psi4.set_num_threads(8)
psi4.set_output_file('output.dat', False)

# Basis sets often used with each area
# SAPT:            aug-cc-pVDZ, aug-cc-pVTZ
# Frequencies:     cc-pVTZ, def2-TZVP
# DFT benchmark:   def2-TZVP, aug-cc-pVTZ, def2-QZVP
```

---

*Document covers Psi4 ≥ 1.4. For the latest options and defaults, consult the official documentation at [psicode.org](https://psicode.org).*
