# PyBEST — Library Reference

PyBEST (Pythonic Black-box Electronic Structure Tool) is an open-source Python platform for electronic structure calculations at the interface between chemistry and physics. Its signature feature is a suite of pCCD-based and post-pCCD strongly correlated methods that are not available in any other quantum chemistry package, combined with a quantum information analysis framework for orbital entanglement and correlation.

The codebase is ~90% Python 3 and ~10% C++ (bound via Pybind11). Integral evaluation relies on LibInt (one- and two-electron integrals) and LibChol (Cholesky-decomposed ERIs). Current stable release: **PyBEST 2.1.0** (2025).

Documentation: <http://fizyka.umk.pl/~pybest/pybest-v2.1.0/>

---

## Table of Contents

1. [Overview of Features](#1-overview-of-features)
2. [Computing the Hamiltonian](#2-computing-the-hamiltonian)
3. [Reading and Writing Data](#3-reading-and-writing-data)
4. [Strongly Correlated Methods](#4-strongly-correlated-methods)
   - 4.1 [Orbital-Optimized pCCD (OO-pCCD) for Bond Breaking](#41-orbital-optimized-pccd-oo-pccd-for-bond-breaking)
   - 4.2 [Frozen-Pair Coupled Cluster (fpCC)](#42-frozen-pair-coupled-cluster-fpcc)
   - 4.3 [DMRG-Tailored Coupled Cluster (DMRG-tCCSD)](#43-dmrg-tailored-coupled-cluster-dmrg-tccsd)
   - 4.4 [Equation-of-Motion pCCD (EOM-pCCD)](#44-equation-of-motion-pccd-eom-pccd)
5. [Orbital Correlation Visualization](#5-orbital-correlation-visualization)
6. [Citation](#6-citation)

---

## 1. Overview of Features

PyBEST is organized as a modular Python library. Each electronic structure module accepts PyBEST objects (integrals, orbitals, occupation models) and returns an `IOData` container holding energies, amplitudes, density matrices, and checkpoint paths.

**Ground-state methods:**

| Module | Class(es) |
|--------|-----------|
| Restricted HF (RHF) | `RHF` |
| MP2 | `RMP2` |
| SAPT(0) | `RSAPT0` |
| pCCD / OO-pCCD | `RpCCD`, `ROOpCCD` |
| Perturbation theory (PT2) | `RpCCDPT2` |
| CCD / CCSD / CCS | `RCCD`, `RCCSD`, `RCCS` |
| Frozen-pair CC | `RfpCCD`, `RfpCCSD` |
| Linearized CC | `RpCCDLCCD`, `RpCCDLCCSD` |
| DMRG-tailored CCSD | `RtCCSD` |
| Restricted CI | `RRCISD`, `RRCIS`, … |
| Ionization-potential CC | `RIPCCSD` |
| Electron-attachment CC | `REACC` |
| Reversed spin-flip CC | `RRSFpCCD` |

**Excited-state methods (EOM):**

| Flavor | Class |
|--------|-------|
| EOM-CCS / EOM-CCD / EOM-CCSD | `REOMCCS`, `REOMCCD`, `REOMCCSD` |
| EOM-pCCD | `REOMpCCD` |
| EOM-pCCD+S | `REOMpCCDS` |
| EOM-pCCD-CCS | `REOMpCCDCCS` |
| EOM-pCCD-LCCD / EOM-pCCD-LCCSD | `REOMpCCDLCCD`, `REOMpCCDLCCSD` |

**Post-processing and analysis:**

- Orbital localization (Pipek–Mezey)
- Orbital entanglement analysis (single-orbital entropy, orbital-pair mutual information)
- Dipole and quadrupole moments
- MO-picture generation and correlation diagrams

**Linear algebra backends:**

- `DenseLinalgFactory` — dense four-index ERIs
- `CholeskyLinalgFactory` — Cholesky-decomposed ERIs (lower memory, recommended for large systems)

---

## 2. Computing the Hamiltonian

### 2.1 Key objects

Before any calculation, three objects must be created:

```python
from pybest.linalg import DenseLinalgFactory          # or CholeskyLinalgFactory
from pybest.occ_model import AufbauOccModel

lf        = DenseLinalgFactory(basis.nbasis)          # linear algebra factory
occ_model = AufbauOccModel(basis, ncore=0)            # Aufbau occupation model
orb_a     = lf.create_orbital(basis.nbasis)           # orbital expansion coefficients
```

Set `ncore` to the number of occupied orbitals to freeze. Leaving it at 0 lets PyBEST determine the frozen core automatically.

### 2.2 Building the molecular Hamiltonian

```python
from pybest.gbasis import (
    get_gobasis,
    compute_overlap,
    compute_kinetic,
    compute_nuclear,
    compute_eri,               # dense ERI
    compute_cholesky_eri,      # Cholesky ERI (use with CholeskyLinalgFactory)
    compute_nuclear_repulsion,
)
from pybest.context import context

fn_xyz = context.get_fn("test/water.xyz")   # or provide your own XYZ path
basis  = get_gobasis("cc-pvdz", fn_xyz)

olp = compute_overlap(basis)
kin = compute_kinetic(basis)
ne  = compute_nuclear(basis)
eri = compute_eri(basis)                    # DenseFourIndex object
nr  = compute_nuclear_repulsion(basis)      # scalar nuclear repulsion energy
```

For the Cholesky variant, replace `compute_eri` with:

```python
eri = compute_cholesky_eri(basis, threshold=1e-8)
lf  = CholeskyLinalgFactory(basis.nbasis)
```

### 2.3 Running Hartree–Fock

All post-HF methods in PyBEST accept the RHF output container as input:

```python
from pybest.wrappers import RHF

hf = RHF(lf, occ_model)
hf_output = hf(kin, ne, eri, nr, olp, orb_a)
# hf_output.e_tot  — total RHF energy
# hf_output.orb_a  — converged MO coefficients
```

### 2.4 Transforming integrals to the MO basis

Required when passing integrals to external codes (e.g., for DMRG-tCC):

```python
from pybest.utility import transform_integrals

mo_ints = transform_integrals(kin, ne, eri, hf_output.orb_a)
one = mo_ints.one[0]   # one-electron integrals in MO basis
two = mo_ints.two[0]   # two-electron integrals in MO basis
```

### 2.5 Active-space Hamiltonian

```python
from pybest.utility import split_core_active

mo_cas = split_core_active(one, two, ncore=4, nactive=6, e_core=nr)
one_cas = mo_cas.one
two_cas = mo_cas.two
```

---

## 3. Reading and Writing Data

PyBEST uses the `IOData` container class for all input/output. Two file formats are supported.

### 3.1 IOData container

```python
from pybest.iodata import IOData

# Create manually and set attributes
data = IOData()
data.kin    = kin
data.ne     = ne
data.eri    = eri
data.e_core = nr

# Or construct from keyword arguments
data = IOData(one=one, two=two, e_core=nr, nelec=10, ms2=0)

# Add/remove attributes on the fly
data.title  = "water molecule"
del data.title
```

Every method returns its result as an `IOData` container and simultaneously writes a checkpoint file to `pybest-results/checkpoint_<METHOD>.h5`.

### 3.2 Internal HDF5 format

Full binary precision; can store integrals in the AO or MO basis.

```python
# Write
data.to_file("hamiltonian_ao.h5")

# Read
data = IOData.from_file("hamiltonian_ao.h5")
print(data.__dict__)
```

The HDF5 layout for the AO Hamiltonian above:

```
HDF5 "hamiltonian_ao.h5" {
  /e_core    (scalar)
  /eri/array
  /kin/array
  /lf
  /ne/array
}
```

### 3.3 FCIDUMP format (ASCII, MO basis)

FCIDUMP is the standard interchange format for active-space codes (Molpro, Budapest DMRG, etc.). PyBEST supports reading and writing; integrals must already be in the MO basis.

**Writing:**

```python
data = IOData(one=one, two=two, e_core=nr, nelec=20, ms2=0)
data.to_file("hamiltonian_mo.FCIDUMP")
```

The file begins with a `&FCI ... &END` header (NORB, NELEC, MS2, ORBSYM, ISYM) followed by two-electron integrals in chemists' notation `(ij|kl)`, then one-electron integrals (last two indices = 0), then the core energy (all indices = 0). PyBEST currently assigns all orbitals symmetry label 1 (C₁ symmetry only).

**Reading:**

```python
data = IOData.from_file("hamiltonian_mo.FCIDUMP")
# data.one   — combined one-electron integrals
# data.two   — two-electron integrals
# data.e_core, data.nelec, data.ms2, data.orb_a, data.olp
```

### 3.4 Checkpoint restart

All CC/EOM modules support warm-start from a checkpoint file:

```python
ccsd = RCCSD(lf, occ_model)
ccsd_output = ccsd(kin, ne, eri, hf_output,
                   restart="pybest-results/checkpoint_RCCSD.h5")

# Or use an h5 path as the initial-guess amplitudes
ccsd_output = ccsd(kin, ne, eri, hf_output,
                   initguess="pybest-results/checkpoint_RCCSD.h5")
```

---

## 4. Strongly Correlated Methods

The pCCD-based hierarchy is PyBEST's primary strength. The key idea is that pair coupled cluster doubles (pCCD) captures static correlation efficiently and cheaply, while post-pCCD corrections recover dynamic correlation. Orbital optimization (OO-pCCD) further improves the reference by rotating orbitals to minimize the pCCD energy, which is essential for bond-breaking and strongly correlated ground states.

### 4.1 Orbital-Optimized pCCD (OO-pCCD) for Bond Breaking

`ROOpCCD` variationally optimizes both the pCCD amplitudes and the orbital rotation parameters simultaneously. This makes it well-suited for describing potential energy curves, bond dissociation, and other processes dominated by static/strong correlation.

```python
from pybest.geminals import ROOpCCD
from pybest.localization import PipekMezey
from pybest.part import get_mulliken_operators

# (Optional but recommended) Localize orbitals before pCCD
# to improve convergence on potential energy surfaces
mulliken = get_mulliken_operators(basis)
loc = PipekMezey(lf, occ_model, mulliken)
loc(orb_a, "occ")
loc(orb_a, "virt")

# Run OO-pCCD
pccd = ROOpCCD(lf, occ_model)
pccd_output = pccd(hf_output, kin, ne, eri)
# pccd_output.e_tot    — total OO-pCCD energy
# pccd_output.t_p      — pair amplitudes
# pccd_output.orb_a    — optimized orbitals
# pccd_output.dm1      — 1-RDM (needed for entanglement analysis)
# pccd_output.dm2      — 2-RDM
```

The `RpCCD` class (without orbital optimization) is available when optimized orbitals are not required, but does not produce response density matrices and therefore cannot be used for entanglement analysis.

**Frozen core** — applies identically to all methods:

```python
occ_model = AufbauOccModel(basis, ncore=1)   # freeze 1 core orbital
```

OO-pCCD then uses CCSD or fpCCSD on top for dynamic correlation (see below).

### 4.2 Frozen-Pair Coupled Cluster (fpCC)

In frozen-pair CC (fpCC), the pCCD pair amplitudes are kept fixed ("frozen") while the remaining single and double amplitudes are optimized. This approach efficiently corrects for dynamic correlation on top of the pCCD reference without re-optimizing the strongly correlated electron pairs.

The `RfpCCSD` (or `RfpCCD` for doubles only) class requires an `ROOpCCD` output as its reference:

```python
from pybest.cc import RfpCCD, RfpCCSD

# First: OO-pCCD
pccd = ROOpCCD(lf, occ_model)
pccd_output = pccd(hf_output, kin, ne, eri, e_core=nr)

# Then: fpCCSD on top of OO-pCCD
fpccsd = RfpCCSD(lf, occ_model)
fpccsd_output = fpccsd(kin, ne, eri, pccd_output)
# fpccsd_output.e_tot     — total fpCCSD energy
# fpccsd_output.e_corr    — fpCCSD correlation energy
# fpccsd_output.t_1, t_2  — single and double amplitudes

# Doubles-only variant
fpccd = RfpCCD(lf, occ_model)
fpccd_output = fpccd(kin, ne, eri, pccd_output)
```

Checkpoint restart:

```python
fpccsd_output = fpccsd(kin, ne, eri, pccd_output,
                        restart="pybest-results/checkpoint_RfpCCSD.h5")
```

Key keyword arguments: `initguess` ("random", "mp2", "const", or path to `.h5`), `maxiter` (default 100), `solver` (default "krylov"), `threshold` (default 1e-8).

> **Reference:** Leszczyk et al., *J. Chem. Theory Comput.* (2022). Cite as `[leszczyk2022]` in PyBEST literature list.

### 4.3 DMRG-Tailored Coupled Cluster (DMRG-tCCSD)

DMRG-tailored CC (tCC) embeds DMRG amplitudes into a CCSD wave function. The strongly correlated active-space amplitudes (from DMRG or any CASCI-type solver) are fixed ("tailored") into the CC cluster operator, while the remaining amplitudes are solved perturbatively. PyBEST interfaces with external DMRG codes (e.g., the QC-DMRG-Budapest program) through the FCIDUMP format.

**Workflow (two steps):**

**Step 0 — Generate integrals:**

```python
from pybest.utility import transform_integrals, split_core_active
from pybest.iodata import IOData

# Transform to MO basis (OO-pCCD or RHF orbitals)
mo_ints = transform_integrals(kin, ne, eri, pccd_output.orb_a)
one = mo_ints.one[0]
two = mo_ints.two[0]

# Write full Hamiltonian (for tCC) and CAS Hamiltonian (for DMRG)
IOData(one=one, two=two, e_core=nr, ms2=0, nelec=14).to_file("all.FCIDUMP")

mo_cas = split_core_active(one, two, ncore=2, nactive=8, e_core=nr)
IOData(one=mo_cas.one, two=mo_cas.two, e_core=nr, ms2=0, nelec=8).to_file("cas.FCIDUMP")
```

Run the external DMRG code on `cas.FCIDUMP`; extract CC amplitudes into a file called `T_DUMP`.

**Step 1 — DMRG-tCCSD in PyBEST:**

```python
from pybest.cc import RtCCSD
from pybest.iodata import IOData

ints  = IOData.from_file("all.FCIDUMP")
pccd  = IOData.from_file("pybest-results/checkpoint_pccd.h5")

tcc = RtCCSD(lf, occ_model)
tcc_output = tcc(
    ints.one,
    ints.two,
    ints.orb_a,
    e_ref=pccd.e_ref,
    external_file="T_DUMP",   # CC amplitudes extracted from DMRG MPS
)
# tcc_output.e_tot, e_corr, t_1, t_2, converged
```

The `external_file` keyword points to the pre-computed CAS amplitudes. The calculation can be performed either in the canonical RHF basis or in the OO-pCCD orbital basis (just substitute `hf_output.orb_a` with `pccd_output.orb_a` during `transform_integrals`).

Key keyword arguments: `external_file`, `e_ref`, `initguess`, `maxiter`, `solver`, `threshold_r` (default 1e-5), `threshold_e` (default 1e-6).

> **Reference:** Leszczyk et al., *J. Chem. Theory Comput.* (2022). Cite as `[leszczyk2022]`.

### 4.4 Equation-of-Motion pCCD (EOM-pCCD)

The EOM formalism in PyBEST targets electronically excited states. The pCCD-based EOM variants are unique to PyBEST and designed for systems with strong static correlation in the ground state.

#### EOM-pCCD — electron-pair excitations

```python
from pybest.eomcc import REOMpCCD

eom = REOMpCCD(lf, occ_model)
eom_output = eom(kin, ne, eri, pccd_output, nroot=3)
# eom_output.e_ee    — array of excitation energies (ground state = 0.0)
# eom_output.civ_ee  — eigenvectors (columns); first column = ground state
# eom_output.t_p     — pCCD pair amplitudes
```

`nroot` is the number of excited states **above** the ground state.

#### EOM-pCCD+S — pair and single excitations

```python
from pybest.eomcc import REOMpCCDS

eom = REOMpCCDS(lf, occ_model)
eom_output = eom(kin, ne, eri, pccd_output, nroot=3)
```

Note: EOM-pCCD+S is not size-intensive; excitation energies must be adjusted relative to the ground-state energy shift. For size-intensive results, use the `REOMpCCDCCS` variant with a pCCD-CCS reference.

#### EOM-pCCD-LCCSD — linearized CC correction on top of pCCD

```python
from pybest.cc import RpCCDLCCSD
from pybest.eomcc import REOMpCCDLCCSD

lccsd = RpCCDLCCSD(lf, occ_model)
lccsd_output = lccsd(kin, ne, eri, pccd_output)

eom = REOMpCCDLCCSD(lf, occ_model)
eom_output = eom(kin, ne, eri, lccsd_output, nroot=3)
```

#### Standard EOM flavors (RHF reference)

For molecules where a single-reference CC is sufficient:

```python
from pybest.cc import RCCSD
from pybest.eomcc import REOMCCSD

ccsd = RCCSD(lf, occ_model)
ccsd_output = ccsd(kin, ne, eri, hf_output)

eom = REOMCCSD(lf, occ_model)
eom_output = eom(kin, ne, eri, ccsd_output, nroot=3)
```

Equivalent classes exist for EOM-CCS (`REOMCCS`) and EOM-CCD (`REOMCCD`).

All EOM modules support frozen cores (set `ncore` in `AufbauOccModel`) and restarts via the `restart` keyword.

> **References:** Boguslawski et al., *Phys. Rev. A* (2016), *J. Chem. Phys.* (2017). Cite as `[boguslawski2016a]` and `[boguslawski2017c]`.

---

## 5. Orbital Correlation Visualization

PyBEST provides a complete workflow for computing and visualizing orbital entanglement, grounded in quantum information theory. The two central quantities are:

- **Single-orbital entropy** $s(1)_i = -\sum_\alpha \omega_{\alpha,i} \ln \omega_{\alpha,i}$ — measures how strongly orbital $i$ is entangled with all other orbitals. Calculated from the eigenvalues $\omega_{\alpha,i}$ of the one-orbital reduced density matrix (one-orbital RDM, dimension 4 for spatial orbitals).

- **Orbital-pair mutual information** $I_{i|j} = s(2)_{i,j} - s(1)_i - s(1)_j$ — quantifies the total correlation between orbital pair $(i,j)$. Calculated from the two-orbital entropy $s(2)_{i,j}$ (eigenvalues of the 16-dimensional two-orbital RDM).

An uncorrelated (single Slater determinant) wave function has zero entropy everywhere. Non-zero values signal static or dynamic correlation.

### 5.1 Computing entanglement from a pCCD wave function

`ROOpCCD` automatically computes the 1- and 2-RDMs. The `OrbitalEntanglementRpCCD` module converts them into orbital entropy and mutual information:

```python
from pybest.orbital_entanglement import OrbitalEntanglementRpCCD

pccd = ROOpCCD(lf, occ_model)
pccd_output = pccd(kin, ne, eri, hf_output)

entanglement = OrbitalEntanglementRpCCD(lf, pccd_output)
entanglement()
# Writes: pybest-results/s1-pccd.dat   (single-orbital entropy)
#         pybest-results/i12-pccd.dat  (orbital-pair mutual information)
```

### 5.2 Computing entanglement from pCCD-LCC wave functions

Requires solving the Λ-equations (`lambda_equations=True`):

```python
from pybest.cc import RpCCDLCCD, RpCCDLCCSD
from pybest.orbital_entanglement import OrbitalEntanglementRpCCDLCC

lccd = RpCCDLCCD(lf, occ_model)
lccd_output = lccd(kin, ne, eri, pccd_output, lambda_equations=True)

lccsd = RpCCDLCCSD(lf, occ_model)
lccsd_output = lccsd(kin, ne, eri, pccd_output, lambda_equations=True)

# Entanglement for pCCD-LCCD
entanglement = OrbitalEntanglementRpCCDLCC(lf, lccd_output)
entanglement()
# Writes: pybest-results/s1-pccd-lcc.dat, i12-pccd-lcc.dat

# Entanglement for pCCD-LCCSD (same class, different IOData input)
entanglement = OrbitalEntanglementRpCCDLCC(lf, lccsd_output)
entanglement()
```

### 5.3 Generating correlation diagrams

PyBEST ships the `pybest-entanglement.py` script (requires matplotlib ≥ 3.7). It reads the `.dat` files and produces two PDF figures: a bar chart of single-orbital entropies and an orbital-pair mutual information diagram, optionally decorated with MO pictures.

```bash
# Basic usage
pybest-entanglement.py \
    --threshold 0.001 \
    --iname pybest-results/i12-pccd-lcc.dat \
    --sname pybest-results/s1-pccd-lcc.dat

# Restrict the plot to the 15 most strongly correlated orbital pairs
pybest-entanglement.py \
    --threshold 0.001 \
    --largest-i 15 \
    --iname pybest-results/i12-pccd-lcc.dat \
    --sname pybest-results/s1-pccd-lcc.dat
```

Key options:

| Option | Description |
|--------|-------------|
| `--threshold` | Lower cutoff for mutual information (orders of magnitude: 0.001, 0.01, …) |
| `--iname` | Path to `i12-*.dat` file |
| `--sname` | Path to `s1-*.dat` file |
| `--indices` | Select specific orbital indices to display |
| `--order` | Reorder orbitals in the plot (1-based) |
| `--zoom` | Scaling factor for embedded MO pictures |
| `--largest-i` | Restrict diagram to at most N orbitals with largest mutual information |

Output: `pybest-results/s1.pdf` (single-orbital entropy) and `pybest-results/i12.pdf` (mutual information diagram). MO pictures must be stored as `mo_[index].png` in the working directory to appear in the plot.

> **References:** Boguslawski et al., *J. Chem. Theory Comput.* (2015), *Phys. Chem. Chem. Phys.* (2016), *J. Chem. Phys.* (2017). Cite as `[boguslawski2015a]`, `[boguslawski2016b]`, `[boguslawski2017b]`.

---

## 6. Citation

When publishing results obtained with PyBEST 2.1.0, please cite **both** of the following papers:

> Katharina Boguslawski, Filip Brzęk, Rahul Chakraborty, Kacper Cieślak, Seyedehdelaram Jahani, Aleksandra Leszczyk, Artur Nowak, Emil Sujkowski, Julian Świerczyński, Somayeh Ahmadkhani, Dariusz Kędziera, Maximilian H. Kriebel, Piotr Szymon Żuchowski, Paweł Tecmer,
> "PyBEST: Improved functionality and enhanced performance,"
> *Computer Physics Communications*, **297**, 109049 (2024).
> <https://doi.org/10.1016/j.cpc.2023.109049>

> Katharina Boguslawski, Aleksandra Leszczyk, Artur Nowak, Filip Brzęk, Piotr Szymon Żuchowski, Dariusz Kędziera, Paweł Tecmer,
> "Pythonic Black-box Electronic Structure Tool (PyBEST). An open-source Python platform for electronic structure calculations at the interface between chemistry and physics,"
> *Computer Physics Communications*, **264**, 107933 (2021).
> <https://doi.org/10.1016/j.cpc.2021.107933>

**Software attribution (version string):**

> Saman Behjou, Filip Brzęk, Rahul Chakraborty, Kacper Cieślak, Antonina Dobrowolska, Seyedehdelaram Jahani, Zahra Karimi, Michał Kopczyński, Aleksandra Leszczyk, Artur Nowak, Ram Dhari Pandey, Emil Sujkowski, Julian Świerczyński, Julia Szczuczko, Lena Szczuczko, Somayeh Ahmadkhani, Katharina Boguslawski, Iulia Brumboiu, Marta Gałyńska, Dariusz Kędziera, Maximilian Kriebel, Paweł Tecmer, Piotr Szymon Żuchowski, *PyBEST 2.1.0*, **2025**.

**Per-module references** (cite in addition to the above when using specific modules):

| Module | Reference key |
|--------|--------------|
| RCC / fpCC / tCC | `[leszczyk2022]` |
| EOM-pCCD, EOM-pCCD+S | `[boguslawski2016a]`, `[boguslawski2017c]` |
| EOM-pCCD-LCC | `[boguslawski2019]` |
| Orbital entanglement (seniority-zero) | `[boguslawski2015a]`, `[boguslawski2016b]`, `[boguslawski2017b]` |
| Orbital entanglement (pCCD-LCC) | `[nowak2021]` |

Full BibTeX entries are available in the PyBEST documentation under [Literature](http://fizyka.umk.pl/~pybest/pybest-v2.1.0/tech_ref_literature.html).
