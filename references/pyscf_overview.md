# PySCF Reference Guide

PySCF (Python-based Simulations of Chemistry Framework) is an open-source quantum chemistry package written in Python with performance-critical components in C. It covers mean-field through multireference methods and periodic systems.

---

## Installation

```bash
pip install pyscf
# Optional extras
pip install pyscf[geomopt]   # geometry optimization
pip install pyscf[dftd3]     # DFT-D3 dispersion
pip install gpu4pyscf        # GPU-accelerated variant (see GPU section)
```

---

## Molecule Construction

```python
from pyscf import gto

mol = gto.Mole()
mol.atom = '''
    O  0.000  0.000  0.117
    H  0.000  0.757 -0.469
    H  0.000 -0.757 -0.469
'''
mol.basis  = 'cc-pVDZ'
mol.charge = 0
mol.spin   = 0          # 2S (number of unpaired electrons)
mol.symmetry = True     # use molecular symmetry
mol.build()
```

### Alternative atom input formats

```python
# Z-matrix string
mol.atom = 'O; H 1 0.96; H 1 0.96 2 104.5'

# List of tuples
mol.atom = [('O', (0, 0, 0.117)), ('H', (0, 0.757, -0.469)), ('H', (0, -0.757, -0.469))]

# From XYZ file
mol = gto.M(atom=open('water.xyz').read(), basis='def2-TZVP')
```

### Common keywords

| Keyword | Description |
|---|---|
| `basis` | Basis set string, dict, or custom dict |
| `charge` | Molecular charge |
| `spin` | 2S (0 = singlet) |
| `symmetry` | Enable point-group symmetry |
| `unit` | `'Angstrom'` (default) or `'Bohr'` |
| `verbose` | Output verbosity (0–9) |
| `output` | Redirect output to file path |
| `max_memory` | Memory limit in MB |

### Custom / mixed basis sets

```python
mol.basis = {
    'O': 'aug-cc-pVTZ',
    'H': 'cc-pVDZ',
}

# Inline basis definition (NWChem format)
mol.basis = {'H': gto.parse('''
    H S
      13.0107010  0.01968
       1.9622572  0.13797
       0.4445680  0.47814
       0.1219492  0.50124
''')}
```

---

## Hartree–Fock: RHF / UHF / ROHF

### Restricted HF (closed-shell)

```python
from pyscf import scf

mf = scf.RHF(mol)
mf.max_cycle = 200
mf.conv_tol  = 1e-10
mf.kernel()                  # run SCF

print('RHF energy:', mf.e_tot)
print('HOMO/LUMO index:', mf.mo_occ)
```

### Unrestricted HF (open-shell, broken-symmetry)

```python
mf = scf.UHF(mol)
mf.kernel()

dm_alpha, dm_beta = mf.make_rdm1()
s2, _ = mf.spin_square()
print('<S^2> =', s2)
```

### Restricted open-shell HF

```python
mf = scf.ROHF(mol)
mf.kernel()
```

### Convergence aids

```python
mf = scf.RHF(mol)
mf = scf.addons.remove_linear_dep_(mf)   # remove near-linear dependencies

# DIIS variants
mf.diis = scf.EDIIS()   # or ADIIS(), CDIIS() (default)

# Fermi smearing (stability for metals)
mf = scf.addons.smearing_(mf, sigma=0.1, method='fermi')

# Initial guess options: 'minao', 'atom', 'hcore', '1e', 'chkfile'
mf.init_guess = 'minao'

# Use checkpoint
mf.chkfile = 'water.chk'
mf.init_guess = 'chkfile'
```

### Stability analysis

```python
mo_new = mf.stability()[0]
mf.kernel(mf.make_rdm1(mo_new, mf.mo_occ))
```

---

## DFT (Density Functional Theory)

```python
from pyscf import dft

mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf.grids.level = 4      # integration grid quality (0–9)
mf.kernel()
```

### UKS and ROKS

```python
mf = dft.UKS(mol)
mf.xc = 'PBE'
mf.kernel()

mf = dft.ROKS(mol)
mf.xc = 'M06-2X'
mf.kernel()
```

### Common functionals

| Category | Examples |
|---|---|
| LDA | `'LDA,VWN'`, `'SVWN'` |
| GGA | `'PBE'`, `'BP86'`, `'BLYP'` |
| meta-GGA | `'TPSS'`, `'M06-L'`, `'SCAN'` |
| Hybrid | `'B3LYP'`, `'PBE0'`, `'M06-2X'`, `'HSE06'` |
| Double hybrid | `'B2PLYP'` (via `dft` + `mp2` correction) |
| Range-separated | `'CAM-B3LYP'`, `'wB97X-D'`, `'LC-wHPBE'` |

PySCF uses libxc under the hood; any libxc functional can be used by its name or numeric ID.

### Dispersion corrections

```python
# DFT-D3 (requires pyscf[dftd3] or simple-dftd3)
from pyscf import dftd3
mf = dftd3.dftd3(dft.RKS(mol))
mf.xc = 'B3LYP'
mf.kernel()
```

### Grid customization

```python
mf.grids.atom_grid = {'O': (99, 590), 'H': (75, 302)}  # (radial, angular)
mf.grids.prune = dft.gen_grid.treutler_prune
```

### NLC (non-local correlation, VV10)

```python
mf.xc = 'B97M-V'   # includes NLC automatically through libxc
mf.nlc = 'VV10'
```

---

## TDDFT (Time-Dependent DFT)

```python
from pyscf import tdscf

# Run ground-state DFT first
mf = dft.RKS(mol).run()

# TDDFT (Casida linear response)
td = tdscf.TDDFT(mf)
td.nstates = 10          # number of excited states
td.kernel()

# Print excitation energies and oscillator strengths
for i, (e, f) in enumerate(zip(td.e, td.oscillator_strength())):
    print(f'State {i+1}: {e*27.211:.3f} eV, f = {f:.4f}')
```

### TDA (Tamm–Dancoff Approximation)

```python
td = tdscf.TDA(mf)
td.nstates = 5
td.kernel()
```

### UHF/UKS TDDFT

```python
mf = dft.UKS(mol).run()
td = tdscf.TDDFT(mf)      # automatically dispatches to UTD class
td.kernel()
```

### Transition dipoles and natural transition orbitals

```python
# Transition dipole moments (already computed during oscillator_strength)
dipoles = td.transition_dipole()

# Natural transition orbitals (NTOs)
weights, ntos = td.get_nto(state=1, verbose=4)
```

### Excited-state gradients

```python
from pyscf.grad import tdrhf as tdgrad

grad = tdgrad.Gradients(td)
grad.state = 1        # 1-indexed
g = grad.kernel()
```

---

## CCSD (Coupled Cluster Singles and Doubles)

```python
from pyscf import cc

# Start from RHF
mf = scf.RHF(mol).run()

mycc = cc.CCSD(mf)
mycc.conv_tol    = 1e-8
mycc.max_cycle   = 300
mycc.kernel()

print('CCSD correlation energy:', mycc.e_corr)
print('CCSD total energy:',       mycc.e_tot)

# Perturbative triples correction
et = mycc.ccsd_t()
print('CCSD(T) total energy:', mycc.e_tot + et)
```

### UCCSD (unrestricted)

```python
mf = scf.UHF(mol).run()
mycc = cc.UCCSD(mf)
mycc.kernel()
et = mycc.ccsd_t()
```

### Density matrices and properties

```python
dm1 = mycc.make_rdm1()           # 1-RDM in MO basis
dm2 = mycc.make_rdm2()           # 2-RDM in MO basis
dm1_ao = mycc.make_rdm1(ao_repr=True)  # 1-RDM in AO basis

# Natural orbital occupation numbers
import numpy as np
from pyscf import lo
noons, natorbs = np.linalg.eigh(dm1)
```

### CCSD gradients

```python
from pyscf.grad import ccsd as ccsd_grad
g = ccsd_grad.Gradients(mycc).kernel()
```

---

## EOM-CCSD (Equation-of-Motion CCSD)

EOM-CCSD gives ionization potentials (IP), electron affinities (EA), and excitation energies (EE) directly from a CCSD reference.

### EOM-IP-CCSD (ionization potentials)

```python
mycc = cc.CCSD(mf).run()

myip = mycc.ipccsd(nroots=5)
# Returns arrays: e (energies in Hartree), v (eigenvectors)
e_ip, v_ip = myip
for i, e in enumerate(e_ip):
    print(f'IP {i+1}: {-e*27.211:.3f} eV')
```

### EOM-EA-CCSD (electron affinities)

```python
e_ea, v_ea = mycc.eaccsd(nroots=3)
for i, e in enumerate(e_ea):
    print(f'EA {i+1}: {e*27.211:.3f} eV')
```

### EOM-EE-CCSD (excitation energies)

```python
e_ee, v_ee = mycc.eeccsd(nroots=5)
for i, e in enumerate(e_ee):
    print(f'EE {i+1}: {e*27.211:.3f} eV')
```

### EOM-CCSD with spin-flip

```python
# Spin-flip: access singlet excited states from triplet reference
from pyscf.cc import eom_rccsd
mycc = cc.CCSD(mf).run()
e_sf, v_sf = mycc.eomee_ccsd_singlet(nroots=3)
```

### EOM with density matrices

```python
# 1-RDM for a specific EOM root
dm1 = mycc.eomip_method().make_rdm1(myip[1][0])
```

---

## CASCI / CASSCF (Multireference Methods)

### CASCI (no orbital optimization)

```python
from pyscf import mcscf

# (n_electrons, n_orbitals) active space
nelec_act = 6
norb_act  = 6

mc = mcscf.CASCI(mf, norb_act, nelec_act)
mc.kernel()
print('CASCI energy:', mc.e_tot)
```

### CASSCF (with orbital optimization)

```python
mc = mcscf.CASSCF(mf, norb_act, nelec_act)
mc.max_cycle_macro = 100
mc.conv_tol = 1e-8
mc.kernel()
print('CASSCF energy:', mc.e_tot)
```

### State-averaged CASSCF

```python
mc = mcscf.CASSCF(mf, 6, 6)
# Average over 3 states with equal weights
mc.state_average_([1/3, 1/3, 1/3])
mc.kernel()
```

### Specifying active space by orbital selection

```python
# Sort orbitals by energy, freeze core, etc.
from pyscf.mcscf import avas

# AVAS: Atomic Valence Active Space
norb, nelec, mo = avas.avas(mf, ['C 2p', 'O 2p'])
mc = mcscf.CASSCF(mf, norb, nelec)
mc.kernel(mo)
```

### NEVPT2 (N-electron valence PT2 on top of CASSCF)

```python
from pyscf import mrpt

mc = mcscf.CASSCF(mf, 6, 6).run()
e_nevpt2 = mrpt.NEVPT2(mc).kernel()
print('NEVPT2 energy:', mc.e_tot + e_nevpt2)
```

### CASSCF gradients

```python
from pyscf.grad import casscf as casscf_grad
g = casscf_grad.Gradients(mc).kernel()
```

### DMRG / FCI solvers (external)

```python
# Plug in Block2 or CheMPS2 as FCI solver for large active spaces
from pyscf import dmrgscf
mc = dmrgscf.DMRGSCF(mf, norb_act, nelec_act)
mc.fcisolver.maxM = 1000   # bond dimension
mc.kernel()
```

---

## Periodic Systems (PBC / Solid-State)

PySCF's `pbc` subpackage handles periodic boundary conditions for molecules, slabs, and 3D crystals.

### Cell construction

```python
from pyscf.pbc import gto as pbcgto, scf as pbcscf, dft as pbcdft

cell = pbcgto.Cell()
cell.atom = '''
    Si  0.000  0.000  0.000
    Si  1.357  1.357  1.357
'''
cell.a = [[2.715, 2.715, 0],   # lattice vectors in Angstrom
           [2.715, 0,     2.715],
           [0,     2.715, 2.715]]
cell.basis  = 'gth-dzvp'       # Goedecker–Teter–Hutter basis
cell.pseudo = 'gth-pade'       # GTH pseudopotential (LDA)
cell.ke_cutoff = 200           # kinetic energy cutoff in Ry
cell.verbose = 5
cell.build()
```

### k-point sampling

```python
from pyscf.pbc.tools import pyscf_ase
import numpy as np

# Monkhorst–Pack grid
kpts = cell.make_kpts([4, 4, 4])

mf = pbcscf.KRHF(cell, kpts)
mf.kernel()
```

### Gamma-point HF/DFT

```python
mf = pbcscf.RHF(cell)    # gamma-point only
mf.kernel()

mf = pbcdft.RKS(cell)
mf.xc = 'PBE'
mf.kernel()
```

### Band structure

```python
from pyscf.pbc.dft import gen_band

# Define k-path (fractional coordinates)
kpath = cell.get_band_kpts([[0,0,0], [0.5,0,0], [0.5,0.5,0]])
e_kn = mf.get_bands(kpath)[0]

# e_kn shape: (nkpts, nbands)
```

### Density of states

```python
from pyscf.pbc.tools import dos

e, dos_vals = dos.get_dos(mf, emin=-10, emax=10, npts=500, sigma=0.05)
```

### Periodic MP2 and CCSD

```python
from pyscf.pbc import mp, cc as pbccc

mf = pbcscf.KRHF(cell, kpts).run()

mymp = mp.KMP2(mf)
mymp.kernel()

mycc = pbccc.KCCSD(mf)
mycc.kernel()
```

---

## Relativistic Hamiltonians

### Scalar relativistic: ZORA / DKH

```python
# Douglas–Kroll–Hess (DKH) via Hamiltonian keyword
mf = scf.RHF(mol)
mf.with_x2c()          # eXact 2-Component (X2C) scalar relativity
mf.kernel()

# Alternatively set on the mol object
mol.basis = 'x2c-TZVPall'   # X2C-contracted basis
mf = scf.RHF(mol).x2c()
mf.kernel()
```

### 4-Component Dirac–Fock (DHF)

```python
from pyscf import scf as scf4c

mol4c = gto.M(
    atom='Au 0 0 0',
    basis='dyall-dz',    # Dyall's relativistic basis sets
    charge=0,
    spin=1,
)

mf = scf4c.DHF(mol4c)
mf.kernel()
```

### X2C for DFT and post-HF

```python
# Any method can be X2C-decorated
mf = dft.RKS(mol).x2c()
mf.xc = 'PBE0'
mf.kernel()

mycc = cc.CCSD(mf).run()   # uses X2C integrals automatically
```

### Spin–orbit coupling (2-component)

```python
from pyscf import x2c

# 2-component GHF with spin–orbit
mf = scf.GHF(mol).x2c1e()
mf.kernel()
```

---

## GPU Support

GPU acceleration is provided via the **gpu4pyscf** package, which reimplements core PySCF routines using CuPy (CUDA).

### Installation

```bash
pip install gpu4pyscf-cuda12x   # for CUDA 12.x
# or
pip install gpu4pyscf-cuda11x   # for CUDA 11.x
```

### Drop-in replacement pattern

```python
import gpu4pyscf

# Convert an existing PySCF object to GPU with .to_gpu()
from pyscf import dft
mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf = mf.to_gpu()     # offload everything to GPU
mf.kernel()
```

### Direct GPU construction

```python
from gpu4pyscf import dft as gpu_dft

mf = gpu_dft.RKS(mol)
mf.xc = 'wB97X-D'
mf.grids.level = 5
mf.kernel()

# Convert result back to CPU for post-processing
e_tot = float(mf.e_tot)
dm = mf.make_rdm1().get()   # .get() copies CuPy array to NumPy
```

### Supported methods (as of 2024–2025)

| Method | GPU support |
|---|---|
| RHF, UHF, GHF | Yes |
| RKS, UKS | Yes |
| MP2 | Yes |
| CCSD, CCSD(T) | Yes (RCCSD) |
| TDDFT/TDA | Yes |
| Gradients, Hessians | Yes (DFT) |
| CASSCF | Partial |
| Periodic (PBC) | Partial |

### Multi-GPU

```python
import cupy as cp
cp.cuda.Device(0).use()   # select GPU device before building

# For multi-GPU: use mpi4py or joblib to distribute k-points / states
```

---

## Exporting Molden and Cube Files

### Molden file (orbitals, basis, geometry)

```python
from pyscf.tools import molden

# Write Molden file from an SCF calculation
with open('water.molden', 'w') as f:
    molden.header(mol, f)
    molden.orbital_coeff(mol, f, mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)

# Convenience one-liner
molden.from_mo(mol, 'water.molden', mf.mo_coeff)
```

### Cube file (volumetric electron density / orbital)

```python
from pyscf.tools import cubegen

# Electron density
cubegen.density(mol, 'density.cube', mf.make_rdm1())

# Specific molecular orbital (0-indexed)
cubegen.orbital(mol, 'homo.cube', mf.mo_coeff[:, homo_idx])

# Electrostatic potential (ESP)
cubegen.mep(mol, 'esp.cube', mf.make_rdm1())
```

### Cube file options

```python
# Control grid resolution
cubegen.density(
    mol, 'density.cube', dm,
    nx=80, ny=80, nz=80,     # number of grid points per axis
    margin=3.0,               # extra margin around molecule (Bohr)
)
```

### Natural orbitals to Molden

```python
import numpy as np

# CASSCF natural orbitals
noons = mc.mo_occ
natorbs = mc.mo_coeff

with open('casscf_natorbs.molden', 'w') as f:
    molden.header(mol, f)
    molden.orbital_coeff(mol, f, natorbs, ene=noons, occ=noons)
```

### Cube from periodic calculation

```python
from pyscf.pbc.tools import cubegen as pbc_cubegen

pbc_cubegen.density(cell, 'density.cube', mf.make_rdm1())
```

---

## Useful Utilities

### Population analysis

```python
# Mulliken charges
mf.mulliken_pop()
mf.mulliken_meta()   # meta-Lowdin

# Natural Population Analysis (NPA)
from pyscf.lo import nbo
nbo.NBO(mf).kernel()
```

### Geometry optimization (via geomeTRIC)

```python
from pyscf.geomopt.geometric_solver import optimize

mol_opt = optimize(mf)   # returns optimized Mole object
```

### Frequency analysis (Hessian)

```python
from pyscf.hessian import rhf as rhf_hess

h = rhf_hess.Hessian(mf)
hess = h.kernel()

from pyscf.prop.freq import rhf as freq
freq.Frequency(mf).kernel()
```

### Chkfile / restart

```python
from pyscf import lib

# Save
mf.chkfile = 'calc.chk'
mf.dump_chk()

# Load MO coefficients
mo_coeff = lib.chkfile.load('calc.chk', 'scf/mo_coeff')
mo_energy = lib.chkfile.load('calc.chk', 'scf/mo_energy')
```

---

## Citations

When publishing work that uses PySCF, cite the following primary references:

### PySCF main papers

```
Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu,
J. D. McClain, E. R. Sayfutyarova, S. Sharma, S. Wouters, G. K.-L. Chan,
"PySCF: the Python-based simulations of chemistry framework",
WIREs Comput. Mol. Sci. 2018, 8, e1340.
DOI: 10.1002/wcms.1340

Q. Sun, X. Zhang, S. Banerjee, P. Bao, M. Barbry, N. S. Blunt, N. A. Bogdanov,
G. H. Booth, J. Chen, Z.-H. Cui, J. J. Eriksen, Y. Gao, S. Guo, J. Hermann,
M. R. Hermes, K. Koh, P. Koval, S. Lehtola, Z. Li, J. Liu, N. Mardirossian,
J. D. McClain, M. Motta, B. Mussard, H. Q. Pham, A. Pulkin, W. Purwanto,
P. J. Robinson, E. Ronca, E. R. Sayfutyarova, M. Scheurer, H. F. Schurkus,
J. E. T. Smith, C. Sun, S.-N. Sun, S. Upadhyay, L. K. Wagner, X. Wang,
A. White, J. D. Whitfield, M. J. Williamson, S. Wouters, J. Yang, J. M. Yu,
T. Zhu, T. C. Berkelbach, S. Sharma, A. Y. Sokolov, G. K.-L. Chan,
"Recent developments in the PySCF program package",
J. Chem. Phys. 2020, 153, 024109.
DOI: 10.1063/5.0006074
```

### gpu4pyscf

```
X. Zhang, R. Zhang, J. Li, A. Bhati, A. Nandi, J. J. Ryu, J. Christodouleas,
"Accelerating quantum chemistry with vectorized and GPU-accelerated
implementations of integral evaluation and mean-field methods on modern hardware",
arXiv:2407.09700 (2024).
```

### Method-specific citations

| Method | Key reference |
|---|---|
| X2C scalar relativity | Liu, W. et al., *J. Chem. Phys.* **2009**, 131, 031104 |
| NEVPT2 | Angeli, C. et al., *J. Chem. Phys.* **2001**, 114, 10252 |
| AVAS | Sayfutyarova, E. R. et al., *J. Chem. Theory Comput.* **2017**, 13, 4063 |
| Density fitting (DF-HF/DFT) | Weigend, F. *Phys. Chem. Chem. Phys.* **2002**, 4, 4285 |
| DFT-D3 | Grimme, S. et al., *J. Chem. Phys.* **2010**, 132, 154104 |
| KMP2 / KCCSD (periodic) | McClain, J. D. et al., *J. Chem. Theory Comput.* **2017**, 13, 1209 |

### BibTeX entries

```bibtex
@article{Sun2018,
  author  = {Q. Sun and T. C. Berkelbach and N. S. Blunt and others},
  title   = {{PySCF}: the Python-based simulations of chemistry framework},
  journal = {WIREs Comput. Mol. Sci.},
  year    = {2018},
  volume  = {8},
  pages   = {e1340},
  doi     = {10.1002/wcms.1340},
}

@article{Sun2020,
  author  = {Q. Sun and X. Zhang and S. Banerjee and others},
  title   = {Recent developments in the {PySCF} program package},
  journal = {J. Chem. Phys.},
  year    = {2020},
  volume  = {153},
  pages   = {024109},
  doi     = {10.1063/5.0006074},
}
```

---

## Quick Reference: Method Hierarchy

```
gto.Mole / pbc.gto.Cell
│
├── scf.RHF / UHF / ROHF / GHF
│   ├── dft.RKS / UKS / ROKS          (xc= functional)
│   │   ├── tdscf.TDDFT / TDA
│   │   └── pbc.dft.KRKS              (periodic)
│   ├── cc.CCSD → .ccsd_t()
│   │   └── .ipccsd / .eaccsd / .eeccsd   (EOM)
│   ├── mcscf.CASCI
│   ├── mcscf.CASSCF → mrpt.NEVPT2
│   └── pbc.scf.KRHF                  (periodic, k-points)
│       └── pbc.cc.KCCSD
│
└── .x2c() / .to_gpu()                (decorators)
```
