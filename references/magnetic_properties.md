# Magnetic Coupling Constants

Magnetic coupling constants (exchange coupling constants J) describe the
magnetic interaction between spin centres in polynuclear metal complexes,
organic diradicals, and molecular magnets.

The Heisenberg-Dirac-Van Vleck (HDVV) Hamiltonian:
    Ĥ = −2J Ŝ_A · Ŝ_B

Positive J → ferromagnetic (parallel spins lower in energy)
Negative J → antiferromagnetic (antiparallel spins lower in energy)

---

## 1. Broken-symmetry DFT (BS-DFT) — standard approach

BS-DFT is the workhorse method for large TM dimers. It gives good J values
at DFT cost. Use a hybrid functional (B3LYP, TPSSh, or M06) with a flexible basis.

```python
from pyscf import gto, scf, dft
import numpy as np

def run_state(atom_str, spin, charge=0):
    """Return total energy for a given spin state."""
    mol = gto.Mole()
    mol.atom   = atom_str
    mol.basis  = "def2-TZVP"
    mol.charge = charge
    mol.spin   = spin        # 2S
    mol.verbose = 3
    mol.build()

    mf = dft.UKS(mol)
    mf.xc = "B3LYP"
    mf.grids.level = 4

    # Broken-symmetry: provide initial guess with localised spins
    if spin == 0:
        # Start from high-spin converged MOs
        dm_hs = dm_highspin   # set this from the high-spin run below
        mf.kernel(dm0=dm_hs)
    else:
        mf.kernel()
    return mf.e_tot, mf

# Step 1: high-spin (HS) state
mol_hs = "Cu 0 0 0; Cu 0 0 2.9; ..."
e_hs, mf_hs = run_state(mol_hs, spin=2)   # S=1 → spin=2S=2

# Step 2: broken-symmetry (BS) singlet
#   Use HS MOs as initial guess, then flip α↔β on site B
dm_alpha, dm_beta = mf_hs.make_rdm1()
dm_highspin = (dm_alpha, dm_beta)

# To flip spin on atom B: exchange α and β density blocks
# This step is system-specific — a convenient approach:
from pyscf.tools import wfn_format
# ... swap density matrices for the second metal site

e_bs, mf_bs = run_state(mol_hs, spin=0)

# Step 3: compute J via Yamaguchi formula (recommended over simple formula)
# <S²>_HS and <S²>_BS from UKS:
s2_hs = mf_hs.spin_square()[0]
s2_bs = mf_bs.spin_square()[0]

# Yamaguchi formula:
# J = (E_BS - E_HS) / (<S²>_HS - <S²>_BS)
J_cm1 = (e_bs - e_hs) / (s2_hs - s2_bs) * 219474.6  # Ha → cm⁻¹
print(f"J = {J_cm1:.1f} cm⁻¹")
print(f"    ({'ferromagnetic' if J_cm1 > 0 else 'antiferromagnetic'})")
```

**Three formulas for J — know which to use**:

| Formula | When | Expression |
|---------|------|-----------|
| Noodleman (weak overlap) | Fully localised spin centres | J = (E_BS − E_HS) / (2S_max²) |
| Yamaguchi (general) | **Default — use this** | J = (E_BS − E_HS) / (⟨S²⟩_HS − ⟨S²⟩_BS) |
| Bencini (strong overlap) | Delocalised cases | J = (E_BS − E_HS) / S_max² |

---

## 2. NEVPT2 on spin-state manifold — higher accuracy

For systems where BS-DFT gives inconsistent results (e.g., very short Cu–Cu
distances, radical ligands):

```python
from pyscf import gto, scf, mcscf, mrpt

mol = gto.Mole(atom="...", basis="def2-TZVP").build()
mf  = scf.RHF(mol).run()

# Active space: d-orbitals + magnetic orbitals of both metal centres
# For a Cu(II) dimer: (2e, 2o) is the minimal magnetic space
mc = mcscf.CASSCF(mf, ncas=2, nelecas=2)
mc.fcisolver.nroots = 3   # singlet, triplet (MS=0,+1,-1)
mc.run()

# NEVPT2 corrections
e_nevpt2 = mrpt.NEVPT2(mc).kernel()

# Extract energies of singlet and triplet
e_singlet = mc.e_states[0] + e_nevpt2[0]
e_triplet = mc.e_states[1] + e_nevpt2[1]

# J from energy splitting: E_T − E_S = −2J (for S_A = S_B = 1/2)
J_ha  = (e_singlet - e_triplet) / 2
J_cm1 = J_ha * 219474.6
print(f"J (NEVPT2) = {J_cm1:.1f} cm⁻¹")
```

---

## 3. Practical broken-symmetry setup with spin-flip

The trickiest part of BS-DFT is creating a proper broken-symmetry initial guess.
PySCF utility:

```python
from pyscf import scf

mf_hs = dft.UKS(mol).run()          # converge high-spin first
dm_a, dm_b = mf_hs.make_rdm1()

# Identify AO indices for the second spin centre (atom indices)
# and swap α↔β density blocks for those AOs:
# (this is geometry-dependent — use mol.aoslice_by_atom())
ao_slices = mol.aoslice_by_atom()
site_B_start = ao_slices[idx_metal_B][2]
site_B_end   = ao_slices[idx_metal_B][3]

dm_a_bs = dm_a.copy()
dm_b_bs = dm_b.copy()
dm_a_bs[site_B_start:site_B_end, :] = dm_b[site_B_start:site_B_end, :]
dm_b_bs[site_B_start:site_B_end, :] = dm_a[site_B_start:site_B_end, :]

mf_bs = dft.UKS(mol)
mf_bs.xc = "B3LYP"
mf_bs.kernel(dm0=(dm_a_bs, dm_b_bs))
```

---

## 4. Functional choice for magnetic coupling

| System | Recommended XC | Notes |
|--------|---------------|-------|
| Cu dimers | TPSSh or B3LYP | 10–25 % HF exchange |
| Mn/Fe dimers | M06 or TPSSh | Avoid pure GGAs |
| Organic diradicals | UB3LYP | Check ⟨S²⟩ ≈ 1 for BS |
| Strong exchange | CASSCF/NEVPT2 | DFT fails for very large J |

---

## 5. Interpretation

- J in cm⁻¹: typical TM dimers range from −500 to +100 cm⁻¹
- Convert: J [meV] = J [cm⁻¹] × 0.12398
- Check ⟨S²⟩_BS ≈ 1.0 (for two S=½ centres): confirms broken symmetry achieved
- If ⟨S²⟩_BS ≈ ⟨S²⟩_HS, the BS guess collapsed to HS — increase initial perturbation
- Large |J| with BS-DFT but small with CASSCF → DFT artefact, use wavefunction method
