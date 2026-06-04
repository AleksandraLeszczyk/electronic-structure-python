# Quantum Chemistry Basics

A reference glossary of core quantum chemistry concepts.

---

## Hamiltonian

The **Hamiltonian** ($\hat{H}$) is the quantum mechanical operator corresponding to the total energy of a system. For a molecular system it contains kinetic energy terms for all electrons and nuclei, plus potential energy terms for electron–nucleus attraction and electron–electron and nucleus–nucleus repulsion:

$$\hat{H} = -\sum_i \frac{\hbar^2}{2m_e}\nabla_i^2 - \sum_A \frac{\hbar^2}{2M_A}\nabla_A^2 - \sum_{i,A}\frac{Z_A e^2}{r_{iA}} + \sum_{i<j}\frac{e^2}{r_{ij}} + \sum_{A<B}\frac{Z_A Z_B e^2}{R_{AB}}$$

Under the **Born–Oppenheimer approximation**, nuclear kinetic energy is neglected and the electronic Hamiltonian $\hat{H}_{\text{el}}$ is solved for fixed nuclear positions. The eigenvalue equation $\hat{H}\Psi = E\Psi$ is the time-independent Schrödinger equation.

---

## Electronic Wave Function

The **electronic wave function** $\Psi(\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_N)$ completely describes the quantum state of $N$ electrons in a molecule (at fixed nuclear geometry). Each coordinate $\mathbf{x}_i = (\mathbf{r}_i, \sigma_i)$ encodes both the spatial position $\mathbf{r}_i$ and the spin $\sigma_i$ of electron $i$.

The square of its modulus, $|\Psi|^2$, gives the probability density for finding electrons at the specified positions with the specified spins. All observable properties of the electronic structure can be derived from $\Psi$.

---

## Slater Determinant

A **Slater determinant** is an antisymmetrized product of $N$ one-electron spin orbitals $\{\chi_i\}$, used as the wave function ansatz in Hartree–Fock theory:

$$\Psi_{\text{SD}} = \frac{1}{\sqrt{N!}} \begin{vmatrix} \chi_1(\mathbf{x}_1) & \chi_2(\mathbf{x}_1) & \cdots & \chi_N(\mathbf{x}_1) \\ \chi_1(\mathbf{x}_2) & \chi_2(\mathbf{x}_2) & \cdots & \chi_N(\mathbf{x}_2) \\ \vdots & & \ddots & \vdots \\ \chi_1(\mathbf{x}_N) & \chi_N(\mathbf{x}_N) & \cdots & \chi_N(\mathbf{x}_N) \end{vmatrix}$$

Swapping any two rows (i.e., exchanging two electrons) changes the sign of the determinant, automatically satisfying the **Pauli antisymmetry principle**. A single Slater determinant is exact only for non-interacting electrons; correlated methods build on multiple determinants.

---

## Molecular Orbitals

**Molecular orbitals (MOs)** are one-electron spatial wave functions $\phi_i(\mathbf{r})$ that extend over an entire molecule. They are typically constructed as linear combinations of atomic orbitals (**LCAO-MO**):

$$\phi_i(\mathbf{r}) = \sum_\mu c_{\mu i}\, \chi_\mu(\mathbf{r})$$

where $\chi_\mu$ are basis functions and $c_{\mu i}$ are the MO coefficients determined by solving the Hartree–Fock (Roothaan–Hall) equations. MOs are classified as bonding, antibonding, or non-bonding based on their effect on molecular stability, and as $\sigma$, $\pi$, or $\delta$ based on their symmetry with respect to the internuclear axis.

---

## Spin Orbitals

A **spin orbital** $\chi_i(\mathbf{x})$ is the product of a spatial molecular orbital and a one-electron spin function:

$$\chi_i(\mathbf{x}) = \phi_i(\mathbf{r})\,\sigma(\omega), \quad \sigma \in \{\alpha, \beta\}$$

where $\alpha$ (spin-up, $m_s = +\tfrac{1}{2}$) and $\beta$ (spin-down, $m_s = -\tfrac{1}{2}$) are the two orthonormal spin eigenfunctions. Spin orbitals are the fundamental building blocks of Slater determinants; $N$ electrons occupy $N$ distinct spin orbitals (Pauli exclusion).

---

## Basis Functions

**Basis functions** $\{\chi_\mu\}$ are a finite set of known analytic functions used to expand molecular orbitals. The two main types are:

- **Slater-type orbitals (STOs)**: $\propto r^{n-1} e^{-\zeta r} Y_l^m(\theta,\phi)$ — physically motivated, but multi-center integrals are costly.
- **Gaussian-type orbitals (GTOs)**: $\propto x^a y^b z^c e^{-\alpha r^2}$ — less accurate individually, but products of Gaussians are Gaussians, making four-center two-electron integrals analytically tractable.

A **basis set** (e.g., STO-3G, 6-31G*, cc-pVDZ) specifies the number and exponents of GTOs. Larger, more flexible basis sets improve accuracy at the cost of computational time.

---

## Density Matrix

The **one-particle density matrix** (1-PDM) $\gamma(\mathbf{r}, \mathbf{r}')$ is defined as:

$$\gamma(\mathbf{r}, \mathbf{r}') = N \int \Psi^*(\mathbf{r}, \mathbf{x}_2,\ldots)\,\Psi(\mathbf{r}', \mathbf{x}_2,\ldots)\,d\mathbf{x}_2\cdots d\mathbf{x}_N$$

Its diagonal $\rho(\mathbf{r}) = \gamma(\mathbf{r},\mathbf{r})$ is the **electron density**. In a basis representation, the density matrix element is $P_{\mu\nu} = \sum_i n_i\, c_{\mu i}^* c_{\nu i}$, where $n_i$ are occupation numbers. The density matrix encodes all one-electron properties (kinetic energy, electron density, dipole moment) and is central to DFT and population analyses.

---

## Potential Energy Surface (PES)

The **potential energy surface (PES)** is the electronic energy $E_{\text{el}}(\{R_A\})$ as a function of all nuclear coordinates under the Born–Oppenheimer approximation. It is a hypersurface in $3N_{\text{nuc}} - 6$ (or $-5$ for linear molecules) internal degrees of freedom.

Key features of a PES include:

- **Minima** — stable molecular geometries (reactants, products, intermediates).
- **Transition states (saddle points)** — first-order saddle points connecting minima along a reaction path.
- **Dissociation channels** — regions where bonds break.

Geometry optimization, reaction mechanism exploration, and molecular dynamics all operate on the PES.

---

## Ground vs. Excited States

The **ground state** is the lowest-energy eigenstate of the Hamiltonian; all higher-energy eigenstates are **excited states**. Formally:

$$E_0 \leq E_1 \leq E_2 \leq \cdots$$

The ground state wave function $\Psi_0$ determines equilibrium properties under normal conditions. Excited states are relevant to spectroscopy (UV/Vis absorption, fluorescence), photochemistry, and non-adiabatic dynamics. Methods targeting excited states include Configuration Interaction Singles (CIS), Time-Dependent DFT (TD-DFT), EOM-CCSD, and CASSCF/CASPT2.

---

## Open-Shell vs. Closed-Shell

- **Closed-shell**: Every occupied molecular orbital contains exactly two electrons with opposite spins ($\alpha$ and $\beta$). The total spin $S = 0$ (singlet). A single restricted Hartree–Fock (RHF) Slater determinant is a good zeroth-order description.

- **Open-shell**: One or more MOs are singly occupied. Arises in radicals (odd number of electrons), excited states, or some transition-metal complexes. Treated with **Unrestricted HF/DFT (UHF/UDFT)**, where $\alpha$ and $\beta$ orbitals have different spatial parts, or with **Restricted Open-shell HF (ROHF)**. Open-shell systems may suffer from spin contamination in unrestricted methods.

---

## Spin Multiplicity

**Spin multiplicity** $M$ is defined as:

$$M = 2S + 1$$

where $S = \sum_i m_{s,i}$ is the total electron spin quantum number ($S = N_\alpha/2 - N_\beta/2$ for $N_\alpha$ spin-up and $N_\beta$ spin-down electrons). Common multiplicities:

| $S$ | $M$ | Name |
|-----|-----|------|
| 0 | 1 | Singlet |
| 1/2 | 2 | Doublet |
| 1 | 3 | Triplet |
| 3/2 | 4 | Quartet |

Spin multiplicity determines the degeneracy of a state and governs selection rules in spectroscopy (e.g., singlet–triplet transitions are spin-forbidden under pure electric-dipole selection rules but can occur via spin–orbit coupling).

---

## Fermionic Antisymmetry

The **fermionic antisymmetry principle** (a consequence of the spin–statistics theorem) states that the electronic wave function must be **antisymmetric** under the exchange of any two electrons:

$$\Psi(\ldots, \mathbf{x}_i, \ldots, \mathbf{x}_j, \ldots) = -\Psi(\ldots, \mathbf{x}_j, \ldots, \mathbf{x}_i, \ldots)$$

This has two immediate consequences:

1. **Pauli exclusion principle**: two electrons cannot occupy the same spin orbital (setting $\mathbf{x}_i = \mathbf{x}_j$ makes $\Psi = 0$).
2. **Exchange correlation**: electrons with the same spin are kept apart by a purely quantum mechanical "exchange hole," reducing the electron–electron repulsion beyond the classical (Coulomb) contribution.

Slater determinants are the simplest many-body wave functions that satisfy antisymmetry exactly. Bosonic wave functions, by contrast, are symmetric under particle exchange.
