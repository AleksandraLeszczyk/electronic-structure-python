# Basis Sets in Quantum Chemistry

A basis set is a collection of mathematical functions (typically Gaussian-type orbitals, GTOs) used to represent molecular orbitals. The choice of basis set is one of the most consequential decisions in any quantum chemical calculation, controlling both the computational cost and the accuracy of results. Below each basis set family is described with its purpose, design philosophy, and original literature references.

---

## 1. STO-nG — Minimal Basis Sets

**Purpose/Best Use:** The simplest possible starting point; suitable for qualitative surveys of large systems, molecular topology studies, or pedagogical use. The accuracy is too low for most quantitative purposes.

**Design:** Each Slater-type orbital (STO) — which has the correct cusp and long-range decay behaviour — is approximated by a fixed linear combination of *n* primitive Gaussian-type functions, whose exponents and coefficients are optimised to minimise the deviation from the true STO. The most common variant is STO-3G (n = 3). Only one set of basis functions per occupied atomic shell is included (minimal basis).

**Key limitations:** No polarisation or diffuse functions; poor description of ionic or anion species; basis set superposition error (BSSE) is large.

**Original references:**
- Hehre, W. J.; Stewart, R. F.; Pople, J. A. "Self-Consistent Molecular-Orbital Methods. I. Use of Gaussian Expansions of Slater-Type Atomic Orbitals." *J. Chem. Phys.* **51**, 2657–2664 (1969). DOI: [10.1063/1.1672392](https://doi.org/10.1063/1.1672392)
- For STO-6G and heavier atoms: Hehre, W. J.; Ditchfield, R.; Stewart, R. F.; Pople, J. A. *J. Chem. Phys.* **52**, 2769 (1970). DOI: [10.1063/1.1673374](https://doi.org/10.1063/1.1673374)

---

## 2. Pople Split-Valence Family

The Pople family is a series of basis sets developed throughout the 1970s–1990s that use a split representation of the valence shell. The notation *X*-*YZG* means: core orbitals are each contracted from *X* primitives; the valence shell is split into two (or three) parts contracted from *Y* (or *Y* and *Z*) primitives. Polarisation and diffuse functions are added as separate designators.

### 2.1 Double-Zeta Valence — 3-21G

**Purpose:** Inexpensive geometry optimisations and frequency calculations for organic/main-group molecules; often used as a starting geometry for subsequent refinement.

**Design:** The valence shell is split into two contractions (inner and outer) drawn from 3-primitive and 2-primitive or 1-primitive sets, respectively. Core orbitals are each contracted from 3 primitives. Exponents and coefficients were optimised at the Hartree–Fock level.

**Original reference:**
- Binkley, J. S.; Pople, J. A.; Hehre, W. J. "Self-Consistent Molecular Orbital Methods. 21. Small Split-Valence Basis Sets for First-Row Elements." *J. Am. Chem. Soc.* **102**, 939–947 (1980). (No DOI assigned; JACS 102(3))

### 2.2 Double-Zeta Valence — 6-31G

**Purpose:** The workhorse Pople double-zeta basis for routine HF/DFT geometry optimisations and single-point energies for organic molecules. Widely benchmarked.

**Design:** Core orbitals contracted from 6 primitives; valence split into a 3-primitive inner and 1-primitive outer contraction. Exponents optimised at HF for ground-state atoms.

**Original reference:**
- Hehre, W. J.; Ditchfield, R.; Pople, J. A. "Self-Consistent Molecular Orbital Methods. XII. Further Extensions of Gaussian-Type Basis Sets for Use in Molecular Orbital Studies of Organic Molecules." *J. Chem. Phys.* **56**, 2257–2261 (1972). DOI: [10.1063/1.1677527](https://doi.org/10.1063/1.1677527)

### 2.3 Triple-Zeta Valence — 6-311G

**Purpose:** Better valence description than 6-31G; useful when accurate energetics or properties are needed without the expense of Dunning cc sets.

**Design:** Valence split into three contractions; exponents were optimised at the MP2 level for H–Ar, providing a basis with improved correlation recovery compared to double-zeta alternatives. Note that despite the name, 6-311G is not a true triple-zeta basis in the Dunning sense and has known deficiencies for properties sensitive to the outer valence region.

**Original reference:**
- Krishnan, R.; Binkley, J. S.; Seeger, R.; Pople, J. A. "Self-Consistent Molecular Orbital Methods. XX. A Basis Set for Correlated Wave Functions." *J. Chem. Phys.* **72**, 650–654 (1980). DOI: [10.1063/1.438955](https://doi.org/10.1063/1.438955)

### 2.4 Polarisation Functions: \* and \*\*

**Purpose:** Adding d-functions on heavy atoms (* notation) and d + p functions on hydrogen (**) captures polarisation of the electron density in bonds and is essential for accurate geometries, dipole moments, and conformational energies.

**Design:** Single sets of d-type Gaussian primitives on first-row heavy atoms (C, N, O, F) and optionally p-type functions on H, optimised by minimising energy corrections to hydrogenation energies. The exponents are a single uncontracted Gaussian, e.g. 6-31G(d) adds one d function on each heavy atom.

**Original reference:**
- Hariharan, P. C.; Pople, J. A. "The Influence of Polarization Functions on Molecular Orbital Hydrogenation Energies." *Theor. Chim. Acta* **28**, 213–222 (1973). DOI: [10.1007/BF00533485](https://doi.org/10.1007/BF00533485)

### 2.5 Diffuse Functions: + and ++

**Purpose:** Diffuse functions extend the basis set into regions far from the nucleus. They are required for anions, excited states, Rydberg states, long-range interactions (van der Waals, hydrogen bonds), and properties such as polarisabilities.

**Design:** A single set of s and p diffuse Gaussians with small exponents is added to each heavy atom (+) and optionally also to hydrogen (++). Exponents were optimised for proton affinities and electron affinities of anions.

**Original reference:**
- Clark, T.; Chandrasekhar, J.; Spitznagel, G. W.; Schleyer, P. v. R. "Efficient Diffuse Function-Augmented Basis Sets for Anion Calculations. III. The 3-21+G Basis Set for First-Row Elements, Li–F." *J. Comput. Chem.* **4**, 294–301 (1983). DOI: [10.1002/jcc.540040303](https://doi.org/10.1002/jcc.540040303)

---

## 3. Karlsruhe / Ahlrichs Family

Developed at the Karlsruhe Institute of Technology, these basis sets emphasise a balanced quality-to-cost ratio, making them the default choice in many DFT-oriented codes (e.g. TURBOMOLE, ORCA).

### 3.1 def-SV, def-SVP, def-TZV, def-TZVP (original Ahlrichs sets)

**Purpose:** The first generation of Karlsruhe basis sets; balanced split-valence (SV/SVP) and triple-zeta valence (TZV/TZVP) for Li–Kr. "P" denotes the inclusion of polarisation functions.

**Design:** Exponents and contractions were fully optimised via energy minimisation at the HF level. The design goal was that basis set errors should be comparable across all elements, rather than achieving high accuracy for specific atoms only.

**Original references:**
- Schäfer, A.; Horn, H.; Ahlrichs, R. "Fully Optimized Contracted Gaussian Basis Sets for Atoms Li to Kr." *J. Chem. Phys.* **97**, 2571–2577 (1992). DOI: [10.1063/1.463096](https://doi.org/10.1063/1.463096)
- Schäfer, A.; Huber, C.; Ahlrichs, R. (1994) extended the series to triple-zeta valence quality: *J. Chem. Phys.* **100**, 5829–5835. DOI: [10.1063/1.467146](https://doi.org/10.1063/1.467146)

### 3.2 def2-SVP, def2-TZVP, def2-TZVPP, def2-QZVP, def2-QZVPP

**Purpose:** The second-generation Karlsruhe family, covering the entire periodic table H–Rn (with relativistic ECPs for Rb–Rn). This is the most commonly recommended starting point for DFT calculations today. def2-SVP is suitable for pre-optimisations; def2-TZVP and def2-TZVPP are the standard for production DFT energetics; def2-QZVP and def2-QZVPP approach the basis set limit.

**Design:** The "def2" (default 2) sets are a systematic redesign for consistent accuracy across the whole periodic table. For 5th-row elements and beyond, small-core relativistic effective core potentials (ECPs) replace inner-shell electrons, while valence basis functions are re-optimised to match properties. The quality descriptors are SVP (split-valence + polarisation), TZVP (triple-zeta valence + polarisation), TZVPP (triple-zeta valence + two sets of polarisation), QZVP/QZVPP (quadruple-zeta).

**Original reference:**
- Weigend, F.; Ahlrichs, R. "Balanced Basis Sets of Split Valence, Triple Zeta Valence and Quadruple Zeta Valence Quality for H to Rn: Design and Assessment of Accuracy." *Phys. Chem. Chem. Phys.* **7**, 3297–3305 (2005). DOI: [10.1039/b508541a](https://doi.org/10.1039/b508541a)

### 3.3 def2-SVPD, def2-TZVPD, def2-TZVPPD, def2-QZVPPD (diffuse-augmented)

**Purpose:** Property-optimised variants of the def2 family with additional diffuse functions, designed for response properties such as polarisabilities, NMR chemical shifts, and excited states. They use fewer diffuse functions than Dunning aug-cc sets, making them cheaper while still describing the outer regions well.

**Design:** One or two diffuse functions per angular momentum per atom are added to each def2 set and optimised for atomic and molecular polarisabilities. The philosophy mirrors Dunning augmentation but applies it to the Ahlrichs segmented-contracted framework.

**Original reference:**
- Rappoport, D.; Furche, F. "Property-Optimized Gaussian Basis Sets for Molecular Response Calculations." *J. Chem. Phys.* **133**, 134105 (2010). DOI: [10.1063/1.3484283](https://doi.org/10.1063/1.3484283)

---

## 4. Dunning Correlation-Consistent Family (cc)

Designed explicitly for systematically convergent correlated calculations. Each step up in the cardinal number *X* (D→T→Q→5→6) adds one more shell of correlating functions per atom, allowing extrapolation to the complete basis set (CBS) limit.

### 4.1 cc-pVXZ (X = D, T, Q, 5, 6)

**Purpose:** The primary tool for wave-function-based correlated calculations (CCSD(T), MRCI, MP2) where systematic CBS extrapolation is needed. Optimised for valence correlation of main-group atoms.

**Design:** For each angular momentum shell, functions are added in groups that recover similar amounts of correlation energy. Starting from an atomic HF reference, CI natural orbitals provide the initial functions; exponents and contractions are then re-optimised for correlated atomic energies. The "cc" designation means the sets are "correlation consistent."

**Original references:**
- Dunning, T. H., Jr. "Gaussian Basis Sets for Use in Correlated Molecular Calculations. I. The Atoms Boron through Neon and Hydrogen." *J. Chem. Phys.* **90**, 1007–1023 (1989). DOI: [10.1063/1.456153](https://doi.org/10.1063/1.456153)
- For Al–Ar: Woon, D. E.; Dunning, T. H., Jr. "Gaussian Basis Sets for Use in Correlated Molecular Calculations. III. The Atoms Aluminum through Argon." *J. Chem. Phys.* **98**, 1358–1371 (1993). DOI: [10.1063/1.464303](https://doi.org/10.1063/1.464303)

### 4.2 aug-cc-pVXZ — Augmented Correlation-Consistent

**Purpose:** Extend the cc-pVXZ sets with diffuse functions to describe anions, excited states, Rydberg series, non-covalent interactions, polarisabilities, and other properties sensitive to the outer density. Essential for accurate electron affinities.

**Design:** One diffuse function of each angular-momentum type present in cc-pVXZ is added per atom, with exponents optimised for the atomic electron affinity. The even-tempered extension keeps the ratio of successive exponents constant.

**Original reference:**
- Kendall, R. A.; Dunning, T. H., Jr.; Harrison, R. J. "Electron Affinities of the First-Row Atoms Revisited. Systematic Basis Sets and Wave Functions." *J. Chem. Phys.* **96**, 6796–6806 (1992). DOI: [10.1063/1.462569](https://doi.org/10.1063/1.462569)

### 4.3 cc-pCVXZ — Core-Valence Correlation

**Purpose:** Core–valence and core–core correlation effects on molecular properties (geometries, harmonic frequencies, bond dissociation energies) that are non-negligible at high accuracy targets (e.g. CCSD(T)/CBS + core). Required when inner-shell electrons must be explicitly correlated.

**Design:** Additional tight (high-exponent) functions are appended to cc-pVXZ to describe the contraction and polarisation of core orbitals under correlation. The extra functions are optimised for core-correlation energies of atoms.

**Original reference:**
- Woon, D. E.; Dunning, T. H., Jr. "Gaussian Basis Sets for Use in Correlated Molecular Calculations. V. Core-Valence Basis Sets for Boron through Neon." *J. Chem. Phys.* **103**, 4572–4585 (1995). DOI: [10.1063/1.470645](https://doi.org/10.1063/1.470645)

### 4.4 cc-pwCVXZ — Weighted Core-Valence

**Purpose:** A more balanced alternative to cc-pCVXZ, weighting core–valence relative to core–core contributions so that valence properties (geometry, binding energies) converge more smoothly with cardinal number. Preferred over cc-pCVXZ for molecular property calculations.

**Design:** The "w" (weighted) sets optimise the additional tight functions against a weighted combination of core–core and core–valence energies, reducing unnecessary expansion in core–core space.

**Original reference:**
- Peterson, K. A.; Dunning, T. H., Jr. "Accurate Correlation Consistent Basis Sets for Molecular Core–Valence Correlation Effects. The Second Row Atoms Al–Ar and the First Row Atoms B–Ne Revisited." *J. Chem. Phys.* **117**, 10548–10560 (2002). DOI: [10.1063/1.1520138](https://doi.org/10.1063/1.1520138)

---

## 5. ANO — Atomic Natural Orbital Basis Sets

### 5.1 Almlöf–Taylor ANO

**Purpose:** Highly accurate all-purpose sets for correlated calculations; well-suited when compact basis sets with efficient correlation recovery are needed.

**Design:** A large primitive Gaussian set is contracted using the eigenvectors (natural orbitals) of the one-particle density matrix averaged over several atomic states (neutral, ions) at the CISD level. The resulting ANOs provide a near-optimal compact contraction: truncating to fewer contracted functions loses less accuracy than with conventional contraction schemes.

**Original reference:**
- Almlöf, J.; Taylor, P. R. "General Contraction of Gaussian Basis Sets. I. Atomic Natural Orbitals for First- and Second-Row Atoms." *J. Chem. Phys.* **86**, 4070–4077 (1987). DOI: [10.1063/1.451917](https://doi.org/10.1063/1.451917)

### 5.2 Roos ANO (Widmark–Malmqvist–Roos)

**Purpose:** Widely used in multireference calculations (CASSCF/CASPT2, MRCI) in MOLCAS/OpenMolcas, where the ability to systematically enlarge the contracted basis without changing primitive exponents is a key advantage.

**Design:** Primitive exponents span a broad range covering both core and valence regions. The contraction coefficients are obtained from a density matrix averaged over the atomic ground state, low-lying excited states, positive and negative ions, and atoms in external electric fields. This averaging guarantees balanced description of ionisation, electron affinity, and polarisability in a single contracted set.

**Original reference:**
- Widmark, P.-O.; Malmqvist, P.-Å.; Roos, B. O. "Density Matrix Averaged Atomic Natural Orbital (ANO) Basis Sets for Correlated Molecular Wave Functions." *Theor. Chim. Acta* **77**, 291–306 (1990). DOI: [10.1007/BF01120130](https://doi.org/10.1007/BF01120130)

---

## 6. Jensen Polarisation-Consistent Family (pc-*n*)

**Purpose:** Optimised for Hartree–Fock and density functional theory calculations where systematic convergence of DFT exchange–correlation energies and properties (rather than MP2/CC correlation) is desired. The complement to cc-pVXZ for DFT-oriented work.

**Design:** Unlike Dunning's sets (where the criterion for adding functions is recovery of correlation energy), the pc-*n* sets add functions by angular-momentum type in decreasing order of importance as estimated from atomic DFT calculations. Each increment in *n* (0 through 4) adds the next most important angular-momentum shell. The primitive exponent and contraction optimisations target DFT exchange–correlation integrals rather than correlated wave-function energies.

**Variants:** aug-pc-*n* adds diffuse functions for anions/excited states; pcS-*n* (property-optimised for NMR shielding); pcJ-*n* (NMR spin–spin coupling).

**Original references:**
- Jensen, F. "Polarization Consistent Basis Sets: Principles." *J. Chem. Phys.* **115**, 9113–9125 (2001). DOI: [10.1063/1.1413524](https://doi.org/10.1063/1.1413524)
- Jensen, F. "Polarization Consistent Basis Sets. III. The Importance of Diffuse Functions." *J. Chem. Phys.* **117**, 9234–9240 (2002). DOI: [10.1063/1.1515484](https://doi.org/10.1063/1.1515484)

---

## 7. Sapporo Basis Sets

**Purpose:** All-electron segmented-contracted basis sets of double-, triple-, and quadruple-zeta quality for elements H–Xe and beyond, with relativistic (Douglas–Kroll–Hess) contracted variants. Useful as an alternative to the Ahlrichs or Dunning families for all-electron relativistic work.

**Design:** Valence-correlating functions are derived from natural orbitals of non-relativistic (or scalar-relativistic DKH2/DKH3) CI calculations. Core and core-valence correlating functions are added by decontracting inner shells. Segmented contraction (each primitive appears in at most one contracted function) is enforced for computational efficiency in conventional integral evaluation codes.

**Original references:**
- Noro, T.; Sekiya, M.; Koga, T. "Segmented Contracted Basis Sets for Atoms H through Xe: Sapporo-(DK)-nZP Sets (n = D, T, Q)." *Theor. Chem. Acc.* **131**, 1124 (2012). DOI: [10.1007/s00214-012-1124-z](https://doi.org/10.1007/s00214-012-1124-z)

---

## 8. IGLO Basis Sets (NMR / Magnetic Properties)

**Purpose:** Specifically constructed for the IGLO (Individual Gauge for Localized Orbitals) method for computing NMR chemical shielding constants, magnetic susceptibilities, and related properties. Used whenever high-quality NMR shielding tensors are needed within an affordable HF or DFT framework.

**Design:** Unlike energy-optimised basis sets, IGLO sets include additional tight p-functions and more diffuse s/p functions whose exponents were optimised to reproduce NMR chemical shifts of reference molecules. The basis is tailored to minimise the gauge-dependence error in the nuclear shielding computation. IGLO-II and IGLO-III are the most commonly used members; IGLO-III approaches the basis set limit for first- and second-row elements.

**Original reference:**
- Kutzelnigg, W.; Fleischer, U.; Schindler, M. "The IGLO-Method: Ab-Initio Calculation and Interpretation of NMR Chemical Shifts and Magnetic Susceptibilities." In *NMR Basic Principles and Progress* vol. **23**, pp. 165–262 (Springer, Berlin, 1990). DOI: [10.1007/978-3-642-75932-1_3](https://doi.org/10.1007/978-3-642-75932-1_3)

---

## 9. Stuttgart–Dresden Effective Core Potentials (Stuttgart ECPs)

**Purpose:** Replace chemically inert core electrons of heavy atoms with a smooth analytical potential, dramatically reducing the number of explicitly treated electrons. This makes calculations on 4d, 5d, 6d transition metals and post-d main-group elements (In–Rn, actinides) computationally tractable. The accompanying valence basis sets are optimised for use with the ECP.

**Design:** Energy-consistent pseudopotentials (ECPs) are derived by adjusting the ECP parameters until the valence orbital energies and spectroscopic constants of the atom reproduce those from an all-electron scalar-relativistic (Douglas–Kroll) calculation. Separate sets of spin–orbit ECPs are available. Small-core ECPs treat more shells explicitly and are more accurate; large-core ECPs are faster but less reliable for properties sensitive to semicore electrons.

**Original references:**
- Dolg, M.; Wedig, U.; Stoll, H.; Preuss, H. "Energy-Adjusted Ab Initio Pseudopotentials for the First Row Transition Elements." *J. Chem. Phys.* **86**, 866–872 (1987). DOI: [10.1063/1.452288](https://doi.org/10.1063/1.452288)
- Bergner, A.; Dolg, M.; Küchle, W.; Stoll, H.; Preuß, H. "Ab Initio Energy-Adjusted Pseudopotentials for Elements of Groups 13–17." *Mol. Phys.* **80**, 1431–1441 (1993). DOI: [10.1080/00268979300103121](https://doi.org/10.1080/00268979300103121)

---

## 10. Dyall All-Electron Relativistic Basis Sets

**Purpose:** Highly accurate all-electron basis sets for four-component relativistic calculations (Dirac–Coulomb Hamiltonian) and scalar-relativistic methods. Designed for the DIRAC program and for high-accuracy work on heavy elements where neither ECPs nor two-component approximations are desired.

**Design:** Primitive exponents and contractions are optimised using four-component Dirac–Fock calculations on the atom. The sets are hierarchical (dyall.cv2z, cv3z, cv4z) and include core-correlating functions. They are generally contracted, meaning each primitive can appear in multiple contracted functions, providing maximum flexibility.

**Original references (selected):**
- Dyall, K. G. "Relativistic and Nonrelativistic Finite Nucleus Optimized Double Zeta Basis Sets for the 4p, 5p and 6p Elements." *Theor. Chem. Acc.* **99**, 366–371 (1998). DOI: [10.1007/s002140050337](https://doi.org/10.1007/s002140050337)
- Dyall, K. G. "Relativistic and Nonrelativistic Finite Nucleus Optimized Triple Zeta Basis Sets for the 4p, 5p and 6p Elements." *Theor. Chem. Acc.* **108**, 335–340 (2002). DOI: [10.1007/s00214-002-0388-0](https://doi.org/10.1007/s00214-002-0388-0)
- Dyall, K. G. "Relativistic Double-Zeta, Triple-Zeta, and Quadruple-Zeta Basis Sets for the 5d Elements Hf–Hg." *Theor. Chem. Acc.* **112**, 403–409 (2004). DOI: [10.1007/s00214-004-0607-y](https://doi.org/10.1007/s00214-004-0607-y)

---

## 11. SARC — Segmented All-Electron Relativistically Contracted Basis Sets

**Purpose:** All-electron basis sets designed for use with the Douglas–Kroll–Hess (DKH2) and ZORA scalar-relativistic Hamiltonians within codes such as ORCA. Intended as a practical replacement for ECPs in routine DFT/CASSCF work on 3d–5d transition metals, lanthanides, and actinides, offering higher accuracy than ECP approaches for core-sensitive properties (EPR, Mössbauer, X-ray spectroscopy).

**Design:** The primitive exponents are adapted from fully relativistic all-electron calculations. The contraction scheme is segmented (each Gaussian belongs to one contracted shell), making the integrals efficient to evaluate in conventional programs. The sets exist in SARC-DZP, SARC-TZP levels, tuned to give reliable results comparable to much larger, generally contracted sets.

**Original references:**
- Pantazis, D. A.; Chen, X.-Y.; Landis, C. R.; Neese, F. "All-Electron Scalar Relativistic Basis Sets for Third-Row Transition Metal Atoms." *J. Chem. Theory Comput.* **4**, 908–919 (2008). DOI: [10.1021/ct800047t](https://doi.org/10.1021/ct800047t)
- Pantazis, D. A.; Neese, F. "All-Electron Scalar Relativistic Basis Sets for the Lanthanides." *J. Chem. Theory Comput.* **5**, 2229–2238 (2009). DOI: [10.1021/ct900090f](https://doi.org/10.1021/ct900090f)
- Pantazis, D. A.; Neese, F. "All-Electron Scalar Relativistic Basis Sets for the Actinides." *J. Chem. Theory Comput.* **7**, 677–684 (2011). DOI: [10.1021/ct100736b](https://doi.org/10.1021/ct100736b)

---

## 12. x2c Basis Sets — Exact Two-Component All-Electron Sets

**Purpose:** Counterparts to the def2 Karlsruhe basis sets for calculations using the exact two-component (X2C) relativistic Hamiltonian, covering the full periodic table H–Rn in a fully all-electron framework. Preferred when an ECP-free treatment is needed together with the efficiency of segmented contraction.

**Design:** Primitive exponents and contractions are taken from the def2 sets and re-optimised at the X2C one-electron level with a finite-size Gaussian nuclear charge distribution. This ensures that the basis set errors are balanced ("error-consistent") across the periodic table and consistent with the X2C approximation. Variants x2c-SVPall, x2c-TZVPall, and x2c-TZVPPall map to SVP, TZVP, and TZVPP quality, respectively. The "-s" variants (x2c-TZVPall-s) add extra tight functions for accurate NMR shielding.

**Original reference:**
- Pollak, P.; Weigend, F. "Segmented Contracted Error-Consistent Basis Sets of Double- and Triple-ζ Valence Quality for One- and Two-Component Relativistic All-Electron Calculations." *J. Chem. Theory Comput.* **13**, 3696–3705 (2017). DOI: [10.1021/acs.jctc.7b00593](https://doi.org/10.1021/acs.jctc.7b00593)

---

## Quick-Reference Summary

| Family | Best Use | Elements | Relativistic? |
|---|---|---|---|
| STO-*n*G | Qualitative / pedagogical | H–Kr | No |
| 3-21G / 6-31G / 6-311G | Routine DFT/HF, organic molecules | H–Kr | No |
| def2-SVP/TZVP/QZVP | Standard DFT (all-purpose) | H–Rn | Via ECP (Rb–Rn) |
| def2-TZVPD / def2-QZVPPD | Response properties, NMR, polarisabilities | H–Rn | Via ECP (Rb–Rn) |
| cc-pVXZ | Correlated WF, CBS extrapolation | H–Ar, heavier via ECP | No (main-group) |
| aug-cc-pVXZ | Anions, excited states, EA/PA | H–Ar | No |
| cc-pCVXZ / cc-pwCVXZ | Core–valence correlation | B–Ar | No |
| ANO (Almlöf–Taylor) | Compact correlated WF | H–Ar | No |
| ANO-Roos | Multireference (CASSCF/CASPT2) | H–Xe | No |
| pc-*n* | DFT, systematic DFT convergence | H–Kr | No |
| aug-pc-*n* / pcS-*n* / pcJ-*n* | DFT NMR, excitations | H–Kr | No |
| Sapporo-(DK)-nZP | All-electron DKH2/DKH3 | H–Xe | Scalar (DKH) |
| IGLO-II/III | NMR shielding (HF/DFT) | H–Se | No |
| Stuttgart ECPs | Heavy elements (4d–6d, lanthanides, actinides) | K–Rn + | Via ECP |
| Dyall cv2z/cv3z/cv4z | 4-component / high-accuracy relativistic | All | 4-component |
| SARC-DZP/TZP | Routine ORCA relativistic DFT, EPR, X-ray spectroscopy | 3d–5d, Ln, An | Scalar (DKH/ZORA) |
| x2c-SVPall / TZVPall | All-electron DFT with X2C | H–Rn | Scalar (X2C) |

---

## Notes on Basis Set Selection

**Completeness vs. cost:** Increasing the cardinal number in a systematic family (cc or def2) always improves accuracy but increases cost steeply (~N⁴ for DFT, ~N⁵–N⁷ for correlated methods). A pragmatic strategy is to optimise geometries with a smaller basis (e.g. def2-SVP or 6-31G(d)) and compute energies with a larger one (def2-TZVP, aug-cc-pVTZ).

**Basis set superposition error (BSSE):** Finite basis sets artificially stabilise complexes because each fragment can "borrow" basis functions from the other. The counterpoise correction (Boys and Bernardi) partially mitigates this; using larger basis sets reduces it intrinsically.

**Matching basis set to method:** cc sets were designed for coupled-cluster and perturbation theory; using them with pure DFT is not wrong but may not be the most efficient choice — pc-*n* or def2 sets often converge DFT energies faster. Conversely, Pople sets were not designed for correlated calculations and should not be used for CBS extrapolation.

**Property-specific sets:** NMR shieldings converge slowly with standard energy-optimised bases; IGLO, pcS-*n*, or def2-TZVPD are strongly preferred. For electron affinities and excited states, augmented sets (aug-cc, def2-SVPD, aug-pc-*n*) are essential.
