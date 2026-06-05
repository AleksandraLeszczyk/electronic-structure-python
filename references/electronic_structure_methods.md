# Electronic Structure Methods: Scaling, Accuracy, and Use Cases

A reference guide to quantum chemistry and condensed-matter electronic structure methods, organized from mean-field to highly correlated approaches.

---

## Notation

| Symbol | Meaning |
|--------|---------|
| N | Number of basis functions (∝ system size) |
| M | Number of active orbitals |
| e | Number of active electrons |
| O(·) | Formal computational scaling |

"Accuracy" refers to the typical error in thermochemical quantities (bond energies, reaction energies) relative to exact non-relativistic results.

---

## Wavefunction-Based Methods

| Method | Full Name | Formal Scaling | Accuracy | Best Use Cases | Original Reference |
|--------|-----------|---------------|----------|----------------|-------------------|
| **HF** | Hartree–Fock | O(N⁴) | ~100–300 kJ/mol | Reference state; qualitative bonding; orbital generation | Hartree, *Proc. Cambridge Phil. Soc.* **24**, 89 (1928); Fock, *Z. Phys.* **61**, 126 (1930) |
| **MP2** | Møller–Plesset 2nd order | O(N⁵) | ~10–40 kJ/mol | Weak interactions, geometries of closed-shell molecules; fast correlation correction | Møller & Plesset, *Phys. Rev.* **46**, 618 (1934) |
| **MP3** | Møller–Plesset 3rd order | O(N⁶) | ~15–40 kJ/mol | Rarely used alone; intermediate step to MP4 | Møller & Plesset, *Phys. Rev.* **46**, 618 (1934) |
| **MP4** | Møller–Plesset 4th order | O(N⁷) | ~5–15 kJ/mol | Higher accuracy single-reference; expensive | Krishnan & Pople, *Int. J. Quantum Chem.* **14**, 91 (1978) |
| **CISD** | Configuration Interaction Singles & Doubles | O(N⁶) | ~10–30 kJ/mol | Benchmark on small systems; not size-consistent | Shavitt, in *Modern Theoretical Chemistry* vol. 3, Plenum (1977) |
| **CCSD** | Coupled Cluster Singles & Doubles | O(N⁶) | ~5–10 kJ/mol | Accurate single-reference ground states; size-consistent | Čížek, *J. Chem. Phys.* **45**, 4256 (1966); Purvis & Bartlett, *J. Chem. Phys.* **76**, 1910 (1982) |
| **pCCD** | pair Coupled Cluster Doubles (= AP1roG) | O(N³) (pair eqs.); O(N⁵) with orbital optimization | Captures static correlation; poor dynamic correlation (similar to CASSCF in quality) | Strongly correlated systems; bond breaking; large active spaces inaccessible to CASSCF; geminal-based reference for post-pCCD corrections | Limacher, Ayers, Johnson, De Baerdemacker, Van Neck & Bultinck, *J. Chem. Theory Comput.* **9**, 1394 (2013); Boguslawski, Tecmer, Ayers, Bultinck, De Baerdemacker & Van Neck, *Phys. Rev. B* **89**, 201106(R) (2014) |
| **fpCCSD** | frozen-pair CCSD | O(N⁶) | Better than pCCD alone; approaches CCSD quality for dynamic correlation on top of a geminal reference | Multireference problems where pCCD orbitals/pair amplitudes are frozen and residual SD correlations recovered; improved over pCCD+linearized corrections | Boguslawski, *J. Chem. Theory Comput.* **15**, 6481 (2019) |
| **CCSD(T)** | CCSD with perturbative triples | O(N⁷) | ~1–5 kJ/mol | "Gold standard"; thermochemistry, spectroscopy; closed-shell or mild open-shell | Raghavachari, Trucks, Pople & Head-Gordon, *Chem. Phys. Lett.* **157**, 479 (1989) |
| **CCSDT** | CC Singles, Doubles & Triples (full) | O(N⁸) | <2 kJ/mol | Sub-kJ/mol thermochemistry; benchmark | Noga & Bartlett, *J. Chem. Phys.* **86**, 7041 (1987) |
| **CCSDTQ** | CC through Quadruples | O(N¹⁰) | <1 kJ/mol | Near-exact for small molecules | Kucharski & Bartlett, *Theor. Chim. Acta* **80**, 387 (1991) |
| **EOM-CC2** | Equation-of-Motion CC2 | O(N⁵) | ~0.2–0.5 eV | Large molecules; fast excited-state screening; valence excitations of closed-shell systems | Christiansen, Koch & Jørgensen, *Chem. Phys. Lett.* **243**, 409 (1995) |
| **EOM-CCSD** | Equation-of-Motion CCSD | O(N⁶) | ~0.1–0.3 eV | Excited states (EE), ionization potentials (IP), electron attachment (EA); workhorse for excited-state CC | Stanton & Bartlett, *J. Chem. Phys.* **98**, 7029 (1993) |
| **SF-EOM-CCSD** | Spin-Flip EOM-CCSD | O(N⁶) | ~0.1–0.3 eV; good for spin-state splittings | Diradicals, triradicals, bond-breaking, open-shell ground states; accesses multi-reference targets from a high-spin reference | Krylov, *Chem. Phys. Lett.* **338**, 375 (2001) |
| **EOM-CC3** | Equation-of-Motion CC3 | O(N⁷) | ~0.05–0.1 eV | High-accuracy excited states and ionization energies; benchmark for valence and Rydberg states | Koch, Christiansen, Jørgensen, Sanchez de Merás & Helgaker, *J. Chem. Phys.* **106**, 1808 (1997) |
| **EOM-CCSDT** | Equation-of-Motion CCSDT | O(N⁸) | ~0.02–0.05 eV | Near-exact excited states; benchmark reference for doubly excited and charge-transfer states | Kucharski, Włoch, Musiał & Bartlett, *J. Chem. Phys.* **115**, 8263 (2001) |
| **ADC(2)** | Algebraic Diagrammatic Construction 2nd order | O(N⁵) | ~0.2–0.5 eV | Excited states; cheaper alternative to EOM-CCSD | Schirmer, *Phys. Rev. A* **26**, 2395 (1982) |
| **ADC(3)** | Algebraic Diagrammatic Construction 3rd order | O(N⁶) | ~0.1–0.2 eV | Accurate excited states; valence and core spectra | Schirmer, Trofimov & Stelter, *J. Chem. Phys.* **109**, 4734 (1998) |
| **CASSCF** | Complete Active Space SCF | O(e! / (e/2)!² × N²) | Qualitative–moderate (no dynamic correlation) | Multireference ground & excited states; bond breaking; transition metals | Roos, Taylor & Siegbahn, *Chem. Phys.* **48**, 157 (1980) |
| **CASPT2** | CASSCF + 2nd-order PT | O(N⁵) per state | ~5–15 kJ/mol | Multireference thermochemistry, spectroscopy, photochemistry | Andersson, Malmqvist, Roos, Sadlej & Wolinski, *J. Phys. Chem.* **94**, 5483 (1990) |
| **NEVPT2** | N-Electron Valence PT2 | O(N⁵) | ~5–15 kJ/mol | Multireference; intruder-state-free alternative to CASPT2 | Angeli, Cimiraglia, Evangelisti, Leininger & Malrieu, *J. Chem. Phys.* **114**, 10252 (2001) |
| **MRCI** | Multireference CI | O(N⁶–N⁸) | ~2–10 kJ/mol | Highly accurate multireference; potential energy surfaces | Werner & Knowles, *J. Chem. Phys.* **89**, 5803 (1988) |
| **MRCI+Q** | MRCI + Davidson size-consistency correction | O(N⁶–N⁸) | ~1–5 kJ/mol | Benchmark PES; improved size-consistency | Davidson & Silver, *Chem. Phys. Lett.* **52**, 403 (1977) |
| **DMRG** | Density Matrix Renormalization Group | O(M³ m³) (m = bond dim.) | High for large active spaces | Large active spaces (>30 orbitals); strongly correlated; 1D-like topology | White, *Phys. Rev. Lett.* **69**, 2863 (1992) |
| **FCI** | Full Configuration Interaction | O(N! / ((N/2)!)²) | Exact (within basis) | Benchmark for small molecules; exact within a given basis set | Knowles & Handy, *Chem. Phys. Lett.* **111**, 315 (1984) |
| **FCIQMC** | FCI Quantum Monte Carlo | Stochastic; sub-exponential | Near-exact | Exact-quality energies for medium molecules (10–30 electrons) | Booth, Thom & Alavi, *J. Chem. Phys.* **131**, 054106 (2009) |
| **VMC** | Variational Monte Carlo | O(N³–N⁴) per sample | Moderate–high | Flexible trial wavefunctions; metallic systems | McMillan, *Phys. Rev.* **138**, A442 (1965) |
| **DMC** | Diffusion Monte Carlo | O(N³–N⁴) per sample | High (fixed-node error ~1–5 kJ/mol) | Large systems; solids; benchmark quality beyond CCSD(T) range | Anderson, *J. Chem. Phys.* **63**, 1499 (1975) |
| **AFQMC** | Auxiliary-Field QMC | O(N³) per step | High | Solids, strongly correlated; controllable sign problem with trial WF | Zhang & Krakauer, *Phys. Rev. Lett.* **90**, 136401 (2003) |

---

## Density Functional Theory (DFT)

| Method | Full Name | Formal Scaling | Accuracy (thermochem.) | Best Use Cases | Original Reference |
|--------|-----------|---------------|------------------------|----------------|-------------------|
| **LDA** | Local Density Approximation | O(N³) | ~50–150 kJ/mol | Metals, solids where error cancellation helps; qualitative | Kohn & Sham, *Phys. Rev.* **140**, A1133 (1965); exchange: Dirac, *Proc. Cambridge Phil. Soc.* **26**, 376 (1930) |
| **GGA (PBE)** | Generalized Gradient Approx. – Perdew–Burke–Ernzerhof | O(N³) | ~20–50 kJ/mol | General-purpose solids, surfaces, large molecules | Perdew, Burke & Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996) |
| **GGA (BLYP)** | Becke 88 exchange + Lee–Yang–Parr correlation | O(N³) | ~20–50 kJ/mol | Organic molecules; geometry optimization | Becke, *Phys. Rev. A* **38**, 3098 (1988); Lee, Yang & Parr, *Phys. Rev. B* **37**, 785 (1988) |
| **meta-GGA (TPSS)** | Tao–Perdew–Staroverov–Scuseria | O(N³) | ~10–30 kJ/mol | Improved over GGA; thermochemistry, lattice constants | Tao, Perdew, Staroverov & Scuseria, *Phys. Rev. Lett.* **91**, 146401 (2003) |
| **meta-GGA (r²SCAN)** | Regularized SCAN | O(N³) | ~10–25 kJ/mol | Balanced accuracy/cost for solids and molecules | Furness, Kaplan, Ning, Perdew & Sun, *J. Phys. Chem. Lett.* **11**, 8208 (2020) |
| **Hybrid (B3LYP)** | Becke 3-parameter + LYP | O(N⁴) | ~10–20 kJ/mol | Organic thermochemistry; most-cited DFT functional | Becke, *J. Chem. Phys.* **98**, 5648 (1993) |
| **Hybrid (PBE0)** | PBE + 25% exact exchange | O(N⁴) | ~10–20 kJ/mol | General-purpose; molecules and solids | Ernzerhof & Scuseria, *J. Chem. Phys.* **110**, 5029 (1999) |
| **Hybrid (HSE06)** | Heyd–Scuseria–Ernzerhof range-separated | O(N⁴) | ~10–20 kJ/mol | Semiconductors and insulators; band gaps | Heyd, Scuseria & Ernzerhof, *J. Chem. Phys.* **118**, 8207 (2003) |
| **Double hybrid (B2-PLYP)** | B2-PLYP | O(N⁵) | ~2–8 kJ/mol | High-accuracy thermochemistry and kinetics | Grimme, *J. Chem. Phys.* **124**, 034108 (2006) |
| **DFT-D3** | DFT + empirical dispersion | O(N³) + O(N²) correction | Significantly improved for non-covalent | Van der Waals complexes, molecular crystals, large biomolecules | Grimme, Antony, Ehrlich & Krieg, *J. Chem. Phys.* **132**, 154104 (2010) |
| **DFT+U** | DFT + Hubbard U correction | O(N³) | Improved for localized d/f electrons | Transition-metal oxides; Mott insulators; f-electron systems | Anisimov, Aryasetiawan & Lichtenstein, *J. Phys.: Condens. Matter* **9**, 767 (1997) |
| **TDDFT** | Time-Dependent DFT | O(N³–N⁴) per state | ~0.2–0.5 eV (excitations) | UV/vis spectra of large molecules; excited-state dynamics | Runge & Gross, *Phys. Rev. Lett.* **52**, 997 (1984) |

---

## Many-Body Perturbation Theory & Green's Function Methods

| Method | Full Name | Formal Scaling | Accuracy | Best Use Cases | Original Reference |
|--------|-----------|---------------|----------|----------------|-------------------|
| **GW** | GW approximation (Hedin's equations) | O(N⁴) (conventional) | ~0.1–0.3 eV (quasiparticle energies) | Band gaps and quasiparticle spectra of solids and molecules | Hedin, *Phys. Rev.* **139**, A796 (1965) |
| **G₀W₀** | Single-shot GW on top of DFT | O(N⁴) | ~0.1–0.3 eV | Band gaps of semiconductors/insulators | Hybertsen & Louie, *Phys. Rev. B* **34**, 5390 (1986) |
| **scGW** | Self-consistent GW | O(N⁴–N⁵) | ~0.1 eV | More robust quasiparticle energies; avoids starting-point dependence | Schilfgaarde, Kotani & Faleev, *Phys. Rev. Lett.* **96**, 226402 (2006) |
| **BSE** | Bethe–Salpeter Equation | O(N⁶) | ~0.1–0.2 eV (optical gaps) | Optical spectra, excitons in solids; charge-transfer excitations | Salpeter & Bethe, *Phys. Rev.* **84**, 1232 (1951); Strinati, *Riv. Nuovo Cimento* **11**, 1 (1988) |
| **DMFT** | Dynamical Mean-Field Theory | O(N³) + impurity solver | Qualitative–high for local correlations | Strongly correlated metals, Mott transitions, heavy-fermion systems | Georges, Kotliar, Krauth & Rozenberg, *Rev. Mod. Phys.* **68**, 13 (1996) |

---

## Semiempirical and Low-Cost Methods

| Method | Full Name | Formal Scaling | Accuracy | Best Use Cases | Original Reference |
|--------|-----------|---------------|----------|----------------|-------------------|
| **Extended Hückel** | Extended Hückel Theory | O(N³) | Qualitative only | Orbital topology, crystal orbital analysis, pedagogical | Hoffmann, *J. Chem. Phys.* **39**, 1397 (1963) |
| **AM1** | Austin Model 1 | O(N³) | ~20–50 kJ/mol | Fast organic thermochemistry; large molecules | Dewar, Zoebisch, Healy & Stewart, *J. Am. Chem. Soc.* **107**, 3902 (1985) |
| **PM3** | Parametric Model 3 | O(N³) | ~15–40 kJ/mol | Organic and some inorganic systems | Stewart, *J. Comput. Chem.* **10**, 209 (1989) |
| **PM7** | Parametric Model 7 | O(N³) | ~10–25 kJ/mol | Large biomolecules, drug-like molecules; best general semiempirical | Stewart, *J. Mol. Model.* **19**, 1 (2013) |
| **GFN2-xTB** | Geometry, Frequency, Non-covalent, extended TB (2nd gen.) | O(N²–N³) | ~15–30 kJ/mol | Conformational sampling, molecular dynamics, large systems (1000s of atoms) | Bannwarth, Ehlert & Grimme, *J. Chem. Theory Comput.* **15**, 1652 (2019) |
| **DFTB** | Density Functional Tight Binding | O(N²–N³) | ~20–50 kJ/mol | Biological systems, carbon allotropes, large-scale MD | Porezag, Frauenheim, Köhler, Seifert & Kaschner, *Phys. Rev. B* **51**, 12947 (1995) |
| **DFTB3** | Third-order DFTB | O(N²–N³) | ~10–25 kJ/mol | Improved for H-bonded systems, charged molecules | Gaus, Cui & Elstner, *J. Chem. Theory Comput.* **7**, 931 (2011) |

---

## Linear-Scaling and Embedding Methods

| Method | Full Name | Formal Scaling | Accuracy | Best Use Cases | Original Reference |
|--------|-----------|---------------|----------|----------------|-------------------|
| **O(N) DFT** | Linear-scaling DFT (density-matrix methods) | O(N) | Same as underlying functional | Large-scale biomolecular DFT (10⁴–10⁶ atoms) | Goedecker, *Rev. Mod. Phys.* **71**, 1085 (1999) |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics | O(N³) for QM region | High in QM region, empirical elsewhere | Enzyme active sites, solvent effects, large biological systems | Warshel & Levitt, *J. Mol. Biol.* **103**, 227 (1976) |
| **ONIOM** | Our own N-layered Integrated Molecular Orbital and Molecular Mechanics | O(N³) for inner layer | Inherits accuracy of high-level method | Multilayer embedding for large molecules | Maseras & Morokuma, *J. Comput. Chem.* **16**, 1170 (1995) |
| **DMET** | Density Matrix Embedding Theory | O(N) outer + O(N³) per fragment | High for fragments | Strongly correlated solids and molecules; fragment-based exact embedding | Knizia & Chan, *Phys. Rev. Lett.* **109**, 186404 (2012) |

---

## Summary: Accuracy vs. Cost Trade-off

```
Increasing accuracy →
│
│  FCI / FCIQMC / CCSDTQ
│  CCSDT
│  CCSD(T)  ←── "gold standard"
│  CCSD / CASPT2 / MRCI
│  MP2 / EOM-CCSD
│  Hybrid DFT (B3LYP, PBE0)
│  GGA DFT (PBE, BLYP)
│  LDA
│  Semiempirical (GFN2-xTB, PM7)
│  Extended Hückel / HF alone
│
└──────────────────────────────── Increasing system size →
```

### Rules of Thumb

- **Small molecules, high accuracy**: CCSD(T)/aug-cc-pVTZ or better.
- **Organic thermochemistry**: B3LYP-D3/6-311+G(d,p) or B2-PLYP-D3.
- **Transition metals, multireference character**: CASPT2 or NEVPT2.
- **Periodic solids, band gaps**: HSE06 or G₀W₀.
- **Strongly correlated**: DMFT, DMET, or DMRG.
- **Large biomolecules / MD**: GFN2-xTB, DFTB3, or QM/MM.
- **Optical spectra**: TDDFT (large systems) or BSE@GW (solids, charge transfer).

---

*Last updated: June 2026*
