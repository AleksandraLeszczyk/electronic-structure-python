# Pitfalls of Electronic Structure Methods

## Spin Contamination

In unrestricted Hartree–Fock (UHF) and unrestricted DFT (UDFT), the wavefunction is not required to be an eigenfunction of the total spin operator $\hat{S}^2$. The computed $\langle S^2 \rangle$ therefore deviates from the ideal value $S(S+1)$, mixing in contributions from higher-multiplicity states. A doublet ($S = 1/2$) should give $\langle S^2 \rangle = 0.75$; values significantly above this indicate contamination from quartet and higher states.

Spin contamination inflates energies, distorts geometries, and corrupts spin-density distributions, making it particularly dangerous for radical reaction barriers, magnetic coupling constants, and broken-symmetry calculations of open-shell singlets. Diagnostics include comparing $\langle S^2 \rangle$ before and after annihilation, and checking whether ROHF or multireference methods give qualitatively different results.

## Basis Set Superposition Error (BSSE)

When two fragments A and B interact, each fragment can "borrow" basis functions from the other, artificially lowering its energy in the complex relative to the isolated monomer calculation. This spurious stabilisation is basis set superposition error. It appears in any property computed as a difference — binding energies, interaction energies, adsorption enthalpies — and is most severe with small or diffuse basis sets.

The standard remedy is the counterpoise (CP) correction of Boys and Bernardi, in which each monomer energy is evaluated in the full dimer basis set with ghost atoms on the partner's positions. CP correction overcorrects when the basis is too small; the right approach is to use a sufficiently large basis (at least triple-zeta with polarisation and diffuse functions) so that BSSE becomes negligible, and to verify convergence by comparing CP-corrected and uncorrected results.

## SCF Convergence

Self-consistent field iterations can fail to converge, oscillate between two states, or converge to a spurious local minimum rather than the ground-state solution. Common causes include near-degeneracy of frontier orbitals, charge-transfer states, heavy elements with large relativistic effects, and poorly chosen initial guesses.

Standard remedies include level shifting (raising virtual orbital energies to damp oscillations), damping of the Fock matrix update, and direct inversion in the iterative subspace (DIIS). For difficult cases, quadratic convergence (QC-SCF) or orbital-optimised methods may be necessary. It is good practice to check the final orbital occupation, inspect the HOMO–LUMO gap, and confirm that the converged energy is lower than alternative solutions obtained from different starting guesses — for example, from a broken-symmetry or mixed initial guess.

## Multireference Character

Single-reference methods (HF, MP2, CCSD, most DFT functionals) assume the ground state is dominated by a single Slater determinant. When two or more determinants contribute comparably — as in bond breaking, diradicals, excited states, transition-metal complexes, and many open-shell systems — single-reference descriptions are qualitatively wrong and quantitatively unreliable.

Common diagnostics include the $T_1$ amplitude in CCSD (values $> 0.02$ for closed-shell, $> 0.05$ for open-shell systems are warning signs), the $D_1$ diagnostic, and the weight of the leading configuration in a CASSCF calculation. When multireference character is significant, the appropriate tools are CASSCF, CASPT2, NEVPT2, MRCI, or multireference coupled-cluster methods. Applying CCSD(T) to a strongly correlated system can produce errors that dwarf those of the method in the single-reference regime.

## DFT Functional Dependence

Kohn–Sham DFT is formally exact, but the exchange–correlation functional is unknown and must be approximated. Different functional families — LDA, GGA (PBE, BLYP), meta-GGA (TPSS, M06-L), hybrid (B3LYP, PBE0), range-separated hybrid (ωB97X-D, CAM-B3LYP), and double hybrid (B2-PLYP) — make different approximations and perform very differently across chemical problems.

Known systematic failures include the self-interaction error in pure functionals (which causes spurious delocalisation, underestimated barriers, and incorrect charge-transfer excitation energies), the absence of long-range dispersion (corrected by empirical D3/D4 or many-body dispersion terms), and poor treatment of strongly correlated systems. A result obtained with a single functional should not be trusted without validation: at minimum, test one additional functional of a different family and compare. For thermochemistry, barrier heights, non-covalent interactions, and electronic spectra, benchmarked functional–basis set combinations exist and should be preferred over ad hoc choices.
