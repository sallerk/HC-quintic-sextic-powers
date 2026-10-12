# HC-quintic-sextic-powers

Code and data supporting the note *The Hodge conjecture for powers of very general split quintic and sextic
hypersurfaces* (Kevin Saller, 2026), DOI [10.5281/zenodo.23287272](https://doi.org/10.5281/zenodo.23287272). Commits up
to 44e0c59 accompanied an earlier draft of the note, titled *Weil classes on split quintic hypersurfaces and their
powers*.

**Result.** The note considers products Y = X_{F^1} × … × X_{F^k} of hypersurfaces
F_0(u_0,v_0) + … + F_m(u_m,v_m) = 0 in P^(2m+1) whose binary forms of degree d are taken, with repetitions, from very
general forms G_1, …, G_r; for instance the powers X_F^k. For d = 5 and d = 6 it proves the Hodge conjecture for all
such products. For d = 5 the main new input is that the Weil classes for Q(ζ_5) on the products J × J' of the
Jacobians of the curves y^5 + G(u,v) = 0 and y^5 + G'(u,v) = 0 are algebraic, for all G and G'. As in work of Schoen,
this follows from the fact that the quotient of the third power of such a curve by a certain finite group is a cyclic
cover of P^3 branched along five hyperplanes, a quotient of the Fermat quintic threefold, and so does not depend on G.
By a degeneration, this can also be deduced from a theorem of Schoen on cyclic covers of curves. For d = 6, the same
argument with the Fermat cubic fourfold shows that the Weil classes for Q(√−3) on the abelian sixfolds A × E × E are
algebraic, where A is the Jacobian of w^3 + G(u,v) = 0 and E is the Fermat cubic curve; the Hodge classes that involve
determinants on four copies of one curve are reduced to these. The note also determines
the identity component of the Zariski closure of the monodromy of the curves y^d + G(u,v) = 0 on their first
cohomology (SL_3 × SL_3 for d = 5, where this is known; SL_4 × Sp_4 on the part with infinite monodromy for d = 6),
and counts the Hodge classes of X_F × X_F for d = 5.

The proofs in the note do not use these computations. [`README.txt`](README.txt) gives the requirements, the rerun
instructions and what each script checks.

**Status:** not yet peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). They were checked with separately written code and in a separate internal review, which were also carried
out with AI assistance. What was not checked is listed in the note and in README.txt.

## Licences

- **Code** (`code/`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files, data and computed outputs in `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
