# HC-quintic-sextic-powers

Code and data supporting the note *Weil classes on split quintic hypersurfaces and their powers* (Kevin Saller, 2026;
draft).

**Result.** The note considers products Y = X_{F^1} × … × X_{F^k} of hypersurfaces
F_0(u_0,v_0) + … + F_m(u_m,v_m) = 0 in P^(2m+1) whose binary forms of degree d are taken, with repetitions, from very
general forms G_1, …, G_r; for instance the powers X_F^k. It determines the identity component of the Zariski closure
of the monodromy of the curves y^d + G(u,v) = 0 on their first cohomology (SL_3 × SL_3 for d = 5; SL_4 × Sp_4 on the
part with infinite monodromy for d = 6), and uses classical invariant theory to describe the Hodge classes of Y. For
d = 5, every Hodge class of Y is algebraic provided that the Weil classes for Q(ζ_5) on the products J_a × J_b of the
Jacobians of the curves y^5 + G_a = 0 and y^5 + G_b = 0 are algebraic (as far as the author knows, it is not known
whether they are). The Hodge conjecture for Y holds unconditionally when at most one of the forms occurs three or more times, for
instance for X_F × X_F, whose Hodge classes the note counts. Conversely, the Hodge conjecture for X_F^3, for the
tenfolds X_(G,G,G,G',G',G') or for X_(G,G,G',G')^2 would imply that these Weil classes are algebraic. For d = 6 the
Hodge conjecture holds for Y when every form occurs at most three times.

The proofs in the note do not use these computations. [`README.txt`](README.txt) gives the requirements, the rerun
instructions and what each script checks.

**Status:** draft, not yet peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). They were checked with separately written code and in a separate internal review, which were also carried
out with AI assistance. What was not checked is listed in the note and in README.txt.

## Licences

- **Code** (`code/`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files, data and computed outputs in `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
