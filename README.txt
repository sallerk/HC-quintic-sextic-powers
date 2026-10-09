Supporting code and data for the note
"Weil classes on split quintic hypersurfaces and their powers"
==============================================================

Status: draft accompanying a draft note; not peer reviewed.

What the note proves: let X_F = {F_0(u_0,v_0) + ... + F_m(u_m,v_m) = 0} in P^(2m+1), with binary forms
F_s of degree d with distinct zeros, and let Y = X_(F^1) x ... x X_(F^k) be a product whose forms are
taken, with repetitions, from very general forms G_1, ..., G_r (a "product over G" in the note);
n_j(Y) is the number of places where G_j occurs.
- Theorem A: the identity component of the Zariski closure of the monodromy of the curves
  y^d + G(u,v) = 0 on their first cohomology is SL_3 x SL_3 for d = 5, and SL_4 x Sp_4 on the part
  with infinite monodromy for d = 6.
- Theorem B (d = 5): if the Weil classes for Q(zeta_5) in W(G_a,G_b), on the product of the
  Jacobians of the curves y^5 + G_a = 0 and y^5 + G_b = 0, are algebraic for all a != b with
  n_a(Y) >= 3 and n_b(Y) >= 3, then the Hodge conjecture holds for Y; in particular it holds if
  n_j(Y) >= 3 for at most one j. Conversely, if the Hodge conjecture holds for Y and Y has a Hodge
  class of Weil type (a,b), then the classes in W(G_a,G_b) are algebraic.
- Corollary C (d = 5, very general forms): the Hodge conjecture holds for X_F x X_F, and the space
  of Hodge classes in H^{2m}(X_F,Q) (x) H^{2m}(X_F,Q) has dimension 293, 4429, 73013 for m = 1, 2, 3
  (closed formula in the note); it holds for X_F^3 if and only if the classes in W(F_s,F_t) are
  algebraic for all s != t; it holds for the tenfold X_(G,G,G,G',G',G'), and for
  X_(G,G,G',G') x X_(G,G,G',G'), if and only if the classes in W(G,G') are algebraic.
- Theorem D (d = 6): the Hodge conjecture holds for Y if n_j(Y) <= 3 for all j.
As far as the author knows, it is not known whether the classes in W(G,G') are algebraic.

IMPORTANT: the proofs in the note do not use any of these computations. The scripts below compute
the dimensions of the spaces of Hodge classes for very general forms in two independent ways
(i2_invariants.py, counts2.py), check Theorem A numerically (monodromy_check.py) and check the closed
formulas exactly (formulas.py), for finitely many cases; they are described in Section 9
("Computations") of the note. monodromy_check.py is numerical evidence for Theorem A, not a proof.

Licence: the code in code/ is under the MIT licence; the data in outputs/ and this text are under
CC BY 4.0 (see the LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3; numpy for counts2.py and monodromy_check.py; i2_invariants.py and formulas.py use only the
standard library. Tested with Python 3.12.9 and numpy 2.1.3. The scripts do not parallelize
themselves; numpy may start several threads, which can be prevented by setting OMP_NUM_THREADS=1,
OPENBLAS_NUM_THREADS=1 and MKL_NUM_THREADS=1 (this was done for the runtimes below). Runtimes
measured on a desktop PC, one script at a time: i2_invariants.py 5 s, counts2.py 277 s (of which
212 s for X_(G,G,G)^2 and 49 s for d = 6, m = 1, k = 4), monodromy_check.py 167 s, formulas.py
0.1 s.


2. How to run
-------------
Each script writes its result file (.json) to the file named as first argument and prints a summary
to the screen. The reference results are in outputs/ (the .out files there are the screen output).
To reproduce them, from the top folder of this repository:

    mkdir rerun
    cd rerun
    python ../code/i2_invariants.py   i2_invariants.json   > i2_invariants.out
    python ../code/counts2.py         counts2.json         > counts2.out
    python ../code/monodromy_check.py monodromy_check.json > monodromy_check.out
    python ../code/formulas.py counts2.json i2_invariants.json formulas.json > formulas.out
                                    (reads the results of the two counting programs, so run it in
                                     the same folder, after them)

counts2.py takes an optional second argument: "all" (the default, used for outputs/), "powers" (only
the powers X_F^k) or "repeated" (only the six products with repeated forms); formulas.py needs the
"powers" entries.

Then compare with the reference results:

    python ../code/compare_outputs.py ../outputs .

It compares every .json file in outputs/ with the file of the same name in the current folder,
ignores the run times (keys "sec"), and compares floating-point numbers with relative tolerance 1e-6
(absolute 1e-12). Expected last line: "ALL EQUAL". The screen output of counts2.py and
monodromy_check.py contains the run times as well. The singular values "gap_nonzero_min" and
"gap_zero_max" in monodromy_check.json are floating-point numbers and may differ in the last digits
with other numpy/BLAS versions; everything else is integer or boolean.

Reproducibility: the files in outputs/ were produced by the commands above in a fresh folder (Windows,
Python and numpy versions as above). They agree with the results of earlier runs of the same code,
which has since changed only in comments: i2_invariants.json and formulas.json byte for byte,
counts2.json up to the run times "sec", and monodromy_check.json up to the run times and to
differences in "gap_nonzero_min" and "gap_zero_max" below 1e-6 relative (or 1e-12 absolute), with the
same invariant dimensions and "agree" values. The counts of i2_invariants.py and counts2.py agree for
all powers X_F^k computed by both (totals and "det").


3. Files
--------
Notation as in the note: d = degree of the binary forms; for G in U_d, P_G(gamma) is the
zeta^gamma-eigenspace of the automorphism (u:v:y) -> (zeta u : zeta v : y) on H^1 of the curve
y^d + G(u,v) = 0. The group H of Theorem A acts, for d = 5, by SL_3 x SL_3 on P(1), P(4) (dual to
each other) and P(2), P(3) (dual to each other); for d = 6, by SL_4 on P(2), P(4) (dual to each other)
and by Sp_4 on P(3), and trivially on P(1), P(5) (finite monodromy). A sector of H^{2m}(X_F) is
c in (Z/d)^(m+1); in the model of the note, V_c is the tensor product of the P(c_s) (P(0) trivial of
dimension d - 1 and type (0,0)) with a line of a fixed Hodge type. Hodge classes are counted in
H^{2m_1}(X_(F^1),Q) (x) ... (x) H^{2m_k}(X_(F^k),Q). "det" (key "hodge_EXC"): classes in sectors whose
invariants involve a determinant (D(c) nonempty in the note).

code/i2_invariants.py  ->  outputs/i2_invariants.json, outputs/i2_invariants.out
  Supports: Section 9 ("Hodge classes", first table: the dimensions for X_F^k) and Corollary C(i)
  for m <= 3. Exact integer arithmetic, standard library only.
  For d = 5, 6, m = 1, 2, 3 and k <= 4 (k <= 3 for m = 3), computes the dimension of the space of
  Hodge classes in H^{2m}(X_F,Q)^{(x)k} for very general F: for each multiset of k sectors (or the
  class h^m), the dimension of the invariants of H, by the Weyl character formula, and their Hodge
  type; a Galois orbit of sector tuples contributes if and only if every member has type (km,km).
  Also checks that the Hodge numbers of the sector model agree with those obtained from the Jacobian
  ring (d = 5, 6, m <= 3), and a few known invariant dimensions of SL_3, SL_4 and Sp_4.
  Expected: "model_vs_jacobian_mismatches": 0; for each "d=..,m=..,k=.." the sum of "hodge_PAIR",
  "hodge_CMFIN" and "hodge_EXC" is the value in the first table of Section 9 (e.g. 17, 293, 5121,
  90801 for d = 5, m = 1, k = 1, ..., 4), and "hodge_EXC" is the "det" column; for k = 1 these are
  the values 1 + R(d,m) of the companion note on very general sums of binary quintics and sextics
  (https://doi.org/10.5281/zenodo.23267886).

code/counts2.py  ->  outputs/counts2.json, outputs/counts2.out
  Supports: both tables of Section 9 ("Hodge classes"), Lemma 4.2 and Lemma 4.3. Written
  independently of i2_invariants.py; numpy, floating point.
  Computes the same dimensions as common kernels of the Lie algebra actions (singular value
  decomposition, tolerance 1e-8), and reads off their Hodge types from a torus acting compatibly with
  the Hodge decompositions, for the powers X_F^k (same range as i2_invariants.py) and for the six
  products with repeated forms of the second table (d = 5: X_(G,G,G,G',G',G'), X_(G,G,G',G')^2,
  X_(G,G,G)^2, X_(G,G,G',G'',G'''); d = 6: X_(G,G,G,G,G',G''), X_(G,G,G,G')). Checks: every space
  of invariants has a single Hodge type ("mixed_type"); the formula of Lemma 4.2 for the Hodge type
  ("formula_fail"); for d = 5, that the determinant content of every Hodge orbit is symmetric
  (Lemma 4.3; "decomp_fail").
  Expected: "mixed_type", "formula_fail" and "decomp_fail" are 0 in every entry; "hodge_total" and
  "hodge_EXC" are the values and the "det" column of the two tables of Section 9 (e.g. 12549 and 4 for
  X_(G,G,G,G',G',G')); for the powers, the same totals and "det" values as i2_invariants.py.

code/monodromy_check.py  ->  outputs/monodromy_check.json, outputs/monodromy_check.out
  Supports: Theorem A (Section 3), as numerical evidence only; numpy, floating point.
  Builds the reduced Burau matrices at t = zeta_d^r (the convention of burau_check.py in the
  repository HC-quintic-sextic) on the quotient by the invariant vector, and the elements
  s_1^k s_i s_1^(-k-1), k in {-1, 0}, i = 2, ..., d-1, of the braid group, which have exponent sum
  zero; their images lie in the group H of Theorem A. For every multiset of
  1 to 6 labels in {1, 2, 3, 4} (d = 5) and of 1 to 5 labels in {2, 3, 4} (d = 6), it compares the
  dimension of the invariants of these elements on the tensor product of the corresponding
  eigenspaces with the dimension of the invariants of the Lie algebra of H, and records the
  numerical gap (smallest singular value treated as nonzero, largest treated as zero, relative to the
  operator norm; threshold 1e-9).
  Expected: 264 entries, all with "agree": true; last screen line "cases 264 disagreements 0
  min nonzero sv 2.58e-07 max zero sv 2.4e-13" (rounded here; the last digits vary between runs),
  i.e. a gap of more than five orders of magnitude.

code/formulas.py  ->  outputs/formulas.json, outputs/formulas.out
  Supports: Lemma 4.3, the formula for N_k and Corollary C(i). Exact integer arithmetic, standard
  library only.
  Checks: (G) for every multiset D of elements of (Z/5)^x of size at most 10,
  sum_{g in D} (2<ug> - 5) = 0 for all units u if and only if g and -g have the same multiplicity in
  D for every g (the equivalence in Lemma 4.3); (N) N_k = #{x in ((Z/5)^x)^k : sum x = 0} equals
  (4^k + 4(-1)^k)/5 for k <= 8; (P) the closed formula of Corollary C(i) against "hodge_total" in
  counts2.json and the totals in i2_invariants.json, for m = 1, 2, 3.
  Expected: "galois_multisets": 1000, "galois_failures": 0, "N_ok": true, and in "XxX" the values
  293, 4429, 73013 (m = 1, 2, 3) with "ok": true.

code/compare_outputs.py
  Compares the .json files of a rerun with those in outputs/ (see Section 2). Standard library only.
  Against the earlier runs mentioned under "Reproducibility" it reports "equal" for every file.


4. What the checks do not show
------------------------------
- The counts cover finitely many cases: d = 5, 6, m <= 3, k <= 4 (k <= 3 for m = 3) for the powers,
  and six products with repeated forms. The statements for all m and all products are proved in the
  note.
- Both counting programs take as input the monodromy groups of Theorem A, Lemma 2.1 (Hodge classes of
  very general members have finite monodromy orbits, by Cattani-Deligne-Kaplan) and the sector model
  (1) of the note; they do not test these inputs, except that i2_invariants.py compares the Hodge
  numbers of the model with those of the Jacobian ring.
- They count Hodge classes; nothing here concerns their algebraicity, or the Weil classes W(G,G').
- monodromy_check.py is numerical evidence for Theorem A, not a proof: it uses floating-point ranks
  with a threshold, finitely many braid group elements and finitely many tensor representations.
  Agreement of the invariant dimensions for these representations is consistent with the elements
  generating a Zariski dense subgroup of H, but does not prove it.
- counts2.py also uses floating-point ranks; i2_invariants.py is exact. They agree in all cases
  computed by both.
