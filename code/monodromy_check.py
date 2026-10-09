"""monodromy_check.py -- monodromy sanity check (numerical; evidence, not proof).

Model: reduced Burau representation at t = zeta_d^r (convention of burau_check.py in the repository HC-quintic-sextic), quotient
E_r = V / C v_t of dimension d - 2; E_r is the monodromy representation on P(r) (up to the identification of r with
gamma, which only permutes labels inside a dual pair).
Gamma_0 = image of the exponent-sum-zero subgroup of the braid group, generated (as a subgroup containing the
generators used here) by g_{k,i} = s_1^k s_i s_1^{-k-1}.  Gamma_0 has determinant 1 on every E_r and preserves the
pairings E_r x E_{-r} -> C, so Gamma_0 lies in the predicted group H (d = 5: SL3 x SL3 on (E1, E4 = E1^*, E2,
E3 = E2^*); d = 6: SL4 x Sp4 on (E2, E4 = E2^*, E3)).  Hence W^H is contained in W^{Gamma_0}, and equality of
dimensions for many tensor representations W is evidence that Gamma_0 is Zariski dense in H.
For each multiset W of labels: dim W^{Gamma_0} (iterated numerical nullspaces) vs dim W^{Lie H} (nullspace of the
Lie algebra action).  Also reports the numerical gap (smallest singular value kept as nonzero / largest treated as 0).
usage: python monodromy_check.py out.json
"""
import itertools, json, sys, time
import numpy as np

TOL = 1e-9


def burau_E(d, r):
    t = np.exp(2j * np.pi * r / d)
    n = d - 1
    Bs = []
    for i in range(n):
        B = np.eye(n, dtype=complex)
        B[i, i] = -t
        if i > 0:
            B[i, i - 1] = t
        if i < n - 1:
            B[i, i + 1] = 1
        Bs.append(B)
    v = np.cumsum([t ** j for j in range(n)])
    for B in Bs:
        assert np.allclose(B @ v, v)
    As = []
    for B in Bs:
        A = np.zeros((n - 1, n - 1), dtype=complex)
        for i in range(1, n):
            for j in range(1, n):
                A[i - 1, j - 1] = B[i, j] - B[0, j] * v[i]
        As.append(A)
    return As


def gamma0_gens(d, labels, K):
    """list of tuples (one matrix per label in `labels`) for g_{k,i}, k in range(-K, K), i = 2..n"""
    mats = {r: burau_E(d, r) for r in labels}
    n = d - 1
    gens = []
    for k in range(-K, K):
        for i in range(1, n):
            g = {}
            for r in labels:
                s1 = mats[r][0]
                P = np.linalg.matrix_power(s1, k) if k >= 0 else np.linalg.matrix_power(np.linalg.inv(s1), -k)
                g[r] = P @ mats[r][i] @ np.linalg.inv(s1) @ np.linalg.inv(P)
            gens.append(g)
    return gens


def kron_all(ms):
    out = np.ones((1, 1), dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def fixed_dim(gens, word):
    D = int(np.prod([gens[0][r].shape[0] for r in word]))
    Q = np.eye(D, dtype=complex)
    gap_lo, gap_hi = np.inf, 0.0
    for g in gens:
        if Q.shape[1] == 0:
            break
        M = kron_all([g[r] for r in word])
        A = M @ Q - Q
        _, s, vh = np.linalg.svd(A, full_matrices=True)
        s_full = np.zeros(Q.shape[1])
        s_full[:len(s)] = s / max(1.0, np.linalg.norm(M, 2))      # relative to the size of the operator
        keep = s_full < TOL
        if (~keep).any():
            gap_lo = min(gap_lo, s_full[~keep].min())
        if keep.any():
            gap_hi = max(gap_hi, s_full[keep].max())
        Q = Q @ vh.conj().T[:, keep]
        Q, _ = np.linalg.qr(Q)
    return Q.shape[1], gap_lo, gap_hi


def sl_basis(n):
    B = []
    for i in range(n):
        for j in range(n):
            if i != j:
                E = np.zeros((n, n)); E[i, j] = 1; B.append(E)
    for i in range(n - 1):
        H = np.zeros((n, n)); H[i, i] = 1; H[i + 1, i + 1] = -1; B.append(H)
    return B


def sp4_basis():
    J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    B = []
    for i in range(4):
        for j in range(i, 4):
            S = np.zeros((4, 4)); S[i, j] = S[j, i] = 1
            B.append(np.linalg.solve(J, S))      # X = J^{-1} S satisfies X^T J + J X = 0
    return B


def lie_reps(d):
    """label -> list over Lie generators of the matrix (0 if the factor does not act)"""
    if d == 5:
        A, Bb = sl_basis(3), sl_basis(3)
        gens = [('A', X) for X in A] + [('B', X) for X in Bb]
        act = {1: ('A', False), 4: ('A', True), 2: ('B', False), 3: ('B', True)}
        dim = 3
    else:
        A, Bb = sl_basis(4), sp4_basis()
        gens = [('A', X) for X in A] + [('B', X) for X in Bb]
        act = {2: ('A', False), 4: ('A', True), 3: ('B', False)}
        dim = 4
    reps = {}
    for r, (grp, dual) in act.items():
        reps[r] = [(X if not dual else -X.T) if g == grp else np.zeros((dim, dim)) for g, X in gens]
    return reps, len(gens), dim


def lie_inv_dim(d, word):
    reps, ng, dim = lie_reps(d)
    D = dim ** len(word)
    I = np.eye(dim)
    tot = np.zeros((D, D))
    for a in range(ng):
        op = np.zeros((D, D))
        for pos, r in enumerate(word):
            X = reps[r][a]
            if not X.any():
                continue
            op += kron_all([X if q == pos else I for q in range(len(word))]).real
        tot += op.T @ op
    w = np.linalg.eigvalsh(tot)
    return int((w < 1e-8).sum())


def main():
    out = []
    for d, labels, maxdeg, K in ((5, (1, 2, 3, 4), 6, 1), (6, (2, 3, 4), 5, 1)):
        gens = gamma0_gens(d, labels, K)
        # sanity: Gamma_0 generators have det 1 on each E_r
        assert all(abs(np.linalg.det(g[r]) - 1) < 1e-9 for g in gens for r in labels)
        for deg in range(1, maxdeg + 1):
            for word in itertools.combinations_with_replacement(labels, deg):
                t0 = time.time()
                fd, glo, ghi = fixed_dim(gens, word)
                ld = lie_inv_dim(d, word)
                res = dict(d=d, word=list(word), gamma0_inv=fd, lie_inv=ld, agree=fd == ld,
                           gap_nonzero_min=float(glo), gap_zero_max=float(ghi), sec=round(time.time() - t0, 2))
                out.append(res)
                print(json.dumps(res), flush=True)
    bad = [r for r in out if not r['agree']]
    print("cases", len(out), "disagreements", len(bad),
          "min nonzero sv", min(r['gap_nonzero_min'] for r in out),
          "max zero sv", max(r['gap_zero_max'] for r in out))
    with open(sys.argv[1], 'w') as fh:
        json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
