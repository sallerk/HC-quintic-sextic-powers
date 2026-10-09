"""counts2.py -- Hodge classes on products of split hypersurfaces with blocks from very general forms G_1, ..., G_r
(d = 5, 6), second implementation, written independently of i2_invariants.py.

Method: numerical linear algebra. For each block label G, the identity component M_G of the monodromy closure
(d=5: SL3 x SL3 on P(1)=std_A, P(4)=dual_A, P(2)=std_B, P(3)=dual_B; d=6: SL4 x Sp4 on P(2)=std, P(4)=dual,
P(3)=std of Sp4; trivial on P(0) (dim d-1), and for d=6 on P(1), P(5)) acts diagonally on all copies of P_G.
Invariants = common kernel of the Lie algebra generators acting on the tensor product (SVD). Their Hodge types are
read off from a Hodge torus acting diagonally in a basis adapted to the Hodge splitting (P(g): first <g>-1 basis
vectors of type (1,0); dual reps: dual basis; Sp4: Lagrangian (1,0)-part), via the ranks of the weight-space
projections of the invariant subspace. The sector model V_c ~ (x)_s P(c_s) (x) L_c, L_c a line of type
(m - sigma, n0 + sigma - 1), and the class h^m (type (m,m)) are as in the notes; a Galois orbit of sector tuples
contributes rational Hodge classes iff every member has type (n,n).

Checks: (1) counts vs i2_invariants.py; (2) single Hodge type of every invariant space; (3) d=5: p - q = (4/5) sum_dets(<g> - 5/2)
for every invariant tuple, and for Hodge ones the det multiset satisfies n_g = n_{-g}; d=6: p - q = (n5 - n1)/3 +
(2/3)(D4 - D2); (4) repeated-block examples.
"""
import itertools, json, sys, time
from functools import lru_cache
import numpy as np

TOL = 1e-8


def sl_basis(n):
    B = []
    for i in range(n):
        for j in range(n):
            if i != j:
                X = np.zeros((n, n)); X[i, j] = 1; B.append(X)
    for i in range(n - 1):
        X = np.zeros((n, n)); X[i, i] = 1; X[i + 1, i + 1] = -1; B.append(X)
    return B


def sp4_basis():
    J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    # X in sp4 iff X^T J + J X = 0; parametrize X = J^{-1} S with S symmetric
    B = []
    for i in range(4):
        for j in range(i, 4):
            S = np.zeros((4, 4)); S[i, j] = 1; S[j, i] = 1
            X = np.linalg.solve(J, S)
            assert np.allclose(X.T @ J + J @ X, 0)
            B.append(X)
    return B


SL3, SL4, SP4 = sl_basis(3), sl_basis(4), sp4_basis()


def rep_table(d):
    """entry -> (group factor 'A'/'B'/None, kind 'std'/'dual'/None, dim, p-weights of basis)."""
    if d == 5:
        return {0: (None, None, 4, [0] * 4),
                1: ('A', 'std', 3, [0, 0, 0]), 4: ('A', 'dual', 3, [1, 1, 1]),
                2: ('B', 'std', 3, [1, 0, 0]), 3: ('B', 'dual', 3, [0, 1, 1])}
    if d == 6:
        return {0: (None, None, 5, [0] * 5), 1: (None, None, 4, [0] * 4), 5: (None, None, 4, [1] * 4),
                2: ('A', 'std', 4, [1, 0, 0, 0]), 4: ('A', 'dual', 4, [0, 1, 1, 1]),
                3: ('B', 'std', 4, [1, 1, 0, 0])}
    raise ValueError


def group_gens(d):
    if d == 5:
        return {'A': SL3, 'B': SL3}
    return {'A': SL4, 'B': SP4}


@lru_cache(maxsize=None)
def inv_info(d, entries):
    """entries: sorted tuple of nontrivial entries (one block label). Returns (dim, {pweight: mult})."""
    T = rep_table(d)
    if not entries:
        return 1, {0: 1}
    dims = [T[e][2] for e in entries]
    N = int(np.prod(dims))
    gens = group_gens(d)
    mats = []
    for fac in ('A', 'B'):
        for X in gens[fac]:
            M = np.zeros((N, N))
            for pos, e in enumerate(entries):
                g, kind, dim, _ = T[e]
                if g != fac:
                    continue
                Y = X if kind == 'std' else -X.T
                ops = [np.eye(dd) for dd in dims]
                ops[pos] = Y
                K = ops[0]
                for o in ops[1:]:
                    K = np.kron(K, o)
                M += K
            if np.any(M):
                mats.append(M)
    if mats:
        A = np.vstack(mats)
        _, sv, Vt = np.linalg.svd(A)
        rank = int((sv > TOL * max(1.0, sv[0])).sum())
        Bnull = Vt[rank:].T  # N x dim
    else:
        Bnull = np.eye(N)
    dim = Bnull.shape[1]
    if dim == 0:
        return 0, {}
    w = np.zeros(N, dtype=int)
    for idx, multi in enumerate(itertools.product(*[range(dd) for dd in dims])):
        w[idx] = sum(T[e][3][i] for e, i in zip(entries, multi))
    mult = {}
    for val in sorted(set(w.tolist())):
        sub = Bnull[w == val, :]
        r = np.linalg.matrix_rank(sub, tol=1e-7) if sub.size else 0
        if r:
            mult[val] = r
    assert sum(mult.values()) == dim, ("torus-unstable invariants", entries, mult, dim)
    return dim, mult


def sectors(d, n):
    return [c for c in itertools.product(range(d), repeat=n) if sum(c) % d == 0]


def units(d):
    return [u for u in range(1, d) if np.gcd(u, d) == 1]


@lru_cache(maxsize=None)
def label_info(d, ent):
    """ent: sorted tuple of all entries of one block label. Returns (dim, weights tuple, dets tuple, n1, n5, D2, D4)."""
    T = rep_table(d)
    dim, w0 = 1, 0
    n1 = n5 = 0
    nontriv = []
    for e in ent:
        g, kind, dd, pw = T[e]
        if g is None:
            dim *= dd
            w0 += pw[0]
            n1 += (d == 6 and e == 1)
            n5 += (d == 6 and e == 5)
        else:
            nontriv.append(e)
    idim, mult = inv_info(d, tuple(nontriv))
    if idim == 0:
        return (0, (), (), 0, 0, 0, 0)
    dets = []
    D2 = D4 = 0
    if d == 5:
        eA = nontriv.count(1) - nontriv.count(4)
        eB = nontriv.count(2) - nontriv.count(3)
        assert eA % 3 == 0 and eB % 3 == 0
        dets += [1] * (eA // 3) if eA > 0 else [4] * (-eA // 3)
        dets += [2] * (eB // 3) if eB > 0 else [3] * (-eB // 3)
    else:
        eA = nontriv.count(2) - nontriv.count(4)
        assert eA % 4 == 0
        D2, D4 = max(eA // 4, 0), max(-eA // 4, 0)
        dets = [2] * D2 + [4] * D4
    weights = tuple(sorted((w0 + k, v) for k, v in mult.items()))
    return (dim * idim, weights, tuple(dets), n1, n5, D2, D4)


def count(d, factors, label=None, max_tuples=None):
    """factors: list of block-label lists. Returns dict of totals."""
    ms = [len(f) - 1 for f in factors]
    n = sum(ms)
    labels = sorted(set(l for f in factors for l in f))
    secs = [[None] + sectors(d, len(f)) for f in factors]  # None = h^m
    # per factor and sector: L-weight and entries per label
    pre = []
    for f, m, sl in zip(factors, ms, secs):
        rows = []
        for c in sl:
            if c is None:
                rows.append((m, tuple(() for _ in labels)))
            else:
                rows.append((m - sum(c) // d, tuple(tuple(c[s] for s, l in enumerate(f) if l == lab) for lab in labels)))
        pre.append(rows)
    # Galois: index maps for u = 2 (d = 5); d = 6 needs none (conjugation)
    umaps = []
    if d == 5:
        for sl in secs:
            pos = {c: i for i, c in enumerate(sl)}
            umaps.append([0 if c is None else pos[tuple((2 * x) % d for x in c)] for c in sl])

    def info(idx):
        dim = 1
        p = 0
        single = True
        dets = []
        n1 = n5 = D2 = D4 = 0
        for i, j in enumerate(idx):
            p += pre[i][j][0]
        for li in range(len(labels)):
            ent = []
            for i, j in enumerate(idx):
                ent.extend(pre[i][j][1][li])
            r = label_info(d, tuple(sorted(ent)))
            if r[0] == 0:
                return 0, None, None, None
            dim *= r[0]
            if len(r[1]) != 1:
                single = False
            p += r[1][0][0]
            dets.extend(r[2])
            n1 += r[3]; n5 += r[4]; D2 += r[5]; D4 += r[6]
        return dim, (p if single else None), dets, (n1, n5, D2, D4)

    res = {"finite_orbit_dim": 0, "hodge_total": 0, "hodge_PAIR": 0, "hodge_CMFIN": 0, "hodge_EXC": 0,
           "mixed_type": 0, "formula_fail": 0, "decomp_fail": 0, "exc_det_multisets": {}}
    cnt = 0
    for idx in itertools.product(*[range(len(sl)) for sl in secs]):
        cnt += 1
        dim, p, dets, aux = info(idx)
        if dim == 0:
            continue
        res["finite_orbit_dim"] += dim
        if p is None:
            res["mixed_type"] += 1
            continue
        pq = 2 * p - 2 * n
        if d == 5:
            if 5 * pq != sum(4 * e - 10 for e in dets):
                res["formula_fail"] += 1
        else:
            n1, n5, D2, D4 = aux
            if 3 * pq != (n5 - n1) + 2 * (D4 - D2):
                res["formula_fail"] += 1
        if p != n:
            continue
        if d == 5:
            _, pu, _, _ = info(tuple(umaps[i][j] for i, j in enumerate(idx)))
            if pu != n:
                continue
        res["hodge_total"] += dim
        if dets:
            res["hodge_EXC"] += dim
            key = str(sorted(dets))
            res["exc_det_multisets"][key] = res["exc_det_multisets"].get(key, 0) + dim
            if d == 5 and (dets.count(1) != dets.count(4) or dets.count(2) != dets.count(3)):
                res["decomp_fail"] += 1
        elif d == 6 and (aux[0] or aux[1]):
            res["hodge_CMFIN"] += dim
        else:
            res["hodge_PAIR"] += dim
    res["tuples"] = cnt
    return res


def main():
    out = {}
    jobs = []
    for m in (1, 2, 3):
        for k in (1, 2, 3, 4):
            if m == 3 and k == 4:
                continue
            jobs.append((5, m, k))
    for m in (1, 2):
        for k in (1, 2, 3, 4):
            jobs.append((6, m, k))
    for k in (1, 2, 3):
        jobs.append((6, 3, k))
    which = sys.argv[2] if len(sys.argv) > 2 else "all"
    if which in ("all", "powers"):
        for d, m, k in jobs:
            t = time.time()
            f = [list(range(m + 1))] * k
            r = count(d, f)
            r["sec"] = round(time.time() - t, 1)
            out[f"d={d},m={m},k={k}"] = r
            print(f"d={d},m={m},k={k}", {kk: vv for kk, vv in r.items() if kk != 'exc_det_multisets'},
                  r["exc_det_multisets"] if len(r["exc_det_multisets"]) < 8 else len(r["exc_det_multisets"]))
            sys.stdout.flush()
    if which in ("all", "repeated"):
        rep = {
            "d=5 X(G,G,G,G',G',G')": (5, [[0, 0, 0, 1, 1, 1]]),
            "d=5 X(G,G,G',G')^2": (5, [[0, 0, 1, 1], [0, 0, 1, 1]]),
            "d=5 X(G,G,G)^2": (5, [[0, 0, 0], [0, 0, 0]]),
            "d=5 X(G,G,G',G'',G''')": (5, [[0, 0, 1, 2, 3]]),
            "d=6 X(G,G,G,G,G',G'')": (6, [[0, 0, 0, 0, 1, 2]]),
            "d=6 X(G,G,G,G')": (6, [[0, 0, 0, 1]]),
        }
        for name, (d, f) in rep.items():
            t = time.time()
            r = count(d, f)
            r["sec"] = round(time.time() - t, 1)
            out[name] = r
            print(name, {kk: vv for kk, vv in r.items() if kk != 'exc_det_multisets'}, r["exc_det_multisets"])
            sys.stdout.flush()
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
